"""
sequencing.py -- VCF, gVCF, BAM and CRAM ingest for DNAInsight v3.0.

WHY THIS FILE EXISTS
--------------------
Everything DNAInsight has read so far is a consumer array export: a few hundred
thousand rows, one genotype per line, GRCh37 plus strand by convention. Whole
genome and whole exome sequencing arrives instead as VCF, gVCF, BAM or CRAM, and
each of those brings a failure mode the array format simply does not have.

Four rules govern this module, and each exists because the alternative produces
a confidently wrong answer rather than an error:

  1. STREAM, NEVER SLURP. ``parsers.parse_dna_file`` calls ``read_text()`` on the
     whole file. That is fine for a 15 MB array export and fatal for a 30x WGS
     VCF, which is gigabytes before compression. Every read path here is a
     generator over an open handle.

  2. COUNT WHAT YOU DROP. A VCF contains indels, structural variants, spanning
     deletions and gVCF reference blocks, none of which DNAInsight can annotate.
     Skipping them is correct. Skipping them silently is not, because "we found
     nothing at that position" and "we never looked at that position" lead to
     opposite clinical conclusions. Every skipped record lands in a counter with
     a named reason.

  3. REFUSE TO GUESS THE BUILD. DNAInsight's bundled reference is GRCh37 plus
     strand. GRCh37 and GRCh38 coordinates disagree by thousands to millions of
     bases on most chromosomes, so a GRCh38 VCF read as GRCh37 does not fail, it
     reports the wrong variants with full confidence. Contig LENGTHS are the only
     signal in a VCF header that cannot be copied wrong by a pipeline author, so
     they decide. When there is no evidence at all the answer is None, never a
     plausible default. This is the same stance ``orientation.py`` takes on
     palindromic SNPs: flag it, do not guess it.

  4. NEVER PASS COORDINATES THROUGH A LIFTOVER THAT DID NOT HAPPEN. There is no
     liftover in the standard library and no chain-file implementation small
     enough to bundle honestly. If the chain file is absent the caller gets the
     documented "unavailable" payload from ``external.unavailable()``, with the
     input coordinates untouched and NOT relabelled.

DEPENDENCIES
------------
Standard library only. No pysam, no cyvcf2, no numpy. BAM and CRAM are binary
formats whose readers are large and fiddly, so they are not reimplemented here;
they go out through the ``external`` adapter to samtools, which is the same
licence and degradation boundary every other v3.0 capability uses.
"""

from __future__ import annotations

import gzip
import re
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

from . import external
from . import orientation
from .parsers import ParseError

__all__ = [
    "BuildMismatch",
    "ParseError",
    "REFERENCE_BUILD",
    "REFERENCE_STRAND",
    "CONTIG_LENGTHS",
    "INFORMATIVE_CONTIGS",
    "SKIP_REASONS",
    "SAMTOOLS_TOOL_ID",
    "SAMTOOLS_CAPABILITY",
    "LIFTOVER_CAPABILITY",
    "VcfHeader",
    "ChainFile",
    "normalize_contig",
    "detect_format",
    "parse_header_lines",
    "read_header",
    "detect_build",
    "assert_build_compatible",
    "split_gt",
    "alleles_from_gt",
    "is_simple_snv",
    "parse_record",
    "select_sample",
    "iter_vcf_records",
    "read_vcf",
    "chain_dir",
    "find_chain_file",
    "parse_chain_file",
    "liftover_available",
    "liftover",
    "write_positions_bed",
    "parse_mpileup_line",
    "call_genotype_from_counts",
    "extract_positions",
    "parse_sequencing_file",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# What the bundled reference is expressed in. Stated once, imported everywhere,
# so a future rebuild changes one line rather than eleven scattered literals.
REFERENCE_BUILD = "GRCh37"
REFERENCE_STRAND = "plus"

# The tool id this module expects in external.TOOLS.
#
# WIRING PASS: samtools is deliberately NOT in the registry yet, and external.py
# is not edited by this wave. Until the entry lands, external.guard() reports
# state "unknown" and every alignment path here degrades to the standard
# unavailable payload, which is the behaviour we want anyway. Adding the entry
# below to external.TOOLS is the ONLY change needed to switch these paths on:
#
#   "samtools": {
#       "id": "samtools", "name": "SAMtools",
#       "purpose": "Targeted pileup extraction from BAM and CRAM alignments",
#       "capability": "sequencing_pileup",
#       "licence": "MIT License", "spdx": "MIT",
#       "composable": True, "commercial_ok": True, "redistributable": True,
#       "homepage": "https://www.htslib.org/",
#       "binaries": ["samtools", "samtools.exe"],
#       "kind": "binary", "verified": "<date read>", "requires": [],
#   }
SAMTOOLS_TOOL_ID = "samtools"
SAMTOOLS_CAPABILITY = "sequencing_pileup"

LIFTOVER_CAPABILITY = "liftover"

# Chromosome lengths, in bases, for the two builds DNAInsight can distinguish.
#
# These are the load-bearing numbers of this module. A ##reference line is free
# text a pipeline author can copy from the wrong template, and an assembly= tag
# is equally a claim rather than a measurement. A contig length comes from the
# FASTA the caller actually used, so it is the one header field that cannot
# quietly disagree with the coordinates in the body of the file.
CONTIG_LENGTHS: dict[str, dict[str, int]] = {
    "GRCh37": {
        "1":  249250621,
        "2":  243199373,
        "3":  198022430,
        "4":  191154276,
        "5":  180915260,
        "X":  155270560,
        "Y":  59373566,
        "MT": 16569,
    },
    "GRCh38": {
        "1":  248956422,
        "2":  242193529,
        "3":  198295559,
        "4":  190214555,
        "5":  181538259,
        "X":  156040895,
        "Y":  57227415,
        "MT": 16569,
    },
}

# MT is excluded from the vote on purpose: both builds carry the identical 16569
# base rCRS, so an MT length is compatible with everything and distinguishes
# nothing. Keeping it in the table but out of the vote records that the omission
# is a decision and not an oversight.
INFORMATIVE_CONTIGS: frozenset[str] = frozenset({"1", "2", "3", "4", "5", "X", "Y"})

# Fixed key set for the skipped-record counters. The shape is frozen so a caller
# can render every reason without inspecting which ones happened to be non-zero.
SKIP_REASONS: tuple[str, ...] = (
    "indel",
    "symbolic",
    "no_id",
    "multiallelic_complex",
    "no_call",
)

_ACGT = frozenset("ACGT")
_RS_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
_GT_SPLIT_RE = re.compile(r"[|/]")

# Tokens that name a build in free text. Matched as whole-ish tokens because a
# bare "38" appears in filenames that have nothing to do with the assembly.
_BUILD_TOKENS: tuple[tuple[str, str], ...] = (
    ("GRCh38", r"grch38|hg38|hs38|b38|gr38"),
    ("GRCh37", r"grch37|hg19|hs37|b37|g1k_v37|human_g1k"),
)

_CONFIDENCE_HIGH = "high"
_CONFIDENCE_MEDIUM = "medium"
_CONFIDENCE_NONE = "none"
_CONFIDENCE_CONFLICT = "conflict"


class BuildMismatch(ParseError):
    """The file's genome build is not the build DNAInsight can annotate.

    Subclasses ParseError on purpose. Callers already wrap ingest in one
    ``except ParseError``, and a build mismatch is an ingest failure, not a
    separate category of problem. Carrying ``detected`` and ``expected`` lets
    the route layer render the two builds without re-parsing the message.
    """

    def __init__(self, detected: str | None, expected: str, message: str = "") -> None:
        super().__init__(message or _build_mismatch_message(detected, expected))
        self.detected = detected
        self.expected = expected


def _build_mismatch_message(detected: str | None, expected: str) -> str:
    if detected is None:
        return (
            f"Genome build could not be determined and DNAInsight requires {expected}. "
            "The header carried no contig lengths, no assembly tag and no usable "
            "##reference line. Coordinates from the wrong build do not fail, they "
            "match the wrong variants, so this file is refused rather than guessed at. "
            "Re-export with contig lines, or state the build explicitly."
        )
    return (
        f"Genome build mismatch: this file is {detected}, DNAInsight's bundled reference "
        f"is {expected} {REFERENCE_STRAND} strand. {detected} and {expected} coordinates "
        "differ by thousands to millions of bases on most chromosomes, so reading one as "
        "the other produces confident findings at the wrong positions. Lift the file over "
        f"to {expected} first, or supply a {expected} export."
    )


# ---------------------------------------------------------------------------
# Contig naming
# ---------------------------------------------------------------------------

def normalize_contig(name: Any) -> str:
    """Normalise a contig name to the convention the rest of DNAInsight uses.

    Consumer exports and the bundled reference use bare "1".."22", "X", "Y" and
    "MT". A VCF may use any of "1", "chr1", "CHR1", "chrM" or "MT". Normalising
    here is what stops a chr-prefixed WGS VCF from matching exactly zero
    reference positions while looking like it parsed perfectly.

    "23" is deliberately NOT translated to "X". PLINK-derived files mean X by 23,
    but a VCF contig literally named "23" is ambiguous, and this module does not
    guess at identities it cannot verify.
    """
    token = str(name or "").strip()
    if not token:
        return ""
    upper = token.upper()
    if upper.startswith("CHR"):
        upper = upper[3:]
    if upper in ("M", "MT"):
        return "MT"
    return upper


def _chrom_sort_key(chrom: str) -> tuple[int, int, str]:
    """Autosomes numerically, then X, Y, MT, then anything else alphabetically."""
    key = normalize_contig(chrom)
    if key.isdigit():
        return (0, int(key), "")
    order = {"X": 1, "Y": 2, "MT": 3}
    if key in order:
        return (1, order[key], "")
    return (2, 0, key)


# ---------------------------------------------------------------------------
# File opening and format detection
# ---------------------------------------------------------------------------

def _is_gzip(path: Path) -> bool:
    """Detect gzip by magic bytes rather than by extension.

    Users rename files, and a BGZF-compressed VCF is a valid gzip stream whatever
    it is called. Trusting the suffix means a ``.vcf`` that is actually gzipped
    fails with a UnicodeDecodeError that tells the user nothing.
    """
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError as exc:
        raise ParseError(f"Cannot read {path}: {exc}") from exc


def _open_text(path: Path):
    """Open a plain or gzipped text file for streaming line iteration."""
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def detect_format(filepath: str | Path) -> str:
    """Identify the container: "vcf", "vcf.gz", "bam", "cram" or "unknown".

    Decided from magic bytes, not the filename, for the same reason ``_is_gzip``
    is. BAM is a BGZF stream whose first four decompressed bytes are "BAM\\1",
    which is how a compressed BAM is told apart from a compressed VCF.
    """
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")
    try:
        with open(path, "rb") as fh:
            magic = fh.read(4)
    except OSError as exc:
        raise ParseError(f"Cannot read {path}: {exc}") from exc

    if magic[:4] == b"CRAM":
        return "cram"
    if magic[:4] == b"BAM\x01":
        return "bam"
    if magic[:2] == b"\x1f\x8b":
        try:
            with gzip.open(path, "rb") as gz:
                head = gz.read(4)
        except (OSError, EOFError, gzip.BadGzipFile) as exc:
            raise ParseError(
                f"{path.name} begins with the gzip magic bytes but could not be "
                f"decompressed: {exc}. A truncated download is the usual cause."
            ) from exc
        if head[:4] == b"BAM\x01":
            return "bam"
        return "vcf.gz"

    try:
        with open(path, "rt", encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
    except OSError as exc:
        raise ParseError(f"Cannot read {path}: {exc}") from exc
    if first.startswith("##") or first.startswith("#CHROM"):
        return "vcf"
    return "unknown"


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

@dataclass
class VcfHeader:
    """The parts of a VCF header this module actually acts on.

    Deliberately not a full header model. Storing every ##INFO and ##FORMAT
    definition would be memory spent on metadata nothing downstream reads.
    """

    fileformat: str = ""
    reference: str = ""
    contigs: list[dict] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    column_count: int = 0
    gvcf_hint: bool = False
    chr_prefixed: bool = False
    line_count: int = 0

    @property
    def sample_count(self) -> int:
        return len(self.samples)

    def contig_length(self, name: str) -> int | None:
        key = normalize_contig(name)
        for contig in self.contigs:
            if normalize_contig(contig.get("id", "")) == key:
                return contig.get("length")
        return None


def _parse_structured(body: str) -> dict[str, str]:
    """Split the ``<ID=1,length=249250621,assembly=b37>`` body into key/value pairs.

    Written by hand rather than with ``str.split(",")`` because descriptions are
    quoted and legitimately contain commas.
    """
    out: dict[str, str] = {}
    key_chars: list[str] = []
    val_chars: list[str] = []
    in_key = True
    in_quote = False
    for ch in body:
        if in_key:
            if ch == "=":
                in_key = False
            else:
                key_chars.append(ch)
            continue
        if ch == '"':
            in_quote = not in_quote
            continue
        if ch == "," and not in_quote:
            out["".join(key_chars).strip()] = "".join(val_chars).strip()
            key_chars, val_chars = [], []
            in_key = True
            continue
        val_chars.append(ch)
    if key_chars:
        out["".join(key_chars).strip()] = "".join(val_chars).strip()
    return out


def parse_header_lines(lines: Iterable[str]) -> VcfHeader:
    """Build a VcfHeader from an iterable of header lines.

    Pure and side effect free so build detection can be tested without a file on
    disk. Lines may or may not carry trailing newlines.
    """
    header = VcfHeader()
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        header.line_count += 1

        if line.startswith("#CHROM"):
            fields = line.lstrip("#").split("\t")
            header.column_count = len(fields)
            # Columns 0..7 are fixed, column 8 is FORMAT, samples follow.
            if len(fields) > 9:
                header.samples = [f.strip() for f in fields[9:] if f.strip()]
            continue

        if not line.startswith("##"):
            continue

        if line.startswith("##fileformat="):
            header.fileformat = line.split("=", 1)[1].strip()
            continue
        if line.startswith("##reference="):
            header.reference = line.split("=", 1)[1].strip()
            continue
        if line.startswith("##contig=<") and line.rstrip().endswith(">"):
            body = line[len("##contig=<"):].rstrip()[:-1]
            attrs = _parse_structured(body)
            contig_id = attrs.get("ID", "")
            if not contig_id:
                continue
            length: int | None = None
            raw_length = attrs.get("length", "")
            if raw_length.isdigit():
                length = int(raw_length)
            if contig_id.lower().startswith("chr"):
                header.chr_prefixed = True
            header.contigs.append({
                "id": contig_id,
                "length": length,
                "assembly": attrs.get("assembly", ""),
            })
            continue
        # GATK writes ##ALT=<ID=NON_REF,...> and ##GVCFBlock lines only into
        # gVCFs, so either one settles the provider before a record is read.
        if "NON_REF" in line or line.startswith("##GVCFBlock"):
            header.gvcf_hint = True
    return header


def _read_header(handle) -> VcfHeader:
    """Consume header lines from an open handle, stopping after #CHROM.

    The handle is left positioned on the first data line so the caller can keep
    streaming. That is the whole point: the header and the body are read in one
    pass over a file too large to hold in memory.
    """
    collected: list[str] = []
    consumed = 0
    saw_chrom = False
    for line in handle:
        consumed += 1
        if not line.strip():
            continue
        if line.startswith("#"):
            collected.append(line)
            if line.startswith("#CHROM"):
                saw_chrom = True
                break
            continue
        raise ParseError(
            "Malformed VCF: a data line appeared before the #CHROM column header. "
            "Verify this is a VCF and not a headerless position list."
        )
    if not saw_chrom:
        raise ParseError(
            "Malformed VCF: no #CHROM header line was found. Without it the sample "
            "columns cannot be identified and no genotype can be read."
        )
    header = parse_header_lines(collected)
    # Overwritten with the raw count so downstream error messages quote the real
    # line number, blank lines included.
    header.line_count = consumed
    return header


def read_header(filepath: str | Path) -> VcfHeader:
    """Read only the header of a VCF or gVCF, then close the file."""
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")
    with _open_text(path) as handle:
        return _read_header(handle)


# ---------------------------------------------------------------------------
# Build detection
# ---------------------------------------------------------------------------

def _build_from_text(text: str) -> str | None:
    """Name the build a free-text string claims, or None when it names neither.

    A string naming both is treated as naming neither. That happens with paths
    like ``/refs/GRCh37/lifted_from_GRCh38.fa``, where believing either half
    would be a coin flip.
    """
    token = str(text or "").lower()
    if not token:
        return None
    hits = [build for build, pattern in _BUILD_TOKENS if re.search(pattern, token)]
    if len(hits) == 1:
        return hits[0]
    return None


def detect_build(header: VcfHeader) -> dict:
    """Infer the genome build from a parsed header.

    Returns ``{"build": "GRCh37" | "GRCh38" | None, "confidence": str,
    "evidence": [...]}``.

    Evidence is ranked, not pooled, because the signals are not equally
    trustworthy:

      contig length   measured from the FASTA that produced the coordinates.
                      Decides on its own, confidence "high".
      assembly= tag   a claim in the header. Used only when no length is present,
                      confidence "medium".
      ##reference     free text, frequently a stale path. Same tier as above.
      chr prefix      recorded and never voted on. "chr1" tells you the file came
                      through UCSC tooling, which distributes both hg19 and hg38.
                      It is evidence about the distributor, not the assembly, and
                      treating it as a build signal is a known way to get this
                      wrong.

    When lengths from both builds appear the answer is None with confidence
    "conflict": a file like that is either concatenated from two sources or
    hand-edited, and neither is safe to annotate.
    """
    evidence: list[dict] = []
    length_votes: dict[str, int] = {"GRCh37": 0, "GRCh38": 0}

    for contig in header.contigs:
        contig_id = normalize_contig(contig.get("id", ""))
        length = contig.get("length")
        if not contig_id or length is None:
            continue
        if contig_id not in INFORMATIVE_CONTIGS:
            if contig_id == "MT" and length == CONTIG_LENGTHS["GRCh37"]["MT"]:
                evidence.append({
                    "signal": "contig_length",
                    "detail": f"MT length {length} is the rCRS, identical in both builds",
                    "supports": None,
                })
            continue
        matches = [b for b in ("GRCh37", "GRCh38") if CONTIG_LENGTHS[b].get(contig_id) == length]
        if len(matches) == 1:
            length_votes[matches[0]] += 1
            evidence.append({
                "signal": "contig_length",
                "detail": f"contig {contig_id} length {length}",
                "supports": matches[0],
            })
        elif not matches:
            evidence.append({
                "signal": "contig_length",
                "detail": f"contig {contig_id} length {length} matches no known build",
                "supports": None,
            })

    for contig in header.contigs:
        assembly = str(contig.get("assembly") or "").strip()
        if not assembly:
            continue
        evidence.append({
            "signal": "assembly_tag",
            "detail": f"contig {contig.get('id')} assembly={assembly}",
            "supports": _build_from_text(assembly),
        })

    if header.reference:
        evidence.append({
            "signal": "reference_header",
            "detail": f"##reference={header.reference}",
            "supports": _build_from_text(header.reference),
        })

    if header.contigs:
        evidence.append({
            "signal": "contig_naming",
            "detail": "chr-prefixed contig names" if header.chr_prefixed else "bare contig names",
            # Always None. Naming style tracks the distributor, not the assembly.
            "supports": None,
        })

    voted = [b for b, n in length_votes.items() if n]
    if len(voted) > 1:
        return {"build": None, "confidence": _CONFIDENCE_CONFLICT, "evidence": evidence}
    if len(voted) == 1:
        chosen = voted[0]
        claimed = {e["supports"] for e in evidence
                   if e["signal"] in ("assembly_tag", "reference_header") and e["supports"]}
        # Lengths outrank text. A disagreement is still worth reporting, so
        # confidence drops rather than the answer changing.
        confidence = _CONFIDENCE_HIGH
        if claimed and claimed != {chosen}:
            confidence = _CONFIDENCE_MEDIUM
        return {"build": chosen, "confidence": confidence, "evidence": evidence}

    claimed = [e["supports"] for e in evidence
               if e["signal"] in ("assembly_tag", "reference_header") and e["supports"]]
    unique = set(claimed)
    if len(unique) == 1:
        return {"build": unique.pop(), "confidence": _CONFIDENCE_MEDIUM, "evidence": evidence}
    if len(unique) > 1:
        return {"build": None, "confidence": _CONFIDENCE_CONFLICT, "evidence": evidence}

    return {"build": None, "confidence": _CONFIDENCE_NONE, "evidence": evidence}


def assert_build_compatible(detected: str | None,
                            expected: str | None = REFERENCE_BUILD,
                            *,
                            allow_unknown: bool = False) -> None:
    """Raise BuildMismatch unless ``detected`` is the build we can annotate.

    ``expected=None`` means the caller has opted out of the check entirely and
    accepts the consequences. ``allow_unknown=True`` downgrades "no evidence" to
    the caller's problem, which is how ``parse_sequencing_file`` reports an
    undetectable build as a warning instead of destroying an otherwise readable
    file. A build that is detected and WRONG always raises, in both modes: that
    is the case where continuing produces findings at the wrong positions.
    """
    if expected is None:
        return
    if detected == expected:
        return
    if detected is None and allow_unknown:
        return
    raise BuildMismatch(detected, expected)


# ---------------------------------------------------------------------------
# Genotype decoding
# ---------------------------------------------------------------------------

def split_gt(gt: str) -> tuple[list[str], bool]:
    """Split a GT value into its allele index tokens and its phasing flag.

    "0|1" is phased, "0/1" is not. Phasing changes nothing about which alleles a
    person carries, so it is recorded and not acted on. Mixed separators, which
    do occur in polyploid and merged records, count as phased only when every
    separator is a pipe.
    """
    token = str(gt or "").strip()
    if not token:
        return [], False
    parts = [p.strip() for p in _GT_SPLIT_RE.split(token)]
    separators = [c for c in token if c in "|/"]
    phased = bool(separators) and all(c == "|" for c in separators)
    return parts, phased


def is_simple_snv(ref: str, alt: str) -> bool:
    """True when REF and ALT are both single unambiguous bases."""
    return (
        len(ref) == 1 and ref.upper() in _ACGT
        and len(alt) == 1 and alt.upper() in _ACGT
    )


def _is_symbolic(allele: str) -> bool:
    """True for <NON_REF>, <DEL>, breakends and the spanning-deletion star."""
    if not allele:
        return True
    if allele.startswith("<") or allele == "*":
        return True
    return "[" in allele or "]" in allele


def alleles_from_gt(gt: str, ref: str, alts: Sequence[str]) -> dict:
    """Resolve a GT against REF and ALT into a two-allele genotype.

    Returns ``{"allele1", "allele2", "phased", "ploidy", "no_call",
    "partial_no_call", "selected", "indexes"}``.

    Three cases the array parsers never had to handle:

      haploid GT   chrX, chrY and chrM in a male sample are called "1", not
                   "1/1". The allele is duplicated, because 23andMe and Ancestry
                   render male X as a single letter which
                   ``parsers._split_genotype`` also duplicates. One convention
                   downstream beats two.
      missing GT   "./." becomes ("N", "N"), matching ``parsers._split_genotype``
                   on "--". The record is still counted as a no-call by the
                   caller rather than emitted.
      multi-allele GT indexes address the full ALT list, so "1/2" on
                   ALT="A,G" is a heterozygote carrying neither reference copy.
    """
    ref = str(ref or "").upper()
    alt_list = [str(a).upper() for a in alts]
    indexes, phased = split_gt(gt)

    if not indexes:
        return {
            "allele1": "N", "allele2": "N", "phased": False, "ploidy": 0,
            "no_call": True, "partial_no_call": False, "selected": [], "indexes": [],
        }

    resolved: list[str | None] = []
    for token in indexes:
        if token == "." or token == "":
            resolved.append(None)
            continue
        if not token.isdigit():
            raise ParseError(f"Unreadable GT allele index {token!r} in genotype {gt!r}.")
        idx = int(token)
        if idx == 0:
            resolved.append(ref)
            continue
        if idx > len(alt_list):
            raise ParseError(
                f"GT {gt!r} references ALT allele {idx} but only {len(alt_list)} "
                "ALT alleles are present on this record."
            )
        resolved.append(alt_list[idx - 1])

    selected = [a for a in resolved if a is not None]
    no_call = not selected
    partial = bool(selected) and len(selected) != len(resolved)

    if len(resolved) == 1:
        first = resolved[0]
        allele1 = allele2 = first if first is not None else "N"
    else:
        allele1 = resolved[0] if resolved[0] is not None else "N"
        allele2 = resolved[1] if resolved[1] is not None else "N"

    return {
        "allele1": allele1,
        "allele2": allele2,
        "phased": phased,
        "ploidy": len(resolved),
        "no_call": no_call,
        "partial_no_call": partial,
        "selected": selected,
        "indexes": indexes,
    }


def _extract_rsid(id_field: str) -> str | None:
    """Return the first rs identifier in the ID column, or None.

    A VCF ID column may carry several semicolon-separated identifiers, and only
    an rs number is usable: the bundled reference, SNPedia and every genoset are
    keyed by rsid. A COSMIC or internal id is not a partial match, it is no match.
    """
    token = str(id_field or "").strip()
    if not token or token == ".":
        return None
    for part in token.split(";"):
        part = part.strip()
        if _RS_RE.match(part):
            return "rs" + part[2:]
    return None


def parse_record(line: str, *, sample_index: int = 0, line_no: int | None = None) -> dict:
    """Parse one VCF data line into a call or a classified skip.

    Returns a dict carrying the raw fields plus:

      ``skip``        None, or one of ``SKIP_REASONS``
      ``is_ref_block`` True for a gVCF <NON_REF>-only block
      ``call``        the ``parsers.py`` shaped dict when ``skip`` is None

    Skip classification runs in a fixed order so every record lands in exactly
    one counter: reference block, then symbolic, then the shape of the SELECTED
    alleles, then no-call, then missing rsid. Classifying on the selected alleles
    rather than the whole ALT list matters: ALT="G,AT" with GT="0/1" is a clean
    A/G heterozygote for this person and there is no reason to discard it.
    """
    where = f" (line {line_no})" if line_no else ""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 8:
        raise ParseError(
            f"Malformed VCF data line{where}: expected at least 8 tab-separated "
            f"columns, found {len(fields)}."
        )

    chrom = normalize_contig(fields[0])
    try:
        position = int(fields[1].strip())
    except ValueError as exc:
        raise ParseError(f"Unreadable POS value {fields[1]!r}{where}.") from exc

    rsid = _extract_rsid(fields[2])
    ref = fields[3].strip().upper()
    alt_field = fields[4].strip().upper()
    filter_field = fields[6].strip()

    alts_full: list[str] = [] if alt_field in (".", "") else [a.strip() for a in alt_field.split(",")]

    record: dict[str, Any] = {
        "chromosome": chrom,
        "position": position,
        "rsid": rsid,
        "ref": ref,
        "alts": alts_full,
        "filter": filter_field,
        "gt": None,
        "phased": False,
        "allele1": "N",
        "allele2": "N",
        "is_ref_block": False,
        "skip": None,
        "call": None,
    }

    # gVCF reference block. Emitting a call here would produce one genotype for a
    # run that can be thousands of bases long, all of it recorded at the block's
    # start position. That is the specific way gVCF ingest goes silently wrong.
    if alts_full and all(a == "<NON_REF>" for a in alts_full):
        record["is_ref_block"] = True
        record["skip"] = "symbolic"
        return record

    if len(fields) < 9:
        raise ParseError(
            f"VCF record{where} has no FORMAT or sample column, so no genotype can "
            "be read. Sites-only VCFs carry no calls for any individual."
        )
    sample_col = 9 + sample_index
    if sample_col >= len(fields):
        raise ParseError(
            f"VCF record{where} has {len(fields) - 9} sample column(s); sample index "
            f"{sample_index} does not exist."
        )

    format_keys = fields[8].strip().split(":")
    sample_values = fields[sample_col].strip().split(":")
    gt_value = ""
    if "GT" in format_keys:
        gt_pos = format_keys.index("GT")
        if gt_pos < len(sample_values):
            gt_value = sample_values[gt_pos].strip()
    record["gt"] = gt_value or None

    decoded = alleles_from_gt(gt_value, ref, alts_full)
    record["allele1"] = decoded["allele1"]
    record["allele2"] = decoded["allele2"]
    record["phased"] = decoded["phased"]

    selected = decoded["selected"]
    # <NON_REF> is a placeholder for "some allele we did not observe", so a GT
    # that lands on it carries no base. It is dropped before shape checks.
    real_alts = [a for a in alts_full if a != "<NON_REF>"]

    if any(_is_symbolic(a) for a in selected):
        record["skip"] = "symbolic"
        return record

    if decoded["ploidy"] > 2:
        record["skip"] = "multiallelic_complex"
        return record

    if not decoded["no_call"]:
        clean = all(is_simple_snv(ref, a) for a in selected)
        if not clean:
            record["skip"] = "multiallelic_complex" if len(real_alts) > 1 else "indel"
            return record

    if decoded["no_call"] or decoded["partial_no_call"]:
        # A half genotype cannot be turned into a comparable pair, so it is a
        # no-call like any other rather than a call with one allele guessed.
        record["skip"] = "no_call"
        return record

    if rsid is None:
        # Every annotation path in DNAInsight is keyed by rsid. A position with
        # ID "." is a real observation with nothing to join it to, which is why
        # it is counted rather than dropped.
        record["skip"] = "no_id"
        return record

    record["call"] = {
        "rsid": rsid,
        "chromosome": chrom,
        "position": position,
        "allele1": record["allele1"],
        "allele2": record["allele2"],
    }
    return record


# ---------------------------------------------------------------------------
# Sample selection
# ---------------------------------------------------------------------------

def select_sample(header: VcfHeader, sample: str | int | None = None) -> tuple[int, str]:
    """Resolve a sample name or index to ``(index, name)``.

    Defaulting to the first sample is a convenience, not a judgement, so the
    caller is handed ``sample_count`` alongside and can refuse the ambiguity
    itself. A trio VCF silently analysed as the proband when it was ordered
    mother-father-child is the failure this exists to make visible.
    """
    if not header.samples:
        raise ParseError(
            "This VCF carries no sample columns, so it contains no genotypes. "
            "A sites-only VCF lists variants observed in a cohort, not a person."
        )
    if sample is None:
        return 0, header.samples[0]
    if isinstance(sample, bool):
        raise ParseError(f"Invalid sample selector: {sample!r}")
    if isinstance(sample, int):
        if sample < 0 or sample >= len(header.samples):
            raise ParseError(
                f"Sample index {sample} is out of range; this file has "
                f"{len(header.samples)} sample(s): {', '.join(header.samples)}."
            )
        return sample, header.samples[sample]
    name = str(sample).strip()
    if name in header.samples:
        return header.samples.index(name), name
    raise ParseError(
        f"Sample {name!r} is not in this file. Available samples: "
        f"{', '.join(header.samples) or 'none'}."
    )


# ---------------------------------------------------------------------------
# Streaming readers
# ---------------------------------------------------------------------------

def iter_vcf_records(filepath: str | Path, *, sample: str | int | None = None) -> Iterator[dict]:
    """Yield one ``parse_record`` dict per data line, streaming.

    Nothing is accumulated. A caller that only needs counts can run this over a
    30x WGS VCF in constant memory.
    """
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")
    with _open_text(path) as handle:
        header = _read_header(handle)
        sample_index, _name = select_sample(header, sample)
        line_no = header.line_count
        for raw in handle:
            line_no += 1
            if not raw.strip() or raw.startswith("#"):
                continue
            yield parse_record(raw, sample_index=sample_index, line_no=line_no)


def read_vcf(filepath: str | Path, *, sample: str | int | None = None) -> dict:
    """Stream a VCF or gVCF and collect the calls DNAInsight can use.

    Returns ``{"header", "sample", "sample_index", "sample_count", "snps",
    "skipped", "ref_blocks", "non_pass", "record_count", "warnings"}``.

    Only accepted calls are retained, which is the same working set
    ``parsers.parse_dna_file`` produces. The file itself is never held in memory.
    """
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")

    skipped = {reason: 0 for reason in SKIP_REASONS}
    snps: list[dict] = []
    ref_blocks = 0
    non_pass = 0
    record_count = 0
    saw_non_ref = False

    with _open_text(path) as handle:
        header = _read_header(handle)
        sample_index, sample_name = select_sample(header, sample)
        line_no = header.line_count
        for raw in handle:
            line_no += 1
            if not raw.strip() or raw.startswith("#"):
                continue
            record_count += 1
            parsed = parse_record(raw, sample_index=sample_index, line_no=line_no)
            if parsed["is_ref_block"]:
                ref_blocks += 1
                saw_non_ref = True
            elif "<NON_REF>" in parsed["alts"]:
                saw_non_ref = True
            if parsed["filter"] and parsed["filter"] not in (".", "PASS"):
                non_pass += 1
            reason = parsed["skip"]
            if reason is not None:
                skipped[reason] += 1
                continue
            snps.append(parsed["call"])

    warnings: list[str] = []
    if header.sample_count > 1 and sample is None:
        warnings.append(
            f"{header.sample_count} samples are present in this file and none was "
            f"requested, so {sample_name!r} (the first) was used. Pass sample= to "
            "choose, or refuse the file if the order is not known."
        )
    if ref_blocks:
        warnings.append(
            f"{ref_blocks} gVCF reference block(s) were skipped rather than expanded. "
            "A block asserts only that no variant was called across its span."
        )
    if non_pass:
        warnings.append(
            f"{non_pass} record(s) carry a FILTER other than PASS and were kept. "
            "DNAInsight does not drop them, because filter thresholds are the "
            "calling pipeline's opinion, but they are lower confidence."
        )
    if not header.fileformat:
        warnings.append("No ##fileformat line was present, so the VCF version is unknown.")

    return {
        "header": header,
        "sample": sample_name,
        "sample_index": sample_index,
        "sample_count": header.sample_count,
        "snps": snps,
        "skipped": skipped,
        "ref_blocks": ref_blocks,
        "non_pass": non_pass,
        "record_count": record_count,
        "gvcf": header.gvcf_hint or saw_non_ref,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Liftover
#
# There is no liftover in the standard library, and the chain files that make
# one possible are UCSC downloads, not something this repository can bundle
# under its own data rules. So the code below is the whole implementation: a
# chain parser small enough to audit, and a refusal when the chain file is not
# there. The refusal is the important half. Relabelling GRCh38 coordinates as
# GRCh37 costs nothing at runtime and produces a report that is wrong in a way
# no user can detect.
# ---------------------------------------------------------------------------

_UCSC_NAMES = {"GRCh37": "hg19", "GRCh38": "hg38"}


def chain_dir() -> Path:
    """Directory where a user-supplied UCSC chain file is looked for."""
    return external.panel_root() / "chains"


def find_chain_file(from_build: str, to_build: str) -> Path | None:
    """Locate a chain file for this build pair, or None.

    Both the UCSC download name and a plain descriptive name are accepted,
    gzipped or not, because a user who has fetched ``hg19ToHg38.over.chain.gz``
    should not have to rename it to satisfy us.
    """
    base = chain_dir()
    if not base.is_dir():
        return None
    src = _UCSC_NAMES.get(from_build, str(from_build).lower())
    dst = _UCSC_NAMES.get(to_build, str(to_build).lower())
    ucsc = f"{src}To{dst[0].upper()}{dst[1:]}"
    candidates = [
        f"{from_build}_to_{to_build}.chain",
        f"{from_build}_to_{to_build}.chain.gz",
        f"{ucsc}.over.chain",
        f"{ucsc}.over.chain.gz",
        f"{ucsc}.chain",
        f"{ucsc}.chain.gz",
    ]
    for name in candidates:
        path = base / name
        if path.is_file():
            return path
    return None


class ChainFile:
    """A parsed UCSC chain file, queried one position at a time.

    Chain format, for the record, because it is easy to get backwards: the "t"
    (target) coordinates belong to the FROM assembly and the "q" (query)
    coordinates to the TO assembly. Both are 0-based half-open. A chain header
    is followed by block lines of ``size dt dq`` and terminated by a bare
    ``size``. ``dt`` and ``dq`` are the gaps that follow the block in each
    assembly, and they are what make liftover a real operation rather than an
    offset.
    """

    def __init__(self, path: Path | None = None) -> None:
        self.path = Path(path) if path else None
        self.chains: dict[str, list[dict]] = {}

    @classmethod
    def load(cls, path: str | Path) -> "ChainFile":
        target = Path(path)
        if not target.is_file():
            raise ParseError(f"Chain file not found: {path}")
        obj = cls(target)
        opener = gzip.open if _is_gzip(target) else open
        with opener(target, "rt", encoding="utf-8", errors="replace") as handle:
            obj._parse(handle)
        return obj

    @classmethod
    def from_text(cls, text: str) -> "ChainFile":
        """Parse chain text directly. Used by tests and by in-memory callers."""
        obj = cls(None)
        obj._parse(text.splitlines())
        return obj

    def _parse(self, lines: Iterable[str]) -> None:
        current: dict | None = None
        t_pos = q_pos = 0
        for raw in lines:
            line = raw.strip()
            if not line or line.startswith("#"):
                current = None
                continue
            if line.startswith("chain"):
                parts = line.split()
                if len(parts) < 12:
                    raise ParseError(
                        f"Malformed chain header, expected at least 12 fields: {line!r}"
                    )
                try:
                    current = {
                        "score": int(parts[1]),
                        "t_name": parts[2], "t_size": int(parts[3]), "t_strand": parts[4],
                        "t_start": int(parts[5]), "t_end": int(parts[6]),
                        "q_name": parts[7], "q_size": int(parts[8]), "q_strand": parts[9],
                        "q_start": int(parts[10]), "q_end": int(parts[11]),
                        "id": parts[12] if len(parts) > 12 else "",
                        "blocks": [],
                    }
                except ValueError as exc:
                    raise ParseError(f"Malformed chain header {line!r}: {exc}") from exc
                t_pos, q_pos = current["t_start"], current["q_start"]
                self.chains.setdefault(normalize_contig(current["t_name"]), []).append(current)
                continue
            if current is None:
                continue
            parts = line.split()
            try:
                size = int(parts[0])
            except ValueError as exc:
                raise ParseError(f"Malformed chain block line {line!r}: {exc}") from exc
            current["blocks"].append((t_pos, q_pos, size))
            if len(parts) >= 3:
                t_pos += size + int(parts[1])
                q_pos += size + int(parts[2])
            else:
                # A bare size terminates the chain. Anything after it belongs to
                # the next header, so dropping the reference here prevents a
                # missing blank line from welding two chains together.
                current = None
        for chain_list in self.chains.values():
            chain_list.sort(key=lambda c: c["score"], reverse=True)

    @property
    def chain_count(self) -> int:
        return sum(len(v) for v in self.chains.values())

    def map_position(self, chrom: str, position: int) -> dict | None:
        """Map a 1-based position, or return None when it is unmappable.

        None means the position falls in a chain gap or on a contig this file
        does not cover. Those are real outcomes: a few percent of GRCh37
        positions have no GRCh38 equivalent at all, and inventing one for them
        is the entire reason this function does not fall back to identity.
        """
        key = normalize_contig(chrom)
        p0 = int(position) - 1
        for chain in self.chains.get(key, []):
            if not (chain["t_start"] <= p0 < chain["t_end"]):
                continue
            for t_block, q_block, size in chain["blocks"]:
                if t_block <= p0 < t_block + size:
                    q0 = q_block + (p0 - t_block)
                    if chain["q_strand"] == "-":
                        q0 = chain["q_size"] - q0 - 1
                    return {
                        "chromosome": normalize_contig(chain["q_name"]),
                        "position": q0 + 1,
                        "strand": chain["q_strand"],
                        "chain_id": chain["id"],
                        "score": chain["score"],
                    }
            return None
        return None


def parse_chain_file(path: str | Path) -> ChainFile:
    """Load and index a UCSC chain file."""
    return ChainFile.load(path)


def liftover_available(from_build: str = "GRCh38", to_build: str = REFERENCE_BUILD) -> bool:
    """True only when a chain file for this pair is actually on disk."""
    return find_chain_file(from_build, to_build) is not None


def _chain_unavailable(from_build: str, to_build: str, *, detail: str = "") -> dict:
    """The degraded payload, shaped exactly like ``external.unavailable()``.

    Reusing that shape is not cosmetic. The frontend already renders one honest
    empty state for it, and a caller can already tell "looked and found nothing"
    from "could not look" without learning a second convention.
    """
    return {
        "available": False,
        "capability": LIFTOVER_CAPABILITY,
        "tool": "UCSC chain file",
        "tool_id": "chain_file",
        "state": "not_installed",
        "reason": detail or (
            f"No {from_build} to {to_build} chain file is installed, so no coordinate "
            "was converted. The input coordinates are returned unchanged and are still "
            f"{from_build}; they have NOT been relabelled as {to_build}."
        ),
        "not_attempted": True,
        "results": [],
        "how_to_enable": {
            "what": f"A UCSC chain file mapping {from_build} to {to_build}.",
            "homepage": "https://hgdownload.soe.ucsc.edu/downloads.html",
            "install_to": str(chain_dir()),
            "expected_files": [
                f"{from_build}_to_{to_build}.chain",
                f"{_UCSC_NAMES.get(from_build, from_build)}To"
                f"{_UCSC_NAMES.get(to_build, to_build).capitalize()}.over.chain.gz",
            ],
            "steps": [
                "Download the over.chain.gz file for this build pair from UCSC.",
                f"Place it in {chain_dir()} (create the folder if needed).",
                "Chain files are UCSC downloads and are not bundled with DNAInsight.",
            ],
        },
        "from_build": from_build,
        "to_build": to_build,
    }


def liftover(records: Iterable[dict], from_build: str, to_build: str) -> dict:
    """Convert record coordinates between builds, or refuse and say why.

    Records are the usual ``{"rsid", "chromosome", "position", "allele1",
    "allele2"}`` dicts. A record whose position lands in a chain gap goes to
    ``unmapped`` and is not carried over with its old coordinate.

    When a chain maps to the minus strand of the destination assembly the
    alleles are complemented through ``orientation.complement_allele`` and the
    record is marked ``strand_flipped``. Moving the coordinate while leaving a
    plus-strand allele alone would produce a genotype that is wrong at a
    position that is right, which is harder to notice than either error alone.
    """
    if from_build == to_build:
        return {
            "available": True,
            "capability": LIFTOVER_CAPABILITY,
            "from_build": from_build,
            "to_build": to_build,
            "chain_file": None,
            "results": [dict(r) for r in records],
            "unmapped": [],
            "mapped_count": 0,
            "unmapped_count": 0,
            "not_attempted": True,
            "reason": "Source and destination builds are the same; nothing to convert.",
        }

    chain_path = find_chain_file(from_build, to_build)
    if chain_path is None:
        return _chain_unavailable(from_build, to_build)

    chain = ChainFile.load(chain_path)
    mapped: list[dict] = []
    unmapped: list[dict] = []
    for record in records:
        try:
            position = int(record.get("position"))
        except (TypeError, ValueError):
            unmapped.append({**dict(record), "reason": "position is not an integer"})
            continue
        hit = chain.map_position(record.get("chromosome", ""), position)
        if hit is None:
            unmapped.append({
                **dict(record),
                "reason": f"no {from_build} to {to_build} chain block covers this position",
            })
            continue
        out = dict(record)
        out["chromosome"] = hit["chromosome"]
        out["position"] = hit["position"]
        out["strand_flipped"] = hit["strand"] == "-"
        if out["strand_flipped"]:
            out["allele1"] = orientation.complement_allele(record.get("allele1"))
            out["allele2"] = orientation.complement_allele(record.get("allele2"))
        mapped.append(out)

    return {
        "available": True,
        "capability": LIFTOVER_CAPABILITY,
        "from_build": from_build,
        "to_build": to_build,
        "chain_file": str(chain_path),
        "results": mapped,
        "unmapped": unmapped,
        "mapped_count": len(mapped),
        "unmapped_count": len(unmapped),
        "not_attempted": False,
    }


# ---------------------------------------------------------------------------
# BAM and CRAM: targeted extraction
#
# DNAInsight annotates around 122 curated positions plus whatever a genoset or
# PRS model needs. Calling variants across three gigabases to reach them would
# be hours of compute to answer a question that a pileup over a few thousand
# BED intervals answers in seconds. The design is therefore targeted extraction,
# never whole-genome calling, and the heavy lifting belongs to samtools.
# ---------------------------------------------------------------------------

def _position_tuple(entry: Any) -> tuple[str, int, str]:
    """Accept a dict, a 2-tuple or a 3-tuple and return (chrom, pos, name)."""
    if isinstance(entry, dict):
        chrom = normalize_contig(entry.get("chromosome") or entry.get("chrom") or "")
        raw_pos = entry.get("position")
        name = str(entry.get("rsid") or entry.get("name") or ".")
    elif isinstance(entry, (list, tuple)) and len(entry) >= 2:
        chrom = normalize_contig(entry[0])
        raw_pos = entry[1]
        name = str(entry[2]) if len(entry) > 2 else "."
    else:
        raise ParseError(f"Unreadable position entry: {entry!r}")
    try:
        pos = int(raw_pos)
    except (TypeError, ValueError) as exc:
        raise ParseError(f"Unreadable position value in entry {entry!r}") from exc
    if not chrom:
        raise ParseError(f"Position entry {entry!r} has no chromosome.")
    if pos < 1:
        raise ParseError(
            f"Position {pos} in entry {entry!r} is not 1-based. VCF and consumer "
            "exports are 1-based; only the BED written here is 0-based."
        )
    return chrom, pos, name


def write_positions_bed(positions: Iterable[Any],
                        dest: str | Path,
                        *,
                        chrom_prefix: str = "") -> Path:
    """Write DNAInsight's positions of interest as a BED file.

    BED is 0-based half-open while everything else in this module is 1-based, so
    a position P becomes the interval [P-1, P). Getting that wrong shifts every
    pileup by one base and yields a plausible genotype for the neighbouring
    position, which is why the conversion lives in exactly one function.

    ``chrom_prefix`` exists because a BED whose contig names do not match the BAM
    header produces an empty pileup and no error at all. The caller that knows
    the BAM's naming style passes "chr" when it needs to.
    """
    seen: set[tuple[str, int]] = set()
    rows: list[tuple[str, int, str]] = []
    for entry in positions:
        chrom, pos, name = _position_tuple(entry)
        if (chrom, pos) in seen:
            continue
        seen.add((chrom, pos))
        rows.append((chrom, pos, name))
    rows.sort(key=lambda r: (_chrom_sort_key(r[0]), r[1]))

    target = Path(dest)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w", encoding="utf-8", newline="\n") as fh:
        for chrom, pos, name in rows:
            fh.write(f"{chrom_prefix}{chrom}\t{pos - 1}\t{pos}\t{name}\n")
    return target


def _count_pileup_bases(bases: str, ref: str) -> dict[str, int]:
    """Count observed bases in an mpileup read-base string.

    The read-base column is a small language, not a list of bases: "^" carries a
    mapping-quality character that must be skipped or it is counted as an allele,
    "+2AG" and "-2AG" embed indel sequences whose bases are not observations at
    this position, and "*", ">" and "<" are placeholders for deleted or skipped
    reference. Naive counting of A/C/G/T over this string overcounts, so the
    walk below is explicit.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    ref_base = str(ref or "").strip().upper()
    text = str(bases or "")
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "^":
            i += 2               # read start, plus its mapping quality char
            continue
        if ch == "$":
            i += 1
            continue
        if ch in "+-":
            j = i + 1
            digits = ""
            while j < n and text[j].isdigit():
                digits += text[j]
                j += 1
            i = j + (int(digits) if digits else 0)
            continue
        if ch in ".,":
            if ref_base in counts:
                counts[ref_base] += 1
            i += 1
            continue
        upper = ch.upper()
        if upper in counts:
            counts[upper] += 1
            i += 1
            continue
        i += 1                   # '*', '>', '<', 'N' and anything else: no allele
    return counts


def parse_mpileup_line(line: str) -> dict | None:
    """Parse one ``samtools mpileup`` output row into base counts.

    Returns None for a blank line. Kept separate from ``extract_positions`` so
    the parsing can be tested without samtools installed, which is the only part
    of the alignment path that can be tested at all in a clean environment.
    """
    text = line.rstrip("\r\n")
    if not text.strip():
        return None
    fields = text.split("\t")
    if len(fields) < 4:
        raise ParseError(
            f"Unreadable mpileup row, expected at least 4 columns: {text[:80]!r}"
        )
    try:
        position = int(fields[1])
        depth = int(fields[3])
    except ValueError as exc:
        raise ParseError(f"Unreadable mpileup row {text[:80]!r}: {exc}") from exc
    ref = fields[2].strip().upper()
    bases = fields[4] if len(fields) > 4 else ""
    counts = _count_pileup_bases(bases, ref)
    return {
        "chromosome": normalize_contig(fields[0]),
        "position": position,
        "ref": ref,
        "depth": depth,
        "counts": counts,
        "observed": sum(counts.values()),
    }


# Thresholds for the two-allele call below. Named constants rather than magic
# numbers because they are the difference between a reported genotype and a
# no-call, and a reviewer needs to see them without reading the arithmetic.
MIN_PILEUP_DEPTH = 8
MIN_ALLELE_FRACTION = 0.20


def call_genotype_from_counts(counts: dict[str, int],
                              *,
                              min_depth: int = MIN_PILEUP_DEPTH,
                              min_fraction: float = MIN_ALLELE_FRACTION) -> tuple[str, str]:
    """Call a diploid genotype from base counts, or return a no-call.

    This is deliberately not a variant caller. There is no error model, no
    mapping-quality weighting and no prior. It answers one narrow question:
    which one or two bases are supported strongly enough at a position we
    already care about. Below ``min_depth`` the answer is ("N", "N"), because a
    genotype from four reads is a guess wearing a result's clothes, and
    ``orientation`` and the scoring engine already treat "N" as "not read".
    """
    observed = {b: int(c) for b, c in (counts or {}).items() if b in _ACGT}
    depth = sum(observed.values())
    if depth < min_depth:
        return "N", "N"
    supported = [b for b, c in observed.items() if c and (c / depth) >= min_fraction]
    if not supported:
        return "N", "N"
    supported.sort(key=lambda b: (-observed[b], b))
    if len(supported) == 1:
        return supported[0], supported[0]
    first, second = sorted(supported[:2])
    return first, second


def _samtools_guard() -> dict | None:
    """Standard guard, with the registry gap explained rather than leaked."""
    blocked = external.guard(SAMTOOLS_TOOL_ID, SAMTOOLS_CAPABILITY)
    if blocked is None:
        return None
    if blocked.get("state") == "unknown":
        # external.py has no samtools entry yet, so status() reports "unknown"
        # and the stock reason reads like an internal error. The shape is left
        # untouched; only the sentence a user sees is replaced. This branch
        # stops firing on its own once the wiring pass adds the entry.
        blocked = dict(blocked)
        blocked["tool"] = "SAMtools"
        blocked["reason"] = (
            "SAMtools is not registered as an external tool in this build, so BAM "
            "and CRAM extraction was not attempted. This is different from finding "
            "nothing in the alignment."
        )
    return blocked


def extract_positions(bam_path: str | Path,
                      positions: Iterable[Any],
                      build: str = REFERENCE_BUILD,
                      *,
                      reference: str | Path | None = None,
                      chrom_prefix: str = "",
                      min_depth: int = MIN_PILEUP_DEPTH,
                      min_fraction: float = MIN_ALLELE_FRACTION,
                      min_mapping_quality: int = 1,
                      min_base_quality: int = 13) -> dict:
    """Pile up an alignment at DNAInsight's positions of interest.

    Returns either the usual degraded payload (samtools absent, or CRAM without
    its reference) or ``{"available": True, "provider", "build", "snps", ...}``
    with ``snps`` in the same shape ``parsers.py`` returns.

    The build check runs BEFORE the tool check on purpose. A GRCh38 BAM piled up
    at GRCh37 positions returns depth and bases at every interval, so it looks
    like a complete success and is wrong everywhere. Refusing costs the user a
    minute; not refusing costs them the report.
    """
    path = Path(bam_path)
    if not path.exists():
        raise ParseError(f"Alignment file not found: {bam_path}")
    assert_build_compatible(build, REFERENCE_BUILD)

    fmt = detect_format(path)
    if fmt not in ("bam", "cram"):
        raise ParseError(
            f"{path.name} is not a BAM or CRAM alignment (detected {fmt!r})."
        )

    requested = [_position_tuple(p) for p in positions]
    if not requested:
        raise ParseError("No positions were requested, so there is nothing to extract.")

    if fmt == "cram" and reference is None:
        # CRAM stores bases as differences from a reference. Without the exact
        # FASTA the file was compressed against, the sequence cannot be
        # reconstructed at all, so this is a refusal and not a degradation of
        # quality.
        return {
            "available": False,
            "capability": SAMTOOLS_CAPABILITY,
            "tool": "SAMtools",
            "tool_id": SAMTOOLS_TOOL_ID,
            "state": "reference_missing",
            "reason": (
                "CRAM is reference-compressed and cannot be decoded without the exact "
                "reference FASTA it was written against. Pass reference= pointing at "
                f"the {REFERENCE_BUILD} FASTA used by the sequencing provider."
            ),
            "not_attempted": True,
            "results": [],
            "how_to_enable": None,
        }

    blocked = _samtools_guard()
    if blocked is not None:
        return blocked

    workdir = Path(tempfile.mkdtemp(prefix="dnainsight-pileup-"))
    try:
        bed = write_positions_bed(requested, workdir / "positions.bed", chrom_prefix=chrom_prefix)
        args: list[str] = [
            "mpileup",
            "-l", str(bed),
            "-q", str(min_mapping_quality),
            "-Q", str(min_base_quality),
            "--no-BAQ",
            "-a",
        ]
        if reference is not None:
            args += ["-f", str(reference)]
        args.append(str(path))
        completed = external.run(SAMTOOLS_TOOL_ID, args)
        rows = [parse_mpileup_line(line) for line in (completed.stdout or "").splitlines()]
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    names = {(chrom, pos): name for chrom, pos, name in requested}
    snps: list[dict] = []
    low_depth = 0
    for row in rows:
        if row is None:
            continue
        key = (row["chromosome"], row["position"])
        rsid = names.get(key)
        if rsid is None or rsid == ".":
            continue
        allele1, allele2 = call_genotype_from_counts(
            row["counts"], min_depth=min_depth, min_fraction=min_fraction
        )
        if allele1 == "N":
            low_depth += 1
            continue
        snps.append({
            "rsid": rsid,
            "chromosome": row["chromosome"],
            "position": row["position"],
            "allele1": allele1,
            "allele2": allele2,
        })

    warnings: list[str] = []
    missing = len(requested) - len(snps) - low_depth
    if low_depth:
        warnings.append(
            f"{low_depth} position(s) had fewer than {min_depth} usable reads and were "
            "reported as no-calls rather than called from thin coverage."
        )
    if missing > 0:
        warnings.append(
            f"{missing} requested position(s) produced no pileup row at all. A contig "
            "naming mismatch between the BED and the alignment header is the usual "
            "cause, and it produces silence rather than an error."
        )
    return {
        "available": True,
        "capability": SAMTOOLS_CAPABILITY,
        "provider": fmt,
        "build": build,
        "positions_requested": len(requested),
        "snp_count": len(snps),
        "snps": snps,
        "no_call_count": low_depth,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def _empty_result(provider: str, fmt: str) -> dict:
    return {
        "provider": provider,
        "format": fmt,
        "build": None,
        "build_confidence": _CONFIDENCE_NONE,
        "build_evidence": [],
        "sample": None,
        "sample_count": 0,
        "snp_count": 0,
        "snps": [],
        "skipped": {reason: 0 for reason in SKIP_REASONS},
        "warnings": [],
    }


def parse_sequencing_file(filepath: str | Path,
                          *,
                          sample: str | int | None = None,
                          expected_build: str | None = REFERENCE_BUILD) -> dict:
    """Parse a VCF, gVCF, BAM or CRAM into DNAInsight's standard SNP shape.

    Returns::

        {"provider", "format", "build", "build_confidence", "build_evidence",
         "sample", "sample_count", "snp_count", "snps", "skipped", "warnings"}

    ``snps`` carries the same dicts ``parsers.parse_dna_file`` produces, so every
    downstream stage is unchanged. ``skipped`` always carries all five reasons,
    zero-valued when nothing was skipped for that reason.

    Raises ParseError, or BuildMismatch (a ParseError subclass) when the file is
    a build DNAInsight cannot annotate. Callers keep one except clause.

    An undetectable build is a warning, not an error: a file with no contig lines
    is still readable, and the caller is handed ``build: None`` with
    ``build_confidence: "none"`` so it can decide. A build that is detected and
    wrong is always fatal, because there is nothing to decide.
    """
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")

    fmt = detect_format(path)

    if fmt in ("bam", "cram"):
        result = _empty_result(fmt, fmt)
        payload = _samtools_guard()
        result["warnings"].append(
            "BAM and CRAM ingest is targeted extraction at DNAInsight's reference "
            "positions, not whole-genome variant calling. Call extract_positions() "
            "with the positions of interest."
        )
        if payload is not None:
            result["warnings"].append(str(payload.get("reason", "")))
        result["warnings"].append(
            "The genome build of an alignment is recorded in its @SQ headers, which "
            "cannot be read without samtools, so no build was determined here."
        )
        return result

    if fmt not in ("vcf", "vcf.gz"):
        raise ParseError(
            f"{path.name} is not a format this module reads (detected {fmt!r}). "
            "Expected VCF, gVCF, BAM or CRAM."
        )

    header = read_header(path)
    detected = detect_build(header)
    build = detected["build"]

    # Checked before the body is streamed. There is no point spending minutes
    # reading a gigabyte of coordinates that are going to be refused.
    assert_build_compatible(build, expected_build, allow_unknown=True)

    parsed = read_vcf(path, sample=sample)

    warnings = list(parsed["warnings"])
    if build is None:
        if detected["confidence"] == _CONFIDENCE_CONFLICT:
            warnings.insert(0, (
                "Genome build evidence in this header is self-contradictory, so no "
                "build was assigned. A file carrying contig lengths from two builds "
                "has been concatenated or hand-edited."
            ))
        else:
            warnings.insert(0, (
                f"Genome build could not be determined, so it was NOT assumed to be "
                f"{expected_build or REFERENCE_BUILD}. Add ##contig lines with lengths, "
                "or confirm the build before trusting any position in this report."
            ))
    else:
        claimed = {e["supports"] for e in detected["evidence"]
                   if e["signal"] in ("assembly_tag", "reference_header") and e["supports"]}
        if claimed and claimed != {build}:
            warnings.append(
                f"Contig lengths say {build} but the header text claims "
                f"{', '.join(sorted(claimed))}. The lengths were believed, because they "
                "come from the reference FASTA and the text does not."
            )

    if parsed["record_count"] == 0:
        raise ParseError(
            f"No variant records were found in {path.name}. The header parsed but the "
            "body is empty, which usually means a truncated download."
        )

    provider = "gvcf" if parsed["gvcf"] else "vcf"
    if not parsed["snps"]:
        # Not an error. Every record was classified and counted, and those counts
        # are more useful to the user than an exception that hides them.
        warnings.append(
            "No annotatable SNVs were produced from this file. The skipped counts "
            "record why every record was set aside."
        )

    return {
        "provider": provider,
        "format": fmt,
        "build": build,
        "build_confidence": detected["confidence"],
        "build_evidence": detected["evidence"],
        "sample": parsed["sample"],
        "sample_count": parsed["sample_count"],
        "snp_count": len(parsed["snps"]),
        "snps": parsed["snps"],
        "skipped": parsed["skipped"],
        "warnings": warnings,
    }
