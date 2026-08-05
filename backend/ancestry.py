"""
ancestry.py -- global and local ancestry, with the panel, the model and the
per-population marker coverage all shipped rather than asserted.

WHY THIS MODULE LOOKS THE WAY IT DOES
-------------------------------------
Every incumbent ships admixture as a black box, and the marketing numbers do
not survive contact with the documentation:

  * 23andMe advertises "4,500+ regions" while its own methodology paper
    documents 78 hierarchical populations plus a separate country-match feature
    built on more than 4,000 recent-ancestor matches. Those are two different
    things presented as one.
  * MyHeritage went from 42 to 79 ethnicities in v2.5 without publishing the
    model, the panel or the validation.
  * GEDmatch runs community calculators whose model files carry no licence at
    all, which is why external.BLOCKED refuses them by name.

None of those users can answer the only question that matters: "could my array
even see this population?" DNAInsight answers it. This module ships the panel
manifest, the model parameters, per-population marker coverage, and an explicit
NOT RESOLVABLE state for populations the user's array cannot distinguish.

That last state is invariant 3 (never confuse "checked and absent" with "never
checked") applied to ancestry. Reporting a population at 0.0 percent when the
array reads 11 of its 900 informative markers is a false statement. The honest
output is that the question was not answerable, and it is a different output.

CONFIDENCE INTERVALS
--------------------
A point estimate of "12.4 percent Iberian" with no interval is a number
pretending to be a fact. Every proportion here carries an interval and, more
importantly, a LABEL saying how the interval was produced, because a bootstrap
over marker subsets and a normal approximation from a marker count are not the
same claim and must not look the same.

LICENCE BOUNDARY
----------------
fastmixture is GPL-3.0, FLARE is Apache-2.0. Neither is imported, linked or
vendored. Every invocation goes through ``external.run``, and THE SUBPROCESS
BOUNDARY IS THE LICENCE BOUNDARY. ADMIXTURE (academic-only) and RFMix
(academic-only) are permanently blocked in external.BLOCKED and are not
reachable from here.

OFFLINE CONTRACT
----------------
No network access, at import or on any read path.

SCOPE
-----
Y and mtDNA haplogroups are NOT here. backend/haplogroups.py owns them. This
module only ever references it through a defensive local import that degrades
when it is absent, so the two can land in either order.
"""

from __future__ import annotations

import gzip
import hashlib
import math
import random
import re
import struct
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import external
from .frequency import POPULATIONS as ARRAY_POPULATIONS

__all__ = [
    "AncestryError",
    "DEFAULT_PANEL", "DEFAULT_MODE", "MODES",
    "MIN_RESOLVABLE_FRACTION", "MIN_RESOLVABLE_MARKERS",
    "CI_LEVEL", "CI_METHODS", "SUPERPOP_COLOURS", "FALLBACK_PALETTE",
    "MANDATORY_CAVEATS", "CAVEAT_PANEL_BIAS", "CAVEAT_MODEL_DEPENDENT",
    "CAVEAT_ARRAY_DENSITY", "CAVEAT_NOT_IDENTITY", "HASH_MAX_BYTES",
    "normal_cdf", "inverse_normal_cdf", "z_for_level",
    "wilson_interval", "percentile_interval",
    "parse_q_file", "parse_fam", "parse_population_map",
    "marker_coverage", "resolvable",
    "global_ancestry", "local_ancestry", "chromosome_painting",
    "panel_manifest", "ancestry_caveats", "haplogroup_note",
    "write_plink", "segments_from_calls", "panel_unavailable",
]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class AncestryError(Exception):
    """Base class for failures raised by this module."""


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_PANEL: str = "onekg_sgdp"

# projection scores the user against a FIXED reference panel. That is the same
# shape as a community admixture calculator, and it is the mode that makes the
# result reproducible and the model publishable. unsupervised re-estimates the
# components from scratch and its output columns are not populations at all,
# which is why it is labelled differently downstream.
MODES: tuple[str, ...] = ("projection", "supervised", "unsupervised")
DEFAULT_MODE: str = "projection"

# RESOLUTION THRESHOLDS, documented here because a threshold in a template is a
# threshold nobody can audit.
#
# A population is reported only when the user's array reads at least
# MIN_RESOLVABLE_FRACTION of the panel's informative markers for it AND at
# least MIN_RESOLVABLE_MARKERS of them in absolute terms. Both are needed: 25
# percent of 40 markers is 10 markers, which cannot separate anything, and 200
# markers out of 90,000 is a rounding error however impressive 200 sounds.
MIN_RESOLVABLE_FRACTION: float = 0.20
MIN_RESOLVABLE_MARKERS: int = 50

CI_LEVEL: float = 0.95

# How an interval was produced. Shown to the user, not just stored.
CI_METHODS: tuple[str, ...] = (
    "tool_reported",                  # the tool supplied its own standard errors
    "bootstrap_marker_subsets",       # we re-ran the tool on resampled markers
    "wilson_marker_count",            # normal approximation from marker counts
    "none",                           # no interval could be produced
)

CI_METHOD_LABELS: dict[str, str] = {
    "tool_reported": (
        "Interval reported by the ancestry tool itself."
    ),
    "bootstrap_marker_subsets": (
        "Interval from a bootstrap over marker subsets: the estimate was "
        "recomputed on resampled sets of your markers and the spread of those "
        "estimates is the interval."
    ),
    "wilson_marker_count": (
        "APPROXIMATE interval, not a bootstrap. It is a Wilson score interval "
        "derived from how many informative markers your array actually reads, "
        "so it reflects sampling of markers only. It does not capture model "
        "error, panel error or reference bias, all of which are larger."
    ),
    "none": "No interval could be produced for this proportion.",
}

# Colour-blind-safe qualitative palette (Wong 2011), which is also what the
# interactive report's colour-blind mode is built around. Keyed by
# superpopulation so a painted chromosome stays legible when several
# populations from the same continent appear.
SUPERPOP_COLOURS: dict[str, str] = {
    "AFR":     "#E69F00",
    "AMR":     "#56B4E9",
    "EAS":     "#009E73",
    "EUR":     "#0072B2",
    "SAS":     "#D55E00",
    "OCE":     "#CC79A7",
    "CAS":     "#F0E442",
    "UNKNOWN": "#999999",
}

# Deterministic fallback for labels with no superpopulation. Assigned by sorted
# label so the same input always paints the same colours; a chromosome painting
# whose colours move between runs is unreadable.
FALLBACK_PALETTE: tuple[str, ...] = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#56B4E9", "#F0E442", "#7A5195",
    "#BC5090", "#003F5C", "#58508D", "#FFA600",
)

# Hashing a 40 GB panel VCF on every manifest request is not acceptable, so
# files above this size record their length instead and say so. A size-only
# record is a weaker claim than a content hash and is labelled as one.
HASH_MAX_BYTES: int = 64 * 1024 * 1024

# GRCh37 chromosome lengths, used only to scale a painted chromosome for
# display. Public domain reference metadata, no licence attached.
CHROMOSOME_LENGTHS: dict[str, int] = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}


# ---------------------------------------------------------------------------
# Mandatory caveats
#
# Same pattern as prs.MANDATORY_CAVEATS: module constants, attached verbatim to
# every result, so no UI change can quietly drop one.
# ---------------------------------------------------------------------------

CAVEAT_PANEL_BIAS = (
    "Reference panels under-represent non-European ancestries. That is true of "
    "every openly licensed panel that exists, not just this one. Populations "
    "with few reference samples are estimated less accurately, and their "
    "ancestry is frequently absorbed into whichever better-sampled population "
    "sits nearest to them in the model."
)
CAVEAT_MODEL_DEPENDENT = (
    "Admixture proportions are a property of the model and the panel, not a "
    "fact about you. Change the reference populations and the percentages "
    "change. Two honest tools can return different numbers for the same person "
    "without either being wrong."
)
CAVEAT_ARRAY_DENSITY = (
    "A consumer array reads a few hundred thousand positions. That is enough to "
    "separate continental ancestries and not enough to separate populations "
    "that only differ at rare or fine-scale markers. Where your array cannot "
    "resolve a population, this report says so instead of returning zero."
)
CAVEAT_NOT_IDENTITY = (
    "Genetic ancestry is not ethnicity, nationality or identity. It describes "
    "which reference samples your DNA most resembles, nothing more."
)

MANDATORY_CAVEATS: tuple[str, ...] = (
    CAVEAT_PANEL_BIAS,
    CAVEAT_MODEL_DEPENDENT,
    CAVEAT_ARRAY_DENSITY,
    CAVEAT_NOT_IDENTITY,
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _as_float(value: Any) -> float | None:
    """Coerce to float, preserving 0.0 and rejecting junk and booleans."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _text(value: Any) -> str:
    return str(value or "").strip()


def _maybe_path(source: Any) -> Path | None:
    """Return ``source`` as an existing Path, or None if it is not one.

    A .Q matrix passed as text is still a str, and handing a multi-line string
    to ``Path.exists()`` raises OSError "File name too long" rather than
    returning False. Checking for a newline first is what keeps "text or a
    path" from being a trap.
    """
    if isinstance(source, Path):
        try:
            return source if source.is_file() else None
        except OSError:
            return None
    # An empty string resolves to Path("."), which exists and is a directory,
    # so is_file() rather than exists() is what keeps "" from being read as a
    # path to the working directory.
    if not isinstance(source, str) or not source or "\n" in source or len(source) > 4096:
        return None
    try:
        candidate = Path(source)
        return candidate if candidate.is_file() else None
    except OSError:
        return None


def _source_lines(source: Any) -> Iterable[Any] | None:
    """Normalise text, a path or an iterable of lines into lines to iterate.

    Every parser below accepts the same three input shapes, so the resolution
    order lives here once rather than three times. The str check must come
    before the Iterable check, because a str IS iterable and iterating one
    yields characters, which would silently turn a whole .Q matrix into
    one-character rows instead of failing.

    Returns None when the input is none of the three shapes. The caller decides
    what that means: parse_q_file raises, because a Q matrix that cannot be read
    must never be reported as an empty ancestry result.
    """
    resolved = _maybe_path(source)
    if resolved is not None:
        source = resolved.read_text(encoding="utf-8", errors="replace")
    if isinstance(source, str):
        return source.splitlines()
    if isinstance(source, Iterable):
        return source
    return None


_ARRAY_POP_BY_CODE: dict[str, dict] = {p["code"]: p for p in ARRAY_POPULATIONS}


def population_label(code: str) -> str:
    """Human label for a population code, from frequency.POPULATIONS if known.

    Reusing that table rather than defining a second one keeps the population
    vocabulary identical to the one the frequency layer already shows the user.
    A code the panel introduced and frequency.py does not know is returned as
    itself rather than blanked, because an unlabelled real population is better
    than a labelled wrong one.
    """
    entry = _ARRAY_POP_BY_CODE.get(_text(code).upper())
    return entry["label"] if entry else _text(code)


def population_superpop(code: str) -> str:
    entry = _ARRAY_POP_BY_CODE.get(_text(code).upper())
    return entry["superpop"] if entry else "UNKNOWN"


# ---------------------------------------------------------------------------
# Statistics, hand-written
#
# Runtime dependencies are Flask and requests. prs.py already hand-writes
# normal_cdf via math.erf rather than pulling in scipy; this follows that
# precedent rather than starting a second one.
# ---------------------------------------------------------------------------

def normal_cdf(z: float) -> float:
    """Standard normal CDF via math.erf. Mirrors prs.normal_cdf deliberately."""
    return 0.5 * (1.0 + math.erf(float(z) / math.sqrt(2.0)))


def inverse_normal_cdf(p: float) -> float:
    """Inverse standard normal CDF, by bisection on :func:`normal_cdf`.

    Bisection rather than a rational approximation because it is ten lines, has
    no magic constants to mistype, and 200 iterations over a bracket of plus or
    minus 40 converges far past the precision any percentage needs. Speed is
    irrelevant: this is called a handful of times per report.
    """
    value = _as_float(p)
    if value is None or not (0.0 < value < 1.0):
        raise AncestryError(f"inverse_normal_cdf needs 0 < p < 1, got {p!r}")
    low, high = -40.0, 40.0
    for _ in range(200):
        mid = (low + high) / 2.0
        if normal_cdf(mid) < value:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def z_for_level(level: float = CI_LEVEL) -> float:
    """Two-sided z multiplier for a confidence level. 0.95 gives about 1.96."""
    value = _as_float(level)
    if value is None or not (0.0 < value < 1.0):
        raise AncestryError(f"confidence level must be between 0 and 1, got {level!r}")
    return inverse_normal_cdf(1.0 - (1.0 - value) / 2.0)


def wilson_interval(proportion: float, n: int,
                    *, level: float = CI_LEVEL) -> tuple[float, float]:
    """Wilson score interval for a proportion observed over n markers.

    Chosen over the normal (Wald) interval because Wald produces intervals that
    run below zero and above one for proportions near the ends, and an ancestry
    proportion of "minus 2 percent" is not a number anybody can be shown.

    ``n`` of zero returns (0.0, 1.0): with no informative markers the estimate
    is unconstrained, and saying so is more honest than returning a tight
    interval around a number that rests on nothing.
    """
    p = _as_float(proportion)
    if p is None:
        return (0.0, 1.0)
    p = min(1.0, max(0.0, p))
    count = int(n or 0)
    if count <= 0:
        return (0.0, 1.0)

    z = z_for_level(level)
    denom = 1.0 + (z * z) / count
    centre = (p + (z * z) / (2.0 * count)) / denom
    spread = (z / denom) * math.sqrt(
        (p * (1.0 - p) / count) + (z * z) / (4.0 * count * count)
    )
    return (round(max(0.0, centre - spread), 6), round(min(1.0, centre + spread), 6))


def percentile_interval(replicates: Sequence[float],
                        *, level: float = CI_LEVEL) -> tuple[float, float] | None:
    """Percentile bootstrap interval from a list of replicate estimates.

    Returns None for fewer than two replicates, because an interval computed
    from one number is not an interval and returning (x, x) would claim
    certainty that does not exist.
    """
    values = sorted(v for v in (_as_float(r) for r in (replicates or []))
                    if v is not None)
    if len(values) < 2:
        return None
    tail = (1.0 - float(level)) / 2.0
    lower_index = max(0, int(math.floor(tail * (len(values) - 1))))
    upper_index = min(len(values) - 1, int(math.ceil((1.0 - tail) * (len(values) - 1))))
    return (round(values[lower_index], 6), round(values[upper_index], 6))


# ---------------------------------------------------------------------------
# Output parsing
#
# fastmixture writes a whitespace-delimited .Q matrix, one row per sample and
# one column per ancestry component, alongside the PLINK .fam that fixes the
# row order. Neither file carries column names: the mapping from column to
# population comes from the reference panel used to build the model, which is
# exactly why this module refuses to invent labels when the panel does not
# supply them.
# ---------------------------------------------------------------------------

def parse_q_file(source: Any) -> list[list[float]]:
    """Parse a .Q ancestry-proportion matrix into rows of floats.

    Accepts text, a path or an iterable of lines. Blank lines are skipped.
    Raises AncestryError for a non-numeric field or a ragged matrix, because a
    Q matrix whose rows have different widths means the run is corrupt and
    guessing which column is which would produce a plausible wrong ancestry
    report, the worst possible outcome.
    """
    lines = _source_lines(source)
    if lines is None:
        raise AncestryError("parse_q_file needs text, a path or an iterable of lines")

    rows: list[list[float]] = []
    width: int | None = None
    for number, line in enumerate(lines, start=1):
        text = _text(line)
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        row: list[float] = []
        for field in fields:
            value = _as_float(field)
            if value is None:
                raise AncestryError(
                    f"non-numeric value {field!r} on line {number} of the .Q matrix"
                )
            row.append(value)
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise AncestryError(
                f"ragged .Q matrix: line {number} has {len(row)} columns, "
                f"expected {width}"
            )
        rows.append(row)
    return rows


def parse_fam(source: Any) -> list[dict]:
    """Parse a PLINK .fam file, which fixes the row order of the .Q matrix.

    Six whitespace-separated columns: family id, individual id, paternal id,
    maternal id, sex, phenotype. Short rows are tolerated and padded, because
    several tools emit four-column .fam files and refusing them would fail a
    run over a field this module does not use.
    """
    lines = _source_lines(source)
    if lines is None:
        return []

    out: list[dict] = []
    for line in lines:
        text = _text(line)
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        while len(fields) < 6:
            fields.append("0")
        out.append({
            "fid": fields[0], "iid": fields[1],
            "father": fields[2], "mother": fields[3],
            "sex": fields[4], "phenotype": fields[5],
        })
    return out


def parse_population_map(source: Any) -> dict:
    """Parse the panel's populations.tsv into a population summary.

    Expected columns, tab or whitespace separated, with an optional header:
    sample id, population code, superpopulation code. Returns::

        {"samples": {sample_id: population},
         "populations": [{"code", "label", "superpop", "samples"}],
         "order": [population codes in first-seen order]}

    First-seen order is preserved because it is the only ordering the panel
    itself asserts, and reordering it would break the column-to-population
    mapping of a .Q matrix built against that panel.
    """
    lines = _source_lines(source)
    if lines is None:
        lines = []

    samples: dict[str, str] = {}
    counts: dict[str, int] = {}
    superpops: dict[str, str] = {}
    order: list[str] = []

    for line in lines:
        text = _text(line)
        if not text or text.startswith("#"):
            continue
        fields = [f for f in re.split(r"[\t,]|\s{1,}", text) if f]
        if len(fields) < 2:
            continue
        if fields[0].lower() in ("sample", "sample_id", "iid", "id"):
            continue  # header row
        sample, code = fields[0], fields[1].upper()
        superpop = fields[2].upper() if len(fields) > 2 else population_superpop(code)
        samples[sample] = code
        if code not in counts:
            counts[code] = 0
            order.append(code)
        counts[code] += 1
        superpops[code] = superpop

    populations = [{
        "code": code,
        "label": population_label(code),
        "superpop": superpops.get(code, "UNKNOWN"),
        "samples": counts[code],
    } for code in order]

    return {"samples": samples, "populations": populations, "order": order}


# ---------------------------------------------------------------------------
# Marker coverage
#
# This is the differentiator made concrete. Nobody else tells the user which
# populations their array is physically capable of resolving.
# ---------------------------------------------------------------------------

def resolvable(markers_read: int, informative_markers: int) -> bool:
    """True when an array reads enough of a population's informative markers.

    Both tests must pass: a fraction and an absolute count. Either alone is
    gameable. See MIN_RESOLVABLE_FRACTION and MIN_RESOLVABLE_MARKERS.
    """
    read = int(markers_read or 0)
    total = int(informative_markers or 0)
    if total <= 0:
        return False
    return read >= MIN_RESOLVABLE_MARKERS and (read / total) >= MIN_RESOLVABLE_FRACTION


def marker_coverage(typed_rsids: Any, informative_markers: Any) -> dict:
    """Per-population marker coverage for one array.

    ``typed_rsids`` is any iterable of rsIDs the array actually reads, or a
    genotype map whose keys are those rsIDs. ``informative_markers`` maps a
    population code to the panel's ancestry-informative markers for it.

    Returns ``{population_code: record}`` where each record carries::

        {"population", "label", "superpop",
         "informative_markers", "markers_read", "coverage",
         "resolvable", "state", "reason"}

    ``state`` is "resolvable" or "not_resolvable". A population that is not
    resolvable must be reported as such and NEVER as 0.0 percent ancestry.
    Zero percent is a measurement: it says we looked and found none. Not
    resolvable says we could not look. They are different claims and this
    project does not conflate them anywhere.
    """
    if isinstance(typed_rsids, dict):
        available = {str(k).strip().lower() for k in typed_rsids}
    elif isinstance(typed_rsids, Iterable) and not isinstance(typed_rsids, (str, bytes)):
        available = {str(k).strip().lower() for k in typed_rsids}
    else:
        available = set()

    out: dict[str, dict] = {}
    for code, markers in (informative_markers or {}).items():
        wanted = [str(m).strip().lower() for m in (markers or [])]
        total = len(wanted)
        read = sum(1 for m in wanted if m in available)
        fraction = (read / total) if total else 0.0
        ok = resolvable(read, total)

        if total == 0:
            reason = (
                "The panel lists no ancestry-informative markers for this "
                "population, so nothing about it can be estimated."
            )
        elif ok:
            reason = ""
        elif read < MIN_RESOLVABLE_MARKERS:
            reason = (
                f"Your array reads {read} of this population's {total} "
                f"informative markers. At least {MIN_RESOLVABLE_MARKERS} are "
                f"needed before an estimate means anything."
            )
        else:
            reason = (
                f"Your array reads {read} of this population's {total} "
                f"informative markers ({fraction * 100:.1f} percent). At least "
                f"{MIN_RESOLVABLE_FRACTION * 100:.0f} percent is needed."
            )

        out[str(code).upper()] = {
            "population":          str(code).upper(),
            "label":               population_label(code),
            "superpop":            population_superpop(code),
            "informative_markers": total,
            "markers_read":        read,
            "coverage":            round(fraction, 4),
            "resolvable":          ok,
            "state":               "resolvable" if ok else "not_resolvable",
            "reason":              reason,
        }
    return out


# ---------------------------------------------------------------------------
# PLINK input writing
#
# fastmixture consumes PLINK bed/bim/fam. Writing the binary .bed here rather
# than shelling out to plink avoids a second external tool dependency for a
# format that is 3 header bytes and 2 bits per genotype.
# ---------------------------------------------------------------------------

_BED_MAGIC = bytes([0x6C, 0x1B, 0x01])   # magic, magic, SNP-major mode

_NOCALL = {"", "N", "-", "--", "0", "D", "I", "."}


def _norm_allele(allele: Any) -> str:
    a = str(allele or "").strip().upper()
    return "" if a in _NOCALL else a


def _bed_code(a1: str, a2: str, allele_a: str, allele_b: str) -> int:
    """Two-bit PLINK genotype code for one sample at one marker.

    PLINK's encoding, which is not the obvious one and is worth stating:
        00 homozygous for the FIRST allele
        01 missing
        10 heterozygous
        11 homozygous for the SECOND allele
    """
    x, y = _norm_allele(a1), _norm_allele(a2)
    if not x or not y:
        return 0b01
    if x == y == allele_a:
        return 0b00
    if x == y == allele_b:
        return 0b11
    if {x, y} == {allele_a, allele_b}:
        return 0b10
    return 0b01


def write_plink(records: Any, prefix: str | Path,
                *, sample: str = "SAMPLE", family: str = "DNAINSIGHT") -> dict:
    """Write a single-sample PLINK trio of .bed, .bim and .fam files.

    ``records`` is an iterable of dicts carrying chromosome, position, rsid,
    allele1 and allele2, or a dict keyed by rsID with those as values. Records
    without coordinates are skipped and counted, never given an invented
    position.

    Returns ``{"prefix", "bed", "bim", "fam", "written", "skipped"}``.
    """
    rows: list[dict] = []
    skipped: list[dict] = []

    if isinstance(records, dict):
        candidates = []
        for rsid, value in records.items():
            if isinstance(value, dict):
                item = dict(value)
                item.setdefault("rsid", rsid)
                candidates.append(item)
            elif isinstance(value, str) and len(value.strip()) == 2:
                candidates.append({"rsid": rsid, "allele1": value.strip()[0],
                                   "allele2": value.strip()[1]})
            else:
                candidates.append({"rsid": rsid})
    elif isinstance(records, Iterable):
        candidates = [dict(r) for r in records if isinstance(r, dict)]
    else:
        candidates = []

    for item in candidates:
        chrom = _text(item.get("chromosome") or item.get("chrom")).upper().replace("CHR", "")
        position = _as_float(item.get("position") or item.get("pos"))
        a1 = _norm_allele(item.get("allele1"))
        a2 = _norm_allele(item.get("allele2"))
        genotype = item.get("genotype")
        if (not a1 or not a2) and isinstance(genotype, str) and len(genotype.strip()) == 2:
            a1, a2 = _norm_allele(genotype.strip()[0]), _norm_allele(genotype.strip()[1])
        if not chrom or position is None:
            skipped.append({"rsid": _text(item.get("rsid")).lower(),
                            "reason": "no chromosome or position"})
            continue
        rows.append({
            "rsid": _text(item.get("rsid")).lower() or ".",
            "chromosome": chrom, "position": int(position),
            "allele1": a1, "allele2": a2,
            "cm": _as_float(item.get("cm")) or 0.0,
        })

    rows.sort(key=lambda r: (_chrom_sort(r["chromosome"]), r["position"]))

    base = Path(prefix)
    base.parent.mkdir(parents=True, exist_ok=True)
    bim_path = Path(str(base) + ".bim")
    fam_path = Path(str(base) + ".fam")
    bed_path = Path(str(base) + ".bed")

    with open(bim_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in rows:
            alleles = sorted({a for a in (r["allele1"], r["allele2"]) if a})
            allele_a = alleles[0] if alleles else "0"
            allele_b = alleles[1] if len(alleles) > 1 else "0"
            r["_a"], r["_b"] = allele_a, allele_b
            fh.write(f"{r['chromosome']}\t{r['rsid']}\t{r['cm']:g}\t"
                     f"{r['position']}\t{allele_a}\t{allele_b}\n")

    with open(fam_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{family} {sample} 0 0 0 -9\n")

    with open(bed_path, "wb") as fh:
        fh.write(_BED_MAGIC)
        for r in rows:
            # One sample means one byte per marker: the genotype in the low two
            # bits and the remaining three sample slots left as missing.
            code = _bed_code(r["allele1"], r["allele2"], r["_a"], r["_b"])
            fh.write(struct.pack("B", code | 0b01010100))

    return {
        "prefix": str(base), "bed": str(bed_path),
        "bim": str(bim_path), "fam": str(fam_path),
        "written": len(rows), "skipped": skipped,
    }


def _chrom_sort(chrom: str) -> tuple:
    text = str(chrom or "").upper()
    if text.isdigit():
        return (0, int(text), "")
    return (1, 0, text)


# ---------------------------------------------------------------------------
# Degraded payloads
# ---------------------------------------------------------------------------

def panel_unavailable(panel_id: str, capability: str, *, detail: str = "") -> dict:
    """Degraded payload for a reference panel that has not been built.

    Kept separate from external.unavailable for the same reason imputation.py
    keeps them separate: "install a tool" and "build a panel" are different
    instructions, and a user sent to the wrong one loses an evening and stops
    believing the error messages.
    """
    st = external.panel_status(panel_id)
    name = st.get("name", panel_id)
    missing = [f for f in st.get("files_expected", [])
               if f not in set(st.get("files_present", []))]

    if st.get("state") == "unknown":
        reason = (
            f"Unknown reference panel {panel_id!r}. Known panels: "
            f"{', '.join(sorted(external.PANELS))}."
        )
    elif st.get("state") == "partial":
        reason = (
            f"The reference panel '{name}' is only partially built. Missing: "
            f"{', '.join(missing)}. Nothing needs to be installed; the panel "
            f"build needs finishing."
        )
    else:
        reason = (
            f"The reference panel '{name}' has not been built. Ancestry "
            f"estimation compares you against reference samples, and with no "
            f"reference samples there is nothing to compare against. This is a "
            f"DATA problem, not a tool problem."
        )

    return {
        "available":     False,
        "capability":    capability,
        "problem":       "panel_missing",
        "state":         f"panel_{st.get('state')}",
        "panel":         st.get("id", panel_id),
        "panel_name":    name,
        "panel_state":   st.get("state"),
        "reason":        detail or reason,
        "not_attempted": True,
        "results":       [],
        "proportions":   [],
        "path":          st.get("path"),
        "files_present": st.get("files_present", []),
        "files_missing": missing,
        "licence":       st.get("licence", ""),
        "ethics_gate":   st.get("ethics_gate", False),
        "how_to_enable": {
            "what":  f"Build the reference panel '{name}'.",
            "where": st.get("path"),
            "files": st.get("files_expected", []),
            "steps": [
                f"Build the panel into {st.get('path')} using the reference "
                f"panel builder under data/. Source terms are recorded in "
                f"data/DATA_SOURCES.md and in external.PANELS.",
                "Panel data is large and lives outside the repository tree, in "
                "the same place as external tools.",
            ],
            "note": st.get("note", ""),
        },
        "caveats": list(MANDATORY_CAVEATS),
    }


def _tool_unavailable(tool_id: str, capability: str, *, extra: dict | None = None) -> dict:
    """external.guard's payload, tagged so callers can branch on the cause."""
    payload = external.unavailable(tool_id, capability)
    payload["problem"] = "tool_missing"
    payload["proportions"] = []
    payload["caveats"] = list(MANDATORY_CAVEATS)
    if extra:
        payload.update(extra)
    return payload


# ---------------------------------------------------------------------------
# Panel transparency
# ---------------------------------------------------------------------------

def _file_digest(path: Path) -> dict:
    """SHA-256 of a panel file, or an honest refusal for a large one."""
    try:
        size = path.stat().st_size
    except OSError:
        return {"present": False, "bytes": None, "sha256": None,
                "hashed": False, "reason": "file could not be read"}

    if size > HASH_MAX_BYTES:
        return {
            "present": True, "bytes": size, "sha256": None, "hashed": False,
            "reason": (
                f"file is larger than the {HASH_MAX_BYTES} byte inline hashing "
                f"limit, so its size is recorded instead. A size is a weaker "
                f"claim than a content hash and is labelled as one."
            ),
        }

    digest = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b""):
                digest.update(block)
    except OSError:
        return {"present": True, "bytes": size, "sha256": None,
                "hashed": False, "reason": "file could not be read"}
    return {"present": True, "bytes": size, "sha256": digest.hexdigest(),
            "hashed": True, "reason": ""}


def panel_manifest(panel: str = DEFAULT_PANEL) -> dict:
    """Everything DNAInsight knows about a reference panel, as data.

    This is the "we ship the model" claim made concrete and testable: the
    populations, the sample count per population, the source, the licence, the
    build, the marker count and a content hash. An incumbent that says
    "4,500 regions" cannot produce this dict for its own product.

    Safe to call when nothing has been built: everything countable comes back
    as None or empty, ``available`` is False, and ``how_to_enable`` explains
    what is missing. None is used rather than 0 throughout, because a panel
    with zero populations and a panel that has not been built are different
    states and this module does not conflate them.
    """
    st = external.panel_status(panel)
    entry = external.PANELS.get(_text(panel).lower(), {})

    # An unknown panel id has no path at all. Falling through would read
    # populations.tsv relative to the working directory, which is how a
    # manifest ends up quietly describing whatever happened to be lying around.
    if st.get("state") == "unknown" or not st.get("path"):
        return {
            "panel": st.get("id", panel), "name": "", "purpose": "",
            "available": False, "state": "unknown", "path": None,
            "source": "", "licence": "", "spdx": None, "commercial_ok": None,
            "licence_verified": "", "ethics_gate": False, "build": "",
            "populations": [], "population_count": None, "sample_count": None,
            "marker_count": None, "files": [], "content_hash": None,
            "excluded": {}, "note": "",
            "reason": (
                f"Unknown reference panel {panel!r}. Known panels: "
                f"{', '.join(sorted(external.PANELS))}."
            ),
            "how_to_enable": None,
            "transparency": (
                "DNAInsight publishes the panel, the model and the "
                "per-population marker coverage behind every ancestry number."
            ),
        }

    base = Path(st["path"])

    files: list[dict] = []
    for name in st.get("files_expected", []):
        candidate = base / name
        if candidate.exists():
            record = _file_digest(candidate)
            record["name"] = name
            files.append(record)
        else:
            files.append({"name": name, "present": False, "bytes": None,
                          "sha256": None, "hashed": False,
                          "reason": "not built"})

    populations: list[dict] = []
    sample_count: int | None = None
    marker_count: int | None = None
    build = ""

    pop_file = base / "populations.tsv"
    if pop_file.exists():
        try:
            parsed = parse_population_map(pop_file)
            populations = parsed["populations"]
            sample_count = sum(p["samples"] for p in populations)
        except OSError:
            populations = []

    # The .map file is one line per marker, which makes a marker count cheap
    # and exact. Counting lines in the VCF instead would mean decompressing
    # tens of gigabytes to answer a question the map already answers.
    map_file = base / "panel.map"
    if map_file.exists():
        try:
            with open(map_file, "r", encoding="utf-8", errors="replace") as fh:
                marker_count = sum(1 for line in fh if _text(line)
                                   and not line.startswith("#"))
        except OSError:
            marker_count = None

    build_file = base / "BUILD.txt"
    if build_file.exists():
        try:
            build = build_file.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            build = ""

    # One hash over the per-file hashes, so a caller can compare two builds
    # with a single string. Files that were too large to hash contribute their
    # size, which is stated in the manifest so the weaker guarantee is visible.
    hashed_parts = [
        f"{f['name']}:{f.get('sha256') or ('size=' + str(f.get('bytes')))}"
        for f in files if f.get("present")
    ]
    content_hash = (hashlib.sha256("|".join(sorted(hashed_parts)).encode("utf-8"))
                    .hexdigest() if hashed_parts else None)

    return {
        "panel":             st.get("id", panel),
        "name":              st.get("name", panel),
        "purpose":           st.get("purpose", ""),
        "available":         bool(st.get("available")),
        "state":             st.get("state"),
        "path":              st.get("path"),
        "source":            entry.get("name", st.get("name", panel)),
        "licence":           st.get("licence", ""),
        "spdx":              entry.get("spdx"),
        "commercial_ok":     st.get("commercial_ok"),
        "licence_verified":  st.get("licence_verified", ""),
        "ethics_gate":       st.get("ethics_gate", False),
        "build":             build,
        "populations":       populations,
        "population_count":  len(populations) if populations else None,
        "sample_count":      sample_count,
        "marker_count":      marker_count,
        "files":             files,
        "content_hash":      content_hash,
        "excluded":          st.get("excluded", {}),
        "note":              st.get("note", ""),
        "how_to_enable":     None if st.get("available") else {
            "what":  f"Build the reference panel '{st.get('name', panel)}'.",
            "where": st.get("path"),
            "files": st.get("files_expected", []),
        },
        "transparency": (
            "DNAInsight publishes the panel, the model and the per-population "
            "marker coverage behind every ancestry number. If a field below is "
            "empty it is because the panel has not been built, not because the "
            "information is withheld."
        ),
    }


# ---------------------------------------------------------------------------
# Global ancestry
# ---------------------------------------------------------------------------

def _proportion_records(row: Sequence[float],
                        labels: Sequence[str],
                        coverage: dict,
                        *,
                        level: float,
                        replicates: dict | None,
                        tool_intervals: dict | None) -> tuple[list[dict], list[dict]]:
    """Build the per-population proportion records and the not-resolvable list.

    A population whose markers the array cannot resolve is REMOVED from the
    reported proportions and returned separately with ``proportion`` None. It
    is not reported as zero. Its estimated mass is reported in
    ``withheld_mass`` by the caller so the arithmetic stays visible rather than
    silently vanishing.
    """
    reported: list[dict] = []
    withheld: list[dict] = []

    for index, label in enumerate(labels):
        value = _as_float(row[index]) if index < len(row) else None
        code = _text(label).upper()
        cov = coverage.get(code, {})
        markers_read = int(cov.get("markers_read") or 0)
        informative = int(cov.get("informative_markers") or 0)
        can_resolve = cov.get("resolvable", True) if cov else True

        record = {
            "population":          code,
            "label":               population_label(code),
            "superpop":            population_superpop(code),
            "proportion":          round(value, 6) if value is not None else None,
            "percent":             round(value * 100.0, 2) if value is not None else None,
            "informative_markers": informative or None,
            "markers_read":        markers_read or None,
            "marker_coverage":     cov.get("coverage"),
            "resolvable":          bool(can_resolve),
            "state":               "resolvable" if can_resolve else "not_resolvable",
        }

        if not can_resolve:
            record["proportion"] = None
            record["percent"] = None
            record["ci_low"] = None
            record["ci_high"] = None
            record["ci_method"] = "none"
            record["ci_method_label"] = CI_METHOD_LABELS["none"]
            record["ci_level"] = level
            record["reason"] = cov.get("reason", "")
            record["raw_estimate"] = round(value, 6) if value is not None else None
            record["note"] = (
                "NOT RESOLVABLE. Your array does not read enough of this "
                "population's informative markers for an estimate to mean "
                "anything. This is not a report of zero percent: zero would "
                "mean we looked and found none."
            )
            withheld.append(record)
            continue

        low = high = None
        method = "none"
        if tool_intervals and code in tool_intervals:
            low, high = tool_intervals[code]
            method = "tool_reported"
        elif replicates and replicates.get(code):
            interval = percentile_interval(replicates[code], level=level)
            if interval is not None:
                low, high = interval
                method = "bootstrap_marker_subsets"
        if method == "none" and value is not None:
            low, high = wilson_interval(value, markers_read, level=level)
            method = "wilson_marker_count"

        record.update({
            "ci_low":           low,
            "ci_high":          high,
            "ci_percent_low":   round(low * 100.0, 2) if low is not None else None,
            "ci_percent_high":  round(high * 100.0, 2) if high is not None else None,
            "ci_method":        method,
            "ci_method_label":  CI_METHOD_LABELS.get(method, ""),
            "ci_level":         level,
            "reason":           "",
        })
        reported.append(record)

    reported.sort(key=lambda r: (-(r["proportion"] or 0.0), r["population"]))
    withheld.sort(key=lambda r: r["population"])
    return reported, withheld


def global_ancestry(genotypes: Any,
                    *,
                    panel: str = DEFAULT_PANEL,
                    mode: str = DEFAULT_MODE,
                    informative_markers: dict | None = None,
                    level: float = CI_LEVEL,
                    bootstrap: int = 0,
                    seed: int = 20260804,
                    workdir: str | Path | None = None,
                    sample: str = "SAMPLE",
                    timeout: int = external.DEFAULT_TIMEOUT) -> dict:
    """Estimate global ancestry proportions against a fixed reference panel.

    Runs fastmixture in projection mode by default: the reference panel is
    fixed and the user is scored against it. That is the same shape as a
    community admixture calculator, minus the unlicensed model file, and it is
    what makes the result reproducible by anybody holding the same panel.

    Degrades rather than raising, with a machine-readable ``problem`` key:

      * ``tool_missing``  fastmixture is absent, its runtime is absent, or its
        licence has not been accepted.
      * ``panel_missing``  the tool is ready but the reference panel is not
        built. Different fix, different message.
      * ``bad_mode``  an unknown mode was requested.
      * ``run_failed``  the tool ran and failed, with its stderr tail.

    On success returns ``proportions`` (resolvable populations only, each with
    an interval and the method that produced it), ``not_resolvable``
    (populations the array cannot see, with ``proportion`` None rather than
    zero), ``marker_coverage``, ``panel_manifest`` and ``caveats``.

    ``bootstrap`` is opt-in and defaults to 0 because each replicate is a full
    re-invocation of the external tool. With 0 replicates the intervals are
    Wilson approximations from the marker counts and are labelled as such;
    presenting an approximation as a bootstrap would be exactly the kind of
    quiet overclaim this module exists to avoid.

    The subprocess boundary in external.run IS the licence boundary.
    fastmixture is GPL-3.0 and is never imported here.
    """
    requested_mode = _text(mode).lower() or DEFAULT_MODE
    if requested_mode not in MODES:
        return {
            "available":     False,
            "capability":    "ancestry_global",
            "problem":       "bad_mode",
            "state":         "bad_mode",
            "reason": (
                f"Unknown ancestry mode {mode!r}. Known modes: "
                f"{', '.join(MODES)}."
            ),
            "not_attempted": True,
            "results": [], "proportions": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    blocked = external.guard("fastmixture", "ancestry_global")
    if blocked is not None:
        return _tool_unavailable("fastmixture", "ancestry_global",
                                 extra={"panel": panel, "mode": requested_mode})

    panel_state = external.panel_status(panel)
    if not panel_state.get("available"):
        payload = panel_unavailable(panel, "ancestry_global")
        payload["mode"] = requested_mode
        return payload

    manifest = panel_manifest(panel)
    base = Path(panel_state["path"])

    # Informative markers come from the panel when the builder wrote them.
    # Without them, marker coverage cannot be computed, and the honest response
    # is to say so rather than to assume every population is resolvable.
    markers = informative_markers
    coverage_known = markers is not None
    if markers is None:
        markers_file = base / "informative_markers.tsv"
        if markers_file.exists():
            markers = _read_informative_markers(markers_file)
            coverage_known = True
        else:
            markers = {}

    coverage = marker_coverage(genotypes, markers) if markers else {}

    owned = workdir is None
    work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="dnainsight_anc_"))
    work.mkdir(parents=True, exist_ok=True)

    written = write_plink(genotypes, work / "study", sample=sample)
    if not written["written"]:
        return {
            "available":     False,
            "capability":    "ancestry_global",
            "problem":       "no_input",
            "state":         "no_usable_calls",
            "reason": (
                "No usable calls were found. Ancestry estimation needs a "
                "chromosome and a position for every marker, and a plain rsID "
                "to genotype map does not carry them."
            ),
            "not_attempted": True,
            "results": [], "proportions": [],
            "skipped": written["skipped"],
            "caveats": list(MANDATORY_CAVEATS),
        }

    args = [
        "--bfile", written["prefix"],
        "--out", str(work / "ancestry"),
        "--seed", str(int(seed)),
    ]
    if requested_mode == "projection":
        args += ["--projection", "--reference", str(base / "panel.vcf.gz")]
    elif requested_mode == "supervised":
        args += ["--supervised", str(base / "populations.tsv")]

    try:
        external.run("fastmixture", args, timeout=timeout, cwd=work)
    except external.ExternalError as exc:
        return {
            "available":     False,
            "capability":    "ancestry_global",
            "problem":       "run_failed",
            "state":         "run_failed",
            "reason":        str(exc),
            "not_attempted": False,
            "results": [], "proportions": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    q_files = sorted(work.glob("ancestry*.Q"))
    if not q_files:
        return {
            "available":     False,
            "capability":    "ancestry_global",
            "problem":       "run_failed",
            "state":         "no_output",
            "reason": (
                "fastmixture reported success but wrote no .Q matrix. Treating "
                "that as a failure rather than as an empty ancestry result."
            ),
            "not_attempted": False,
            "results": [], "proportions": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    rows = parse_q_file(q_files[0])
    fam = parse_fam(written["fam"])
    row_index = 0
    if fam:
        ids = [f["iid"] for f in fam]
        if sample in ids:
            row_index = ids.index(sample)
    if row_index >= len(rows):
        row_index = 0
    row = rows[row_index] if rows else []

    labels, labelled = _column_labels(base, len(row))

    replicates = None
    if bootstrap and bootstrap > 0:
        replicates = _bootstrap_replicates(
            genotypes, base, work, labels,
            replicates=int(bootstrap), seed=seed, mode=requested_mode,
            timeout=timeout, sample=sample,
        )

    reported, withheld = _proportion_records(
        row, labels, coverage, level=level,
        replicates=replicates, tool_intervals=None,
    )

    total = sum(r["proportion"] or 0.0 for r in reported)
    withheld_mass = sum(r.get("raw_estimate") or 0.0 for r in withheld)

    caveats = ancestry_caveats(panel=panel)
    if not labelled:
        caveats.append(
            "The reference panel did not supply population labels for the "
            "model's components, so they are reported as numbered components. "
            "A numbered component is not a population and must not be read as "
            "one."
        )
    if not coverage_known:
        caveats.append(
            "This panel build ships no ancestry-informative marker list, so "
            "per-population marker coverage could not be computed and no "
            "population could be declared not resolvable. Treat every "
            "proportion below as unverified for resolvability."
        )

    result = {
        "available":       True,
        "capability":      "ancestry_global",
        "problem":         None,
        "state":           "ready",
        "tool":            "fastmixture",
        "tool_licence":    "GNU General Public License v3.0",
        "panel":           panel,
        "panel_name":      panel_state.get("name", panel),
        "panel_manifest":  manifest,
        "mode":            requested_mode,
        "components":      len(row),
        "components_labelled": labelled,
        "sample":          sample,
        "proportions":     reported,
        "results":         reported,
        "not_resolvable":  withheld,
        "marker_coverage": coverage,
        "marker_coverage_available": coverage_known,
        "reported_total":  round(total, 6),
        "withheld_mass":   round(withheld_mass, 6),
        "ci_level":        level,
        "ci_method":       ("bootstrap_marker_subsets" if replicates
                            else "wilson_marker_count"),
        "bootstrap_replicates": int(bootstrap or 0),
        "markers_used":    written["written"],
        "caveats":         caveats,
        "not_attempted":   False,
    }
    if owned:
        result["workdir"] = ""
    return result


def _read_informative_markers(path: Path) -> dict:
    """Read a panel's ancestry-informative marker list.

    Two columns, population code then rsID, one pair per line. Kept trivial on
    purpose: this file is produced by the panel builder in this repository, so
    a rich format would be complexity for its own sake.
    """
    out: dict[str, list[str]] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                text = _text(line)
                if not text or text.startswith("#"):
                    continue
                fields = text.split()
                if len(fields) < 2:
                    continue
                out.setdefault(fields[0].upper(), []).append(fields[1].strip().lower())
    except OSError:
        return {}
    return out


def _column_labels(base: Path, width: int) -> tuple[list[str], bool]:
    """Map .Q columns to population codes, or refuse to.

    A .Q matrix carries no column names. The mapping comes from the panel that
    built the model. When the panel does not state it, the columns are labelled
    ``component_1``, ``component_2`` and so on, and ``labelled`` comes back
    False so the caller can say plainly that these are model components rather
    than populations. Guessing the mapping would produce a confident, wrong
    ancestry report, which is worse than no report.
    """
    order_file = base / "q_columns.tsv"
    if order_file.exists():
        try:
            codes = [_text(line).split()[0].upper()
                     for line in order_file.read_text(encoding="utf-8",
                                                      errors="replace").splitlines()
                     if _text(line) and not line.startswith("#")]
            if len(codes) == width:
                return codes, True
        except (OSError, IndexError):
            pass

    pop_file = base / "populations.tsv"
    if pop_file.exists():
        try:
            order = parse_population_map(pop_file)["order"]
            if len(order) == width:
                return order, True
        except OSError:
            pass

    return [f"component_{i + 1}" for i in range(width)], False


def _bootstrap_replicates(genotypes: Any, base: Path, work: Path,
                          labels: Sequence[str], *, replicates: int,
                          seed: int, mode: str, timeout: int,
                          sample: str) -> dict:
    """Re-run the estimator on resampled marker subsets.

    Each replicate is a full external invocation, which is why this is opt-in.
    Failed replicates are dropped rather than retried: a bootstrap missing a
    few draws is still a bootstrap, and hanging a report on a flaky subprocess
    is not.
    """
    rng = random.Random(seed)
    records = _as_record_list(genotypes)
    if not records:
        return {}

    draws: dict[str, list[float]] = {str(code).upper(): [] for code in labels}
    for i in range(int(replicates)):
        subset = [records[rng.randrange(len(records))] for _ in range(len(records))]
        prefix = work / f"boot_{i}"
        written = write_plink(subset, prefix, sample=sample)
        if not written["written"]:
            continue
        args = ["--bfile", written["prefix"], "--out", str(work / f"boot_{i}_out"),
                "--seed", str(seed + i)]
        if mode == "projection":
            args += ["--projection", "--reference", str(base / "panel.vcf.gz")]
        try:
            external.run("fastmixture", args, timeout=timeout, cwd=work)
        except external.ExternalError:
            continue
        found = sorted(work.glob(f"boot_{i}_out*.Q"))
        if not found:
            continue
        try:
            rows = parse_q_file(found[0])
        except AncestryError:
            continue
        if not rows:
            continue
        for index, code in enumerate(labels):
            if index < len(rows[0]):
                draws[str(code).upper()].append(rows[0][index])
    return draws


def _as_record_list(genotypes: Any) -> list[dict]:
    """Coerce any accepted genotype shape into a list of record dicts."""
    if isinstance(genotypes, dict):
        out = []
        for rsid, value in genotypes.items():
            if isinstance(value, dict):
                item = dict(value)
                item.setdefault("rsid", rsid)
                out.append(item)
        return out
    if isinstance(genotypes, Iterable) and not isinstance(genotypes, (str, bytes)):
        return [dict(r) for r in genotypes if isinstance(r, dict)]
    return []


# ---------------------------------------------------------------------------
# Local ancestry
# ---------------------------------------------------------------------------

_PHASED_GT = re.compile(r"^[0-9.]+\|[0-9.]+")


def is_phased_vcf(source: Any, *, sample_lines: int = 200) -> bool:
    """True when a VCF's genotypes use the phased separator.

    Local ancestry assigns each of your two haplotypes independently, so it
    needs to know which allele sits on which copy. Unphased input produces
    output that looks like a painted chromosome and is not one. Checking is
    cheap; discovering it from a wrong result is not.
    """
    path = _maybe_path(source)
    if path is not None:
        try:
            opener = (gzip.open(path, "rt", encoding="utf-8", errors="replace")
                      if path.suffix == ".gz"
                      else open(path, "r", encoding="utf-8", errors="replace"))
            with opener as fh:
                return _scan_phase(fh, sample_lines)
        except OSError:
            return False
    if isinstance(source, str):
        return _scan_phase(source.splitlines(), sample_lines)
    if isinstance(source, Iterable):
        return _scan_phase(source, sample_lines)
    return False


def _scan_phase(lines: Iterable[Any], limit: int) -> bool:
    seen = 0
    for line in lines:
        text = str(line or "").rstrip("\r\n")
        if not text or text.startswith("#"):
            continue
        fields = text.split("\t")
        if len(fields) < 10:
            fields = text.split()
        if len(fields) < 10:
            continue
        seen += 1
        if not _PHASED_GT.match(fields[9]):
            return False
        if seen >= limit:
            break
    return seen > 0


def local_ancestry(phased_vcf: Any,
                   *,
                   panel: str = DEFAULT_PANEL,
                   map_path: str | Path | None = None,
                   workdir: str | Path | None = None,
                   timeout: int = external.DEFAULT_TIMEOUT) -> dict:
    """Infer local (per-segment) ancestry with FLARE.

    Requires PHASED input. FLARE assigns ancestry to each haplotype separately,
    so unphased input silently produces a picture that looks like a painted
    chromosome without being one. That failure is caught here and reported as
    ``problem="input_not_phased"`` with the fix (phase first, with Beagle or
    SHAPEIT5), rather than being discovered later by a user who trusted the
    picture.

    Degrades with ``problem`` set to ``tool_missing``, ``panel_missing``,
    ``input_not_phased``, ``no_input`` or ``run_failed``.

    On success returns per-chromosome painted segments, ready for
    :func:`chromosome_painting`.

    FLARE is Apache-2.0, so no copyleft attaches, but it is still invoked
    across the subprocess boundary like everything else so there is one rule to
    understand rather than two.
    """
    blocked = external.guard("flare", "ancestry_local")
    if blocked is not None:
        return _tool_unavailable("flare", "ancestry_local",
                                 extra={"panel": panel, "segments": []})

    panel_state = external.panel_status(panel)
    if not panel_state.get("available"):
        payload = panel_unavailable(panel, "ancestry_local")
        payload["segments"] = []
        return payload

    study = _maybe_path(phased_vcf)
    if study is None:
        return {
            "available":     False,
            "capability":    "ancestry_local",
            "problem":       "no_input",
            "state":         "no_input",
            "reason": (
                "Local ancestry needs a phased study VCF on disk. No readable "
                "input file was supplied."
            ),
            "not_attempted": True,
            "results": [], "segments": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    if not is_phased_vcf(study):
        return {
            "available":     False,
            "capability":    "ancestry_local",
            "problem":       "input_not_phased",
            "state":         "input_not_phased",
            "reason": (
                "The supplied VCF is not phased. Local ancestry paints each of "
                "your two haplotypes separately, which is impossible without "
                "knowing which allele sits on which copy. Phase the file first "
                "with Beagle or SHAPEIT5, then run this again. Running it "
                "anyway would produce a picture that looks correct and is not."
            ),
            "not_attempted": True,
            "results": [], "segments": [],
            "how_to_enable": {
                "what": "Phase the study VCF before painting it.",
                "steps": [
                    "Run the phasing capability (SHAPEIT5, MIT, or Beagle) "
                    "against the same reference panel.",
                    "Feed the phased output back into local_ancestry.",
                ],
            },
            "caveats": list(MANDATORY_CAVEATS),
        }

    base = Path(panel_state["path"])
    owned = workdir is None
    work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="dnainsight_flare_"))
    work.mkdir(parents=True, exist_ok=True)

    out_prefix = work / "painted"
    args = [
        f"ref={base / 'panel.vcf.gz'}",
        f"ref-panel={base / 'populations.tsv'}",
        f"gt={study}",
        f"map={map_path or (base / 'panel.map')}",
        f"out={out_prefix}",
    ]

    try:
        external.run("flare", args, timeout=timeout, cwd=work)
    except external.ExternalError as exc:
        return {
            "available":     False,
            "capability":    "ancestry_local",
            "problem":       "run_failed",
            "state":         "run_failed",
            "reason":        str(exc),
            "not_attempted": False,
            "results": [], "segments": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    produced = None
    for candidate in (Path(str(out_prefix) + ".anc.vcf.gz"),
                      Path(str(out_prefix) + ".anc.vcf"),
                      Path(str(out_prefix) + ".vcf.gz")):
        if candidate.exists():
            produced = candidate
            break
    if produced is None:
        return {
            "available":     False,
            "capability":    "ancestry_local",
            "problem":       "run_failed",
            "state":         "no_output",
            "reason": (
                f"FLARE reported success but produced no ancestry VCF at "
                f"{out_prefix}.anc.vcf.gz."
            ),
            "not_attempted": False,
            "results": [], "segments": [],
            "caveats": list(MANDATORY_CAVEATS),
        }

    labels = parse_population_map(base / "populations.tsv")["order"] \
        if (base / "populations.tsv").exists() else []
    calls = _parse_flare_vcf(produced, labels)
    segments = segments_from_calls(calls)

    by_chromosome: dict[str, list[dict]] = {}
    for seg in segments:
        by_chromosome.setdefault(seg["chromosome"], []).append(seg)

    return {
        "available":      True,
        "capability":     "ancestry_local",
        "problem":        None,
        "state":          "ready",
        "tool":           "FLARE",
        "tool_licence":   "Apache License 2.0",
        "panel":          panel,
        "panel_name":     panel_state.get("name", panel),
        "phased":         True,
        "segments":       segments,
        "results":        segments,
        "chromosomes":    [
            {"chromosome": chrom, "segments": segs,
             "length": CHROMOSOME_LENGTHS.get(chrom)}
            for chrom, segs in sorted(by_chromosome.items(),
                                      key=lambda kv: _chrom_sort(kv[0]))
        ],
        "markers":        len(calls),
        "ancestries":     sorted({s["ancestry"] for s in segments}),
        "caveats":        ancestry_caveats(panel=panel),
        "not_attempted":  False,
        "workdir":        "" if owned else str(work),
    }


def _parse_flare_vcf(path: Path, labels: Sequence[str]) -> list[dict]:
    """Extract per-marker haplotype ancestry assignments from FLARE output.

    ASSUMPTION, FLAGGED: FLARE writes the per-haplotype ancestry as the AN1 and
    AN2 FORMAT fields, as integer indices into the reference panel file's
    population order. This has NOT been verified against a real FLARE run in
    this environment because the tool is not installed here. If the field names
    differ, this function is the only place that needs changing.
    """
    opener = (gzip.open(path, "rt", encoding="utf-8", errors="replace")
              if path.suffix == ".gz"
              else open(path, "r", encoding="utf-8", errors="replace"))

    calls: list[dict] = []
    with opener as fh:
        for line in fh:
            text = line.rstrip("\r\n")
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0].replace("chr", "").replace("CHR", "").upper()
            position = _as_float(fields[1])
            keys = fields[8].split(":")
            values = fields[9].split(":")
            record = dict(zip(keys, values))
            an1 = record.get("AN1")
            an2 = record.get("AN2")
            if position is None or an1 is None or an2 is None:
                continue
            calls.append({
                "chromosome": chrom,
                "position": int(position),
                "haplotype1": _label_for_index(an1, labels),
                "haplotype2": _label_for_index(an2, labels),
            })
    return calls


def _label_for_index(value: Any, labels: Sequence[str]) -> str:
    """Translate an ancestry index into a population code, or keep the index."""
    index = _as_float(value)
    if index is None:
        return _text(value) or "UNKNOWN"
    i = int(index)
    if 0 <= i < len(labels):
        return str(labels[i]).upper()
    return f"component_{i + 1}"


def segments_from_calls(calls: Any) -> list[dict]:
    """Merge consecutive per-marker ancestry calls into contiguous segments.

    A run of markers assigned the same ancestry on the same haplotype is one
    segment. Merging happens per haplotype, because the two copies of a
    chromosome have independent ancestry and averaging them would invent a
    third thing that is neither.

    Segment boundaries are the positions of the first and last marker in the
    run, not the midpoints between runs. The true recombination breakpoint sits
    somewhere between two markers and this module does not know where; using
    marker positions understates each segment slightly and never claims
    precision the data does not support.
    """
    ordered = sorted(
        (c for c in (calls or []) if isinstance(c, dict)),
        key=lambda c: (_chrom_sort(str(c.get("chromosome") or "")),
                       int(c.get("position") or 0)),
    )

    segments: list[dict] = []
    for hap_key, hap_number in (("haplotype1", 1), ("haplotype2", 2)):
        current: dict | None = None
        for call in ordered:
            chrom = _text(call.get("chromosome")).upper()
            position = int(call.get("position") or 0)
            ancestry = _text(call.get(hap_key)).upper() or "UNKNOWN"

            if (current is not None and current["chromosome"] == chrom
                    and current["ancestry"] == ancestry):
                current["end"] = position
                current["markers"] += 1
                continue
            if current is not None:
                segments.append(current)
            current = {
                "chromosome": chrom, "haplotype": hap_number,
                "start": position, "end": position,
                "ancestry": ancestry, "markers": 1,
            }
        if current is not None:
            segments.append(current)

    for seg in segments:
        seg["span"] = max(0, seg["end"] - seg["start"])
    segments.sort(key=lambda s: (_chrom_sort(s["chromosome"]), s["haplotype"],
                                 s["start"]))
    return segments


# ---------------------------------------------------------------------------
# Chromosome painting
#
# Shape note: backend/interactive_report.py embeds its payload as a JSON island
# and draws charts as hand-written SVG with no charting library. So everything
# below is plain JSON-serialisable primitives, every segment carries the
# fractional offsets an SVG rect needs so the template does no arithmetic, and
# the colour key is data rather than CSS. This is intended to travel in the
# report's ``extras`` dict.
# ---------------------------------------------------------------------------

def _colour_for(ancestry: str, assigned: dict[str, str]) -> str:
    """Stable colour for an ancestry label.

    Superpopulation first, so two European populations do not come out as two
    unrelated colours. Otherwise a deterministic slot in FALLBACK_PALETTE keyed
    by insertion order of the sorted labels, so the same input always paints
    the same picture. A painting whose colours move between runs cannot be
    compared with the one printed last month.
    """
    code = _text(ancestry).upper()
    if code in assigned:
        return assigned[code]
    superpop = population_superpop(code)
    if superpop in SUPERPOP_COLOURS and superpop != "UNKNOWN":
        colour = SUPERPOP_COLOURS[superpop]
    elif code in SUPERPOP_COLOURS:
        colour = SUPERPOP_COLOURS[code]
    else:
        colour = FALLBACK_PALETTE[len(assigned) % len(FALLBACK_PALETTE)]
    assigned[code] = colour
    return colour


def chromosome_painting(local_result: Any) -> dict:
    """Turn a :func:`local_ancestry` result into SVG-ready segment data.

    Returns::

        {"available": bool,
         "chromosomes": [{"chromosome", "length", "segments": [...]}],
         "segments": [...],
         "colour_key": [{"ancestry", "label", "colour"}],
         "scale": {"units", "max_length", "lengths_source"},
         "totals": [{"ancestry", "span", "fraction"}],
         "caveats": [...]}

    Every segment carries ``x1`` and ``x2``, its start and end as fractions of
    the chromosome length, so a template can draw a rect without doing
    arithmetic on genomic coordinates. A chromosome whose length is unknown
    gets ``x1``/``x2`` of None rather than a fabricated scale, because a bar
    drawn to an invented length is a lie in picture form.

    A degraded local_result passes straight through with ``available`` False,
    so a caller can render this unconditionally.
    """
    if not isinstance(local_result, dict):
        return {"available": False, "chromosomes": [], "segments": [],
                "colour_key": [], "scale": {}, "totals": [],
                "reason": "no local ancestry result supplied",
                "caveats": list(MANDATORY_CAVEATS)}

    if not local_result.get("available"):
        return {
            "available": False,
            "chromosomes": [], "segments": [], "colour_key": [],
            "scale": {}, "totals": [],
            "problem": local_result.get("problem"),
            "reason": local_result.get("reason", ""),
            "caveats": local_result.get("caveats") or list(MANDATORY_CAVEATS),
        }

    segments = [s for s in (local_result.get("segments") or [])
                if isinstance(s, dict)]
    assigned: dict[str, str] = {}
    for code in sorted({_text(s.get("ancestry")).upper() for s in segments}):
        _colour_for(code, assigned)

    painted: list[dict] = []
    for seg in segments:
        chrom = _text(seg.get("chromosome")).upper()
        length = CHROMOSOME_LENGTHS.get(chrom)
        start = int(seg.get("start") or 0)
        end = int(seg.get("end") or 0)
        ancestry = _text(seg.get("ancestry")).upper() or "UNKNOWN"
        painted.append({
            "chromosome": chrom,
            "haplotype":  int(seg.get("haplotype") or 1),
            "start":      start,
            "end":        end,
            "span":       max(0, end - start),
            "markers":    int(seg.get("markers") or 0),
            "ancestry":   ancestry,
            "label":      population_label(ancestry),
            "colour":     _colour_for(ancestry, assigned),
            "x1":         round(start / length, 6) if length else None,
            "x2":         round(min(end, length) / length, 6) if length else None,
        })

    by_chromosome: dict[str, list[dict]] = {}
    for seg in painted:
        by_chromosome.setdefault(seg["chromosome"], []).append(seg)

    totals: dict[str, int] = {}
    for seg in painted:
        totals[seg["ancestry"]] = totals.get(seg["ancestry"], 0) + seg["span"]
    grand = sum(totals.values())

    return {
        "available": True,
        "chromosomes": [
            {"chromosome": chrom,
             "length": CHROMOSOME_LENGTHS.get(chrom),
             "segments": sorted(segs, key=lambda s: (s["haplotype"], s["start"]))}
            for chrom, segs in sorted(by_chromosome.items(),
                                      key=lambda kv: _chrom_sort(kv[0]))
        ],
        "segments": painted,
        "colour_key": [
            {"ancestry": code, "label": population_label(code), "colour": colour}
            for code, colour in sorted(assigned.items())
        ],
        "scale": {
            "units": "bp",
            "build": "GRCh37",
            "max_length": max(CHROMOSOME_LENGTHS.values()),
            "lengths_source": "GRCh37 chromosome lengths, public domain reference metadata",
        },
        "totals": [
            {"ancestry": code, "label": population_label(code),
             "span": span,
             "fraction": round(span / grand, 6) if grand else None}
            for code, span in sorted(totals.items(), key=lambda kv: -kv[1])
        ],
        "caveats": (local_result.get("caveats") or list(MANDATORY_CAVEATS)) + [
            "Segment boundaries are drawn at the outermost marker in each run. "
            "The real crossover sits somewhere in the gap between two markers, "
            "and on a consumer array that gap can be tens of thousands of bases."
        ],
    }


# ---------------------------------------------------------------------------
# Caveats
# ---------------------------------------------------------------------------

def ancestry_caveats(*, panel: str = DEFAULT_PANEL) -> list[str]:
    """Return every caveat that must accompany an ancestry result.

    The four MANDATORY_CAVEATS always come first and always in full, following
    prs.build_caveats. The panel's own recorded note is appended when the panel
    has one, because external.PANELS states the ancestry-bias limitation in the
    panel's own words and losing that would be losing a licence-adjacent
    statement somebody wrote deliberately.
    """
    out = list(MANDATORY_CAVEATS)
    note = _text(external.panel_status(panel).get("note"))
    if note and note not in out:
        out.append(note)
    return out


# ---------------------------------------------------------------------------
# Haplogroups live elsewhere
# ---------------------------------------------------------------------------

def haplogroup_note() -> dict:
    """Report whether Y and mtDNA haplogroup calling is available.

    Haplogroups are backend/haplogroups.py's job, not this module's. The import
    is local and defensive so the two modules can land in either order and so a
    missing haplogroups.py degrades to a documented absence rather than an
    ImportError at scan time.

    Haplogroups are also a genuinely different claim from admixture: a
    haplogroup traces one unbroken line (all-male or all-female) and says
    nothing about the rest of a family tree. Keeping them in a separate module
    keeps that distinction structural rather than editorial.
    """
    try:
        from . import haplogroups  # noqa: F401  (probe import, deliberate)
    except ImportError:
        return {
            "available": False,
            "reason": (
                "Y and mitochondrial haplogroup calling is provided by "
                "backend/haplogroups.py, which is not present in this build."
            ),
            "note": (
                "A haplogroup traces a single unbroken line of descent and is "
                "not a summary of your ancestry. It is reported separately for "
                "that reason."
            ),
        }
    return {
        "available": True,
        "module": "backend.haplogroups",
        "note": (
            "A haplogroup traces a single unbroken line of descent and is not "
            "a summary of your ancestry."
        ),
    }
