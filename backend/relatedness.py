"""
relatedness.py -- IBD segments, relationship prediction and household analysis.

THE FRAMING, WHICH DECIDES WHAT THIS MODULE IS
----------------------------------------------
GEDmatch's moat is 1.5 million uploaded profiles. That is a network effect, not
a software moat. A tool that runs on one person's laptop over the kits that
person already has cannot copy it and should not pretend to try. Building a
worse matching database would be a worse product AND a worse privacy position.

So DNAInsight does not try. This module is PRIVATE FAMILY GENOMICS over the
kits the user has loaded: their own files, their parents' files, their
children's files. Nothing is uploaded, nothing is compared against strangers,
and no database of other people exists to be searched.

That is also the one form of DNA matching that is structurally safe from the
opt-out circumvention documented at GEDmatch in 2023 and MyHeritage in 2025. An
opt-out only protects you while the operator honours it and while the database
stays in the hands you agreed to. A comparison that never leaves your machine
has no operator to honour anything.

WHAT IS APPROXIMATE HERE, STATED UP FRONT
-----------------------------------------
Everything in the pure-Python path is an approximation and every payload says
so in a returned field rather than in documentation nobody reads:

  1. Segment detection is unphased IBS-based. Two people who share at least one
     allele at every marker across a stretch look identical by descent, and
     sometimes are not. The classic false positive is a long run of common
     homozygous genotypes in a low-diversity region. Real IBD callers use
     population allele frequencies, phasing or both. This one does not, so
     ``approximate: True`` is on every result and short segments are not
     trustworthy at all.
  2. Centimorgans are interpolated from a real genetic map when the user has
     built one, and estimated from a per-chromosome average rate when they have
     not. An estimated cM is flagged ``cm_estimated: True`` and must never be
     rendered as a measured one, because the recombination rate varies by more
     than an order of magnitude along a chromosome.
  3. Relationship prediction returns a RANGE of possibilities, never one
     confident answer. Shared cM alone genuinely cannot separate a half sibling
     from a grandparent from an aunt or uncle: all three sit on the same part
     of the distribution. Any tool that returns one of those three as the answer
     is guessing and not telling you.

EXTERNAL TOOLS
--------------
IBIS is used when the user has installed it. IBIS was chosen for this wave for
one reason: it is phase-free, so it runs on raw unphased array data with no
imputation step, which is exactly what a consumer file is. It is GPL-3.0.
Nothing here imports, links or vendors it. It runs through ``external.run``,
and the subprocess boundary IS the licence boundary: a separate process running
a program the user installed into their own home directory on explicit consent
does not relicense this MIT tree. With IBIS absent the pure-Python path runs
and the payload carries ``method: "builtin_ibs"``.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

try:  # the backend package is the normal import path
    from backend import external as _external
    from backend import merge as _merge
except ImportError:  # pragma: no cover - direct-module import fallback
    import external as _external  # type: ignore
    import merge as _merge  # type: ignore

external = _external
merge = _merge


__all__ = [
    "CHROM_LENGTHS_GRCH37", "AVERAGE_CM_PER_MB", "GENOME_AVERAGE_CM_PER_MB",
    "RELATIONSHIP_BANDS", "SHARED_CM_CAVEAT", "IBD_CAVEAT",
    "MIN_SNPS", "MIN_CM", "MAX_OPPOSITE_HOMOZYGOTES",
    "ibs_state", "load_genetic_map", "centimorgans",
    "shared_segments", "total_shared_cm", "longest_segment_cm",
    "predict_relationship", "detect_ibd", "parse_ibis_segments",
    "write_plink_map", "write_plink_ped",
    "phase_by_parents", "analyse_household", "chromosome_browser_data",
    "expected_for_roles",
]


# ---------------------------------------------------------------------------
# Defaults
#
# 300 SNPs and 7 cM are the conventional consumer thresholds. Below roughly
# 7 cM the false positive rate on unphased array data climbs fast enough that a
# segment stops being evidence of anything, which is why the floor exists and
# why it is a named constant rather than a magic number at a call site.
# ---------------------------------------------------------------------------

MIN_SNPS = 300
MIN_CM = 7.0
MAX_OPPOSITE_HOMOZYGOTES = 1

_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "D", "I", "?", "."}


# ---------------------------------------------------------------------------
# Chromosome lengths, GRCh37 / hg19. Public domain reference constants, used
# only to scale the chromosome browser bars.
# ---------------------------------------------------------------------------

CHROM_LENGTHS_GRCH37: dict[str, int] = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    "MT": 16569,
}


# ---------------------------------------------------------------------------
# Fallback recombination rates, cM per Mb, sex averaged.
#
# These are whole-chromosome averages: published total genetic length divided
# by physical length. They are ONLY used when the user has no genetic map
# installed, and any cM computed from them is flagged cm_estimated=True.
#
# Why that flag matters more than the numbers do: recombination is wildly
# non-uniform. Hotspots run an order of magnitude above the chromosome average
# and centromeric regions run far below it, so a 10 Mb stretch can be 2 cM or
# 25 cM. An average rate gives a defensible order of magnitude and nothing more.
#
# The data builder should confirm these against a published map before release.
# ---------------------------------------------------------------------------

AVERAGE_CM_PER_MB: dict[str, float] = {
    "1": 1.14, "2": 1.09, "3": 1.11, "4": 1.09, "5": 1.13, "6": 1.10,
    "7": 1.19, "8": 1.16, "9": 1.29, "10": 1.30, "11": 1.20, "12": 1.28,
    "13": 1.15, "14": 1.28, "15": 1.42, "16": 1.60, "17": 1.68, "18": 1.28,
    "19": 2.05, "20": 1.61, "21": 1.87, "22": 2.12, "X": 1.02,
}

GENOME_AVERAGE_CM_PER_MB = 1.16

IBD_CAVEAT = (
    "Segments were detected from unphased genotypes by identity by state. Two "
    "people who share one allele at every marker across a stretch look "
    "identical by descent here and are sometimes not, particularly across long "
    "runs of common homozygous genotypes. Treat short segments as noise."
)

SHARED_CM_CAVEAT = (
    "Shared centimorgans alone cannot tell a half sibling from a grandparent "
    "from an aunt or uncle. All three share roughly a quarter of their DNA and "
    "sit on the same part of the distribution. Ages, the number and length of "
    "the segments, and X chromosome sharing are what separate them, and none of "
    "those is a total. Any tool that hands you one confident answer from a "
    "total is guessing without telling you."
)


# ---------------------------------------------------------------------------
# Genotype access
# ---------------------------------------------------------------------------

def _pair(value: Any) -> tuple[str, str] | None:
    """Normalise a genotype to a sorted 2-tuple, or None for a no-call.

    Accepts the project's ``(allele1, allele2)`` tuple, a 2-character string
    and the merged-genotype dict shape, because household analysis walks all
    three at once: pooled primary entries are dicts, comparison rows are dicts
    with no coordinates, and a caller-supplied map is usually tuples.
    """
    if value is None:
        return None
    if isinstance(value, dict):
        raw: Sequence[Any] = (value.get("allele1"), value.get("allele2"))
    elif isinstance(value, (list, tuple)):
        raw = value
    else:
        text = str(value).strip().upper()
        for sep in ("/", "|", ";", " ", "\t"):
            text = text.replace(sep, "")
        raw = list(text)
    alleles = [str(a if a is not None else "").strip().upper() for a in raw]
    if len(alleles) != 2:
        return None
    if any(a in _NOCALL_ALLELES for a in alleles):
        return None
    return (alleles[0], alleles[1]) if alleles[0] <= alleles[1] else (alleles[1], alleles[0])


def ibs_state(gt_a: Any, gt_b: Any) -> int | None:
    """Identity by state between two genotypes: 2, 1, 0, or None for a no-call.

    At a biallelic SNP an IBS of 0 means opposite homozygotes, AA against GG.
    That cannot happen between two people who share that segment by descent, so
    it is the signal a run is over. IBS 0 from more than two alleles across the
    pair is a data problem rather than biology and is counted the same way.
    """
    left, right = _pair(gt_a), _pair(gt_b)
    if left is None or right is None:
        return None
    if left == right:
        return 2
    return 1 if set(left) & set(right) else 0


def _norm_chrom(value: Any) -> str:
    chrom = str(value if value is not None else "").strip().upper()
    if chrom.startswith("CHR"):
        chrom = chrom[3:]
    if chrom in ("23",):
        return "X"
    if chrom in ("24",):
        return "Y"
    if chrom in ("25",):
        return "X"
    if chrom in ("26", "M"):
        return "MT"
    return chrom


def _chrom_sort_key(chrom: str) -> tuple:
    try:
        return (0, int(chrom), chrom)
    except (TypeError, ValueError):
        return (1, 0, chrom)


def _coords(genotypes: dict, positions: dict | None) -> dict[str, tuple[str, int]]:
    """Build {rsid: (chromosome, position)} from whatever the caller supplied.

    A merged-genotype map carries its own coordinates. A comparison row does
    not, because ``merge.py`` deliberately stores only the call for a relative.
    So an explicit ``positions`` index wins, and the genotype map is used to
    fill any gap.
    """
    out: dict[str, tuple[str, int]] = {}
    for source in (genotypes or {}, ):
        for rsid, value in source.items():
            if not isinstance(value, dict):
                continue
            chrom = _norm_chrom(value.get("chromosome"))
            try:
                pos = int(value.get("position") or 0)
            except (TypeError, ValueError):
                pos = 0
            if chrom and pos > 0:
                out[str(rsid).strip().lower()] = (chrom, pos)
    for rsid, value in (positions or {}).items():
        key = str(rsid).strip().lower()
        if isinstance(value, dict):
            chrom = _norm_chrom(value.get("chromosome"))
            try:
                pos = int(value.get("position") or 0)
            except (TypeError, ValueError):
                pos = 0
        elif isinstance(value, (list, tuple)) and len(value) == 2:
            chrom = _norm_chrom(value[0])
            try:
                pos = int(value[1])
            except (TypeError, ValueError):
                pos = 0
        else:
            continue
        if chrom and pos > 0:
            out[key] = (chrom, pos)
    return out


# ---------------------------------------------------------------------------
# Genetic maps
# ---------------------------------------------------------------------------

_map_cache: dict[str, list[tuple[int, float]] | None] = {}


def _maps_root() -> Path:
    return external.panel_root() / "maps"


def reset_map_cache() -> None:
    """Drop cached genetic maps. Tests and the panel builder call this."""
    _map_cache.clear()


def load_genetic_map(chrom: Any, root: str | Path | None = None
                     ) -> list[tuple[int, float]] | None:
    """Load a PLINK-format cM map for one chromosome, or None if absent.

    Looks under ``external.panel_root()/maps/``, which is outside the
    repository tree for the same reason every other external artefact is: the
    map files are not ours to redistribute and the licence gate lives on the
    user's own installation, not in this repo.

    A PLINK map row is ``chromosome  marker  genetic_position_cM  base_pair``.
    Rows that do not parse are skipped rather than raising, because one bad
    line in a 3 million line map must not take out the whole analysis.
    """
    chrom = _norm_chrom(chrom)
    if root is None and chrom in _map_cache:
        return _map_cache[chrom]

    base = Path(root) if root is not None else _maps_root()
    candidates: list[Path] = []
    if base.is_dir():
        for pattern in (f"chr{chrom}.map", f"{chrom}.map",
                        f"*chr{chrom}.*map", f"*chr{chrom}_*map"):
            candidates.extend(sorted(base.glob(pattern)))

    points: list[tuple[int, float]] = []
    for candidate in candidates:
        try:
            text = candidate.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                cm = float(parts[2])
                bp = int(parts[3])
            except (TypeError, ValueError):
                continue
            points.append((bp, cm))
        if points:
            break

    points.sort()
    result = points or None
    if root is None:
        _map_cache[chrom] = result
    return result


def _interpolate(points: Sequence[tuple[int, float]], bp: int) -> float:
    """Cumulative cM at ``bp``, linearly interpolated, clamped at both ends."""
    if not points:
        return 0.0
    if bp <= points[0][0]:
        return points[0][1]
    if bp >= points[-1][0]:
        return points[-1][1]
    low, high = 0, len(points) - 1
    while low + 1 < high:
        mid = (low + high) // 2
        if points[mid][0] <= bp:
            low = mid
        else:
            high = mid
    left_bp, left_cm = points[low]
    right_bp, right_cm = points[high]
    if right_bp == left_bp:
        return left_cm
    frac = (bp - left_bp) / (right_bp - left_bp)
    return left_cm + frac * (right_cm - left_cm)


def centimorgans(chrom: Any, start_bp: int, end_bp: int,
                 genetic_map: Any = None) -> dict:
    """Genetic length of a physical interval.

    With a real map the answer is interpolated between the two bracketing map
    points and ``cm_estimated`` is False. Without one the per-chromosome average
    rate is used and ``cm_estimated`` is True.

    That flag is not decoration. A cM figure drives every relationship
    prediction downstream, and presenting a rate-derived estimate as a measured
    genetic length would launder a guess into a number the user trusts.
    """
    chrom = _norm_chrom(chrom)
    try:
        start_bp = int(start_bp)
        end_bp = int(end_bp)
    except (TypeError, ValueError):
        start_bp = end_bp = 0
    if end_bp < start_bp:
        start_bp, end_bp = end_bp, start_bp
    span_bp = max(0, end_bp - start_bp)

    points: Sequence[tuple[int, float]] | None = None
    if isinstance(genetic_map, dict):
        raw = genetic_map.get(chrom) or genetic_map.get(_norm_chrom(chrom))
        points = sorted(raw) if raw else None
    elif isinstance(genetic_map, (list, tuple)) and genetic_map:
        points = sorted(genetic_map)
    elif genetic_map is None:
        points = load_genetic_map(chrom)

    if points:
        cm = _interpolate(points, end_bp) - _interpolate(points, start_bp)
        return {
            "chromosome": chrom,
            "start_bp": start_bp,
            "end_bp": end_bp,
            "span_bp": span_bp,
            "cm": round(max(0.0, cm), 4),
            "cm_estimated": False,
            "source": "genetic_map",
            "rate_cm_per_mb": None,
            "note": "Interpolated from an installed genetic map.",
        }

    rate = AVERAGE_CM_PER_MB.get(chrom, GENOME_AVERAGE_CM_PER_MB)
    cm = (span_bp / 1_000_000.0) * rate
    return {
        "chromosome": chrom,
        "start_bp": start_bp,
        "end_bp": end_bp,
        "span_bp": span_bp,
        "cm": round(cm, 4),
        "cm_estimated": True,
        "source": "average_rate",
        "rate_cm_per_mb": rate,
        "note": (
            "Estimated from the whole-chromosome average recombination rate "
            "because no genetic map is installed. Recombination is very "
            "uneven along a chromosome, so this is an order of magnitude, not "
            "a measurement."
        ),
    }


# ---------------------------------------------------------------------------
# Segment detection
# ---------------------------------------------------------------------------

def _shared_markers(genotypes_a: dict, genotypes_b: dict,
                    coords: dict[str, tuple[str, int]]
                    ) -> list[tuple[str, int, str]]:
    """Markers both samples call, with coordinates, in chromosome order."""
    genotypes_a = genotypes_a or {}
    genotypes_b = genotypes_b or {}
    smaller, larger = ((genotypes_a, genotypes_b)
                       if len(genotypes_a) <= len(genotypes_b)
                       else (genotypes_b, genotypes_a))
    out: list[tuple[str, int, str]] = []
    for rsid in smaller:
        key = str(rsid).strip().lower()
        if key not in larger and rsid not in larger:
            continue
        where = coords.get(key)
        if not where:
            continue
        out.append((where[0], where[1], key))
    out.sort(key=lambda row: (_chrom_sort_key(row[0]), row[1], row[2]))
    return out


def shared_segments(genotypes_a: dict, genotypes_b: dict, *,
                    positions: dict | None = None,
                    min_snps: int = MIN_SNPS,
                    min_cm: float = MIN_CM,
                    max_opposite_homozygotes: int = MAX_OPPOSITE_HOMOZYGOTES,
                    genetic_map: Any = None,
                    chromosomes: Iterable[str] | None = None) -> dict:
    """Detect shared segments between two samples from unphased genotypes.

    The walk goes along each chromosome in coordinate order. A run opens at the
    first marker where the two samples share at least one allele, IBS >= 1, and
    stays open while that holds. An opposite homozygote, IBS 0, is impossible
    inside a genuinely shared segment, so it closes the run; up to
    ``max_opposite_homozygotes`` of them are tolerated first because array
    genotyping error is real and one bad probe should not chop a 60 cM segment
    into two 30 cM segments. When the allowance is exceeded the run is closed at
    the marker before the offending one and a new run may open after it.

    This is the classic unphased IBS approach. It is APPROXIMATE. It uses no
    population allele frequencies, so it cannot tell a genuinely shared segment
    from a long stretch where both people happen to carry the common genotype,
    which is why ``approximate`` is True in the payload and why the 300 SNP and
    7 cM floors exist.
    """
    coords = _coords(genotypes_a, positions)
    coords.update(_coords(genotypes_b, None))
    if positions:
        coords.update(_coords({}, positions))
    markers = _shared_markers(genotypes_a, genotypes_b, coords)

    wanted = None
    if chromosomes is not None:
        wanted = {_norm_chrom(c) for c in chromosomes}

    segments: list[dict] = []
    compared = 0
    ibs0_total = 0
    by_chrom: dict[str, list[tuple[int, str]]] = {}
    for chrom, pos, rsid in markers:
        if wanted is not None and chrom not in wanted:
            continue
        by_chrom.setdefault(chrom, []).append((pos, rsid))

    estimated_any = False
    for chrom in sorted(by_chrom, key=_chrom_sort_key):
        rows = by_chrom[chrom]
        states: list[int] = []
        usable: list[tuple[int, str]] = []
        for pos, rsid in rows:
            state = ibs_state(genotypes_a.get(rsid), genotypes_b.get(rsid))
            if state is None:
                continue
            usable.append((pos, rsid))
            states.append(state)
        compared += len(usable)
        ibs0_total += sum(1 for s in states if s == 0)

        run_start: int | None = None
        run_oh: list[int] = []
        for index, state in enumerate(states):
            if run_start is None:
                if state >= 1:
                    run_start, run_oh = index, []
                continue
            if state >= 1:
                continue
            if len(run_oh) < max_opposite_homozygotes:
                run_oh.append(index)
                continue
            seg = _close_run(chrom, usable, run_start, index - 1, run_oh,
                             genetic_map, min_snps, min_cm)
            if seg is not None:
                segments.append(seg)
                estimated_any = estimated_any or seg["cm_estimated"]
            run_start, run_oh = None, []
        if run_start is not None:
            seg = _close_run(chrom, usable, run_start, len(states) - 1, run_oh,
                             genetic_map, min_snps, min_cm)
            if seg is not None:
                segments.append(seg)
                estimated_any = estimated_any or seg["cm_estimated"]

    segments.sort(key=lambda s: (_chrom_sort_key(s["chromosome"]), s["start_bp"]))
    total = total_shared_cm(segments)
    longest = longest_segment_cm(segments)
    return {
        "method": "builtin_ibs",
        "approximate": True,
        "segments": segments,
        "segment_count": len(segments),
        "total_cm": total,
        "longest_cm": longest,
        "cm_estimated": estimated_any or not segments,
        "compared_snps": compared,
        "opposite_homozygotes": ibs0_total,
        "chromosomes": sorted(by_chrom, key=_chrom_sort_key),
        "thresholds": {
            "min_snps": min_snps,
            "min_cm": min_cm,
            "max_opposite_homozygotes": max_opposite_homozygotes,
        },
        "caveats": [IBD_CAVEAT],
        "note": (
            "Comparison is limited to the kits loaded on this machine. There is "
            "no database of other people here and nothing was uploaded."
        ),
    }


def _close_run(chrom: str, usable: Sequence[tuple[int, str]],
               start_index: int, end_index: int, oh_indices: Sequence[int],
               genetic_map: Any, min_snps: int, min_cm: float) -> dict | None:
    """Turn an open run into a segment, or None when it fails the thresholds."""
    oh = set(oh_indices)
    while end_index > start_index and end_index in oh:
        end_index -= 1
    while start_index < end_index and start_index in oh:
        start_index += 1
    if end_index < start_index:
        return None
    snps = end_index - start_index + 1
    start_bp = usable[start_index][0]
    end_bp = usable[end_index][0]
    length = centimorgans(chrom, start_bp, end_bp, genetic_map)
    if snps < min_snps or length["cm"] < min_cm:
        return None
    return {
        "chromosome": chrom,
        "start_bp": start_bp,
        "end_bp": end_bp,
        "start_rsid": usable[start_index][1],
        "end_rsid": usable[end_index][1],
        "snps": snps,
        "cm": length["cm"],
        "cm_estimated": length["cm_estimated"],
        "cm_source": length["source"],
        "opposite_homozygotes": len([i for i in oh if start_index <= i <= end_index]),
    }


def total_shared_cm(segments: Any) -> float:
    """Sum of segment lengths in cM. Accepts a segment list or a full payload."""
    return round(sum(float(s.get("cm") or 0.0) for s in _segments_of(segments)), 4)


def longest_segment_cm(segments: Any) -> float:
    """Longest single segment in cM, 0.0 when there are none."""
    values = [float(s.get("cm") or 0.0) for s in _segments_of(segments)]
    return round(max(values), 4) if values else 0.0


def _segments_of(value: Any) -> list[dict]:
    if isinstance(value, dict):
        return list(value.get("segments") or [])
    return [s for s in (value or []) if isinstance(s, dict)]


# ---------------------------------------------------------------------------
# Relationship prediction
#
# LICENCE NOTE. These ranges are DNAInsight's own approximate table, written to
# the shape of the published shared-cM literature. They are NOT a copy of the
# Shared cM Project dataset or of any other licensed compilation, because this
# repository bundles only CC0 and US public domain data. The data builder must
# confirm the numbers against a public-domain source before release.
#
# The bands OVERLAP on purpose. The overlap is the finding: it is the reason a
# single confident answer cannot be given.
# ---------------------------------------------------------------------------

RELATIONSHIP_BANDS: tuple[dict, ...] = (
    {"relationship": "identical twin, or the same person tested twice",
     "low": 3300, "high": 3800, "degree": 0},
    {"relationship": "parent or child", "low": 3200, "high": 3800, "degree": 1},
    {"relationship": "full sibling", "low": 2200, "high": 3400, "degree": 2},
    {"relationship": "half sibling", "low": 1160, "high": 2440, "degree": 2},
    {"relationship": "grandparent or grandchild", "low": 1150, "high": 2320, "degree": 2},
    {"relationship": "aunt, uncle, niece or nephew", "low": 1200, "high": 2280, "degree": 3},
    {"relationship": "great-grandparent or great-grandchild",
     "low": 460, "high": 1490, "degree": 3},
    {"relationship": "half aunt, half uncle, half niece or half nephew",
     "low": 500, "high": 1450, "degree": 4},
    {"relationship": "first cousin", "low": 390, "high": 1400, "degree": 4},
    {"relationship": "first cousin once removed", "low": 100, "high": 980, "degree": 5},
    {"relationship": "half first cousin", "low": 155, "high": 980, "degree": 5},
    {"relationship": "second cousin", "low": 40, "high": 600, "degree": 6},
    {"relationship": "second cousin once removed", "low": 14, "high": 355, "degree": 7},
    {"relationship": "third cousin", "low": 0, "high": 235, "degree": 8},
    {"relationship": "fourth cousin", "low": 0, "high": 140, "degree": 10},
)

# Below this total the shared cM is indistinguishable from background sharing
# between people with no genealogical relationship anyone could trace.
UNRELATED_CEILING_CM = 20.0


def predict_relationship(total_cm: Any, *, longest_cm: Any = None,
                         segment_count: Any = None) -> dict:
    """Plausible relationships for a shared total, as a RANGE and never one answer.

    Every band whose published range contains ``total_cm`` is returned. Bands
    overlap heavily by design: a half sibling, a grandparent and an aunt all
    share about a quarter of their DNA and there is no total that separates
    them. Returning one of the three would be inventing certainty.

    ``longest_cm`` and ``segment_count`` are carried through when supplied
    because they are what actually discriminates within a band, but no rule is
    applied to them here: doing that properly needs the segment size
    distribution, which is Wave 4 work.
    """
    try:
        value = float(total_cm)
    except (TypeError, ValueError):
        value = 0.0

    candidates = [
        {"relationship": band["relationship"], "low": band["low"],
         "high": band["high"], "degree": band["degree"]}
        for band in RELATIONSHIP_BANDS
        if band["low"] <= value <= band["high"]
    ]
    candidates.sort(key=lambda band: (band["degree"], band["relationship"]))

    if value > 3800:
        state = "above_range"
        summary = (
            "This total is higher than a parent and child share. The usual "
            "cause is the same person loaded twice, not an unusual family."
        )
    elif value < UNRELATED_CEILING_CM:
        # The distant-cousin bands run down to zero, so they technically match
        # a total of 5 cM. Returning them would dress noise up as a finding:
        # two unrelated people routinely share this much. The candidates are
        # cleared rather than listed.
        candidates = []
        state = "unrelated_or_distant"
        summary = (
            "This total is at or below the level unrelated people share by "
            "chance. No relationship can be claimed from it."
        )
    elif candidates:
        state = "range"
        names = [c["relationship"] for c in candidates]
        if len(names) == 1:
            summary = f"Consistent with {names[0]}."
        else:
            summary = (
                "Consistent with any of: " + ", ".join(names[:-1])
                + f", or {names[-1]}."
            )
    else:
        state = "no_band"
        summary = (
            "This total falls between the published bands, which usually means "
            "the segments themselves need looking at rather than the total."
        )

    degrees = sorted({c["degree"] for c in candidates})
    return {
        "total_cm": round(value, 4),
        "longest_cm": (round(float(longest_cm), 4)
                       if longest_cm not in (None, "") else None),
        "segment_count": segment_count,
        "state": state,
        "candidates": candidates,
        "relationships": [c["relationship"] for c in candidates],
        "degree_range": [degrees[0], degrees[-1]] if degrees else None,
        "single_answer": None,
        "summary": summary,
        "caveat": SHARED_CM_CAVEAT,
    }


# ---------------------------------------------------------------------------
# PLINK writers for the external path
# ---------------------------------------------------------------------------

def _workdir(workdir: str | Path | None) -> Path:
    if workdir is not None:
        path = Path(workdir)
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix="dnainsight_ibd_"))


def write_plink_map(markers: Sequence[tuple[str, int, str]],
                    path: str | Path) -> str:
    """Write a PLINK .map file: chromosome, marker, genetic position, base pair.

    The genetic position column is written as 0. IBIS reads the cM column when
    it is populated and falls back to a supplied map file when it is not, so
    writing a fabricated genetic position here would be worse than writing
    nothing: it would silently override the user's real map.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{chrom}\t{rsid}\t0\t{pos}" for chrom, pos, rsid in markers]
    target.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return str(target)


def write_plink_ped(samples: Sequence[dict],
                    markers: Sequence[tuple[str, int, str]],
                    path: str | Path) -> str:
    """Write a PLINK .ped file, one row per sample.

    Each sample is ``{"id": str, "genotypes": {rsid: (a1, a2)}}`` with optional
    ``family``, ``father``, ``mother``, ``sex`` and ``phenotype``. A no-call is
    written as ``0 0``, which is what PLINK means by missing, so a failed probe
    never enters an analysis as a fabricated homozygote.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for sample in samples:
        genotypes = sample.get("genotypes") or {}
        cells = [
            str(sample.get("family") or sample.get("id") or "FAM"),
            str(sample.get("id") or "SAMPLE"),
            str(sample.get("father") or 0),
            str(sample.get("mother") or 0),
            str(sample.get("sex") or 0),
            str(sample.get("phenotype") if sample.get("phenotype") is not None else -9),
        ]
        for _chrom, _pos, rsid in markers:
            pair = _pair(genotypes.get(rsid))
            cells.extend(["0", "0"] if pair is None else [pair[0], pair[1]])
        lines.append(" ".join(cells))
    target.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return str(target)


def parse_ibis_segments(text: str) -> list[dict]:
    """Parse an IBIS .seg table into the project's segment shape.

    IBIS writes whitespace separated rows carrying the pair, the chromosome,
    the physical start and end, the IBD type, the genetic start, end and
    length, and the marker count. The parser is positional but defensive:
    a row it cannot read is skipped rather than raising, because a partially
    readable result is still worth showing and a crash is not.
    """
    out: list[dict] = []
    for line in str(text or "").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            chrom = _norm_chrom(parts[2])
            start_bp = int(float(parts[3]))
            end_bp = int(float(parts[4]))
            cm = float(parts[8])
        except (TypeError, ValueError):
            continue
        snps = 0
        if len(parts) >= 10:
            try:
                snps = int(float(parts[9]))
            except (TypeError, ValueError):
                snps = 0
        out.append({
            "chromosome": chrom,
            "start_bp": start_bp,
            "end_bp": end_bp,
            "start_rsid": None,
            "end_rsid": None,
            "snps": snps,
            "cm": round(cm, 4),
            "cm_estimated": False,
            "cm_source": "ibis",
            "ibd_type": parts[5],
            "opposite_homozygotes": None,
        })
    return out


# IBIS flags have moved between releases. Collected here so the wiring pass has
# one place to correct them.
IBIS_ARGS = ("{bed}", "{bim}", "{fam}", "-f", "{outprefix}", "-ibd2")


def detect_ibd(genotypes_a: dict, genotypes_b: dict, *,
               positions: dict | None = None,
               label_a: str = "A", label_b: str = "B",
               workdir: str | Path | None = None,
               **kwargs: Any) -> dict:
    """IBD between two samples, using IBIS when it is installed.

    IBIS is phase-free, which is the whole reason it was chosen for this wave:
    it runs directly on unphased array data with no imputation step, so
    household IBD lands in Wave 3 rather than waiting for a phasing pipeline.

    THE SUBPROCESS BOUNDARY IS THE LICENCE BOUNDARY. IBIS is GPL-3.0 and is
    never imported, linked or vendored; ``external.run`` starts a separate
    process running the copy the user installed themselves. With IBIS absent
    this returns the pure-Python result carrying ``method: "builtin_ibs"`` plus
    the standard degraded payload explaining what was not attempted.
    """
    builtin = shared_segments(genotypes_a, genotypes_b, positions=positions, **kwargs)

    blocked = external.guard("ibis", "ibd_unphased")
    if blocked is not None:
        payload = dict(builtin)
        payload["external"] = blocked
        payload["external_available"] = False
        payload["note"] = (
            "IBIS was not run, so these segments come from DNAInsight's own "
            "unphased identity-by-state walk. IBIS would give a better answer "
            "on the same unphased data. " + str(blocked.get("reason", ""))
        ).strip()
        return payload

    try:
        work = _workdir(workdir)
        coords = _coords(genotypes_a, positions)
        coords.update(_coords(genotypes_b, None))
        markers = _shared_markers(genotypes_a, genotypes_b, coords)
        prefix = work / "household"
        write_plink_map(markers, work / "household.map")
        write_plink_ped(
            [{"id": label_a, "genotypes": genotypes_a},
             {"id": label_b, "genotypes": genotypes_b}],
            markers, work / "household.ped",
        )
        completed = external.run("ibis", _format_args(
            IBIS_ARGS,
            bed=str(work / "household.bed"), bim=str(work / "household.bim"),
            fam=str(work / "household.fam"), outprefix=str(prefix),
        ))
        text = ""
        seg_file = Path(str(prefix) + ".seg")
        if seg_file.exists():
            text = seg_file.read_text(encoding="utf-8", errors="replace")
        if not text.strip():
            text = getattr(completed, "stdout", "") or ""
        segments = parse_ibis_segments(text)
    except external.ExternalError as exc:
        payload = dict(builtin)
        payload["external"] = {
            "available": False, "capability": "ibd_unphased", "tool": "IBIS",
            "tool_id": "ibis", "state": "failed", "reason": str(exc),
            "not_attempted": False, "results": [],
        }
        payload["external_available"] = False
        payload["note"] = (
            "IBIS is installed but failed, so these segments come from "
            "DNAInsight's own unphased walk instead."
        )
        return payload

    if not segments:
        payload = dict(builtin)
        payload["external"] = {
            "available": False, "capability": "ibd_unphased", "tool": "IBIS",
            "tool_id": "ibis", "state": "empty",
            "reason": "IBIS ran and reported no segments above its thresholds.",
            "not_attempted": False, "results": [],
        }
        payload["external_available"] = True
        payload["note"] = (
            "IBIS ran and found no segments. The builtin walk's result is "
            "shown alongside so the two can be compared."
        )
        return payload

    return {
        "method": "ibis",
        "approximate": False,
        "segments": segments,
        "segment_count": len(segments),
        "total_cm": total_shared_cm(segments),
        "longest_cm": longest_segment_cm(segments),
        "cm_estimated": False,
        "compared_snps": builtin["compared_snps"],
        "opposite_homozygotes": builtin["opposite_homozygotes"],
        "chromosomes": sorted({s["chromosome"] for s in segments}, key=_chrom_sort_key),
        "thresholds": builtin["thresholds"],
        "external_available": True,
        "builtin": builtin,
        "caveats": [
            "Segments came from IBIS on unphased data. The builtin walk's "
            "result is kept alongside for comparison; where the two disagree, "
            "neither is silently preferred."
        ],
        "note": "Comparison is limited to the kits loaded on this machine.",
    }


def _format_args(template: Sequence[str], **values: str) -> list[str]:
    return [str(part).format(**values) for part in template]


# ---------------------------------------------------------------------------
# Parental phasing
# ---------------------------------------------------------------------------

def _role_rows(merged: dict, role: str) -> dict:
    """Every comparison row for a role, merged across sources with that role.

    ``merge.py`` keys comparison sets by role where the role is unambiguous and
    by label when the same role was uploaded twice, so the key alone is not
    enough. The row's own ``role`` field is authoritative and that is what is
    matched here.
    """
    wanted = merge.normalize_role(role)
    out: dict[str, dict] = {}
    for key, rows in ((merged or {}).get("comparison") or {}).items():
        for rsid, row in (rows or {}).items():
            row_role = merge.normalize_role(row.get("role") or key)
            if row_role != wanted:
                continue
            out.setdefault(str(rsid).strip().lower(), row)
    return out


def phase_by_parents(merged: dict) -> dict:
    """Determine which allele came from which parent, where it is determinable.

    The determinable case is simple and needs no statistics: where the child is
    heterozygous and one parent is homozygous, that parent could only have
    transmitted its one allele, so the other allele came from the other parent.
    Where both parents are heterozygous at a heterozygous child position,
    nothing is determined and the position is reported ambiguous rather than
    guessed.

    Returns ``{rsid: {"maternal": a, "paternal": b}}`` for resolvable positions
    and ``{rsid: None}`` for the rest, so a caller can tell "we worked it out"
    from "we could not" without a second lookup. Homozygous child positions are
    trivially resolved and counted separately, because they carry no phasing
    information about anything else.
    """
    genotypes = (merged or {}).get("genotypes") or {}
    mother = _role_rows(merged, "mother")
    father = _role_rows(merged, "father")

    phased: dict[str, dict | None] = {}
    counts = {
        "child_positions": 0,
        "resolvable": 0,
        "ambiguous": 0,
        "trivial_homozygous": 0,
        "informative_heterozygous": 0,
        "inconsistent": 0,
        "no_parent_data": 0,
    }
    inconsistent: list[str] = []

    for rsid in sorted(genotypes):
        child = _pair(genotypes[rsid])
        if child is None:
            continue
        counts["child_positions"] += 1
        mrow = mother.get(rsid)
        frow = father.get(rsid)
        mgt = _pair(mrow) if mrow else None
        fgt = _pair(frow) if frow else None

        if mgt is None and fgt is None:
            counts["no_parent_data"] += 1
            counts["ambiguous"] += 1
            phased[rsid] = None
            continue

        # A parent that cannot have transmitted either child allele is a
        # Mendelian inconsistency. merge.TRIO_NOTE explains why that is almost
        # always a strand or chip mismatch rather than non-paternity, so it is
        # recorded and skipped, never acted on.
        bad = False
        if mgt is not None and not (set(child) & set(mgt)):
            bad = True
        if fgt is not None and not (set(child) & set(fgt)):
            bad = True
        if bad:
            counts["inconsistent"] += 1
            counts["ambiguous"] += 1
            inconsistent.append(rsid)
            phased[rsid] = None
            continue

        if child[0] == child[1]:
            counts["trivial_homozygous"] += 1
            counts["resolvable"] += 1
            phased[rsid] = {"maternal": child[0], "paternal": child[0],
                            "basis": "child_homozygous"}
            continue

        maternal = paternal = None
        basis = ""
        if mgt is not None and mgt[0] == mgt[1] and mgt[0] in child:
            maternal = mgt[0]
            paternal = child[1] if child[0] == maternal else child[0]
            basis = "mother_homozygous"
        elif fgt is not None and fgt[0] == fgt[1] and fgt[0] in child:
            paternal = fgt[0]
            maternal = child[1] if child[0] == paternal else child[0]
            basis = "father_homozygous"
        elif mgt is not None and fgt is not None:
            # Both parents heterozygous, or one parent absent from this
            # position. Either way the transmitted allele is not determined.
            only_mother = [a for a in set(child) if a in set(mgt) and a not in set(fgt)]
            only_father = [a for a in set(child) if a in set(fgt) and a not in set(mgt)]
            if len(only_mother) == 1 and len(only_father) == 1:
                maternal, paternal = only_mother[0], only_father[0]
                basis = "exclusive_alleles"

        if maternal is None or paternal is None:
            counts["ambiguous"] += 1
            phased[rsid] = None
            continue

        counts["resolvable"] += 1
        counts["informative_heterozygous"] += 1
        phased[rsid] = {"maternal": maternal, "paternal": paternal, "basis": basis}

    resolvable = counts["resolvable"]
    considered = counts["child_positions"]
    return {
        "available": bool(mother or father),
        "phased": phased,
        "counts": counts,
        "resolvable": resolvable,
        "ambiguous": counts["ambiguous"],
        "rate": round(resolvable / considered, 4) if considered else 0.0,
        "inconsistent_rsids": sorted(inconsistent),
        "parents_present": {
            "mother": bool(mother),
            "father": bool(father),
        },
        "note": (
            "Only positions where the child is heterozygous and a parent is "
            "homozygous are determined outright. Everything else is reported "
            "ambiguous rather than guessed."
        ),
        "trio_note": merge.TRIO_NOTE,
    }


# ---------------------------------------------------------------------------
# Household analysis
# ---------------------------------------------------------------------------

# Expected total cM for a declared pair of roles. Roles are relative to the
# primary person, which is what makes mother-versus-child a grandparent pair.
# A None expectation means the roles do not pin down a relationship and no
# agreement check is possible; that is reported rather than assumed.
_ROLE_EXPECTATIONS: dict[tuple[str, str], str] = {
    ("self", "mother"): "parent or child",
    ("self", "father"): "parent or child",
    ("self", "child"): "parent or child",
    ("self", "sibling"): "full sibling",
    ("self", "mate"): "unrelated",
    ("father", "mother"): "unrelated",
    ("mate", "mother"): "unrelated",
    ("mate", "father"): "unrelated",
    ("mate", "sibling"): "unrelated",
    ("child", "mother"): "grandparent or grandchild",
    ("child", "father"): "grandparent or grandchild",
    ("father", "sibling"): "parent or child",
    ("mother", "sibling"): "parent or child",
    ("child", "sibling"): "aunt, uncle, niece or nephew",
    ("sibling", "sibling"): "full sibling",
}

_UNRELATED_EXPECTATION = {
    "relationship": "unrelated", "low": 0.0, "high": 90.0,
    "note": "Two people with no common ancestor share essentially nothing.",
}


def expected_for_roles(role_a: Any, role_b: Any) -> dict | None:
    """Expected shared cM range for a declared pair of roles, or None.

    None means the declared roles do not determine a relationship. A mate is
    not necessarily a child's other parent, for instance, so that pair gets no
    expectation and therefore never raises a false disagreement.
    """
    first = merge.normalize_role(role_a)
    second = merge.normalize_role(role_b)
    # Both orders are tried because a pair has no inherent direction and the
    # caller should not have to know which way round the table was written.
    name = (_ROLE_EXPECTATIONS.get((first, second))
            or _ROLE_EXPECTATIONS.get((second, first)))
    if name is None:
        return None
    if name == "unrelated":
        return dict(_UNRELATED_EXPECTATION)
    for band in RELATIONSHIP_BANDS:
        if band["relationship"] == name:
            return {"relationship": name, "low": float(band["low"]),
                    "high": float(band["high"]), "note": ""}
    return None


def _household_kits(merged: dict) -> list[dict]:
    """Every loaded kit, primary first, then each comparison role."""
    merged = merged or {}
    kits: list[dict] = []
    genotypes = merged.get("genotypes") or {}
    if genotypes:
        labels = merged.get("primary_labels") or []
        kits.append({
            "key": "self",
            "label": ", ".join(labels) if labels else "self",
            "role": "self",
            "genotypes": genotypes,
        })
    for key in sorted((merged.get("comparison") or {}).keys()):
        rows = merged["comparison"][key] or {}
        if not rows:
            continue
        first = next(iter(rows.values()), {})
        kits.append({
            "key": key,
            "label": first.get("label") or key,
            "role": merge.normalize_role(first.get("role") or key),
            "genotypes": rows,
        })
    return kits


def analyse_household(merged: dict, *, min_snps: int = MIN_SNPS,
                      min_cm: float = MIN_CM,
                      max_opposite_homozygotes: int = MAX_OPPOSITE_HOMOZYGOTES,
                      genetic_map: Any = None,
                      use_external: bool = True) -> dict:
    """Pairwise IBD across every loaded kit, checked against the declared roles.

    For each pair this reports the shared total, the plausible relationships and
    whether the DNA agrees with what the user said the relationship was. When it
    does not, the pair is flagged and the message says so plainly: a declared
    sibling sharing about 3,400 cM is a full sibling, one sharing about 1,700 cM
    is a half sibling, and the person deserves to be told that rather than shown
    a number and left to work it out.

    Comparison is over loaded kits only. Nothing is uploaded and there is no
    database of other people to search, which is the whole design position of
    this module.
    """
    merged = merged or {}
    kits = _household_kits(merged)
    positions = _coords(merged.get("genotypes") or {}, None)

    pairs: list[dict] = []
    disagreements: list[dict] = []
    for i in range(len(kits)):
        for j in range(i + 1, len(kits)):
            left, right = kits[i], kits[j]
            if use_external:
                ibd = detect_ibd(
                    left["genotypes"], right["genotypes"],
                    positions=positions,
                    label_a=left["key"], label_b=right["key"],
                    min_snps=min_snps, min_cm=min_cm,
                    max_opposite_homozygotes=max_opposite_homozygotes,
                    genetic_map=genetic_map,
                )
            else:
                ibd = shared_segments(
                    left["genotypes"], right["genotypes"],
                    positions=positions,
                    min_snps=min_snps, min_cm=min_cm,
                    max_opposite_homozygotes=max_opposite_homozygotes,
                    genetic_map=genetic_map,
                )
            total = ibd["total_cm"]
            prediction = predict_relationship(
                total, longest_cm=ibd["longest_cm"],
                segment_count=ibd["segment_count"],
            )
            expectation = expected_for_roles(left["role"], right["role"])
            agreement = _check_agreement(left, right, total, expectation, prediction)
            pair = {
                "a": {"key": left["key"], "label": left["label"], "role": left["role"]},
                "b": {"key": right["key"], "label": right["label"], "role": right["role"]},
                "declared": expectation["relationship"] if expectation else None,
                "expected_range": ([expectation["low"], expectation["high"]]
                                   if expectation else None),
                "total_cm": total,
                "longest_cm": ibd["longest_cm"],
                "segment_count": ibd["segment_count"],
                "cm_estimated": ibd["cm_estimated"],
                "method": ibd["method"],
                "ibd": ibd,
                "prediction": prediction,
                "agreement": agreement,
            }
            pairs.append(pair)
            if agreement["disagrees"]:
                disagreements.append({
                    "a": pair["a"], "b": pair["b"],
                    "declared": pair["declared"],
                    "total_cm": total,
                    "message": agreement["message"],
                })

    return {
        "kits": [{"key": k["key"], "label": k["label"], "role": k["role"],
                  "positions": len(k["genotypes"])} for k in kits],
        "pairs": pairs,
        "pair_count": len(pairs),
        "disagreements": disagreements,
        "disagreement_count": len(disagreements),
        "thresholds": {"min_snps": min_snps, "min_cm": min_cm,
                       "max_opposite_homozygotes": max_opposite_homozygotes},
        "caveats": [IBD_CAVEAT, SHARED_CM_CAVEAT],
        "scope": (
            "Private family genomics over the kits loaded on this machine. "
            "Nothing was uploaded, no database of other people was searched, "
            "and no third party can see or opt into these comparisons."
        ),
    }


def _check_agreement(left: dict, right: dict, total_cm: float,
                     expectation: dict | None, prediction: dict) -> dict:
    """Does the DNA match what the user declared the relationship to be."""
    if expectation is None:
        return {
            "checkable": False, "disagrees": False, "expected": None,
            "message": (
                f"The declared roles {left['role']} and {right['role']} do not "
                f"pin down an expected amount of shared DNA, so there is "
                f"nothing to check."
            ),
        }
    low, high = expectation["low"], expectation["high"]
    inside = low <= total_cm <= high
    if inside:
        return {
            "checkable": True, "disagrees": False,
            "expected": expectation["relationship"],
            "expected_range": [low, high],
            "message": (
                f"{total_cm:,.0f} cM is consistent with the declared "
                f"relationship, {expectation['relationship']}."
            ),
        }
    alternatives = prediction.get("relationships") or []
    if alternatives:
        instead = "Consistent with " + ", ".join(alternatives) + " instead."
    else:
        instead = "It is not consistent with any relationship in the table."
    return {
        "checkable": True, "disagrees": True,
        "expected": expectation["relationship"],
        "expected_range": [low, high],
        "message": (
            f"The declared relationship does not match the DNA. "
            f"{left['label']} and {right['label']} are recorded as "
            f"{expectation['relationship']}, which shares roughly {low:,.0f} to "
            f"{high:,.0f} cM, but they share {total_cm:,.0f} cM. {instead} "
            f"Before drawing any family conclusion from this, check that both "
            f"files are from the same chip generation and the same strand "
            f"orientation, which is the far more common explanation."
        ),
    }


# ---------------------------------------------------------------------------
# Chromosome browser
# ---------------------------------------------------------------------------

def chromosome_browser_data(pair: Any) -> dict:
    """Per-chromosome segment layout, ready for the interactive report's SVG.

    Fractional start and end are included alongside the base pair coordinates
    because the renderer draws its own SVG with no charting library, so it needs
    a 0 to 1 position rather than a pixel. Chromosomes with no shared segment
    are still returned, empty: an absent bar reads as "not analysed" and a
    present empty bar reads as "analysed, nothing shared", and those are
    different statements.
    """
    if isinstance(pair, dict) and "ibd" in pair:
        payload = pair["ibd"]
        labels = {"a": (pair.get("a") or {}).get("label"),
                  "b": (pair.get("b") or {}).get("label")}
    else:
        payload = pair
        labels = {"a": None, "b": None}
    segments = _segments_of(payload)
    method = (payload or {}).get("method") if isinstance(payload, dict) else None

    by_chrom: dict[str, list[dict]] = {}
    for segment in segments:
        by_chrom.setdefault(_norm_chrom(segment.get("chromosome")), []).append(segment)

    analysed = (payload or {}).get("chromosomes") if isinstance(payload, dict) else None
    order = [c for c in CHROM_LENGTHS_GRCH37 if c not in ("Y", "MT")]
    for chrom in sorted(by_chrom, key=_chrom_sort_key):
        if chrom not in order:
            order.append(chrom)

    rows: list[dict] = []
    for chrom in order:
        length = CHROM_LENGTHS_GRCH37.get(chrom, 0)
        items = []
        for segment in sorted(by_chrom.get(chrom, []),
                              key=lambda s: int(s.get("start_bp") or 0)):
            start = int(segment.get("start_bp") or 0)
            end = int(segment.get("end_bp") or 0)
            items.append({
                "start_bp": start,
                "end_bp": end,
                "cm": float(segment.get("cm") or 0.0),
                "snps": int(segment.get("snps") or 0),
                "cm_estimated": bool(segment.get("cm_estimated")),
                "start_fraction": round(start / length, 6) if length else 0.0,
                "end_fraction": round(min(end, length) / length, 6) if length else 0.0,
            })
        rows.append({
            "chromosome": chrom,
            "length_bp": length,
            "segments": items,
            "segment_count": len(items),
            "cm": round(sum(i["cm"] for i in items), 4),
            "analysed": (chrom in analysed) if analysed is not None else None,
        })

    return {
        "chromosomes": rows,
        "labels": labels,
        "method": method,
        "total_cm": total_shared_cm(segments),
        "longest_cm": longest_segment_cm(segments),
        "segment_count": len(segments),
        "max_length_bp": max(CHROM_LENGTHS_GRCH37[c] for c in order
                             if c in CHROM_LENGTHS_GRCH37) if order else 0,
        "build": "GRCh37",
        "cm_estimated": any(s.get("cm_estimated") for s in segments),
        "caveats": [IBD_CAVEAT],
    }
