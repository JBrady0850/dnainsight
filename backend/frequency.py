"""
frequency.py -- population genotype-frequency layer for DNAInsight.

Purpose
-------
A variant can be globally rare and still be perfectly ordinary inside the
user's own ancestry. Reporting only a global number therefore produces
needless alarm. This module answers the narrower and more useful question:
what percentage of people in a chosen reference population carry the user's
EXACT genotype?

Units
-----
Anything named "genotype frequency" is a PERCENTAGE from 0 to 100.
Anything named "allele frequency" (including ``gmaf``) is a FRACTION from
0 to 1. GMAF is the single allele-level number exposed here and it is the
global minor allele frequency.

Zero is not the same as unknown
-------------------------------
A stored frequency of 0.0 means "this genotype was not observed in that
sample". It is a real measurement. ``None`` means "we have no data for this
rsID, or none for this population". The two are never conflated: every
lookup in this module tests for key presence rather than truthiness, and
callers must do the same.

This distinction matters because the sample sizes behind the bundled data
are small. The HapMap CEU panel is roughly 120 individuals, so it cannot
resolve any frequency below about 1.7 percent (1 chromosome in 240). A 0.0
coming out of a panel that size is a statement about the sample, not a
claim that the genotype does not exist in the population.

Derived values
--------------
When observed genotype counts are unavailable, the genotype frequency is
derived from allele frequencies under Hardy-Weinberg equilibrium
(p squared, 2pq, q squared). Derived results are flagged ``derived=True``
with ``method="hardy_weinberg"`` so the UI can hedge its wording. Observed
genotype tables are reported as ``method="observed"``.

The builder may pre-compute the genotype table to keep read-time work down.
When it does so from allele frequencies rather than real counts it sets
``"genotypes_derived": true`` on the entry, and this module still reports
those values as ``hardy_weinberg``. A pre-computed table is never silently
promoted to "observed".

Data file: data/frequencies.json, produced by data/build_frequencies.py.
A missing file is not an error. Every accessor degrades to None or empty so
the app still runs on a fresh checkout.
"""

import json
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
_BASE = Path(__file__).parent.parent
FREQUENCY_FILE = _BASE / "data" / "frequencies.json"

DEFAULT_POPULATION = "CEU"

# MAX / AVG / MIN aggregate across the populations present in the data.
# GLOBAL ignores the population table and derives from GMAF instead.
AGGREGATE_MODES: tuple[str, ...] = ("MAX", "AVG", "MIN", "GLOBAL")

# Allele tokens an array emits when it could not make a call. Mirrors
# scanner._NOCALL so the two modules agree on what "no data" looks like.
_NOCALL = {"", "N", "-", "--", "0", "D", "I"}


# ---------------------------------------------------------------------------
# Reference populations (1000 Genomes / HapMap panel codes)
# ---------------------------------------------------------------------------
POPULATIONS: list[dict] = [
    {"code": "CEU", "label": "Utah residents of Northern and Western European ancestry",
     "brief": "Northern European (Utah)",           "superpop": "EUR"},
    {"code": "TSI", "label": "Toscani in Italy",
     "brief": "Italian (Tuscany)",                  "superpop": "EUR"},
    {"code": "FIN", "label": "Finnish in Finland",
     "brief": "Finnish",                            "superpop": "EUR"},
    {"code": "GBR", "label": "British in England and Scotland",
     "brief": "British",                            "superpop": "EUR"},
    {"code": "IBS", "label": "Iberian populations in Spain",
     "brief": "Iberian (Spain)",                    "superpop": "EUR"},
    {"code": "YRI", "label": "Yoruba in Ibadan, Nigeria",
     "brief": "Yoruba (Nigeria)",                   "superpop": "AFR"},
    {"code": "LWK", "label": "Luhya in Webuye, Kenya",
     "brief": "Luhya (Kenya)",                      "superpop": "AFR"},
    {"code": "ASW", "label": "African ancestry in Southwest USA",
     "brief": "African American (Southwest USA)",   "superpop": "AFR"},
    {"code": "CHB", "label": "Han Chinese in Beijing",
     "brief": "Han Chinese (Beijing)",              "superpop": "EAS"},
    {"code": "JPT", "label": "Japanese in Tokyo",
     "brief": "Japanese",                           "superpop": "EAS"},
    {"code": "CHS", "label": "Southern Han Chinese",
     "brief": "Southern Han Chinese",               "superpop": "EAS"},
    {"code": "GIH", "label": "Gujarati Indians in Houston",
     "brief": "Gujarati Indian",                    "superpop": "SAS"},
    {"code": "PJL", "label": "Punjabi in Lahore",
     "brief": "Punjabi (Lahore)",                   "superpop": "SAS"},
    {"code": "MXL", "label": "Mexican ancestry in Los Angeles",
     "brief": "Mexican ancestry (Los Angeles)",     "superpop": "AMR"},
    {"code": "PUR", "label": "Puerto Ricans in Puerto Rico",
     "brief": "Puerto Rican",                       "superpop": "AMR"},
    {"code": "CLM", "label": "Colombians in Medellin",
     "brief": "Colombian (Medellin)",               "superpop": "AMR"},
]

POPULATION_CODES: tuple[str, ...] = tuple(p["code"] for p in POPULATIONS)

_POP_BY_CODE: dict[str, dict] = {p["code"]: p for p in POPULATIONS}


# ---------------------------------------------------------------------------
# Rarity colour ramp
#
# Ordered from most common to most rare. The first stop whose lower bound the
# percentage meets or exceeds wins. Near white for everyday genotypes, then
# progressively more saturated red as the genotype gets rarer, so a long scan
# page reads at a glance without the user having to parse numbers.
#
#   >= 50 %    #FFFFFF   white         majority genotype in that population
#   >= 20 %    #FFF3EE   blush         common
#   >= 10 %    #FFE0D4   pale salmon   uncommon, upper half
#   >=  5 %    #FFC5B2   salmon        uncommon, lower half
#   >=  1 %    #FA9A80   coral         rare
#   >= 0.1 %   #EF6A4C   orange red    very rare
#   >=  0 %    #D8351B   deep red      vanishingly rare or not observed
#
# None (no data at all) gets a light neutral grey instead, so "unknown" never
# looks like "alarming".
# ---------------------------------------------------------------------------
RARITY_RAMP: tuple[tuple[float, str], ...] = (
    (50.0, "#FFFFFF"),
    (20.0, "#FFF3EE"),
    (10.0, "#FFE0D4"),
    (5.0,  "#FFC5B2"),
    (1.0,  "#FA9A80"),
    (0.1,  "#EF6A4C"),
    (0.0,  "#D8351B"),
)
UNKNOWN_COLOR = "#EFEFEF"


# ---------------------------------------------------------------------------
# Module-level cache
# ---------------------------------------------------------------------------
_freq_cache: dict = {}
_freq_meta: dict = {}
_freq_loaded: bool = False
_freq_path: str = ""
_avail_cache: list[dict] | None = None
_freq_override: str = ""


def reset_source() -> None:
    """Forget any explicit path passed to load_frequencies and drop the cache.

    Restores the module to reading the bundled data/frequencies.json. Mainly
    for tests, which otherwise leak a fixture path into every later lookup.
    """
    global _freq_cache, _freq_meta, _freq_loaded, _freq_path
    global _avail_cache, _freq_override
    _freq_cache, _freq_meta = {}, {}
    _freq_loaded, _freq_path = False, ""
    _avail_cache, _freq_override = None, ""


def load_frequencies(path: str | Path | None = None) -> dict:
    """Load and cache data/frequencies.json, returning the frequency table.

    Accepts both the versioned ``{"_meta": ..., "frequencies": ...}`` shape
    and a flat ``{rsid: entry}`` shape. Returns ``{}`` when the file is
    missing or unreadable, so a fresh checkout with no built data still runs.
    Repeat calls are served from the module cache; passing a different
    ``path`` forces a reload.

    An explicit ``path`` is REMEMBERED for the rest of the process. Every
    accessor in this module calls ``load_frequencies()`` with no argument, so
    without that the override would be discarded on the very next lookup and a
    caller who pointed the module at a test fixture would silently be reading
    the bundled file again. Call :func:`reset_source` to go back to the default.
    """
    global _freq_cache, _freq_meta, _freq_loaded, _freq_path, _avail_cache
    global _freq_override

    if path:
        _freq_override = str(Path(path))
    target = _freq_override or str(FREQUENCY_FILE)
    if _freq_loaded and target == _freq_path:
        return _freq_cache

    _freq_cache, _freq_meta, _avail_cache = {}, {}, None
    _freq_path, _freq_loaded = target, True

    p = Path(target)
    if not p.exists():
        return _freq_cache
    try:
        with open(p, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except (OSError, ValueError):
        return _freq_cache

    if isinstance(raw, dict) and "frequencies" in raw:
        _freq_meta = raw.get("_meta", {}) or {}
        _freq_cache = raw.get("frequencies", {}) or {}
    elif isinstance(raw, dict):
        _freq_cache = raw
    return _freq_cache


def get_metadata() -> dict:
    """Return the ``_meta`` block of frequencies.json, or ``{}`` when absent."""
    load_frequencies()
    return _freq_meta


def available_populations() -> list[dict]:
    """Return the population dicts actually represented in the loaded data.

    Sorted in POPULATIONS order. A population counts as available when at
    least one rsID carries an allele table, a genotype table or a sample
    count for it.
    """
    global _avail_cache
    data = load_frequencies()
    if _avail_cache is not None:
        return _avail_cache

    present: set[str] = set()
    for entry in data.values():
        if not isinstance(entry, dict):
            continue
        for block in ("alleles", "genotypes", "n"):
            table = entry.get(block)
            if isinstance(table, dict):
                present.update(table.keys())

    _avail_cache = [p for p in POPULATIONS if p["code"] in present]
    return _avail_cache


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _norm_allele(allele: Any) -> str:
    """Normalise one reported allele to an uppercase base, or '' for a no-call."""
    a = str(allele or "").strip().upper()
    return "" if a in _NOCALL else a


def _entry(rsid: str) -> dict | None:
    """Return the raw frequency entry for an rsID, or None when absent."""
    data = load_frequencies()
    e = data.get(str(rsid or "").strip().lower())
    return e if isinstance(e, dict) else None


def _table(rsid: str, block: str, population: str) -> dict | None:
    """Return one population's sub-table from an entry block, or None."""
    e = _entry(rsid)
    if e is None:
        return None
    outer = e.get(block)
    if not isinstance(outer, dict):
        return None
    inner = outer.get(population)
    return inner if isinstance(inner, dict) else None


def _as_float(value: Any) -> float | None:
    """Coerce a stored value to float, preserving 0.0 and rejecting junk."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Strand reconciliation
#
# This layer exists because of a real and easily missed data mismatch.
#
# The frequency table is built from Ensembl / 1000 Genomes, which reports
# alleles in dbSNP orientation. Consumer arrays (23andMe, AncestryDNA and the
# rest) report the PLUS strand of GRCh37 regardless of how dbSNP defined the
# variant. For any SNP dbSNP defined on the minus strand the two disagree by a
# complement.
#
# rs1801133 (MTHFR C677T) is the canonical example: a 23andMe file calls it
# C/T, Ensembl stores it as G/A. Querying the table with the raw array alleles
# returns nothing at all, which silently reads as "no frequency data" when the
# data is in fact present. Every lookup below therefore resolves the queried
# alleles against the alleles the table actually stores, complementing when
# that is the only interpretation that matches.
#
# Palindromic sites (an A/T or a C/G heterozygote) cannot be resolved this way,
# because complementing simply yields the other observed allele. Those are
# reported with ambiguous True and the unflipped reading is kept, so the caller
# can warn instead of quietly guessing.
# ---------------------------------------------------------------------------

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}
_REAL_BASES = frozenset("ACGT")


def _complement(allele: str) -> str:
    """Complement a single base, returning '' for anything that is not ACGT."""
    return _COMPLEMENT.get(_norm_allele(allele), "")


def observed_alleles(rsid: str, population: str | None = None) -> set:
    """Return the ACGT alleles the frequency table actually stores for an rsID.

    Unions across every population when ``population`` is None, which is what
    strand resolution wants: a variant sits on one strand for the whole table,
    so evidence from any panel settles the question.
    """
    e = _entry(rsid)
    if e is None:
        return set()
    tables = e.get("alleles")
    if not isinstance(tables, dict):
        return set()
    found = set()
    if population and isinstance(tables.get(population), dict):
        found.update(str(k).strip().upper() for k in tables[population])
    else:
        for t in tables.values():
            if isinstance(t, dict):
                found.update(str(k).strip().upper() for k in t)
    return {a for a in found if a in _REAL_BASES}


def is_palindromic(a1: str, a2: str) -> bool:
    """True for an A/T or C/G heterozygote, which no metadata can disambiguate."""
    x, y = _norm_allele(a1), _norm_allele(a2)
    if not x or not y or x == y:
        return False
    return _COMPLEMENT.get(x) == y


def resolve_strand(rsid: str, a1: str, a2: str,
                   population: str | None = None) -> dict:
    """Reconcile array-reported alleles against the stored table's strand.

    Returns a dict with keys allele1, allele2, flipped, ambiguous, resolved
    and observed. Falls through unchanged with resolved False when there is no
    data for the rsID, so callers always get a well-formed dict.
    """
    x, y = _norm_allele(a1), _norm_allele(a2)
    out = {
        "allele1": x, "allele2": y,
        "flipped": False, "ambiguous": is_palindromic(x, y),
        "resolved": False, "observed": [],
    }
    if not x or not y:
        return out

    observed = observed_alleles(rsid, population)
    out["observed"] = sorted(observed)
    if not observed:
        return out

    cx, cy = _complement(x), _complement(y)
    direct = x in observed and y in observed
    flipped_ok = bool(cx) and bool(cy) and cx in observed and cy in observed

    if direct:
        # Palindromic pairs match both ways. Keep the unflipped reading and
        # let the caller decide how loudly to warn.
        out["resolved"] = True
        return out

    if flipped_ok:
        out.update({"allele1": cx, "allele2": cy, "flipped": True, "resolved": True})
        return out

    # Partial match: one allele lands, the other does not. Prefer whichever
    # reading places more alleles inside the observed set.
    direct_hits = (1 if x in observed else 0) + (1 if y in observed else 0)
    flip_hits = (1 if cx and cx in observed else 0) + (1 if cy and cy in observed else 0)
    if flip_hits > direct_hits:
        out.update({"allele1": cx, "allele2": cy, "flipped": True})
    return out


# ---------------------------------------------------------------------------
# Allele-level lookups
# ---------------------------------------------------------------------------

def allele_frequency(rsid: str, allele: str,
                     population: str = DEFAULT_POPULATION,
                     strand_tolerant: bool = True) -> float | None:
    """Return the frequency of one allele in one population as a fraction.

    Range 0 to 1. A stored 0.0 is returned as 0.0 (the allele was not seen in
    that sample); ``None`` means the rsID, the population or the allele is not
    in the data at all.

    With ``strand_tolerant`` left on, an allele that is absent from the table
    is retried as its complement, because the table is in dbSNP orientation
    while consumer arrays report the GRCh37 plus strand. Pass False to query
    the table literally.
    """
    a = _norm_allele(allele)
    if not a:
        return None
    table = _table(rsid, "alleles", population)
    if table is None:
        return None
    if a in table:
        return _as_float(table[a])
    # Strand fallback: the table is in dbSNP orientation, the caller is most
    # likely holding array plus-strand alleles. Try the complement before
    # concluding the allele is absent.
    if strand_tolerant:
        c = _complement(a)
        if c and c in table:
            return _as_float(table[c])
    return None


def gmaf(rsid: str) -> float | None:
    """Return the global minor allele frequency as a fraction 0 to 1, or None."""
    e = _entry(rsid)
    if e is None:
        return None
    return _as_float(e.get("gmaf"))


def minor_allele(rsid: str) -> str:
    """Return the global minor allele for an rsID, or '' when unknown."""
    e = _entry(rsid)
    if e is None:
        return ""
    return _norm_allele(e.get("minor_allele"))


def sample_size(rsid: str, population: str = DEFAULT_POPULATION) -> int | None:
    """Return the number of individuals sampled for a population, or None.

    Useful to the UI because it bounds the resolution of any percentage shown
    next to it: a panel of n individuals cannot resolve below 1 / (2n).
    """
    e = _entry(rsid)
    if e is None:
        return None
    counts = e.get("n")
    if not isinstance(counts, dict) or population not in counts:
        return None
    v = _as_float(counts[population])
    return int(v) if v is not None else None


# ---------------------------------------------------------------------------
# Genotype-level lookups
# ---------------------------------------------------------------------------

def genotype_frequency_detail(rsid: str, a1: str, a2: str,
                              population: str = DEFAULT_POPULATION) -> dict:
    """Return the genotype frequency for one population plus its provenance.

    Shape::

        {"frequency": float | None,   # percentage 0 to 100
         "population": str,           # the population code that was queried
         "derived": bool,             # True when Hardy-Weinberg was used
         "method": "observed" | "hardy_weinberg" | "unavailable",
         "n": int | None}             # individuals in that panel, if known

    A stored genotype table is used when present, otherwise the value is
    derived from the allele table under Hardy-Weinberg. A stored table that
    the builder itself derived (entry flag ``genotypes_derived``) is still
    reported as ``hardy_weinberg``, never as ``observed``.

    A frequency of 0.0 with ``method="observed"`` means the genotype was
    genuinely absent from that sample. ``frequency=None`` with
    ``method="unavailable"`` means unknown. These are different facts.
    """
    out: dict = {
        "frequency":  None,
        "population": population,
        "derived":    False,
        "method":     "unavailable",
        "n":          sample_size(rsid, population),
        "flipped":    False,
        "ambiguous":  False,
        "queried":    "",
    }

    x, y = _norm_allele(a1), _norm_allele(a2)
    if not x or not y:
        return out

    entry = _entry(rsid)
    if entry is None:
        return out

    # Reconcile the caller's array-strand alleles with the table's dbSNP
    # orientation before any lookup. Without this, every minus-strand variant
    # reports "unavailable" while its data sits right there in the table.
    strand = resolve_strand(rsid, x, y, population)
    x, y = strand["allele1"], strand["allele2"]
    out["flipped"] = strand["flipped"]
    out["ambiguous"] = strand["ambiguous"]
    out["queried"] = f"{x}{y}"

    # 1. Stored genotype table, if the builder wrote one for this population.
    stored_derived = bool(entry.get("genotypes_derived"))
    stored = _table(rsid, "genotypes", population)
    if stored is not None:
        for key in (x + y, y + x):
            if key in stored:
                value = _as_float(stored[key])
                if value is not None:
                    out["frequency"] = round(value, 2)
                    if stored_derived:
                        out["derived"] = True
                        out["method"] = "hardy_weinberg"
                    else:
                        out["method"] = "observed"
                    return out

    # 2. Fall back to deriving from the allele frequencies.
    p = allele_frequency(rsid, x, population)
    q = allele_frequency(rsid, y, population)
    if p is None or q is None:
        return out

    freq = (p * p) if x == y else (2.0 * p * q)
    out["frequency"] = round(freq * 100.0, 2)
    out["derived"] = True
    out["method"] = "hardy_weinberg"
    return out


def genotype_frequency(rsid: str, a1: str, a2: str,
                       population: str = DEFAULT_POPULATION) -> float | None:
    """Return the percentage of a population carrying this exact genotype.

    Range 0 to 100, or ``None`` when unknown. Call
    :func:`genotype_frequency_detail` when you need to know whether the value
    was observed or derived under Hardy-Weinberg.
    """
    return genotype_frequency_detail(rsid, a1, a2, population)["frequency"]


def population_series(rsid: str, a1: str, a2: str,
                      population: str = DEFAULT_POPULATION) -> list[dict]:
    """Return this genotype's frequency across every available population.

    One entry per available population, in POPULATIONS order::

        {"code": str, "label": str, "brief": str,
         "frequency": float | None, "yours": bool}

    Exactly one entry is flagged ``yours`` whenever at least one population is
    available: the requested population when the data covers it, otherwise the
    first available population as a fallback so the UI always has an anchor.
    """
    avail = available_populations()
    if not avail:
        return []

    codes = [p["code"] for p in avail]
    yours = population if population in codes else codes[0]

    series: list[dict] = []
    for pop in avail:
        series.append({
            "code":      pop["code"],
            "label":     pop["label"],
            "brief":     pop["brief"],
            "frequency": genotype_frequency(rsid, a1, a2, pop["code"]),
            "yours":     pop["code"] == yours,
        })
    return series


def _global_genotype_frequency(rsid: str, a1: str, a2: str) -> float | None:
    """Derive a genotype percentage from GMAF alone under Hardy-Weinberg.

    GMAF gives the minor allele fraction q, so the major allele is 1 - q. This
    is the only estimate available when no per-population table exists, and it
    is an average across all of humanity rather than any one group.
    """
    g = gmaf(rsid)
    minor = minor_allele(rsid)
    if g is None or not minor:
        return None

    x, y = _norm_allele(a1), _norm_allele(a2)
    if not x or not y:
        return None

    # Reconcile strand before comparing against the stored minor allele.
    # minor_allele comes from the same dbSNP-oriented record as the allele
    # tables, so array plus-strand input must be complemented here too.
    # Without this, GLOBAL mode returned None for every minus-strand variant
    # while MAX, AVG and MIN returned real numbers, which is an inconsistency
    # that looks like missing data rather than a bug.
    strand = resolve_strand(rsid, x, y)
    if strand["flipped"]:
        x, y = strand["allele1"], strand["allele2"]
    elif not strand["resolved"]:
        # No allele table to settle it. Fall back to matching the minor allele
        # directly, trying the complement before giving up.
        if minor not in (x, y):
            cx, cy = _complement(x), _complement(y)
            if minor in (cx, cy):
                x, y = cx, cy

    # Biallelic assumption: anything that is not the minor allele is major.
    if x != y and x != minor and y != minor:
        return None

    fx = g if x == minor else (1.0 - g)
    fy = g if y == minor else (1.0 - g)
    freq = (fx * fx) if x == y else (2.0 * fx * fy)
    return round(freq * 100.0, 2)


def aggregate_frequency(rsid: str, a1: str, a2: str,
                        mode: str = "MAX") -> float | None:
    """Aggregate this genotype's frequency across populations.

    ``mode`` is one of AGGREGATE_MODES. MAX, AVG and MIN reduce over every
    available population that has a value for this genotype; populations with
    no value are excluded rather than treated as zero, so a sparse row cannot
    drag the minimum down to a fake zero. GLOBAL ignores the population table
    and derives from GMAF under Hardy-Weinberg. Returns ``None`` for an
    unknown mode or when nothing is available.
    """
    key = str(mode or "").strip().upper()
    if key == "GLOBAL":
        return _global_genotype_frequency(rsid, a1, a2)
    if key not in ("MAX", "AVG", "MIN"):
        return None

    values = [
        f for f in (
            genotype_frequency(rsid, a1, a2, p["code"])
            for p in available_populations()
        ) if f is not None
    ]
    if not values:
        return None
    if key == "MAX":
        return round(max(values), 2)
    if key == "MIN":
        return round(min(values), 2)
    return round(sum(values) / len(values), 2)


# ---------------------------------------------------------------------------
# Presentation helpers
# ---------------------------------------------------------------------------

def rarity_band(freq: float | None) -> str:
    """Classify a genotype percentage into a coarse rarity band.

    Bands: ``very_rare`` under 1 percent, ``rare`` from 1 up to 5,
    ``uncommon`` from 5 up to 20, ``common`` from 20 through 50,
    ``majority`` above 50, and ``unknown`` for ``None``.
    """
    if freq is None:
        return "unknown"
    f = _as_float(freq)
    if f is None:
        return "unknown"
    if f < 1.0:
        return "very_rare"
    if f < 5.0:
        return "rare"
    if f < 20.0:
        return "uncommon"
    if f <= 50.0:
        return "common"
    return "majority"


def rarity_color(freq: float | None) -> str:
    """Return the 7-character hex colour for a genotype percentage.

    Walks RARITY_RAMP from common to rare and returns the first stop the value
    meets. ``None`` returns the light neutral UNKNOWN_COLOR, because absent
    data must never be rendered as though it were an alarming finding.
    """
    if freq is None:
        return UNKNOWN_COLOR
    f = _as_float(freq)
    if f is None:
        return UNKNOWN_COLOR
    for lower, color in RARITY_RAMP:
        if f >= lower:
            return color
    return RARITY_RAMP[-1][1]


# ---------------------------------------------------------------------------
# Finding annotation
# ---------------------------------------------------------------------------

# Every key annotate() guarantees. Callers may index these unconditionally.
ANNOTATION_KEYS: tuple[str, ...] = (
    "freq",
    "freq_population",
    "freq_derived",
    "freq_method",
    "freq_band",
    "freq_color",
    "gmaf",
    "minor_allele",
    "population_series",
    "freq_flipped",
    "freq_ambiguous",
    "freq_queried",
)


def annotate(finding: dict, population: str = DEFAULT_POPULATION) -> dict:
    """Attach population-frequency fields to a finding dict, in place.

    Reads ``finding["rsid"]``, ``finding["allele1"]`` and
    ``finding["allele2"]``, all through ``.get`` so a bare or partial dict is
    safe. Sets every key in ANNOTATION_KEYS, using ``None``, ``"unknown"`` or
    ``[]`` when there is no data, so downstream template and report code can
    never raise KeyError. Returns the same dict it was given.
    """
    if not isinstance(finding, dict):
        return finding

    rsid = finding.get("rsid") or ""
    a1 = finding.get("allele1") or ""
    a2 = finding.get("allele2") or ""

    detail = genotype_frequency_detail(rsid, a1, a2, population)
    freq = detail["frequency"]

    finding["freq"] = freq
    finding["freq_population"] = detail["population"]
    finding["freq_derived"] = bool(detail["derived"])
    finding["freq_method"] = detail["method"]
    finding["freq_band"] = rarity_band(freq)
    finding["freq_color"] = rarity_color(freq)
    finding["gmaf"] = gmaf(rsid)
    finding["minor_allele"] = minor_allele(rsid)
    finding["population_series"] = population_series(rsid, a1, a2, population)

    # Strand provenance. The UI needs these to be honest about two things:
    # that a frequency was read after complementing the reported alleles, and
    # that a palindromic site (A/T or C/G heterozygote) cannot be strand
    # verified at all, so its frequency may belong to the other reading.
    finding["freq_flipped"] = bool(detail.get("flipped"))
    finding["freq_ambiguous"] = bool(detail.get("ambiguous"))
    finding["freq_queried"] = detail.get("queried", "")
    return finding


def annotate_all(findings: list[dict],
                 population: str = DEFAULT_POPULATION) -> list[dict]:
    """Annotate every finding in a list in place and return the same list."""
    for f in findings or []:
        annotate(f, population)
    return findings


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def coverage_report() -> dict:
    """Summarise what the bundled frequency data actually covers.

    Returns ``{"rsids", "populations", "with_gmaf", "source", "built_at"}``.
    Safe to call with no data file present: counts come back as zero and the
    strings as empty.
    """
    data = load_frequencies()
    meta = get_metadata()

    with_gmaf = sum(
        1 for e in data.values()
        if isinstance(e, dict) and _as_float(e.get("gmaf")) is not None
    )
    return {
        "rsids":       len(data),
        "populations": len(available_populations()),
        "with_gmaf":   with_gmaf,
        "source":      meta.get("source", ""),
        "built_at":    meta.get("built_at", ""),
    }
