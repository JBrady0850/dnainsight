"""
prs.py -- additive polygenic score layer for DNAInsight.

What this does
--------------
Computes simple additive polygenic scores entirely locally, from weight files
bundled in data/prs_models.json. No network, no upload, no third-party
scoring service. A score is the sum over its variants of (copies of the
effect allele) times (the published per-allele effect size), and the effect
sizes are natural-log odds ratios.

What this deliberately does not do
----------------------------------
It does not diagnose, and it does not pretend the numbers are better than
they are. Three limitations are structural rather than fixable here:

1. A consumer array genotypes a few hundred thousand positions. Published
   scores often use hundreds of thousands to millions of variants. Coverage
   is reported on every result and a score below 90 percent coverage is
   flagged ``reliable=False``.
2. Missing variants are not silently dropped, because dropping them biases
   the score downward. Where the model supplies an effect allele frequency
   the population mean dosage (2f) is imputed instead, and the result is
   flagged ``mean_imputed=True``. Where it does not, the variant is skipped
   and coverage falls accordingly.
3. The reference mean and standard deviation behind every percentile came
   from a European-ancestry panel. A percentile computed against the wrong
   panel can be badly wrong. This is stated verbatim in the caveats on every
   single result, not buried in documentation.

Percentiles come from the standard normal CDF evaluated with ``math.erf``.
There is no numpy or scipy dependency; the project ships flask and requests
only.

A missing data/prs_models.json is not an error. ``load_models`` returns {}
and ``compute_all`` returns [], so the app still runs on a fresh checkout.
"""

import json
import math
from pathlib import Path
from typing import Any, Sequence


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
_BASE = Path(__file__).parent.parent
MODELS_FILE = _BASE / "data" / "prs_models.json"

# Allele tokens an array emits when it could not make a call. Mirrors
# scanner._NOCALL so the modules agree on what "no data" looks like.
_NOCALL = {"", "N", "-", "--", "0", "D", "I"}

# Below this fraction of directly genotyped variants a score is not trusted.
RELIABLE_COVERAGE = 0.90

# At most this many missing rsIDs are echoed back; a 1M-variant model would
# otherwise return a megabyte of rsIDs no one is going to read.
MAX_MISSING_REPORTED = 50

BANDS: tuple[str, ...] = (
    "low", "below_average", "average", "above_average", "high", "unknown",
)

# A PRS band is a position in a distribution, not a severity. The magnitudes
# here only drive UI sort order: extremes first, unremarkable last.
BAND_MAGNITUDE: dict[str, float] = {
    "high":          3.0,
    "low":           3.0,
    "above_average": 2.0,
    "below_average": 2.0,
    "average":       1.0,
    "unknown":       0.0,
}

BAND_LABELS: dict[str, str] = {
    "low":           "in the lowest decile",
    "below_average": "below average",
    "average":       "average",
    "above_average": "above average",
    "high":          "in the highest decile",
    "unknown":       "not determined",
}


# ---------------------------------------------------------------------------
# Mandatory caveats
#
# These three strings are attached verbatim to every result, plus a fourth
# carrying the real coverage number. They are constants rather than template
# text so that a UI change cannot quietly drop them.
# ---------------------------------------------------------------------------
CAVEAT_NOT_DIAGNOSTIC = (
    "A polygenic score is a statistical predictor, not a diagnostic test."
)
CAVEAT_BOTH_DIRECTIONS = (
    "People with high scores may never develop the trait, and people with "
    "low scores sometimes do."
)
CAVEAT_ANCESTRY = (
    "Scores were developed largely in European-ancestry cohorts and transfer "
    "poorly to other ancestries. Your percentile may be materially wrong if "
    "your ancestry differs from the reference panel."
)
CAVEAT_COVERAGE_TEMPLATE = (
    "Consumer arrays genotype only a fraction of the variants in most "
    "published scores. Coverage for this score was {pct} percent."
)

MANDATORY_CAVEATS: tuple[str, ...] = (
    CAVEAT_NOT_DIAGNOSTIC,
    CAVEAT_BOTH_DIRECTIONS,
    CAVEAT_ANCESTRY,
)


# ---------------------------------------------------------------------------
# Model loading (module-level cache)
# ---------------------------------------------------------------------------
_models_cache: dict = {}
_models_meta: dict = {}
_models_loaded: bool = False
_models_path: str = ""


def load_models(path: str | Path | None = None) -> dict:
    """Load and cache data/prs_models.json, returning ``{model_id: model}``.

    Accepts the versioned ``{"_meta": ..., "models": ...}`` shape and a flat
    ``{model_id: model}`` shape. Returns ``{}`` when the file is missing or
    unreadable. Passing a different ``path`` forces a reload.
    """
    global _models_cache, _models_meta, _models_loaded, _models_path

    target = str(Path(path)) if path else str(MODELS_FILE)
    if _models_loaded and target == _models_path:
        return _models_cache

    _models_cache, _models_meta = {}, {}
    _models_path, _models_loaded = target, True

    p = Path(target)
    if not p.exists():
        return _models_cache
    try:
        with open(p, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except (OSError, ValueError):
        return _models_cache

    if isinstance(raw, dict) and "models" in raw:
        _models_meta = raw.get("_meta", {}) or {}
        block = raw.get("models")
    elif isinstance(raw, dict):
        block = raw
    else:
        block = None

    if isinstance(block, dict):
        _models_cache = {str(k): v for k, v in block.items()
                         if isinstance(v, dict)}
    elif isinstance(block, list):
        # Tolerate a list of models by keying on each model's own id.
        _models_cache = {str(m.get("id")): m for m in block
                         if isinstance(m, dict) and m.get("id")}
    return _models_cache


def get_metadata() -> dict:
    """Return the ``_meta`` block of prs_models.json, or ``{}`` when absent."""
    load_models()
    return _models_meta


def model_ids() -> list[str]:
    """Return the ids of all bundled models, in file order."""
    return list(load_models().keys())


def get_model(model_id: str) -> dict | None:
    """Return one bundled model dict including its weights, or None."""
    return load_models().get(str(model_id))


def list_models() -> list[dict]:
    """Return one summary dict per bundled model, with no weights attached.

    Each entry carries ``id``, ``trait``, ``efo``, ``variant_count``,
    ``source``, ``citation``, ``license``, ``build`` and ``description``. Use
    this for menus and index pages; use :func:`get_model` when the weights
    themselves are needed.
    """
    out: list[dict] = []
    for model_id, model in load_models().items():
        variants = model.get("variants") or []
        out.append({
            "id":            model.get("id", model_id),
            "trait":         model.get("trait", ""),
            "efo":           model.get("efo", ""),
            "variant_count": int(model.get("variant_count") or len(variants)),
            "source":        model.get("source", ""),
            "citation":      model.get("citation", ""),
            "license":       model.get("license", ""),
            "build":         model.get("build", ""),
            "description":   model.get("description", ""),
        })
    return out


# ---------------------------------------------------------------------------
# Genotype handling
# ---------------------------------------------------------------------------

def _norm(allele: Any) -> str:
    """Normalise one allele to an uppercase base, or '' for a no-call."""
    a = str(allele or "").strip().upper()
    return "" if a in _NOCALL else a


def _split_genotype(value: Any) -> tuple[str, str] | None:
    """Coerce a genotype value into an ``(allele1, allele2)`` pair.

    Accepts a 2-element sequence or a 2-character string, both being shapes
    the parsers already produce. Returns None for anything else.
    """
    if isinstance(value, str):
        text = value.strip()
        return (text[0], text[1]) if len(text) == 2 else None
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        items = list(value)
        if len(items) == 2:
            return (str(items[0]), str(items[1]))
    return None


def _normalise_genotypes(genotypes: Any) -> dict[str, tuple[str, str]]:
    """Index a caller's genotype map by lowercase rsID with validated pairs."""
    out: dict[str, tuple[str, str]] = {}
    if not isinstance(genotypes, dict):
        return out
    for rsid, value in genotypes.items():
        pair = _split_genotype(value)
        if pair is not None:
            out[str(rsid).strip().lower()] = pair
    return out


def dosage(a1: str, a2: str, effect_allele: str) -> int | None:
    """Return the number of copies of the effect allele: 0, 1 or 2.

    Returns ``None`` for a no-call on either side, or when the effect allele
    itself is missing. Comparison is strand-naive: the caller must supply
    alleles already oriented to the same strand the model was built on, which
    is what the ORIENTED genotype map passed to :func:`compute_model` means.
    """
    ea = _norm(effect_allele)
    x, y = _norm(a1), _norm(a2)
    if not ea or not x or not y:
        return None
    return int(x == ea) + int(y == ea)


def _as_float(value: Any) -> float | None:
    """Coerce a stored value to float, preserving 0.0 and rejecting junk."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def normal_cdf(z: float) -> float:
    """Standard normal cumulative distribution function, via ``math.erf``.

    Phi(z) = 0.5 * (1 + erf(z / sqrt(2))). Kept explicit so that the absence
    of scipy is a deliberate choice rather than an accident.
    """
    return 0.5 * (1.0 + math.erf(float(z) / math.sqrt(2.0)))


def percentile_from_z(z: float | None) -> float | None:
    """Convert a z score to a percentile from 0 to 100, or None."""
    if z is None:
        return None
    return round(normal_cdf(z) * 100.0, 1)


def band_for_percentile(percentile: float | None) -> str:
    """Bucket a percentile into one of the BANDS.

    Under the 10th percentile is ``low``, 10 to 30 ``below_average``, 30 to 70
    ``average``, 70 to 90 ``above_average`` and above 90 ``high``. ``None``
    gives ``unknown``, which is what happens when a model ships no reference
    distribution.
    """
    if percentile is None:
        return "unknown"
    p = _as_float(percentile)
    if p is None:
        return "unknown"
    if p < 10.0:
        return "low"
    if p < 30.0:
        return "below_average"
    if p <= 70.0:
        return "average"
    if p <= 90.0:
        return "above_average"
    return "high"


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def coverage_caveat(coverage: float) -> str:
    """Return the coverage caveat with the real percentage filled in."""
    return CAVEAT_COVERAGE_TEMPLATE.format(pct=f"{coverage * 100.0:.1f}")


def build_caveats(coverage: float) -> list[str]:
    """Return the four caveats that must accompany every PRS result."""
    return list(MANDATORY_CAVEATS) + [coverage_caveat(coverage)]


def compute_model(model_id: str, genotypes: dict,
                  model: dict | None = None) -> dict:
    """Compute one additive polygenic score from a genotype map.

    ``genotypes`` maps lowercase rsID to either a 2-tuple or a 2-character
    string of ORIENTED alleles, meaning alleles already on the same strand the
    model's effect alleles are stated on.

    ``model`` may be supplied directly to score against a model that is not
    bundled, which is what the PGS Catalog import path in
    data/build_prs.py uses. Otherwise the id is looked up in the bundled file
    and a KeyError is raised when it is not there.

    Scoring: raw_score is the sum over matched variants of dosage times
    weight. A variant with no call is not silently skipped, because dropping
    variants biases the total downward and makes two people's scores
    incomparable. Where the model gives an effect allele frequency f, the
    population mean dosage 2f is imputed and ``mean_imputed`` is set. Where it
    does not, the variant is dropped and coverage falls.

    ``coverage`` counts only directly genotyped variants, never imputed ones,
    so an imputed score cannot flatter itself into looking well covered.
    """
    if model is None:
        model = get_model(model_id)
        if model is None:
            known = ", ".join(model_ids()) or "(no models bundled)"
            raise KeyError(
                f"unknown PRS model id {model_id!r}. Available ids: {known}")

    variants = model.get("variants") or []
    total = len(variants)
    lookup = _normalise_genotypes(genotypes)

    raw = 0.0
    used = 0
    missing: list[str] = []
    mean_imputed = False

    for variant in variants:
        if not isinstance(variant, dict):
            missing.append("?")
            continue
        rsid = str(variant.get("rsid") or "").strip().lower()
        weight = _as_float(variant.get("weight"))
        if not rsid or weight is None:
            missing.append(rsid or "?")
            continue

        pair = lookup.get(rsid)
        copies = dosage(pair[0], pair[1], variant.get("effect_allele", "")) \
            if pair else None

        if copies is not None:
            raw += copies * weight
            used += 1
            continue

        missing.append(rsid)
        eaf = _as_float(variant.get("effect_allele_frequency"))
        if eaf is not None:
            raw += 2.0 * eaf * weight
            mean_imputed = True

    coverage = (used / total) if total else 0.0

    reference = model.get("reference") or {}
    ref_mean = _as_float(reference.get("mean"))
    ref_sd = _as_float(reference.get("sd"))
    if ref_mean is not None and ref_sd is not None and ref_sd > 0.0:
        z = (raw - ref_mean) / ref_sd
    else:
        z = None
    percentile = percentile_from_z(z)

    return {
        "id":                    model.get("id", model_id),
        "trait":                 model.get("trait", ""),
        "raw_score":             round(raw, 6),
        "variants_used":         used,
        "variants_total":        total,
        "coverage":              round(coverage, 4),
        "missing_rsids":         missing[:MAX_MISSING_REPORTED],
        "mean_imputed":          mean_imputed,
        "percentile":            percentile,
        "band":                  band_for_percentile(percentile),
        "z":                     round(z, 4) if z is not None else None,
        "population_reference":  reference.get("population", "unknown"),
        "license":               model.get("license", ""),
        "citation":              model.get("citation", ""),
        "caveats":               build_caveats(coverage),
        "reliable":              coverage >= RELIABLE_COVERAGE,
    }


def compute_all(genotypes: dict) -> list[dict]:
    """Compute every bundled model against one genotype map.

    Returns a list of result dicts in bundled order. A model that fails to
    score is skipped rather than taking the whole report down with it.
    """
    results: list[dict] = []
    for model_id, model in load_models().items():
        try:
            results.append(compute_model(model_id, genotypes, model=model))
        except Exception:
            continue
    return results


# ---------------------------------------------------------------------------
# Conversion to findings
# ---------------------------------------------------------------------------

def _interpretation(result: dict) -> str:
    """Assemble the user-facing sentence for one PRS result."""
    trait = result.get("trait") or result.get("id") or "this trait"
    band = result.get("band", "unknown")
    percentile = result.get("percentile")
    coverage = _as_float(result.get("coverage")) or 0.0

    parts = [f"Polygenic score for {trait}: {BAND_LABELS.get(band, band)}."]
    if percentile is None:
        parts.append("No reference distribution was available for this model, "
                     "so no percentile can be given.")
    else:
        parts.append(
            f"The raw score falls at roughly percentile {percentile:g} of the "
            f"{result.get('population_reference', 'reference')} reference "
            f"distribution.")
    parts.append(
        f"Computed from {result.get('variants_used', 0)} of "
        f"{result.get('variants_total', 0)} model variants "
        f"({coverage * 100.0:.0f} percent coverage).")
    if result.get("mean_imputed"):
        parts.append("Variants your file did not cover were filled in at the "
                     "population average, which pulls the score toward the "
                     "middle of the distribution.")
    if not result.get("reliable", False):
        parts.append("Coverage is below 90 percent, so treat this as "
                     "indicative only.")
    return " ".join(parts)


def to_findings(results: list[dict]) -> list[dict]:
    """Convert PRS results into finding dicts the rest of the app understands.

    Each finding gets ``entity_type`` ``"prs"``, ``rsid`` set to the model id
    (a PRS has no single position), ``silo`` ``"informational"`` and
    ``magnitude`` derived from the band. ``repute`` is always the empty string:
    a polygenic score is a position in a distribution, not a good or bad
    variant, and forcing it into the Repute vocabulary would misrepresent it.

    The full result dict is carried through under the ``prs`` key so a detail
    view can show the caveats and the missing-variant list.
    """
    findings: list[dict] = []
    for result in results or []:
        if not isinstance(result, dict):
            continue
        band = result.get("band", "unknown")
        findings.append({
            "rsid":           result.get("id", ""),
            "entity_type":    "prs",
            "gene":           "",
            "genotype":       "",
            "zygosity":       "",
            "clinical_sig":   "",
            "conditions":     result.get("trait", ""),
            "interpretation": _interpretation(result),
            "category":       "PRS",
            "silo":           "informational",
            "magnitude":      BAND_MAGNITUDE.get(band, 0.0),
            "repute":         "",
            "sources":        ["prs"],
            "prs":            result,
        })
    return findings


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def coverage_report(genotypes: dict) -> dict:
    """Summarise how well one genotype map covers the bundled models."""
    results = compute_all(genotypes)
    return {
        "models":    len(results),
        "reliable":  sum(1 for r in results if r.get("reliable")),
        "mean_coverage": (
            round(sum(_as_float(r.get("coverage")) or 0.0
                      for r in results) / len(results), 4)
            if results else 0.0
        ),
        "version":   get_metadata().get("version", ""),
        "built_at":  get_metadata().get("built_at", ""),
    }
