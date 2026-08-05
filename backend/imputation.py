"""
imputation.py -- genotype imputation adapter, with per-variant quality as a
first-class field rather than a hidden implementation detail.

WHY THIS MODULE LOOKS THE WAY IT DOES
-------------------------------------
Imputation is the single largest number a consumer product can put on a box.
SelfDecode imputes roughly 650,000 typed positions up to a claimed 83 million
and advertises 99.7 percent block-prediction accuracy. Two things are true at
once about that number: the pipeline is real, and the number is close to
meaningless at the level a user actually cares about, which is one variant on
one report card.

The reason is that imputation accuracy is not a property of a run, it is a
property of a variant. Beagle already computes it per variant and writes it to
the INFO field as DR2, the estimated squared correlation between the imputed
dosage and the true dosage. SelfDecode's own reviewers note the two places it
falls apart: variants below 1 percent minor allele frequency, and ancestries
under-represented in the reference panel. Both are exactly the cases where DR2
drops. The information exists. It is simply not shown.

DNAInsight already refuses to guess a strand it cannot verify (invariant 2), and
already separates "checked and absent" from "never checked" (invariant 3).
An imputed call with an unknown or poor DR2 is the same class of problem, so it
gets the same treatment:

  1. DR2 is parsed, stored, displayed and filterable. Never dropped.
  2. An imputed call is capped so it can never reach magnitude parity with a
     directly typed call, however good its DR2. A prediction and a measurement
     are not the same kind of evidence and must not sort as though they were.
  3. A low-DR2 imputed call is hard-capped at IMPUTED_MAGNITUDE_CAP and the cap
     appears in the finding's audit trail as a named step, because an
     unexplained score change is precisely what this project refuses to do.
  4. An imputed pathogenic call with no quality evidence is treated as a
     violation of invariant 1 (do not alarm someone about a variant they may
     not carry) and is detectable in code, not just in review.

LICENCE BOUNDARY
----------------
Beagle is GPL-3.0-or-later. DNAInsight is MIT. This module never imports,
links or vendors it. Every invocation goes through ``external.run``, and THE
SUBPROCESS BOUNDARY IS THE LICENCE BOUNDARY: the copyleft attaches to the copy
of Beagle the user installed into their own home directory, not to this file.
See backend/external.py for the full argument.

OFFLINE CONTRACT
----------------
Nothing here touches the network, at import time or on any read path. The
reference panel is built by the user, out of tree, on explicit action.

TWO DIFFERENT FAILURES, TWO DIFFERENT MESSAGES
----------------------------------------------
"Beagle is not installed" and "the reference panel has not been built" are
distinct problems with distinct fixes. Collapsing them into one "imputation
unavailable" string sends the user to install software they may already have.
They are reported separately, and the payloads carry a machine-readable
``problem`` key so the UI never has to string-match a sentence.

VCF PARSING
-----------
The INFO parser below is deliberately minimal and local. There is no
backend/sequencing.py in the tree at the time of writing. If one lands, this
parser should be deleted and its callers pointed at it; a second VCF parser is
a second place for the same bug to live.
"""

from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

from . import external

__all__ = [
    "ImputationError", "ImputedQualityViolation",
    "IMPUTED_MAGNITUDE_CAP", "TYPED_MAGNITUDE_CEILING", "PARITY_MARGIN",
    "IMPUTED_MAGNITUDE_CEILING", "DEFAULT_DR2_THRESHOLD", "DEFAULT_PANEL",
    "QUALITY_BANDS", "BAND_CUT_POINTS", "DR2_KEYS", "RARE_MAF",
    "CAP_STEP_PREFIX", "MANDATORY_CAVEATS",
    "CAVEAT_PREDICTION_NOT_MEASUREMENT", "CAVEAT_PANEL_BIAS",
    "CAVEAT_RARE_VARIANTS", "CAVEAT_NOT_CONFIRMATORY",
    "parse_info", "dr2_from_info", "af_from_info", "parse_vcf_line",
    "parse_vcf", "read_vcf_text", "quality_band", "max_imputed_magnitude",
    "apply_imputation_cap", "apply_imputation_cap_all",
    "coverage_report", "coverage_caveat", "build_caveats",
    "filter_tokens", "filter_flag_pattern", "filter_op_pattern",
    "assert_no_imputed_pathogenic_without_quality",
    "write_study_vcf", "panel_unavailable", "impute",
]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class ImputationError(Exception):
    """Base class for failures raised by this module."""


class ImputedQualityViolation(ImputationError):
    """An imputed pathogenic call carries no usable quality evidence.

    This is invariant 1 (do not alarm a user about a variant they may not
    carry) failing at scale, because imputation multiplies the number of calls
    by roughly a hundred. Carries the offending findings so a caller can log
    exactly which ones failed rather than just that something did.
    """

    def __init__(self, violations: list[dict]) -> None:
        super().__init__(
            f"{len(violations)} imputed pathogenic call(s) carry no usable "
            f"imputation quality. Refusing to present them."
        )
        self.violations = violations


# ---------------------------------------------------------------------------
# Magnitude policy
#
# The parity rule is structural, not cosmetic. A typed call may reach the full
# scoring range. An imputed call may not, ever, at any DR2. PARITY_MARGIN is
# what enforces "never reaches parity": with a margin of 0.5 a perfectly imputed
# variant tops out at 9.5 while a typed one can reach 10.0, so sorting by
# magnitude always places the measurement above the prediction when the
# underlying evidence is equal.
# ---------------------------------------------------------------------------

# Kept in step with scoring.MAGNITUDE_MAX. Imported defensively rather than at
# module scope so that a future circular import in the wiring pass cannot break
# imputation, and so this module is testable on its own.
def _typed_ceiling() -> float:
    try:
        from .scoring import MAGNITUDE_MAX  # noqa: WPS433 (local import is deliberate)
        return float(MAGNITUDE_MAX)
    except Exception:
        return 10.0


TYPED_MAGNITUDE_CEILING: float = _typed_ceiling()

# The gap that guarantees no imputed call ever ties with a typed one.
PARITY_MARGIN: float = 0.5

# Best magnitude any imputed call may reach, however high its DR2.
IMPUTED_MAGNITUDE_CEILING: float = round(TYPED_MAGNITUDE_CEILING - PARITY_MARGIN, 2)

# Hard cap for an imputed call whose DR2 is below the working threshold, or
# whose DR2 is missing entirely. 3.0 is the same shelf the scoring module uses
# for "interesting enough to show, not strong enough to act on".
IMPUTED_MAGNITUDE_CAP: float = 3.0

# Beagle's own documentation and the imputation literature converge on 0.8 as
# the usual working floor for DR2. It is a default, not a law, which is why it
# is a parameter on every entry point.
DEFAULT_DR2_THRESHOLD: float = 0.8

DEFAULT_PANEL: str = "onekg_sgdp"

# Below this minor allele frequency, imputation accuracy degrades sharply even
# on a well-matched panel. This is the number SelfDecode's reviewers name and
# the one the rare-variant caveat is built on.
RARE_MAF: float = 0.01


# ---------------------------------------------------------------------------
# Quality bands
#
# Cut points are documented here, in one place, because a band shown in the UI
# whose boundary lives in a template is a band nobody can audit.
#
#   high      DR2 >= 0.90   dosage correlation good enough to treat the call as
#                           a working genotype, still not a measurement
#   moderate  0.80 to 0.90  usable, worth showing, worth hedging
#   low       0.30 to 0.80  below the usual working floor, capped hard
#   unusable  DR2 <  0.30   the call carries essentially no information
#   unknown   DR2 absent    NOT the same as unusable. We could not look. This
#                           mirrors frequency.py, where 0.0 and None are kept
#                           strictly apart, and genosets.py, where "not
#                           testable" is its own state.
# ---------------------------------------------------------------------------

QUALITY_BANDS: tuple[str, ...] = ("high", "moderate", "low", "unusable", "unknown")

BAND_CUT_POINTS: tuple[tuple[float, str], ...] = (
    (0.90, "high"),
    (0.80, "moderate"),
    (0.30, "low"),
    (0.00, "unusable"),
)

# Bands that mean "do not let this drive anything". Used by the safety check and
# by the dubious flag.
UNTRUSTED_BANDS: frozenset[str] = frozenset({"low", "unusable", "unknown"})

# INFO keys that may carry an imputation r-squared, most specific first. Beagle
# writes DR2. Minimac and IMPUTE write R2 and INFO respectively, and a user who
# imputed elsewhere and handed us the VCF should not silently lose their quality
# column just because a different tool named it differently.
DR2_KEYS: tuple[str, ...] = ("DR2", "R2", "INFO")

# The named audit-trail step. Tests assert on this prefix, and so should any
# UI that wants to render the cap specially.
CAP_STEP_PREFIX: str = "imputation cap"


# ---------------------------------------------------------------------------
# Mandatory caveats
#
# Same pattern as prs.MANDATORY_CAVEATS: constants, not template text, so a UI
# refactor cannot quietly drop one. The panel-bias caveat and the rare-variant
# caveat are required by policy and are asserted in the tests.
# ---------------------------------------------------------------------------

CAVEAT_PREDICTION_NOT_MEASUREMENT = (
    "An imputed genotype is a statistical prediction made from the typed "
    "markers around it. It was never measured on your array."
)
CAVEAT_PANEL_BIAS = (
    "Imputation is only as good as the reference panel behind it. Every openly "
    "licensed panel under-represents non-European ancestries, so imputed calls "
    "are less accurate for people whose ancestry is not well represented, and "
    "the DR2 figures themselves are less trustworthy for those people."
)
CAVEAT_RARE_VARIANTS = (
    "Accuracy falls sharply for rare variants. Below about 1 percent minor "
    "allele frequency an imputed call is frequently wrong even when the "
    "overall run looks excellent, and rare variants are exactly the ones a "
    "clinical report cares about."
)
CAVEAT_NOT_CONFIRMATORY = (
    "No imputed call is confirmatory. Anything that would change a medical "
    "decision must be re-tested with a clinically validated assay that "
    "directly genotypes the position."
)

MANDATORY_CAVEATS: tuple[str, ...] = (
    CAVEAT_PREDICTION_NOT_MEASUREMENT,
    CAVEAT_PANEL_BIAS,
    CAVEAT_RARE_VARIANTS,
    CAVEAT_NOT_CONFIRMATORY,
)

CAVEAT_COVERAGE_TEMPLATE = (
    "Of {total} imputed variants, {above} ({pct} percent) reached DR2 "
    "{threshold} or better. The remainder are shown with their quality band "
    "and are capped so they cannot outrank a directly typed call."
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


# ---------------------------------------------------------------------------
# VCF INFO parsing
#
# Local and minimal on purpose. Consolidate into backend/sequencing.py when
# that module exists; two VCF parsers is two places for the same bug.
# ---------------------------------------------------------------------------

def parse_info(info: Any) -> dict[str, Any]:
    """Parse a VCF INFO column into a dict.

    ``KEY=VALUE`` pairs become string values. Bare flags such as Beagle's
    ``IMP`` become ``True``, because a flag's presence IS its value and coercing
    it to an empty string would make ``if info.get("IMP")`` silently false.
    A "." or an empty column yields ``{}``.
    """
    text = _text(info)
    if not text or text == ".":
        return {}
    out: dict[str, Any] = {}
    for field in text.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, _, value = field.partition("=")
            key = key.strip()
            if key:
                out[key] = value.strip()
        else:
            out[field] = True
    return out


def _first_number(value: Any) -> float | None:
    """First comma-separated number from an INFO value.

    Multiallelic sites carry one value per ALT. Taking the first is the
    conservative reading for a biallelic array pipeline and is recorded here so
    nobody later assumes the field was scalar.
    """
    if value is None or value is True:
        return None
    head = str(value).split(",")[0].strip()
    return _as_float(head)


def dr2_from_info(info: Any) -> float | None:
    """Return the imputation r-squared from an INFO string or dict, or None.

    Tries DR2, then R2, then INFO, in that order. ``None`` means the field was
    absent, which is NOT the same as a DR2 of zero: absent means the producing
    tool told us nothing about quality, and that is treated as untrusted rather
    than as a measured failure.
    """
    fields = info if isinstance(info, dict) else parse_info(info)
    for key in DR2_KEYS:
        if key in fields:
            value = _first_number(fields[key])
            if value is not None:
                return value
    return None


def af_from_info(info: Any) -> float | None:
    """Return the alternate allele frequency from INFO, or None.

    Beagle writes AF for the imputed ALT allele. Used only to compute the minor
    allele frequency for the rare-variant warning, never to alter a call.
    """
    fields = info if isinstance(info, dict) else parse_info(info)
    for key in ("AF", "MAF", "AF1"):
        if key in fields:
            value = _first_number(fields[key])
            if value is not None:
                return value
    return None


def minor_allele_frequency(af: float | None) -> float | None:
    """Fold an alternate allele frequency onto the minor allele, 0 to 0.5."""
    if af is None:
        return None
    return round(min(af, 1.0 - af), 6)


def parse_vcf_line(line: Any) -> dict | None:
    """Parse one VCF data line into a variant record, or None.

    Returns ``None`` for headers, blanks and malformed rows rather than
    raising. A single corrupt line in a multi-million line VCF must not take the
    whole report down, and a row we cannot read is reported as a skip, not
    silently invented.
    """
    text = str(line or "").rstrip("\r\n")
    if not text or text.startswith("#"):
        return None
    parts = text.split("\t")
    if len(parts) < 8:
        # Tolerate space-delimited output from tools that ignore the spec.
        parts = text.split()
    if len(parts) < 8:
        return None

    chrom, pos, rsid, ref, alt, qual, filt, info = parts[:8]
    position = _as_float(pos)
    fields = parse_info(info)
    dr2 = dr2_from_info(fields)
    af = af_from_info(fields)

    # Beagle flags an imputed record with the bare IMP key. A record without it
    # in a Beagle output VCF was genotyped in the study file and merely carried
    # through, so it is typed, not imputed. That distinction is the whole point
    # of this module and must not be inferred from DR2 alone.
    imputed = bool(fields.get("IMP")) or bool(fields.get("IMPUTED"))

    return {
        "chromosome": _text(chrom).replace("chr", "").replace("CHR", ""),
        "position": int(position) if position is not None else None,
        "rsid": _text(rsid).lower() if _text(rsid) not in ("", ".") else "",
        "ref": _text(ref).upper(),
        "alt": _text(alt).upper(),
        "qual": _text(qual),
        "filter": _text(filt),
        "info": fields,
        "dr2": dr2,
        "af": af,
        "maf": minor_allele_frequency(af),
        "imputed": imputed,
        "typed": not imputed,
        "quality_band": quality_band(dr2) if imputed else "typed",
        "genotype_field": _text(parts[9]) if len(parts) > 9 else "",
    }


def _maybe_path(source: Any) -> Path | None:
    """Return ``source`` as an existing Path, or None if it is not one.

    A multi-megabyte VCF passed as text is still a str, and handing it to
    ``Path.exists()`` raises OSError "File name too long" rather than returning
    False. Checking for a newline first is what keeps "text or a path" from
    being a trap.
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


def read_vcf_text(path: str | Path) -> str:
    """Read a VCF as text, transparently handling gzip.

    Beagle writes ``<out>.vcf.gz``. gzip is standard library, so this costs no
    dependency.
    """
    p = Path(path)
    if p.suffix == ".gz":
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    with open(p, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def parse_vcf(source: Any) -> list[dict]:
    """Parse a VCF given as text, an iterable of lines, or a path.

    Returns one record per readable data line. Unreadable lines are skipped,
    never guessed at.
    """
    path = _maybe_path(source)
    if path is not None:
        source = read_vcf_text(path)
    if isinstance(source, str):
        lines: Iterable[Any] = source.splitlines()
    elif isinstance(source, Iterable):
        lines = source
    else:
        return []

    out: list[dict] = []
    for line in lines:
        record = parse_vcf_line(line)
        if record is not None:
            out.append(record)
    return out


# ---------------------------------------------------------------------------
# Quality bands
# ---------------------------------------------------------------------------

def quality_band(dr2: Any) -> str:
    """Classify a DR2 into one of QUALITY_BANDS.

    Cut points, inclusive at the lower bound:

        DR2 >= 0.90                 "high"
        0.80 <= DR2 < 0.90          "moderate"
        0.30 <= DR2 < 0.80          "low"
        0.00 <= DR2 < 0.30          "unusable"
        DR2 is None or unparseable  "unknown"

    ``None`` maps to "unknown" and not to "unusable" on purpose. A missing DR2
    means the producing tool reported no quality at all; that is a different
    fact from a measured near-zero correlation, and this project does not
    conflate absent with zero anywhere else either. A negative DR2 (which some
    tools emit as a rounding artefact) is clamped into "unusable" rather than
    rejected, because it is still a measurement.
    """
    value = _as_float(dr2)
    if value is None:
        return "unknown"
    if value < 0.0:
        return "unusable"
    for lower, band in BAND_CUT_POINTS:
        if value >= lower:
            return band
    return "unusable"


def max_imputed_magnitude(dr2: Any,
                          *, threshold: float = DEFAULT_DR2_THRESHOLD) -> float:
    """Highest magnitude an imputed call with this DR2 may ever reach.

    Two tiers, and both matter:

      * DR2 at or above the threshold: IMPUTED_MAGNITUDE_CEILING. Strictly
        below TYPED_MAGNITUDE_CEILING, so a prediction never ties a
        measurement, no matter how confident the prediction is.
      * DR2 below the threshold, or missing: IMPUTED_MAGNITUDE_CAP. A hard
        shelf, because below the working floor the call is a lead, not a
        finding.

    The return value is always strictly less than TYPED_MAGNITUDE_CEILING. That
    property is the parity guarantee and is asserted in the tests.
    """
    value = _as_float(dr2)
    limit = _as_float(threshold)
    if limit is None:
        limit = DEFAULT_DR2_THRESHOLD
    if value is None or value < limit:
        return IMPUTED_MAGNITUDE_CAP
    return IMPUTED_MAGNITUDE_CEILING


# ---------------------------------------------------------------------------
# The cap
# ---------------------------------------------------------------------------

def apply_imputation_cap(finding: Any,
                         *,
                         dr2: Any = None,
                         threshold: float = DEFAULT_DR2_THRESHOLD) -> Any:
    """Apply the imputation magnitude policy to one finding, in place.

    Does nothing at all to a typed call. A finding counts as imputed when it
    carries ``imputed`` True, or when a DR2 is supplied here or already present
    on the finding. Anything else is a direct measurement and this function must
    not touch it, because silently shaving a typed call would be the same class
    of dishonesty this module exists to prevent.

    For an imputed call it sets ``imputed``, ``dr2``,
    ``imputation_quality_band``, ``imputation_capped`` and
    ``magnitude_ceiling``, clamps ``magnitude`` to the ceiling from
    :func:`max_imputed_magnitude`, and appends a step to ``magnitude_factors``
    whose text begins with CAP_STEP_PREFIX.

    The step is appended even when the magnitude did not move. "The imputation
    ceiling was considered and did not bind" is information; a silent no-op
    would leave the reader unable to tell that the rule ran at all. The scoring
    module's audit trail exists for exactly this reason.

    Returns the same dict it was given.
    """
    if not isinstance(finding, dict):
        return finding

    # An explicit ``imputed: False`` is the caller stating this call was
    # measured, and it wins over everything else. A merged VCF can carry a DR2
    # on a typed row (Beagle emits one for the markers it phased), and treating
    # that as evidence of imputation would shave a real measurement.
    if finding.get("imputed") is False:
        return finding

    supplied = _as_float(dr2)
    existing = _as_float(finding.get("dr2"))
    value = supplied if supplied is not None else existing

    declared = finding.get("imputed") is True
    has_dr2_key = supplied is not None or "dr2" in finding
    if not declared and not has_dr2_key:
        return finding

    band = quality_band(value)
    limit = _as_float(threshold)
    if limit is None:
        limit = DEFAULT_DR2_THRESHOLD
    ceiling = max_imputed_magnitude(value, threshold=limit)

    before = _as_float(finding.get("magnitude"))
    if before is None:
        before = 0.0
    after = min(before, ceiling)

    finding["imputed"] = True
    finding["dr2"] = value
    finding["imputation_quality_band"] = band
    finding["imputation_dr2_threshold"] = round(limit, 4)
    finding["magnitude_ceiling"] = ceiling
    finding["magnitude"] = round(after, 2)
    finding["imputation_capped"] = after < before

    # An imputed call whose quality we cannot vouch for is dubious in exactly
    # the sense scoring.py already means it: shown, badged, and not trusted.
    if band in UNTRUSTED_BANDS:
        finding["dubious"] = True

    steps = finding.get("magnitude_factors")
    if not isinstance(steps, list):
        steps = []
        finding["magnitude_factors"] = steps

    dr2_text = "not reported by the imputation tool" if value is None else f"{value:.3f}"
    if after < before:
        steps.append(
            f"{CAP_STEP_PREFIX}: imputed call, DR2 {dr2_text} ({band}), "
            f"capped from {before:.2f} to {after:.2f} "
            f"(ceiling {ceiling:.2f}, threshold {limit:.2f})"
        )
    else:
        steps.append(
            f"{CAP_STEP_PREFIX}: imputed call, DR2 {dr2_text} ({band}), "
            f"ceiling {ceiling:.2f} did not bind, magnitude {after:.2f} unchanged"
        )
    return finding


def apply_imputation_cap_all(findings: Any,
                             *, threshold: float = DEFAULT_DR2_THRESHOLD) -> Any:
    """Apply :func:`apply_imputation_cap` across a list, in place."""
    for f in findings or []:
        apply_imputation_cap(f, threshold=threshold)
    return findings


# ---------------------------------------------------------------------------
# Safety invariant
# ---------------------------------------------------------------------------

# ClinVar significance codes the scoring module already uses. 5 is pathogenic,
# 4 likely pathogenic.
_PATHOGENIC_CODES: frozenset[int] = frozenset({4, 5})

_PATHOGENIC_TERMS: tuple[str, ...] = ("pathogenic", "risk factor", "drug response")


def _is_pathogenic(finding: dict) -> bool:
    """True when a finding would be presented to the user as clinically serious.

    Deliberately generous. A false positive here costs one extra check; a false
    negative lets an unverified imputed pathogenic call reach a person's screen.
    """
    code = finding.get("clinvar_sig_code")
    if isinstance(code, int) and code in _PATHOGENIC_CODES:
        return True
    sig = _text(finding.get("clinical_sig")).lower()
    if "benign" in sig and "pathogenic" not in sig:
        return False
    return any(term in sig for term in _PATHOGENIC_TERMS)


def assert_no_imputed_pathogenic_without_quality(findings: Any,
                                                 *,
                                                 strict: bool = True,
                                                 threshold: float = DEFAULT_DR2_THRESHOLD
                                                 ) -> list[dict]:
    """Guarantee that no imputed pathogenic call is presented without quality.

    The worst thing this whole wave can do is show somebody a pathogenic
    variant they do not actually carry. Imputation makes that risk roughly a
    hundred times larger by volume, so the invariant is enforced in code and
    tested, not left to review.

    A finding is a violation when it is imputed AND pathogenic AND any of:

      * ``dr2`` is missing, so quality was never established;
      * its quality band is ``low``, ``unusable`` or ``unknown``, so quality
        was established and is not good enough;
      * the imputation cap never ran, detectable by the absence of
        ``imputation_quality_band`` or of the named cap step in
        ``magnitude_factors``. An uncapped imputed pathogenic call can sort
        above a typed one, which is the parity failure this module exists to
        prevent.

    Returns the list of violations, each a dict with ``rsid``, ``reason``,
    ``dr2``, ``band`` and ``magnitude``. With ``strict`` left True it also
    raises :class:`ImputedQualityViolation` when the list is non-empty. Pass
    ``strict=False`` to audit without raising, which is what a diagnostics
    endpoint wants.
    """
    violations: list[dict] = []
    for f in findings or []:
        if not isinstance(f, dict):
            continue
        if f.get("imputed") is not True:
            continue
        if not _is_pathogenic(f):
            continue

        value = _as_float(f.get("dr2"))
        band = _text(f.get("imputation_quality_band"))
        steps = f.get("magnitude_factors")
        capped = isinstance(steps, list) and any(
            str(s).startswith(CAP_STEP_PREFIX) for s in steps
        )

        reasons: list[str] = []
        if value is None:
            reasons.append("no DR2 was reported, so imputation quality is unknown")
        elif value < (_as_float(threshold) or DEFAULT_DR2_THRESHOLD):
            reasons.append(
                f"DR2 {value:.3f} is below the {threshold} threshold"
            )
        if not band:
            reasons.append("no imputation_quality_band was set")
        elif band in UNTRUSTED_BANDS:
            reasons.append(f"quality band is {band}")
        if not capped:
            reasons.append(
                "the imputation cap never ran, so this call can outrank a typed one"
            )

        if not reasons:
            continue

        violations.append({
            "rsid": f.get("rsid", ""),
            "gene": f.get("gene", ""),
            "clinical_sig": f.get("clinical_sig", ""),
            "dr2": value,
            "band": band or "unknown",
            "magnitude": _as_float(f.get("magnitude")),
            "reason": "; ".join(reasons),
        })

    if violations and strict:
        raise ImputedQualityViolation(violations)
    return violations


# ---------------------------------------------------------------------------
# Coverage reporting
#
# Same shape as frequency.coverage_report and prs.coverage_report: a flat dict
# of plain counts, safe to call when nothing has been built, no exceptions.
# ---------------------------------------------------------------------------

def coverage_report(variants: Any,
                    *,
                    dr2_threshold: float = DEFAULT_DR2_THRESHOLD,
                    panel: str = DEFAULT_PANEL,
                    typed_count: int | None = None) -> dict:
    """Summarise the quality of an imputation run.

    ``variants`` is the record list from :func:`parse_vcf`. Counts are honest
    about the three-way split that matters: typed, imputed with a usable DR2,
    and imputed without one. A run that imputes 80 million variants of which
    half are unusable is not an 80 million variant run, and this report is what
    stops the headline number from being reported alone.
    """
    records = [v for v in (variants or []) if isinstance(v, dict)]

    typed = [v for v in records if not v.get("imputed")]
    imputed = [v for v in records if v.get("imputed")]

    bands = {band: 0 for band in QUALITY_BANDS}
    above = 0
    below = 0
    unknown = 0
    rare = 0
    rare_above = 0
    values: list[float] = []

    limit = _as_float(dr2_threshold)
    if limit is None:
        limit = DEFAULT_DR2_THRESHOLD

    for v in imputed:
        value = _as_float(v.get("dr2"))
        band = v.get("quality_band") or quality_band(value)
        if band not in bands:
            band = "unknown"
        bands[band] += 1

        if value is None:
            unknown += 1
        else:
            values.append(value)
            if value >= limit:
                above += 1
            else:
                below += 1

        maf = _as_float(v.get("maf"))
        if maf is not None and maf < RARE_MAF:
            rare += 1
            if value is not None and value >= limit:
                rare_above += 1

    total_imputed = len(imputed)
    mean_dr2 = round(sum(values) / len(values), 4) if values else None

    if values:
        ordered = sorted(values)
        middle = len(ordered) // 2
        median_dr2 = (ordered[middle] if len(ordered) % 2
                      else round((ordered[middle - 1] + ordered[middle]) / 2.0, 4))
    else:
        median_dr2 = None

    panel_state = external.panel_status(panel)

    return {
        "panel":             panel,
        "panel_available":   bool(panel_state.get("available")),
        "dr2_threshold":     round(limit, 4),
        "typed":             len(typed) if typed_count is None else int(typed_count),
        "imputed":           total_imputed,
        "total":             (len(records) if typed_count is None
                              else int(typed_count) + total_imputed),
        "above_threshold":   above,
        "below_threshold":   below,
        "unknown_dr2":       unknown,
        "usable_fraction":   round(above / total_imputed, 4) if total_imputed else 0.0,
        "bands":             bands,
        "mean_dr2":          mean_dr2,
        "median_dr2":        median_dr2,
        "rare_variants":     rare,
        "rare_above_threshold": rare_above,
        "rare_maf_cutoff":   RARE_MAF,
    }


def coverage_caveat(coverage: dict | None) -> str:
    """Render the coverage caveat with the real numbers filled in."""
    coverage = coverage or {}
    total = int(coverage.get("imputed") or 0)
    above = int(coverage.get("above_threshold") or 0)
    threshold = coverage.get("dr2_threshold", DEFAULT_DR2_THRESHOLD)
    pct = f"{(above / total * 100.0):.1f}" if total else "0.0"
    return CAVEAT_COVERAGE_TEMPLATE.format(
        total=total, above=above, pct=pct, threshold=threshold
    )


def build_caveats(coverage: dict | None = None,
                  *, panel: str = DEFAULT_PANEL) -> list[str]:
    """Return every caveat that must accompany an imputation result.

    The four MANDATORY_CAVEATS always come first and always in full, followed
    by the coverage sentence with the run's real numbers, followed by the
    panel's own recorded note when there is one. The panel note is where the
    ancestry-bias statement in external.PANELS lands, and it is additive: the
    panel-bias caveat above is present whether or not the panel is built,
    because the bias is a property of the field, not of this run.
    """
    out = list(MANDATORY_CAVEATS)
    out.append(coverage_caveat(coverage))

    note = _text(external.panel_status(panel).get("note"))
    if note and note not in out:
        out.append(note)
    return out


# ---------------------------------------------------------------------------
# Query grammar contribution
#
# filters.py owns the grammar. This module owns the meaning of its own tokens,
# so it publishes them here and the wiring pass registers them. Editing
# filters.py from this module would put the grammar in two files.
# ---------------------------------------------------------------------------

def _flag_imputed(finding: dict) -> bool:
    return isinstance(finding, dict) and finding.get("imputed") is True


def _flag_typed(finding: dict) -> bool:
    """Typed means directly genotyped on the array.

    Note this is not simply "not imputed": a finding that never went through
    imputation at all has no ``imputed`` key, and it is typed. Absence of the
    flag is the common case and must resolve to True.
    """
    return isinstance(finding, dict) and finding.get("imputed") is not True


def _op_dr2(finding: dict, operator: str, value: Any) -> bool:
    """Compare a finding's DR2 against a threshold.

    A finding with no DR2 never matches a /dr2 comparison, in either direction.
    Treating a missing DR2 as 0 would make ``/dr2<0.5`` quietly return every
    typed variant in the file, which is the opposite of what the user asked.
    """
    if not isinstance(finding, dict):
        return False
    actual = _as_float(finding.get("dr2"))
    target = _as_float(value)
    if actual is None or target is None:
        return False
    if operator == ">=":
        return actual >= target
    if operator == "<=":
        return actual <= target
    if operator == ">":
        return actual > target
    if operator == "<":
        return actual < target
    if operator == "=":
        return abs(actual - target) < 1e-9
    return False


# Every key a token dict is guaranteed to carry. The wiring pass may rely on
# these unconditionally.
FILTER_TOKEN_KEYS: tuple[str, ...] = (
    "token", "name", "kind", "field", "pattern", "operators",
    "description", "example", "predicate",
)


def filter_tokens() -> list[dict]:
    """Return the query-grammar tokens this module contributes.

    Three tokens: ``/imputed``, ``/typed`` and ``/dr2>=N``.

    Each dict carries ``kind``, which is ``"flag"`` or ``"operator"`` matching
    the two families filters.py already has, a ``pattern`` fragment to fold into
    that module's _FLAG_RE or _OP_RE alternation, and a ``predicate``:

      * a flag predicate takes ``(finding)`` and returns bool;
      * an operator predicate takes ``(finding, operator, value)``.

    Publishing the predicate rather than only the field name matters because
    ``/typed`` is not the negation of a truthy field lookup, and ``/dr2`` must
    not treat a missing DR2 as zero. Both rules would be lost if filters.py had
    to reinvent them from a field name.
    """
    return [
        {
            "token": "/imputed",
            "name": "imputed",
            "kind": "flag",
            "field": "imputed",
            "pattern": "imputed",
            "operators": (),
            "description": "Only calls produced by imputation, never measured on the array.",
            "example": "/imputed",
            "predicate": _flag_imputed,
        },
        {
            "token": "/typed",
            "name": "typed",
            "kind": "flag",
            "field": "imputed",
            "pattern": "typed",
            "operators": (),
            "description": "Only calls the array directly genotyped.",
            "example": "/typed",
            "predicate": _flag_typed,
        },
        {
            "token": "/dr2>=N",
            "name": "dr2",
            "kind": "operator",
            "field": "dr2",
            "pattern": "dr2",
            "operators": (">=", "<=", ">", "<", "="),
            "description": (
                "Filter by imputation quality. DR2 is the estimated squared "
                "correlation between the imputed dosage and the true dosage. "
                "Calls with no DR2 never match."
            ),
            "example": "/dr2>=0.8",
            "predicate": _op_dr2,
        },
    ]


def filter_flag_pattern() -> str:
    """Regex alternation fragment for this module's flags, for _FLAG_RE."""
    return "|".join(t["pattern"] for t in filter_tokens() if t["kind"] == "flag")


def filter_op_pattern() -> str:
    """Regex alternation fragment for this module's operators, for _OP_RE."""
    return "|".join(t["pattern"] for t in filter_tokens() if t["kind"] == "operator")


# ---------------------------------------------------------------------------
# Study VCF writing
# ---------------------------------------------------------------------------

_NOCALL = {"", "N", "-", "--", "0", "D", "I", "."}

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##source=DNAInsight\n"
    '##INFO=<ID=DNAI_REF_ASSUMED,Number=0,Type=Flag,'
    'Description="REF allele inferred from the observed array alleles, '
    'not read from a reference sequence">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
)


def _norm_allele(allele: Any) -> str:
    a = str(allele or "").strip().upper()
    return "" if a in _NOCALL else a


def _call_records(genotypes: Any) -> tuple[list[dict], list[dict]]:
    """Normalise caller input into VCF-writable records.

    Accepts a list of dicts, or a dict keyed by rsID whose values are either
    record dicts or two-character genotype strings. Returns
    ``(records, skipped)``.

    A bare ``{rsid: "AG"}`` map carries no coordinates, and a VCF row without a
    chromosome and position cannot be placed against a reference panel. Those
    entries are returned in ``skipped`` with a reason rather than being dropped
    silently or given a made-up position, because a silently dropped variant
    reads downstream as "your array did not cover it", which is a different and
    false claim.
    """
    raw: list[dict] = []
    if isinstance(genotypes, dict):
        for rsid, value in genotypes.items():
            if isinstance(value, dict):
                record = dict(value)
                record.setdefault("rsid", rsid)
                raw.append(record)
            elif isinstance(value, str) and len(value.strip()) == 2:
                text = value.strip()
                raw.append({"rsid": rsid, "allele1": text[0], "allele2": text[1]})
            elif isinstance(value, Sequence) and len(value) == 2:
                raw.append({"rsid": rsid,
                            "allele1": str(value[0]), "allele2": str(value[1])})
            else:
                raw.append({"rsid": rsid})
    elif isinstance(genotypes, Iterable):
        raw = [dict(r) for r in genotypes if isinstance(r, dict)]

    records: list[dict] = []
    skipped: list[dict] = []
    for record in raw:
        rsid = _text(record.get("rsid")).lower()
        chrom = _text(record.get("chromosome") or record.get("chrom")).upper()
        chrom = chrom.replace("CHR", "")
        position = _as_float(record.get("position") or record.get("pos"))

        genotype = record.get("genotype")
        a1 = _norm_allele(record.get("allele1"))
        a2 = _norm_allele(record.get("allele2"))
        if (not a1 or not a2) and isinstance(genotype, str) and len(genotype.strip()) == 2:
            a1, a2 = _norm_allele(genotype.strip()[0]), _norm_allele(genotype.strip()[1])

        if not chrom or position is None:
            skipped.append({"rsid": rsid, "reason": "no chromosome or position"})
            continue
        if not a1 or not a2:
            skipped.append({"rsid": rsid, "reason": "no call"})
            continue

        records.append({
            "rsid": rsid, "chromosome": chrom, "position": int(position),
            "allele1": a1, "allele2": a2,
            "ref": _norm_allele(record.get("ref") or record.get("ref_allele")),
        })

    records.sort(key=lambda r: (_chrom_sort(r["chromosome"]), r["position"]))
    return records, skipped


def _chrom_sort(chrom: str) -> tuple:
    """Sort chromosomes numerically where possible, then X, Y, MT."""
    text = str(chrom or "").upper()
    if text.isdigit():
        return (0, int(text), "")
    return (1, 0, text)


def write_study_vcf(genotypes: Any,
                    path: str | Path,
                    *,
                    sample: str = "SAMPLE") -> dict:
    """Write array calls out as a minimal single-sample VCF.

    Returns ``{"path", "written", "skipped", "ref_assumed", "warnings"}``.

    KNOWN LIMITATION, STATED RATHER THAN HIDDEN: a consumer array file does not
    say which allele is the reference. Without a reference FASTA the REF and ALT
    columns have to be inferred from the two observed alleles, and the inference
    is wrong whenever the person is homozygous for the alternate allele. Each
    such row is flagged with the DNAI_REF_ASSUMED INFO flag and counted in
    ``ref_assumed``, and a warning is returned. Supply ``ref`` on a record to
    skip the guess. The wiring pass should feed real REF alleles from the
    bundled reference table where it has them.
    """
    records, skipped = _call_records(genotypes)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)

    ref_assumed = 0
    with open(target, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(_VCF_HEADER)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 f"{sample}\n")
        for r in records:
            a1, a2 = r["allele1"], r["allele2"]
            ref = r["ref"]
            assumed = False
            if not ref:
                # Alphabetical is arbitrary but deterministic, and determinism
                # is what makes the flag auditable.
                ref = min(a1, a2)
                assumed = True
            alts = sorted({a for a in (a1, a2) if a != ref})
            alt = ",".join(alts) if alts else "."

            allele_index = {ref: 0}
            for i, a in enumerate(alts, start=1):
                allele_index[a] = i
            gt = f"{allele_index.get(a1, 0)}/{allele_index.get(a2, 0)}"

            info = "DNAI_REF_ASSUMED" if assumed else "."
            if assumed:
                ref_assumed += 1
            fh.write(
                f"{r['chromosome']}\t{r['position']}\t{r['rsid'] or '.'}\t"
                f"{ref}\t{alt}\t.\tPASS\t{info}\tGT\t{gt}\n"
            )

    warnings: list[str] = []
    if ref_assumed:
        warnings.append(
            f"{ref_assumed} of {len(records)} variants had no reference allele "
            f"available, so REF was inferred from the observed alleles. Rows "
            f"where the person is homozygous for the alternate allele will be "
            f"inverted. Those rows carry the DNAI_REF_ASSUMED flag."
        )
    if skipped:
        warnings.append(
            f"{len(skipped)} calls were not written because they carried no "
            f"coordinates or no genotype. They are absent from the study VCF "
            f"and were not imputed."
        )

    return {
        "path": str(target),
        "written": len(records),
        "skipped": skipped,
        "ref_assumed": ref_assumed,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Degraded payloads
#
# Two failures, two payloads. The ``problem`` key exists so no caller ever has
# to string-match a sentence to tell them apart.
# ---------------------------------------------------------------------------

def panel_unavailable(panel_id: str, capability: str, *, detail: str = "") -> dict:
    """The degraded payload for a reference panel that has not been built.

    Deliberately not external.unavailable(). That function describes a missing
    TOOL and its ``how_to_enable`` sends the user to a download page. Sending
    somebody to reinstall Beagle when Beagle is already working and the panel is
    what is missing wastes their evening and teaches them the error messages are
    not to be trusted.
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
            f"The reference panel '{name}' is only partially built. "
            f"Missing: {', '.join(missing)}. This is a DATA problem, not a "
            f"tool problem: nothing needs to be installed. Finish or rerun the "
            f"panel build."
        )
    else:
        reason = (
            f"The reference panel '{name}' has not been built. Imputation "
            f"without a reference panel is not possible, and no amount of "
            f"reinstalling the imputation tool will change that. This is a "
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
        "variants":      [],
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
                "Panel data is large and is deliberately kept outside the "
                "repository tree, in the same place as external tools.",
            ],
            "note": st.get("note", ""),
        },
        "caveats": list(MANDATORY_CAVEATS),
    }


def _tool_unavailable(tool_id: str, capability: str) -> dict:
    """external.guard's payload, tagged so callers can branch on the cause."""
    payload = external.unavailable(tool_id, capability)
    payload["problem"] = "tool_missing"
    payload["variants"] = []
    payload["caveats"] = list(MANDATORY_CAVEATS)
    return payload


# ---------------------------------------------------------------------------
# The entry point
# ---------------------------------------------------------------------------

def impute(genotypes_or_vcf_path: Any,
           *,
           panel: str = DEFAULT_PANEL,
           dr2_threshold: float = DEFAULT_DR2_THRESHOLD,
           workdir: str | Path | None = None,
           sample: str = "SAMPLE",
           chromosome: str | None = None,
           nthreads: int | None = None,
           timeout: int = external.DEFAULT_TIMEOUT,
           keep_workdir: bool = False) -> dict:
    """Impute a genotype set against a reference panel using Beagle.

    ``genotypes_or_vcf_path`` is either a path to an existing VCF, or a
    genotype collection in any of the shapes :func:`write_study_vcf` accepts.

    Degrades rather than raising, in three distinguishable ways, each carrying
    a ``problem`` key:

      * ``tool_missing``  Beagle is absent, its runtime is absent, or its
        licence has not been accepted. external.guard decides which.
      * ``panel_missing``  Beagle is ready but the reference panel is not
        built. Different fix, different message.
      * ``run_failed``  Beagle ran and failed. The stderr tail is carried
        through, because "imputation failed" with no detail is unactionable.

    On success returns ``available`` True plus ``variants`` (each with dr2,
    quality_band, af and maf), ``coverage`` from :func:`coverage_report`, and
    ``caveats`` from :func:`build_caveats`. ``results`` aliases ``variants`` so
    that every payload this module produces, degraded or not, answers to the
    same key.

    The subprocess boundary in external.run IS the licence boundary. Beagle is
    GPL-3.0-or-later and is never imported here.
    """
    blocked = external.guard("beagle", "imputation")
    if blocked is not None:
        return _tool_unavailable("beagle", "imputation")

    panel_state = external.panel_status(panel)
    if not panel_state.get("available"):
        return panel_unavailable(panel, "imputation")

    limit = _as_float(dr2_threshold)
    if limit is None:
        limit = DEFAULT_DR2_THRESHOLD

    owned_workdir = workdir is None
    work = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="dnainsight_impute_"))
    work.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    write_result: dict = {}

    try:
        source = genotypes_or_vcf_path
        supplied_vcf = _maybe_path(source)
        if supplied_vcf is not None:
            study_vcf = supplied_vcf
        else:
            study_vcf = work / "study.vcf"
            write_result = write_study_vcf(source, study_vcf, sample=sample)
            warnings.extend(write_result.get("warnings", []))
            if not write_result.get("written"):
                return {
                    "available":     False,
                    "capability":    "imputation",
                    "problem":       "no_input",
                    "state":         "no_usable_calls",
                    "reason": (
                        "No usable calls were found in the supplied genotypes. "
                        "Imputation needs a chromosome and a position for every "
                        "variant, and a plain rsID to genotype map does not "
                        "carry them."
                    ),
                    "not_attempted": True,
                    "results": [], "variants": [],
                    "skipped": write_result.get("skipped", []),
                    "caveats": list(MANDATORY_CAVEATS),
                }

        base = Path(panel_state["path"])
        out_prefix = work / "imputed"
        args = [
            f"gt={study_vcf}",
            f"ref={base / 'panel.vcf.gz'}",
            f"map={base / 'panel.map'}",
            f"out={out_prefix}",
        ]
        if chromosome:
            args.append(f"chrom={chromosome}")
        if nthreads:
            args.append(f"nthreads={int(nthreads)}")

        try:
            external.run("beagle", args, timeout=timeout, cwd=work)
        except external.ExternalError as exc:
            return {
                "available":     False,
                "capability":    "imputation",
                "problem":       "run_failed",
                "state":         "run_failed",
                "tool":          "Beagle 5.5",
                "panel":         panel,
                "reason":        str(exc),
                "not_attempted": False,
                "results": [], "variants": [],
                "caveats": list(MANDATORY_CAVEATS),
            }

        # Beagle writes <out>.vcf.gz. Accept the uncompressed form too, because
        # a user who post-processed the output should not hit a false negative.
        produced = None
        for candidate in (out_prefix.with_suffix(".vcf.gz"),
                          Path(str(out_prefix) + ".vcf.gz"),
                          Path(str(out_prefix) + ".vcf")):
            if candidate.exists():
                produced = candidate
                break
        if produced is None:
            return {
                "available":     False,
                "capability":    "imputation",
                "problem":       "run_failed",
                "state":         "no_output",
                "reason": (
                    "Beagle reported success but produced no output VCF at "
                    f"{out_prefix}.vcf.gz. Treating that as a failure rather "
                    f"than as an empty result."
                ),
                "not_attempted": False,
                "results": [], "variants": [],
                "caveats": list(MANDATORY_CAVEATS),
            }

        variants = parse_vcf(produced)
        coverage = coverage_report(variants, dr2_threshold=limit, panel=panel)

        return {
            "available":      True,
            "capability":     "imputation",
            "problem":        None,
            "state":          "ready",
            "tool":           "Beagle 5.5",
            "tool_licence":   "GNU General Public License v3.0 or later",
            "panel":          panel,
            "panel_name":     panel_state.get("name", panel),
            "dr2_threshold":  round(limit, 4),
            "output_vcf":     str(produced) if keep_workdir else "",
            "study_vcf":      write_result.get("path", str(study_vcf))
                              if keep_workdir else "",
            "variants":       variants,
            "results":        variants,
            "coverage":       coverage,
            "caveats":        build_caveats(coverage, panel=panel),
            "warnings":       warnings,
            "not_attempted":  False,
        }
    finally:
        if owned_workdir and not keep_workdir:
            shutil.rmtree(work, ignore_errors=True)
