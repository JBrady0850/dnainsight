"""
scoring.py -- DNAInsight magnitude, repute and confidence.

WHY THIS MODULE EXISTS
----------------------
Promethease reports are legible mainly because of two SNPedia fields:

    Magnitude   a hand-typed 0 to 10 "how interesting is this" number
    Repute      a hand-typed Good / Bad / blank direction flag

Both are curated by wiki editors and both are licensed CC-BY-NC-SA-3.0-US.
They cannot be redistributed inside this repository, and there is no open-data
field that means the same thing. So DNAInsight computes its own from CC0 and
public-domain evidence: ClinVar review status, CPIC assignment level, FDA label
tier, replicated GWAS support, publication counts and population frequency.

The output deliberately uses the same 0 to 10 shape, because anyone who has read
a Promethease report can read this one. It is NOT the same number and must never
be presented as though it were. Call it DNAInsight magnitude in the UI.

WHAT MAKES THIS BETTER THAN A NAIVE REPORT
------------------------------------------
Three things, in order of how much harm they prevent:

1. Carrier awareness. A ClinVar classification describes an ALLELE, not a
   position. A report that shows "pathogenic" against someone who carries two
   reference copies is simply wrong, and it is the single most common way these
   tools frighten people for no reason. Non-carriers are multiplied down to a
   quarter and their repute is cleared.
2. Strand honesty. A palindromic site (an A/T or C/G heterozygote) cannot be
   strand-verified, so its genotype might be the complement of what was read.
   Those are capped and flagged dubious rather than allowed to top the list.
3. No-calls score zero. A failed probe is not a finding.

Every score carries a ``magnitude_factors`` list recording each step, so a user
or a reviewer can see exactly why a number came out the way it did. An opaque
interest score would be worse than none at all.
"""

from __future__ import annotations

import math
from typing import Any

__all__ = [
    "REVIEW_STATUS_STARS", "CLINVAR_SIG_CODES", "CPIC_LEVELS",
    "FDA_LABEL_TIERS", "BASE_SCORES", "MAGNITUDE_MIN", "MAGNITUDE_MAX",
    "UNSCORED_SORT_VALUE", "RISK_TERMS", "PROTECTIVE_TERMS", "NEUTRAL_TERMS",
    "review_stars", "clinvar_sig_code", "normalize_cpic_level",
    "evidence_label", "base_magnitude", "compute_magnitude", "compute_repute",
    "compute_confidence", "score_finding", "score_all",
]


# ---------------------------------------------------------------------------
# ClinVar review status to star rating
#
# These are the strings that appear in the DATA, which are not identical to the
# strings in the published documentation. The docs say "no classification for
# the individual variant"; variant_summary.txt.gz emits "no classification for
# the single variant". Two more appear only in the data and in no docs table:
# "no classifications from unflagged records" and a bare "-".
#
# Somatic aggregate records use "criteria provided, multiple submitters" with no
# ", no conflicts" suffix, because consensus is not assessed for somatic
# assertions. Submitted (SCV) records add "flagged submission" and can never be
# 2 stars. All three vocabularies are folded in here deliberately: conflating
# them is a known bug source, but a lookup that accepts all of them and returns
# the right star count for each is safe.
# ---------------------------------------------------------------------------
REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no assertion provided": 0,
    "no classification for the single variant": 0,
    "no classification for the individual variant": 0,
    "no classifications from unflagged records": 0,
    "flagged submission": 0,
    "-": 0,
    "": 0,
}

# ClinVar CLNSIG numeric codes, the same set the reference product filtered on.
CLINVAR_SIG_CODES: dict[str, int] = {
    "pathogenic": 5,
    "likely pathogenic": 4,
    "pathogenic/likely pathogenic": 5,
    "pathogenic, low penetrance": 5,
    "likely pathogenic, low penetrance": 4,
    "likely benign": 3,
    "benign": 2,
    "benign/likely benign": 2,
    "uncertain significance": 1,
    "uncertain risk allele": 1,
    "drug response": 6,
    "histocompatibility": 7,
    "risk factor": 255,
    "protective": 255,
    "association": 255,
    "affects": 255,
    "confers sensitivity": 255,
    "other": 255,
    "not provided": 255,
}

# CPIC levels. Eight values, not four: there are three split levels and a
# non-letter "Retired". A four-value enum silently drops 76 real pairs.
CPIC_LEVELS: tuple[str, ...] = ("A", "A/B", "B", "B/C", "C", "C/D", "D", "Retired")

# FDA drug-label pharmacogenomic tiers, title-cased as published. Six tiers,
# not four: "No Clinical PGx" and "Criteria Not Met" were added in 2024.
FDA_LABEL_TIERS: tuple[str, ...] = (
    "Testing Required", "Testing Recommended", "Actionable PGx",
    "Informative PGx", "No Clinical PGx", "Criteria Not Met",
)

# Base scores by strongest available evidence. See docs/API_V2.md section 6.
BASE_SCORES: dict[str, float] = {
    "cpic_a":            6.0,
    "clinvar_path_3star": 6.0,
    "cpic_b":            4.5,
    "clinvar_path_2star": 4.5,
    "fda_testing":       4.0,
    "clinvar_lp_2star":  3.5,
    "gwas_replicated":   2.5,
    "clinvar_single":    1.5,
    "default":           1.0,
}

MAGNITUDE_MIN = 0.0
MAGNITUDE_MAX = 10.0

# A blank magnitude sorts as 1, matching the documented convention that an
# unscored variant sits just above "common and boring" and below "interesting".
UNSCORED_SORT_VALUE = 1.0

# Direction-of-effect vocabularies. Ordered longest-first inside each tuple so a
# specific phrase wins over a substring of itself.
RISK_TERMS: tuple[str, ...] = (
    "increased toxicity", "severe toxicity", "life-threatening",
    "loss of function", "no function", "null function", "no activity",
    "decreased function", "reduced function", "poor function",
    "poor metabolizer", "ultrarapid metabolizer", "slow metabolizer",
    "contraindicated", "do not use", "avoid",
    "increased risk", "higher risk", "elevated risk", "susceptibility",
    "treatment failure", "reduced efficacy", "decreased response",
    "no response", "adverse reaction", "malignant hyperthermia",
    "disease causing", "pathogenic", "risk factor", "toxicity",
    "deficiency", "overload", "thrombophilia", "impaired",
)
PROTECTIVE_TERMS: tuple[str, ...] = (
    "normal metabolizer", "extensive metabolizer", "normal function",
    "standard response", "standard dose", "favourable response",
    "favorable response", "not at increased", "no increased",
    "decreased risk", "reduced risk", "lower risk", "protective",
    "wildtype", "wild type", "reference genotype", "sustained virologic",
)
NEUTRAL_TERMS: tuple[str, ...] = (
    "eye colour", "eye color", "hair colour", "hair color", "earwax",
    "chronotype", "ancestry", "haplogroup", "taste", "trait",
    "uncertain significance", "conflicting", "not provided",
    "no clinical", "unknown significance",
)


# ---------------------------------------------------------------------------
# Normalisers
# ---------------------------------------------------------------------------

def _text(value: Any) -> str:
    """Lowercase, stripped string for vocabulary matching."""
    return str(value or "").strip().lower()


def _num(value: Any, default: float | None = None) -> float | None:
    """Coerce to float, rejecting bools and junk, preserving 0.0."""
    if value is None or isinstance(value, bool):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def review_stars(review_status: Any) -> int:
    """Map a ClinVar review-status string to its 0 to 4 star rating.

    Unknown strings return 0 rather than raising, because ClinVar adds new
    status strings between releases and an unrecognised one should degrade to
    "no stated confidence", never crash a scan.
    """
    return REVIEW_STATUS_STARS.get(_text(review_status), 0)


def clinvar_sig_code(clinical_sig: Any) -> int | None:
    """Map a ClinVar germline classification to its numeric CLNSIG code."""
    key = _text(clinical_sig)
    if not key:
        return None
    if key in CLINVAR_SIG_CODES:
        return CLINVAR_SIG_CODES[key]

    # A CONFLICTING record is its own tier and must never be resolved to one of
    # the components it conflicts over. This check has to run BEFORE the
    # compound fallback below, because "conflicting classifications of
    # pathogenicity" contains the substring "pathogenic" inside the word
    # "pathogenicity". Without this guard that string scored 5, which swept a
    # genuinely disputed record into the default pathogenic-only whitelist and
    # then coloured it Bad, producing exactly the false alarm this module exists
    # to prevent. ClinVar renamed the string in January 2024, so both the
    # current "classifications" spelling and the legacy "interpretations" one
    # are caught by the same substring.
    if "conflicting" in key:
        return 255

    # Compound strings such as "pathogenic/likely pathogenic, risk factor"
    # take the strongest component present.
    for name in ("pathogenic/likely pathogenic", "likely pathogenic",
                 "pathogenic", "likely benign", "benign",
                 "uncertain significance", "drug response",
                 "histocompatibility", "risk factor", "protective"):
        if name in key:
            return CLINVAR_SIG_CODES[name]
    return 255


def normalize_cpic_level(level: Any) -> str:
    """Return a canonical CPIC level string, or '' when unrecognised."""
    raw = str(level or "").strip()
    if not raw:
        return ""
    upper = raw.upper().replace(" ", "")
    for known in CPIC_LEVELS:
        if upper == known.upper().replace(" ", ""):
            return known
    if upper == "RETIRED":
        return "Retired"
    return ""


def _is_pathogenic(code: int | None) -> bool:
    return code in (4, 5)


def _gwas_replicated(finding: dict) -> bool:
    """True when the finding rests on a replicated GWAS association."""
    studies = _num(finding.get("gwas_studies"), 0) or 0
    if studies >= 2:
        return True
    sources = finding.get("sources") or []
    if isinstance(sources, str):
        sources = [sources]
    return any("gwas" in _text(s) for s in sources)


# ---------------------------------------------------------------------------
# Base magnitude
# ---------------------------------------------------------------------------

def base_magnitude(finding: dict) -> tuple[float, str]:
    """Return (base score, the evidence key that produced it).

    Evaluated strongest first so a variant with both a CPIC A assignment and a
    weak ClinVar record scores on the CPIC assignment.
    """
    cpic = normalize_cpic_level(finding.get("cpic_level"))
    code = finding.get("clinvar_sig_code")
    if code is None:
        code = clinvar_sig_code(finding.get("clinical_sig"))
    stars = finding.get("review_stars")
    stars = int(stars) if isinstance(stars, int) else review_stars(finding.get("review_status"))
    fda = str(finding.get("pgx_level") or "").strip()

    if cpic == "A":
        return BASE_SCORES["cpic_a"], "cpic_a"
    if _is_pathogenic(code) and stars >= 3:
        return BASE_SCORES["clinvar_path_3star"], "clinvar_path_3star"
    if cpic in ("A/B", "B"):
        return BASE_SCORES["cpic_b"], "cpic_b"
    if code == 5 and stars == 2:
        return BASE_SCORES["clinvar_path_2star"], "clinvar_path_2star"
    if fda in ("Testing Required", "Testing Recommended"):
        return BASE_SCORES["fda_testing"], "fda_testing"
    if code == 4 and stars >= 2:
        return BASE_SCORES["clinvar_lp_2star"], "clinvar_lp_2star"
    if _gwas_replicated(finding):
        return BASE_SCORES["gwas_replicated"], "gwas_replicated"
    if code is not None and stars <= 1 and code in (1, 4, 5, 255):
        return BASE_SCORES["clinvar_single"], "clinvar_single"
    return BASE_SCORES["default"], "default"


def evidence_label(finding: dict) -> str:
    """Short human label for the strongest evidence behind a finding."""
    cpic = normalize_cpic_level(finding.get("cpic_level"))
    if cpic and cpic != "Retired":
        return f"CPIC Level {cpic}"
    stars = finding.get("review_stars")
    stars = int(stars) if isinstance(stars, int) else review_stars(finding.get("review_status"))
    code = finding.get("clinvar_sig_code")
    if code is None:
        code = clinvar_sig_code(finding.get("clinical_sig"))
    if code is not None and stars >= 1:
        star_text = "star" if stars == 1 else "stars"
        sig = str(finding.get("clinical_sig") or "ClinVar record").strip()
        return f"ClinVar {sig}, {stars} {star_text}"
    fda = str(finding.get("pgx_level") or "").strip()
    if fda:
        return f"FDA label: {fda}"
    if _gwas_replicated(finding):
        n = int(_num(finding.get("gwas_studies"), 0) or 0)
        return f"GWAS, {n} studies" if n else "GWAS association"
    if str(finding.get("clinical_sig") or "").strip():
        return str(finding["clinical_sig"]).strip()
    return ""


# ---------------------------------------------------------------------------
# Magnitude
# ---------------------------------------------------------------------------

def compute_magnitude(finding: dict) -> dict:
    """Compute the DNAInsight magnitude for one finding.

    Returns::

        {"magnitude": float,
         "base": float,
         "evidence_key": str,
         "factors": list of str,
         "dubious": bool}

    ``factors`` is the audit trail. Every multiplier and addend that fired
    appears in it, in the order applied, so the number is explainable.
    """
    factors: list[str] = []
    dubious = bool(finding.get("dubious"))

    base, key = base_magnitude(finding)
    score = base
    factors.append(f"base {base:.2f} from {key}")

    # 1. No-call. A failed probe is not a finding.
    if _text(finding.get("zygosity")) == "no_call":
        factors.append("no-call, forced to 0")
        return {"magnitude": 0.0, "base": base, "evidence_key": key,
                "factors": factors, "dubious": True}

    # 2. Carrier status. The largest honesty correction available.
    copies = finding.get("variant_copies")
    if isinstance(copies, int):
        if copies == 0:
            score *= 0.25
            factors.append("not a carrier, x0.25")
        elif copies == 2:
            score *= 1.3
            factors.append("homozygous for the variant, x1.3")

    # 3. Rarity. Rare genotypes are more interesting, majority ones less.
    band = _text(finding.get("freq_band"))
    if band == "very_rare":
        score += 0.5
        factors.append("very rare genotype, +0.5")
    elif band == "rare":
        score += 0.25
        factors.append("rare genotype, +0.25")
    elif band == "majority":
        score -= 0.5
        factors.append("majority genotype, -0.5")

    # 4. Literature weight, saturating so a heavily studied common variant
    #    cannot outrank a well-evidenced actionable one on citation count alone.
    pubs = int(_num(finding.get("publications"), 0) or 0)
    if pubs > 0:
        bump = min(1.0, math.log10(1 + pubs) / 2.0)
        score += bump
        factors.append(f"{pubs} publications, +{bump:.2f}")

    # 5. Strand ambiguity. An unverifiable call must not outrank a verifiable
    #    one, however strong its evidence looks.
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        if score > 2.0:
            factors.append(f"palindromic site, capped from {score:.2f} to 2.00")
            score = 2.0
        else:
            factors.append("palindromic site, strand not verifiable")
        dubious = True

    score = max(MAGNITUDE_MIN, min(MAGNITUDE_MAX, score))
    return {"magnitude": round(score, 2), "base": base, "evidence_key": key,
            "factors": factors, "dubious": dubious}


# ---------------------------------------------------------------------------
# Repute
# ---------------------------------------------------------------------------

def _direction(text: str) -> str:
    """Classify free text as 'risk', 'protective', 'neutral' or ''."""
    if not text:
        return ""
    for term in NEUTRAL_TERMS:
        if term in text:
            return "neutral"
    risk_hit = next((t for t in RISK_TERMS if t in text), "")
    prot_hit = next((t for t in PROTECTIVE_TERMS if t in text), "")
    if risk_hit and prot_hit:
        # Longer phrase wins; a tie is genuinely mixed, so stay neutral.
        if len(risk_hit) > len(prot_hit):
            return "risk"
        if len(prot_hit) > len(risk_hit):
            return "protective"
        return "neutral"
    if risk_hit:
        return "risk"
    if prot_hit:
        return "protective"
    return ""


def compute_repute(finding: dict) -> str:
    """Return 'Good', 'Bad' or '' for a finding.

    Rules, in order:

    1. Traits and polygenic scores are ALWAYS ''. A trait is not good or bad,
       and colouring one green or red is editorialising about a person.
    2. A no-call is ''. Nothing is known.
    3. A confirmed non-carrier is '', because the classification describes an
       allele they do not have.
    4. Otherwise the direction of effect decides, from the genotype-specific
       text first and the position-level text second.
    """
    entity = _text(finding.get("entity_type"))
    if entity in ("trait", "prs"):
        return ""
    if _text(finding.get("zygosity")) == "no_call":
        return ""
    if finding.get("variant_copies") == 0:
        return ""

    code = finding.get("clinvar_sig_code")
    if code is None:
        code = clinvar_sig_code(finding.get("clinical_sig"))

    # Genotype-specific text is the better signal, so it is checked first.
    for field in ("summary", "interpretation", "conditions", "clinical_sig"):
        direction = _direction(_text(finding.get(field)))
        if direction == "risk":
            return "Bad"
        if direction == "protective":
            return "Good"
        if direction == "neutral":
            return ""

    if code in (4, 5):
        return "Bad"
    if code in (2, 3):
        return "Good"
    return ""


# ---------------------------------------------------------------------------
# Confidence
# ---------------------------------------------------------------------------

def compute_confidence(finding: dict) -> str:
    """Return 'high', 'moderate', 'low' or 'none' for a finding."""
    if _text(finding.get("zygosity")) == "no_call":
        return "none"
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        return "none"

    cpic = normalize_cpic_level(finding.get("cpic_level"))
    stars = finding.get("review_stars")
    stars = int(stars) if isinstance(stars, int) else review_stars(finding.get("review_status"))

    if cpic == "A" or stars >= 3:
        return "high"
    if cpic in ("A/B", "B") or stars == 2:
        return "moderate"
    if stars == 1 or _gwas_replicated(finding):
        return "low"
    if str(finding.get("clinical_sig") or "").strip():
        return "low"
    return "none"


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def score_finding(finding: dict, prefer_snpedia: bool = True) -> dict:
    """Attach magnitude, repute, confidence and evidence to a finding, in place.

    Sets: magnitude, magnitude_source, magnitude_base, magnitude_factors,
    repute, confidence, evidence, review_stars, clinvar_sig_code, dubious.

    When ``prefer_snpedia`` is true and the finding already carries values from
    a local SNPedia cache (``snpedia_magnitude`` or ``snpedia_repute``), those
    win and ``magnitude_source`` is set to ``snpedia``. That keeps the curated
    human judgement when the user has opted in to fetching it, while never
    requiring it.
    """
    if not isinstance(finding, dict):
        return finding

    finding.setdefault("entity_type", "snp")
    finding["review_stars"] = review_stars(finding.get("review_status")) \
        if not isinstance(finding.get("review_stars"), int) else finding["review_stars"]
    if finding.get("clinvar_sig_code") is None:
        finding["clinvar_sig_code"] = clinvar_sig_code(finding.get("clinical_sig"))

    result = compute_magnitude(finding)
    computed_mag = result["magnitude"]
    computed_repute = compute_repute(finding)

    snp_mag = _num(finding.get("snpedia_magnitude"))
    snp_repute = finding.get("snpedia_repute")

    if prefer_snpedia and snp_mag is not None:
        finding["magnitude"] = round(max(MAGNITUDE_MIN, min(MAGNITUDE_MAX, snp_mag)), 2)
        finding["magnitude_source"] = "snpedia"
        result["factors"].append(
            f"local SNPedia cache supplied magnitude {snp_mag}, computed was {computed_mag}"
        )
    else:
        finding["magnitude"] = computed_mag
        finding["magnitude_source"] = "computed"

    if prefer_snpedia and snp_repute in ("Good", "Bad"):
        entity = _text(finding.get("entity_type"))
        finding["repute"] = "" if entity in ("trait", "prs") else snp_repute
    else:
        finding["repute"] = computed_repute

    finding["magnitude_base"] = result["base"]
    finding["magnitude_factors"] = result["factors"]
    finding["dubious"] = bool(result["dubious"])
    finding["confidence"] = compute_confidence(finding)
    finding["evidence"] = finding.get("evidence") or evidence_label(finding)
    return finding


def score_all(findings: list, prefer_snpedia: bool = True) -> list:
    """Score every finding in a list in place and return the same list."""
    for f in findings or []:
        score_finding(f, prefer_snpedia=prefer_snpedia)
    return findings


def sort_key_magnitude(finding: dict) -> float:
    """Sort key that treats an unscored magnitude as 1, per the convention."""
    value = _num(finding.get("magnitude"))
    return UNSCORED_SORT_VALUE if value is None else value
