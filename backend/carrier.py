"""
carrier.py -- carrier screening, residual risk and ACMG secondary-findings coverage.

WHY THIS MODULE EXISTS
----------------------
The arithmetic in :func:`residual_risk` is the piece nobody else in this market
ships, and it is the piece that decides whether a carrier result means anything.

  CFTR has more than two thousand known pathogenic variants. 23andMe reports
  roughly thirty of them.
  Xcode suspended its carrier report in April 2023 behind a notice promising a
  return "within the next few weeks". That notice was still live in August 2026.
  SelfDecode gates carrier status behind its most expensive kit and refuses
  uploaded files.

None of them states residual risk. That omission is not a rounding error. For a
couple planning a pregnancy, the difference between "no pathogenic variant was
detected in the thirty positions we looked at" and "you are not a carrier" is
the entire answer, and the second sentence is one nobody is entitled to say.

This module is the project's invariant 3 in its most consequential form: NOT
PRESENT and NEVER CHECKED are different states, and the wording that
distinguishes them is enforced in code rather than left to a report template.

THREE THINGS ENFORCED HERE RATHER THAN DOCUMENTED
--------------------------------------------------
1. The bare phrase "not a carrier" is never emitted. Every negative statement is
   scoped: "not a carrier for the N variants tested". :func:`audit_wording`
   exists so a test can grep every produced string, and it does.

2. Residual risk returns None with a stated reason when the detection rate or
   the baseline carrier frequency is unknown for that population. Substituting a
   plausible-looking number would be the single most harmful thing this file
   could do, because it would look like the honest feature while being a guess.

3. Every published figure carries ``verified`` True or False. Nothing in the
   frequency or detection-rate tables was invented; the ones marked False were
   recalled from the carrier-screening literature and were NOT re-verified at
   source in this build. They are usable, they are flagged, and
   :func:`unverified_figures` produces the list for the documentation.

DETECTION RATE IS NOT DNAINSIGHT'S DETECTION RATE
--------------------------------------------------
This is the subtlety that makes the whole calculation honest or dishonest. The
published detection rates below belong to CLINICAL PANELS, for example the
23-variant ACMG/ACOG CFTR panel. DNAInsight tests whatever handful of positions
a consumer array happens to carry, which is fewer. So a residual risk computed
from a published panel detection rate is a LOWER BOUND: the real residual risk
for an array user is higher, not lower. Every result carries
``is_lower_bound`` and says so in words.

OFFLINE CONTRACT
----------------
No network access on any path, at import or at call time.
"""

from __future__ import annotations

from typing import Any, Iterable

__all__ = [
    "CARRIER_PANEL", "PANEL_GENES", "CARRIER_POPULATIONS", "POPULATION_ALIASES",
    "STATUS_CARRIER", "STATUS_NEGATIVE_FOR_TESTED", "STATUS_UNTESTABLE",
    "FORBIDDEN_PHRASES", "DISCLAIMER", "UNCERTAINTY_FACTOR",
    "ACMG_SF_GENES", "ACMG_SF_VERSION", "ACMG_SF_LIST_VERIFIED",
    "ACMG_ARRAY_PROBES", "ACMG_KNOWN_RELEVANT_POSITIONS",
    "carrier_status", "carrier_report", "residual_risk", "residual_risk_value",
    "joint_reproductive_risk", "acmg_coverage_report",
    "negative_wording", "has_forbidden_phrasing", "audit_wording",
    "unverified_figures", "as_fraction", "normalise_population",
]


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------

STATUS_CARRIER = "carrier"
STATUS_NEGATIVE_FOR_TESTED = "not_carrier_for_tested_variants"
STATUS_UNTESTABLE = "untestable"

_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "?", ".", "NN"}

# Carrier-screening literature reports frequencies by broad ancestry group, not
# by the 1000 Genomes panels that backend/frequency.py uses. Inventing a mapping
# from CEU to "European" and pretending the numbers transfer cleanly would be
# false precision, so this module keeps its own codes and publishes the alias
# map rather than hiding it inside a function.
CARRIER_POPULATIONS: tuple[str, ...] = (
    "EUROPEAN", "ASHKENAZI", "AFRICAN", "HISPANIC", "EAST_ASIAN",
    "SOUTH_ASIAN", "MIDDLE_EASTERN", "GENERAL",
)

POPULATION_ALIASES: dict[str, str] = {
    # 1000 Genomes codes used elsewhere in this project.
    "CEU": "EUROPEAN", "GBR": "EUROPEAN", "TSI": "EUROPEAN",
    "IBS": "EUROPEAN", "FIN": "EUROPEAN", "EUR": "EUROPEAN", "NFE": "EUROPEAN",
    "YRI": "AFRICAN", "LWK": "AFRICAN", "ASW": "AFRICAN", "AFR": "AFRICAN",
    "CHB": "EAST_ASIAN", "JPT": "EAST_ASIAN", "CHS": "EAST_ASIAN",
    "EAS": "EAST_ASIAN",
    "GIH": "SOUTH_ASIAN", "PJL": "SOUTH_ASIAN", "SAS": "SOUTH_ASIAN",
    "MXL": "HISPANIC", "PUR": "HISPANIC", "CLM": "HISPANIC",
    "AMR": "HISPANIC",
    # Words a user might type.
    "EUROPEAN": "EUROPEAN", "WHITE": "EUROPEAN", "CAUCASIAN": "EUROPEAN",
    "ASHKENAZI": "ASHKENAZI", "ASHKENAZI JEWISH": "ASHKENAZI", "AJ": "ASHKENAZI",
    "AFRICAN": "AFRICAN", "AFRICAN AMERICAN": "AFRICAN", "BLACK": "AFRICAN",
    "HISPANIC": "HISPANIC", "LATINO": "HISPANIC",
    "EAST ASIAN": "EAST_ASIAN", "EAST_ASIAN": "EAST_ASIAN", "ASIAN": "EAST_ASIAN",
    "SOUTH ASIAN": "SOUTH_ASIAN", "SOUTH_ASIAN": "SOUTH_ASIAN",
    "MIDDLE EASTERN": "MIDDLE_EASTERN", "MIDDLE_EASTERN": "MIDDLE_EASTERN",
    "GENERAL": "GENERAL", "": "GENERAL",
}

DISCLAIMER = (
    "DNAInsight is not a medical device and does not perform clinical carrier "
    "screening. A "
    "consumer DNA array reads a small, fixed set of positions chosen for other "
    "reasons, so it can find a variant it happens to carry a probe for and can "
    "never exclude the ones it does not. A negative result here reduces your "
    "chance of being a carrier; it does not remove it. Anyone making a "
    "reproductive decision needs clinical carrier screening ordered through a "
    "genetic counsellor or clinician, which is a different test with a stated "
    "detection rate. The American Board of Genetic Counseling maintains a "
    "directory at findageneticcounselor.com ."
)

# Phrases that must never appear in any string this module produces. Matched
# case-insensitively as substrings, with one exception handled in
# has_forbidden_phrasing: "not a carrier" is allowed only when it is immediately
# scoped by "for the".
FORBIDDEN_PHRASES: tuple[str, ...] = (
    "not a carrier",
    "non-carrier",
    "noncarrier",
    "you are clear",
    "no risk",
    "rules out",
    "ruled out",
    "negative for cystic fibrosis",
)

# Both inputs to a joint reproductive risk are population averages built on
# figures this build could not verify, so a single number would be false
# precision. A factor of two either way is the honest resolution of the
# underlying data, and it is stated rather than implied.
UNCERTAINTY_FACTOR = 2.0


# ---------------------------------------------------------------------------
# The carrier panel.
#
# Per gene:
#   tested_variants   what DNAInsight could look for on a consumer array. Each
#                     carries array_testable, because an indel or a repeat is
#                     not readable even when the position "exists".
#   total_known_pathogenic  the approximate TOTAL number of known pathogenic
#                     variants in the gene. This is the denominator that makes
#                     a negative result interpretable, and it is the number
#                     every competitor omits.
#   carrier_frequency per population, as a probability.
#   detection_rate    per population, for a named CLINICAL panel. Not for the
#                     handful of variants DNAInsight can read. See the module
#                     docstring.
#
# Every numeric figure carries verified True or False and a basis string. None
# of them were invented. The ones marked False were recalled from the carrier
# screening literature and not re-verified at source during this build.
# ---------------------------------------------------------------------------

def _fig(value: float | None, verified: bool, basis: str) -> dict:
    return {"value": value, "verified": verified, "basis": basis}


_ACMG_CFTR_BASIS = (
    "Detection rates commonly published for the 23-variant ACMG/ACOG CFTR "
    "panel. Recalled, not re-verified at source in this build."
)

CARRIER_PANEL: dict[str, dict] = {

    "CFTR": {
        "gene": "CFTR",
        "condition": "Cystic fibrosis",
        "inheritance": "autosomal recessive",
        "chromosome": "7",
        "strand": "+",
        "tested_variants": (
            {"rsid": "rs113993960", "name": "F508del (c.1521_1523delCTT)",
             "allele": "D", "array_testable": False, "verified": False,
             "note": "The single most common CF-causing variant, roughly 70 "
                     "percent of CF alleles in European ancestry. It is a "
                     "three-base deletion. Arrays report substitutions; vendors "
                     "that carry this position report it through proprietary "
                     "i-numbered probes with inconsistent D/I tokens, so this "
                     "module treats it as not readable rather than risk the "
                     "failure mode backend/traits.py documents for rs8176719, "
                     "where a failed probe read as a deletion produced a "
                     "confident wrong answer."},
            {"rsid": "rs113993959", "name": "G542X (c.1624G>T)",
             "allele": "T", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping; recalled, not "
                     "corroborated in-tree."},
            {"rsid": "rs75527207", "name": "G551D (c.1652G>A)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping; recalled, not "
                     "corroborated in-tree."},
            {"rsid": "rs77010898", "name": "W1282X (c.3846G>A)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. Common Ashkenazi "
                     "founder allele."},
            {"rsid": "rs78655421", "name": "R117H (c.350G>A)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. Consequence depends "
                     "on the poly-T tract in cis, which an array cannot read, "
                     "so even a positive result here is not interpretable "
                     "without clinical testing."},
        ),
        "total_known_pathogenic": _fig(
            2100, False,
            "Order of magnitude for CFTR variants catalogued in CFTR1, of which "
            "roughly 400 are classified CF-causing in CFTR2. Recalled, not "
            "re-verified at source."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 25, False, "Widely published as roughly 1 in 25 "
                                            "in non-Hispanic White populations."),
            "ASHKENAZI": _fig(1 / 24, False, "Widely published as roughly 1 in 24."),
            "HISPANIC": _fig(1 / 58, False, "Widely published as roughly 1 in 58."),
            "AFRICAN": _fig(1 / 61, False, "Widely published as roughly 1 in 61."),
            "EAST_ASIAN": _fig(1 / 94, False, "Widely published as roughly 1 in 94."),
        },
        "detection_rate": {
            "EUROPEAN": _fig(0.88, False, _ACMG_CFTR_BASIS),
            "ASHKENAZI": _fig(0.94, False, _ACMG_CFTR_BASIS),
            "HISPANIC": _fig(0.72, False, _ACMG_CFTR_BASIS),
            "AFRICAN": _fig(0.65, False, _ACMG_CFTR_BASIS),
            "EAST_ASIAN": _fig(0.49, False, _ACMG_CFTR_BASIS),
        },
        "detection_rate_panel": "23-variant ACMG/ACOG CFTR panel",
        "detection_rate_panel_size": 23,
        "notes": (
            "CFTR is the clearest example of why residual risk matters. Even a "
            "full 23-variant clinical panel leaves roughly one in eight European "
            "carriers undetected, and the array positions this module can read "
            "are fewer than that panel.",
        ),
    },

    "HEXA": {
        "gene": "HEXA",
        "condition": "Tay-Sachs disease",
        "inheritance": "autosomal recessive",
        "chromosome": "15",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs387906309", "name": "c.1274_1277dupTATC",
             "allele": "I", "array_testable": False, "verified": False,
             "note": "UNVERIFIED rsID. A four-base duplication and the most "
                     "common Ashkenazi allele. Not readable on an array."},
            {"rsid": "rs147324677", "name": "c.1421+1G>C",
             "allele": "C", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping."},
            {"rsid": "rs121907954", "name": "p.Gly269Ser (adult onset)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base. "
                     "HEXA is on the minus strand."},
        ),
        "total_known_pathogenic": _fig(
            130, False,
            "Order of magnitude for HEXA pathogenic variants reported in the "
            "literature. Recalled, not re-verified at source."),
        "carrier_frequency": {
            "ASHKENAZI": _fig(1 / 27, False, "Widely published as roughly 1 in 27."),
            "EUROPEAN": _fig(1 / 300, False, "Widely published as roughly 1 in 300 "
                                             "in the general population."),
            "GENERAL": _fig(1 / 300, False, "Widely published as roughly 1 in 300."),
        },
        "detection_rate": {
            "ASHKENAZI": _fig(0.94, False,
                              "Detection rate commonly published for the "
                              "three-variant Ashkenazi molecular panel. Recalled, "
                              "not re-verified."),
        },
        "detection_rate_panel": "Three-variant Ashkenazi HEXA molecular panel",
        "detection_rate_panel_size": 3,
        "notes": (
            "Molecular HEXA panels perform well in Ashkenazi ancestry and poorly "
            "elsewhere, which is why enzyme assay remains the standard for "
            "non-Ashkenazi carrier screening. No molecular detection rate is "
            "recorded here for other populations, so residual risk for them "
            "returns None rather than a borrowed number.",
        ),
    },

    "SMN1": {
        "gene": "SMN1",
        "condition": "Spinal muscular atrophy",
        "inheritance": "autosomal recessive",
        "chromosome": "5",
        "strand": "-",
        # Deliberately empty. This is the honest answer, not an oversight.
        "tested_variants": (),
        "total_known_pathogenic": _fig(
            None, False,
            "Not expressible as a variant count. Roughly 95 percent of SMA is "
            "caused by homozygous deletion of SMN1 exon 7, which is a copy "
            "number state rather than a sequence variant."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 47, False, "Widely published as roughly 1 in 47."),
            "ASHKENAZI": _fig(1 / 67, False, "Widely published as roughly 1 in 67."),
            "AFRICAN": _fig(1 / 72, False, "Widely published as roughly 1 in 72."),
            "EAST_ASIAN": _fig(1 / 59, False, "Widely published as roughly 1 in 59."),
            "HISPANIC": _fig(1 / 68, False, "Widely published as roughly 1 in 68."),
        },
        "detection_rate": {
            "EUROPEAN": _fig(0.0, True,
                             "Zero, and this one IS verified, because it is a "
                             "property of the assay rather than of any "
                             "population. A SNP array measures bases, not gene "
                             "copy number, so SMN1 carrier status is undetectable "
                             "from array data by construction."),
            "GENERAL": _fig(0.0, True,
                            "Zero. SMN1 carrier status is a copy number state "
                            "and a SNP array cannot measure copy number."),
        },
        "detection_rate_panel": "None. Array data cannot test SMN1 at all.",
        "detection_rate_panel_size": 0,
        "notes": (
            "SMN1 COPY NUMBER IS UNDETECTABLE ON A CONSUMER ARRAY. SMA carrier "
            "screening counts SMN1 copies by quantitative PCR or MLPA. Nothing "
            "in a raw array file speaks to it. Residual risk after an array "
            "'test' for SMN1 is therefore identical to the population carrier "
            "frequency: the test carries no information at all.",
            "Even a proper copy number assay misses the roughly two percent of "
            "carriers who have two SMN1 copies on one chromosome and none on the "
            "other, so a clinical SMA negative also leaves residual risk.",
        ),
    },

    "HBB": {
        "gene": "HBB",
        "condition": "Sickle cell disease and beta thalassemia",
        "inheritance": "autosomal recessive",
        "chromosome": "11",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs334", "name": "HbS, p.Glu6Val (sickle)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED plus-strand base. HBB is transcribed from the "
                     "minus strand, so the c.20A>T of the literature is not the "
                     "base an array reports. A flip here would invert every "
                     "sickle carrier call in the report, which is why it is "
                     "flagged rather than assumed."},
            {"rsid": "rs33930165", "name": "HbC, p.Glu6Lys",
             "allele": "T", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base."},
        ),
        "total_known_pathogenic": _fig(
            400, False,
            "Order of magnitude for beta thalassemia alleles catalogued in "
            "HbVar, alongside several hundred structural haemoglobin variants. "
            "Recalled, not re-verified at source."),
        "carrier_frequency": {
            "AFRICAN": _fig(1 / 13, False,
                            "Sickle cell trait, widely published as roughly 1 in "
                            "13 in African American populations."),
            "MIDDLE_EASTERN": _fig(1 / 30, False,
                                   "Beta thalassemia trait, order of magnitude "
                                   "for Mediterranean and Middle Eastern "
                                   "populations. Recalled, not verified."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. Clinical haemoglobinopathy screening starts with a "
            "full blood count and haemoglobin electrophoresis, not a variant "
            "panel, so a molecular detection rate is not the right measure and "
            "none is invented here."),
        "detection_rate_panel_size": 0,
        "notes": (
            "Testing two point variants says nothing about the several hundred "
            "beta thalassemia alleles, and nothing at all about alpha "
            "thalassemia, which is caused by HBA1 and HBA2 deletions that an "
            "array cannot see.",
        ),
    },

    "GJB2": {
        "gene": "GJB2",
        "condition": "Nonsyndromic hearing loss (DFNB1)",
        "inheritance": "autosomal recessive",
        "chromosome": "13",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs80338939", "name": "c.35delG",
             "allele": "D", "array_testable": False, "verified": False,
             "note": "UNVERIFIED rsID. A single-base deletion accounting for "
                     "most GJB2 alleles in European ancestry, and not readable "
                     "on an array."},
            {"rsid": "rs80338943", "name": "c.167delT",
             "allele": "D", "array_testable": False, "verified": False,
             "note": "UNVERIFIED rsID. Ashkenazi founder deletion, not "
                     "array-readable."},
            {"rsid": "rs72474224", "name": "c.109G>A (p.Val37Ile)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. Mild and common in "
                     "East Asian ancestry, with incomplete penetrance."},
        ),
        "total_known_pathogenic": _fig(
            100, False,
            "Order of magnitude for GJB2 pathogenic variants. Recalled, not "
            "re-verified at source."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 33, False, "Order of magnitude commonly "
                                            "published for GJB2 carriers."),
            "ASHKENAZI": _fig(1 / 25, False, "Driven by c.167delT. Recalled, not "
                                             "verified."),
            "EAST_ASIAN": _fig(1 / 50, False, "Order of magnitude. Recalled, not "
                                              "verified."),
        },
        "detection_rate": {
            "EUROPEAN": _fig(0.0, True,
                             "Zero for array data. The dominant European allele "
                             "c.35delG is a deletion, and the array cannot read "
                             "it, so the array tests none of the common alleles "
                             "in this population."),
        },
        "detection_rate_panel": (
            "None for array data. Clinical GJB2 screening sequences the gene."),
        "detection_rate_panel_size": 0,
        "notes": (
            "GJB2 is the cleanest demonstration in this panel that a gene can be "
            "famous, common and completely untestable on a consumer array, "
            "because the alleles that matter are deletions.",
        ),
    },

    "PAH": {
        "gene": "PAH",
        "condition": "Phenylketonuria",
        "inheritance": "autosomal recessive",
        "chromosome": "12",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs5030858", "name": "p.Arg408Trp",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base."},
            {"rsid": "rs75193786", "name": "p.Arg261Gln",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base."},
        ),
        "total_known_pathogenic": _fig(
            1000, False,
            "Order of magnitude for PAH variants catalogued in PAHvdb. "
            "Recalled, not re-verified at source."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 50, False, "Order of magnitude commonly "
                                            "published for PKU carriers in "
                                            "European ancestry."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. PAH is highly allelically heterogeneous and clinical "
            "testing sequences the gene rather than using a variant panel, so no "
            "panel detection rate is recorded and none is invented."),
        "detection_rate_panel_size": 0,
        "notes": (
            "PKU is detected by newborn screening in most countries, which is a "
            "biochemical test and does not depend on knowing the genotype.",
        ),
    },

    "ATP7B": {
        "gene": "ATP7B",
        "condition": "Wilson disease",
        "inheritance": "autosomal recessive",
        "chromosome": "13",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs76151636", "name": "p.His1069Gln",
             "allele": "T", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base. "
                     "The most common European ATP7B allele, but it still "
                     "accounts for well under half of European alleles."},
        ),
        "total_known_pathogenic": _fig(
            800, False,
            "Order of magnitude for ATP7B pathogenic variants reported in the "
            "literature. Recalled, not re-verified at source."),
        "carrier_frequency": {
            "GENERAL": _fig(1 / 90, False, "Order of magnitude commonly published."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. One variant is not a panel and no published "
            "single-variant detection rate is recorded here."),
        "detection_rate_panel_size": 0,
        "notes": (
            "Wilson disease is treatable, which makes it one of the few carrier "
            "panel genes where a finding in an adult can matter to that adult "
            "and not only to a future pregnancy.",
        ),
    },

    "GALT": {
        "gene": "GALT",
        "condition": "Classic galactosemia",
        "inheritance": "autosomal recessive",
        "chromosome": "9",
        "strand": "+",
        "tested_variants": (
            {"rsid": "rs75391579", "name": "p.Gln188Arg",
             "allele": "G", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base."},
            {"rsid": "rs2070074", "name": "p.Asn314Asp (Duarte)",
             "allele": "G", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. The Duarte variant is "
                     "NOT classic galactosemia and reporting it as a carrier "
                     "finding without that distinction would be misleading."},
        ),
        "total_known_pathogenic": _fig(
            330, False,
            "Order of magnitude for GALT variants reported in the literature. "
            "Recalled, not re-verified at source."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 107, False, "Order of magnitude commonly published."),
        },
        "detection_rate": {},
        "detection_rate_panel": "None recorded.",
        "detection_rate_panel_size": 0,
        "notes": (
            "Classic galactosemia is on newborn screening panels in most "
            "countries.",
        ),
    },

    "ACADM": {
        "gene": "ACADM",
        "condition": "Medium-chain acyl-CoA dehydrogenase deficiency",
        "inheritance": "autosomal recessive",
        "chromosome": "1",
        "strand": "+",
        "tested_variants": (
            {"rsid": "rs77931234", "name": "c.985A>G (p.Lys304Glu)",
             "allele": "G", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. Accounts for the large "
                     "majority of MCAD alleles in Northern European ancestry, "
                     "which is unusual and makes single-variant testing more "
                     "informative here than for most genes in this panel."},
        ),
        "total_known_pathogenic": _fig(
            100, False,
            "Order of magnitude for ACADM pathogenic variants. Recalled, not "
            "re-verified at source."),
        "carrier_frequency": {
            "EUROPEAN": _fig(1 / 65, False, "Order of magnitude commonly "
                                            "published for Northern European "
                                            "ancestry."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. The c.985A>G share of alleles is often quoted "
            "around 80 percent in Northern European ancestry, but that figure "
            "was not verified in this build and is therefore not used as a "
            "detection rate."),
        "detection_rate_panel_size": 0,
        "notes": (
            "MCAD deficiency is on newborn screening panels in most countries.",
        ),
    },

    "BTD": {
        "gene": "BTD",
        "condition": "Biotinidase deficiency",
        "inheritance": "autosomal recessive",
        "chromosome": "3",
        "strand": "+",
        "tested_variants": (
            {"rsid": "rs13078881", "name": "p.Asp444His",
             "allele": "C", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping. D444H alone causes "
                     "only partial deficiency, and its consequence depends on "
                     "what is on the other chromosome, which unphased array "
                     "data cannot establish."},
        ),
        "total_known_pathogenic": _fig(
            250, False,
            "Order of magnitude for BTD pathogenic variants. Recalled, not "
            "re-verified at source."),
        "carrier_frequency": {
            "GENERAL": _fig(1 / 120, False, "Order of magnitude commonly "
                                            "published for profound deficiency."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. Biotinidase deficiency is diagnosed by enzyme "
            "activity, not by genotype."),
        "detection_rate_panel_size": 0,
        "notes": (
            "Biotinidase deficiency is treatable with an inexpensive vitamin and "
            "is on newborn screening panels in most countries.",
        ),
    },

    "G6PD": {
        "gene": "G6PD",
        "condition": "Glucose-6-phosphate dehydrogenase deficiency",
        "inheritance": "X-linked",
        "chromosome": "X",
        "strand": "-",
        "tested_variants": (
            {"rsid": "rs1050828", "name": "G6PD A- (p.Val68Met)",
             "allele": "T", "array_testable": True, "verified": True,
             "note": "Plus-strand base T is corroborated by "
                     "data/evidence_overlay.py, which records rs1050828 with "
                     "risk allele T at CPIC level A. This is the only variant in "
                     "the whole carrier panel with in-tree corroboration."},
            {"rsid": "rs1050829", "name": "G6PD A (p.Asn126Asp)",
             "allele": "C", "array_testable": True, "verified": False,
             "note": "UNVERIFIED plus-strand base. Needed in cis with rs1050828 "
                     "to define the A- allele, and unphased array data cannot "
                     "establish cis."},
            {"rsid": "rs5030868", "name": "G6PD Mediterranean (p.Ser188Phe)",
             "allele": "A", "array_testable": True, "verified": False,
             "note": "UNVERIFIED rsID-to-variant mapping and plus-strand base. "
                     "The clinically severe variant in Mediterranean, Middle "
                     "Eastern and South Asian ancestry."},
        ),
        "total_known_pathogenic": _fig(
            200, False,
            "Order of magnitude for G6PD variants reported in the literature. "
            "Recalled, not re-verified at source."),
        "carrier_frequency": {
            "AFRICAN": _fig(0.11, False,
                            "Approximate G6PD A- allele frequency in African "
                            "ancestry. Because G6PD is X-linked this is an "
                            "allele frequency, and the proportion of "
                            "heterozygous females is higher than the proportion "
                            "of hemizygous males. Recalled, not verified."),
            "MIDDLE_EASTERN": _fig(0.05, False,
                                   "Approximate Mediterranean allele frequency. "
                                   "Recalled, not verified."),
        },
        "detection_rate": {},
        "detection_rate_panel": (
            "None recorded. G6PD deficiency is confirmed by enzyme activity, and "
            "the enzyme assay can be falsely normal shortly after a haemolytic "
            "episode."),
        "detection_rate_panel_size": 0,
        "notes": (
            "G6PD is X-linked, so carrier arithmetic differs from the rest of "
            "this panel and joint_reproductive_risk handles it separately.",
            "A heterozygous female can still have clinically significant "
            "deficiency because of X inactivation, so 'carrier' understates the "
            "situation for G6PD in a way it does not for the recessive genes.",
        ),
    },
}

PANEL_GENES: tuple[str, ...] = tuple(CARRIER_PANEL)


# ---------------------------------------------------------------------------
# Wording enforcement
# ---------------------------------------------------------------------------

def negative_wording(count: int, gene: str = "") -> str:
    """The ONLY sanctioned way to phrase a negative carrier result.

    "not a carrier" on its own is a claim about a gene. "not a carrier for the N
    variants tested" is a claim about a test, and only the second one is true.
    Every negative statement in this module is built from this function so the
    scoping cannot be dropped by an edit somewhere else.
    """
    n = max(0, int(count))
    noun = "variant" if n == 1 else "variants"
    where = f" in {gene}" if gene else ""
    return f"not a carrier for the {n} {noun} tested{where}"


def has_forbidden_phrasing(text: Any) -> list[str]:
    """Return the forbidden phrases present in ``text``. Empty means clean.

    "not a carrier" is permitted only when immediately followed by "for the",
    which is what :func:`negative_wording` always produces.
    """
    lowered = str(text or "").lower()
    hits: list[str] = []
    for phrase in FORBIDDEN_PHRASES:
        start = 0
        while True:
            idx = lowered.find(phrase, start)
            if idx < 0:
                break
            start = idx + 1
            if phrase == "not a carrier":
                tail = lowered[idx + len(phrase):].lstrip()
                if tail.startswith("for the"):
                    continue
            hits.append(phrase)
            break
    return hits


def audit_wording(payload: Any) -> list[dict]:
    """Walk any nested payload and report strings using forbidden phrasing.

    Lives here rather than in the test file so the phrase list and the walker
    cannot drift apart, the same arrangement backend/diplotype.py uses for
    prescriptive language.
    """
    hits: list[dict] = []

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                walk(value, f"{path}.{key}" if path else str(key))
        elif isinstance(node, (list, tuple)):
            for i, value in enumerate(node):
                walk(value, f"{path}[{i}]")
        elif isinstance(node, str):
            found = has_forbidden_phrasing(node)
            if found:
                hits.append({"path": path, "phrases": found, "text": node})

    walk(payload, "")
    return hits


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def normalise_population(population: Any) -> str:
    """Map a population label or 1000 Genomes code onto a carrier-panel group."""
    raw = str(population or "").strip().upper().replace("-", " ")
    if raw in POPULATION_ALIASES:
        return POPULATION_ALIASES[raw]
    collapsed = raw.replace(" ", "_")
    if collapsed in POPULATION_ALIASES:
        return POPULATION_ALIASES[collapsed]
    if collapsed in CARRIER_POPULATIONS:
        return collapsed
    return ""


def as_fraction(probability: float | None, *, max_denominator: int = 1000000) -> str:
    """Render a probability as '1 in N', which is how people actually read risk.

    Returns '' for None. A percentage with four decimal places is technically the
    same number and is much harder to reason about when comparing 1 in 200 with
    1 in 4000.
    """
    if probability is None:
        return ""
    try:
        value = float(probability)
    except (TypeError, ValueError):
        return ""
    if value <= 0:
        return "0"
    if value >= 1:
        return "1 in 1"
    denominator = round(1.0 / value)
    if denominator > max_denominator:
        return f"below 1 in {max_denominator}"
    return f"1 in {denominator}"


def _index(genotypes: Any) -> dict:
    if not isinstance(genotypes, dict):
        return {}
    return {str(k).strip().lower(): v for k, v in genotypes.items()}


def _read(genotypes: dict, rsid: str) -> str | None:
    key = str(rsid or "").strip().lower()
    if key not in genotypes:
        return None
    value = genotypes[key]
    if isinstance(value, (list, tuple)) and len(value) == 2:
        first, second = str(value[0]).strip().upper(), str(value[1]).strip().upper()
    else:
        raw = str(value if value is not None else "").strip().upper()
        if raw in _NOCALL_ALLELES or len(raw) != 2:
            return None
        first, second = raw[0], raw[1]
    if first in _NOCALL_ALLELES or second in _NOCALL_ALLELES:
        return None
    ordered = sorted((first, second))
    return f"{ordered[0]}{ordered[1]}"


def _copies(genotypes: dict, rsid: str, allele: str) -> int | None:
    key = _read(genotypes, rsid)
    if key is None:
        return None
    wanted = str(allele or "").strip().upper()
    if not wanted:
        return None
    return sum(1 for base in key if base == wanted)


def _canonical_gene(gene: Any) -> str:
    raw = str(gene or "").strip().upper().replace(" ", "")
    for known in CARRIER_PANEL:
        if known.upper() == raw:
            return known
    return ""


# ---------------------------------------------------------------------------
# Carrier status
# ---------------------------------------------------------------------------

def carrier_status(gene: str, genotypes: dict) -> dict:
    """Determine carrier status for one panel gene from array genotypes.

    Returns exactly one of three statuses, and the third one is the reason this
    function exists:

        carrier                          at least one tested variant detected
        not_carrier_for_tested_variants  every variant that COULD be read was
                                         read and none was detected
        untestable                       nothing in this gene could be read, so
                                         nothing was found and nothing was
                                         excluded

    The negative statement is always scoped to the number of variants actually
    read, never to the gene, and never to the panel. If the panel lists five
    variants and this file could read two of them, the answer is "not a carrier
    for the 2 variants tested", not "not a carrier for the 5 variants tested"
    and certainly not "not a carrier".
    """
    canonical = _canonical_gene(gene)
    if not canonical:
        return {
            "gene": str(gene or ""), "known_gene": False,
            "status": STATUS_UNTESTABLE, "copies": None,
            "tested_count": 0, "detected": [], "not_read": [],
            "not_array_testable": [],
            "statement": f"{gene!r} is not a gene in this carrier panel.",
            "caveats": [], "disclaimer": DISCLAIMER,
        }

    entry = CARRIER_PANEL[canonical]
    g = _index(genotypes)

    detected: list[dict] = []
    read_negative: list[dict] = []
    not_read: list[dict] = []
    not_array_testable: list[dict] = []
    unverified_used: list[str] = []
    total_copies = 0

    for variant in entry["tested_variants"]:
        summary = {"rsid": variant["rsid"], "name": variant["name"],
                   "verified": bool(variant.get("verified"))}
        if not variant.get("array_testable"):
            not_array_testable.append({**summary, "reason": variant["note"]})
            continue
        copies = _copies(g, variant["rsid"], variant["allele"])
        if copies is None:
            not_read.append({**summary, "reason": "This position is not in the "
                                                  "file, or the probe returned no "
                                                  "call."})
            continue
        if not variant.get("verified"):
            unverified_used.append(variant["rsid"])
        if copies > 0:
            total_copies += copies
            detected.append({**summary, "copies": copies})
        else:
            read_negative.append(summary)

    tested_count = len(detected) + len(read_negative)
    total_known = entry["total_known_pathogenic"]
    caveats: list[str] = list(entry.get("notes", ()))

    if detected:
        status = STATUS_CARRIER
        names = ", ".join(d["name"] for d in detected)
        if total_copies >= 2:
            statement = (
                f"Two copies of a tested pathogenic {canonical} variant were "
                f"detected ({names}). Two copies is a different situation from "
                f"carrying one, and it needs clinical confirmation before it is "
                f"used for anything, because array calls are not clinical-grade "
                f"and because unphased data cannot show whether two different "
                f"variants sit on the same chromosome or on opposite ones."
            )
        else:
            statement = (
                f"One copy of a tested pathogenic {canonical} variant was "
                f"detected ({names}). That is a carrier result for that variant "
                f"and needs confirmation by a clinical-grade test before it is "
                f"used for anything."
            )
    elif tested_count > 0:
        status = STATUS_NEGATIVE_FOR_TESTED
        statement = (
            f"No variant was detected at the {tested_count} {canonical} "
            f"position(s) this file could read, so this result is "
            f"{negative_wording(tested_count)}."
        )
        if total_known["value"]:
            statement += (
                f" {canonical} has approximately {total_known['value']} known "
                f"pathogenic variants, so this leaves most of the gene "
                f"unexamined."
            )
        else:
            statement += (
                f" The number of pathogenic {canonical} variants is not "
                f"expressible as a count here, so the share examined cannot be "
                f"quoted."
            )
    else:
        status = STATUS_UNTESTABLE
        statement = (
            f"No {canonical} variant in this panel could be read from this file, "
            f"so carrier status was not tested. Nothing was found and nothing "
            f"was excluded."
        )

    if not_read:
        caveats.append(
            "Not read from this file: "
            + ", ".join(f"{v['name']} ({v['rsid']})" for v in not_read)
        )
    if not_array_testable:
        caveats.append(
            "Not readable on any SNP array, regardless of vendor: "
            + ", ".join(f"{v['name']} ({v['rsid']})" for v in not_array_testable)
        )
    if unverified_used:
        caveats.append(
            "This result relies on rsID-to-variant mappings that were not "
            "corroborated inside this repository: " + ", ".join(unverified_used)
            + ". See the note on each variant in backend/carrier.py."
        )

    panel_size = len(entry["tested_variants"])
    coverage = (tested_count / total_known["value"]
                if total_known["value"] else None)

    return {
        "gene": canonical,
        "known_gene": True,
        "condition": entry["condition"],
        "inheritance": entry["inheritance"],
        "status": status,
        "copies": total_copies if (detected or tested_count) else None,
        "tested_count": tested_count,
        "panel_size": panel_size,
        "detected": detected,
        "negative": read_negative,
        "not_read": not_read,
        "not_array_testable": not_array_testable,
        "unverified_variants_used": unverified_used,
        "total_known_pathogenic": total_known["value"],
        "total_known_pathogenic_verified": total_known["verified"],
        "total_known_pathogenic_basis": total_known["basis"],
        "fraction_of_known_variants_tested": (
            round(coverage, 6) if coverage is not None else None
        ),
        "statement": statement,
        "caveats": caveats,
        "disclaimer": DISCLAIMER,
    }


def carrier_report(genotypes: dict, *, population: Any = "") -> dict:
    """Carrier status for every gene in the panel, plus residual risk where known."""
    pop = normalise_population(population)
    results: list[dict] = []
    for gene in PANEL_GENES:
        status = carrier_status(gene, genotypes)
        if status["status"] == STATUS_NEGATIVE_FOR_TESTED:
            status["residual_risk"] = residual_risk(gene, True, pop or population)
        elif status["status"] == STATUS_UNTESTABLE:
            status["residual_risk"] = residual_risk(gene, 0, pop or population)
        else:
            status["residual_risk"] = residual_risk(gene, False, pop or population)
        results.append(status)
    counts = {
        STATUS_CARRIER: sum(1 for r in results if r["status"] == STATUS_CARRIER),
        STATUS_NEGATIVE_FOR_TESTED: sum(
            1 for r in results if r["status"] == STATUS_NEGATIVE_FOR_TESTED),
        STATUS_UNTESTABLE: sum(
            1 for r in results if r["status"] == STATUS_UNTESTABLE),
    }
    return {
        "population": pop,
        "population_as_supplied": str(population or ""),
        "results": results,
        "counts": counts,
        "summary": (
            f"{counts[STATUS_CARRIER]} gene(s) with a detected variant, "
            f"{counts[STATUS_NEGATIVE_FOR_TESTED]} gene(s) where every readable "
            f"position was read and nothing was detected, and "
            f"{counts[STATUS_UNTESTABLE]} gene(s) that could not be tested at "
            f"all from this file."
        ),
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------
# Residual risk
# ---------------------------------------------------------------------------

def residual_risk(gene: str,
                  tested_variants_negative: Any,
                  population: Any,
                  *,
                  detection_rate: float | None = None,
                  carrier_frequency: float | None = None) -> dict:
    """Post-test residual carrier risk after a negative panel result.

    DERIVATION
    ----------
    Let f be the baseline carrier frequency in this population and DR the
    detection rate of the panel that was run, meaning the proportion of carriers
    in that population whose variant the panel would find.

        P(carrier)                     = f
        P(negative result | carrier)   = 1 - DR
        P(negative result | carrier)   = 1        for a non-carrier

    By Bayes:

        P(carrier | negative) = P(neg | carrier) P(carrier)
                               -------------------------------------------
                               P(neg | carrier) P(carrier)
                                 + P(neg | non-carrier) P(non-carrier)

                             =        f (1 - DR)
                               ---------------------------
                               f (1 - DR) + (1 - f)

    and since f(1 - DR) + 1 - f = 1 - f DR, this simplifies to

        residual risk = f (1 - DR) / (1 - f DR)

    which is what this function computes.

    WHAT IT REFUSES TO DO
    ---------------------
    When either f or DR is unknown for this gene and population, the returned
    ``residual_risk`` is None and ``reason`` says which input was missing.
    Substituting a plausible-looking number here would produce a figure that
    looks like the honest feature while being a guess, which is worse than
    printing nothing.

    WHY THE ANSWER IS A LOWER BOUND
    -------------------------------
    Published detection rates belong to clinical panels. DNAInsight reads
    whatever positions a consumer array happens to carry, which is fewer, so the
    real residual risk for an array user is HIGHER than the number returned.
    ``is_lower_bound`` is set and the reason says so.

    ``tested_variants_negative`` may be a bool, an int count of variants that
    came back negative, or a list of them. A falsy bool means there was no
    negative result to condition on and the function says so rather than
    computing anything.
    """
    canonical = _canonical_gene(gene)
    pop = normalise_population(population)

    out: dict = {
        "gene": canonical or str(gene or ""),
        "population": pop,
        "population_as_supplied": str(population or ""),
        "residual_risk": None,
        "residual_risk_as_fraction": "",
        "prior_risk": None,
        "prior_risk_as_fraction": "",
        "carrier_frequency": None,
        "detection_rate": None,
        "verified": False,
        "is_lower_bound": True,
        "formula": "residual = f * (1 - DR) / (1 - f * DR)",
        "reason": "",
        "basis": {},
        "disclaimer": DISCLAIMER,
    }

    if not canonical:
        out["reason"] = f"{gene!r} is not a gene in this carrier panel."
        return out

    entry = CARRIER_PANEL[canonical]

    if isinstance(tested_variants_negative, bool):
        negative_count = None if tested_variants_negative else 0
        had_negative = tested_variants_negative
    elif isinstance(tested_variants_negative, (list, tuple, set)):
        negative_count = len(tested_variants_negative)
        had_negative = negative_count > 0
    else:
        try:
            negative_count = int(tested_variants_negative)
        except (TypeError, ValueError):
            negative_count = 0
        had_negative = negative_count > 0

    if not pop:
        out["reason"] = (
            f"No population was supplied, or {str(population)!r} does not map to "
            f"one of the carrier-screening ancestry groups this module holds "
            f"figures for. Carrier frequency and detection rate both depend on "
            f"population, so neither can be chosen without it."
        )
        return out

    freq_entry = entry["carrier_frequency"].get(pop)
    if carrier_frequency is not None:
        f = float(carrier_frequency)
        freq_verified = False
        freq_basis = "Supplied by the caller, not from this module's table."
    elif freq_entry is None or freq_entry["value"] is None:
        out["reason"] = (
            f"No baseline carrier frequency is recorded for {canonical} in the "
            f"{pop} population, so residual risk cannot be computed. Borrowing a "
            f"frequency from another population would produce a number that "
            f"looks precise and is not about this person."
        )
        out["basis"] = {"carrier_frequency": "unknown"}
        return out
    else:
        f = float(freq_entry["value"])
        freq_verified = bool(freq_entry["verified"])
        freq_basis = freq_entry["basis"]

    dr_entry = entry["detection_rate"].get(pop)
    if detection_rate is not None:
        dr = float(detection_rate)
        dr_verified = False
        dr_basis = "Supplied by the caller, not from this module's table."
        dr_panel = "caller supplied"
    elif dr_entry is None or dr_entry["value"] is None:
        out["carrier_frequency"] = f
        out["prior_risk"] = f
        out["prior_risk_as_fraction"] = as_fraction(f)
        out["reason"] = (
            f"No detection rate is recorded for {canonical} in the {pop} "
            f"population, so residual risk cannot be computed. The baseline "
            f"carrier frequency is known and is shown as the prior risk, but "
            f"without a detection rate there is no way to say how much a "
            f"negative result lowered it. {entry['detection_rate_panel']}"
        )
        out["basis"] = {"carrier_frequency": freq_basis, "detection_rate": "unknown"}
        return out
    else:
        dr = float(dr_entry["value"])
        dr_verified = bool(dr_entry["verified"])
        dr_basis = dr_entry["basis"]
        dr_panel = entry["detection_rate_panel"]

    out["carrier_frequency"] = f
    out["detection_rate"] = dr
    out["prior_risk"] = f
    out["prior_risk_as_fraction"] = as_fraction(f)
    out["verified"] = bool(freq_verified and dr_verified)
    out["basis"] = {"carrier_frequency": freq_basis, "detection_rate": dr_basis,
                    "panel": dr_panel}

    if isinstance(tested_variants_negative, bool) and not had_negative:
        out["reason"] = (
            "Residual risk is defined only after a negative result. No negative "
            "result was supplied for this gene, so there is nothing to update."
        )
        return out

    # A count of zero negatives means nothing was actually tested on this file,
    # whatever detection rate the published panel achieves. Applying the panel's
    # DR here would credit this person with a test they did not have, which is
    # the exact error this module exists to prevent, so the effective detection
    # rate drops to zero and the residual risk stays at the population baseline.
    nothing_tested = (not isinstance(tested_variants_negative, bool)
                      and negative_count == 0)
    if nothing_tested:
        dr = 0.0
        out["detection_rate"] = 0.0
        dr_panel = ("none, because no variant in this gene was actually read "
                    "from this file")

    denominator = 1.0 - (f * dr)
    if denominator <= 0:
        out["reason"] = (
            "The carrier frequency and detection rate supplied give a "
            "degenerate denominator, so no residual risk is reported."
        )
        return out

    residual = (f * (1.0 - dr)) / denominator
    out["residual_risk"] = residual
    out["residual_risk_as_fraction"] = as_fraction(residual)

    if nothing_tested:
        out["reason"] = (
            f"No {canonical} variant was read from this file, so nothing was "
            f"tested and nothing can be subtracted. The residual carrier risk is "
            f"identical to the {pop} population carrier frequency of "
            f"{as_fraction(f)}."
        )
    elif dr <= 0:
        out["reason"] = (
            f"The detection rate for {canonical} in this population is zero, so "
            f"a negative result carries no information and the residual risk is "
            f"identical to the population carrier frequency. "
            f"{entry['detection_rate_panel']}"
        )
    else:
        out["reason"] = (
            f"Prior carrier risk {as_fraction(f)} in the {pop} population, "
            f"updated by a negative result on a panel with a detection rate of "
            f"{dr:.0%}, gives a residual carrier risk of {as_fraction(residual)}."
        )

    tested_here = sum(1 for v in entry["tested_variants"] if v.get("array_testable"))
    out["reason"] += (
        f" This is a LOWER BOUND. The detection rate above belongs to "
        f"{dr_panel}, while this module can read at most {tested_here} "
        f"{canonical} position(s) from a consumer array, so the real residual "
        f"risk is higher than the figure shown."
    )
    if not out["verified"]:
        out["reason"] += (
            " At least one of the two inputs is recorded in this build as "
            "unverified, so treat the figure as an order of magnitude."
        )
    return out


def residual_risk_value(gene: str,
                        tested_variants_negative: Any,
                        population: Any,
                        **kwargs: Any) -> float | None:
    """Just the number from :func:`residual_risk`, or None when it is unknown.

    Provided because callers that only want to render a figure should not have
    to remember which key holds it, and because "returns None when unknown" is a
    contract worth being able to state in one line.
    """
    return residual_risk(gene, tested_variants_negative, population,
                         **kwargs)["residual_risk"]


# ---------------------------------------------------------------------------
# Joint reproductive risk
# ---------------------------------------------------------------------------

def _risk_input(value: Any) -> tuple[float | None, bool, str]:
    """Accept a float, a residual_risk payload, or None. Return (risk, verified, label)."""
    if value is None:
        return None, False, "not supplied"
    if isinstance(value, dict):
        risk = value.get("residual_risk")
        if risk is None:
            risk = value.get("prior_risk")
            if risk is None:
                return None, False, value.get("reason", "unknown")
            return float(risk), bool(value.get("verified")), "prior carrier risk"
        return float(risk), bool(value.get("verified")), "residual carrier risk"
    try:
        return float(value), False, "supplied directly"
    except (TypeError, ValueError):
        return None, False, "not a number"


def joint_reproductive_risk(gene: str,
                            risk_a: Any,
                            risk_b: Any,
                            inheritance: str = "") -> dict:
    """Chance that a pregnancy between two people is affected, given each carrier risk.

    AUTOSOMAL RECESSIVE
        Both parents must carry, and two carriers have a one in four chance per
        pregnancy:

            P(affected) = risk_a * risk_b * 0.25

    X-LINKED RECESSIVE
        Different arithmetic and a different question. Only the parent with two
        X chromosomes can pass a recessive X allele to an affected son:

            P(affected son, per pregnancy) = risk_a * 0.5 * 0.5

        The 0.5s are the chance the pregnancy is male and the chance the
        affected X is the one transmitted. ``risk_b`` is not used, and the
        payload says so explicitly rather than silently ignoring it, because
        silently ignoring an argument is how a caller ends up believing the
        other parent was accounted for.

    A RANGE, NOT A NUMBER
        Both inputs are built on population figures this build could not verify,
        so a single value would be false precision. When either input is
        unverified the result is a range spanning a factor of
        ``UNCERTAINTY_FACTOR`` either side of the point estimate, clamped at the
        arithmetic maximum, and ``range_basis`` states why.
    """
    canonical = _canonical_gene(gene)
    entry = CARRIER_PANEL.get(canonical)
    mode = str(inheritance or (entry or {}).get("inheritance") or "").strip().lower()

    a_value, a_verified, a_label = _risk_input(risk_a)
    b_value, b_verified, b_label = _risk_input(risk_b)

    out: dict = {
        "gene": canonical or str(gene or ""),
        "condition": (entry or {}).get("condition", ""),
        "inheritance": mode or "unknown",
        "risk_a": a_value,
        "risk_b": b_value,
        "risk_a_basis": a_label,
        "risk_b_basis": b_label,
        "point": None,
        "low": None,
        "high": None,
        "point_as_fraction": "",
        "range_as_fraction": "",
        "verified": bool(a_verified and b_verified),
        "range_basis": "",
        "assumptions": [],
        "reason": "",
        "caveats": [],
        "disclaimer": DISCLAIMER,
    }

    if "recessive" not in mode and "x-linked" not in mode and "x linked" not in mode:
        out["reason"] = (
            f"Joint reproductive risk is implemented for autosomal recessive and "
            f"X-linked inheritance only. {mode or 'unknown'!r} is neither, so no "
            f"number is produced."
        )
        return out

    x_linked = "x-linked" in mode or "x linked" in mode

    if a_value is None:
        out["reason"] = (
            "The first carrier risk is unknown, so no joint risk can be "
            f"computed ({a_label}). An unknown input cannot be replaced by an "
            "average without inventing the answer."
        )
        return out

    if x_linked:
        point = a_value * 0.5 * 0.5
        ceiling = 0.25
        out["assumptions"] = [
            "risk_a is the carrier risk of the parent with two X chromosomes.",
            "risk_b is not used for X-linked recessive inheritance: the other "
            "parent's carrier status does not put a son at risk. An affected "
            "father has daughters who are obligate carriers and sons who are "
            "unaffected, which is a different question from this one.",
            "The figure is the chance per pregnancy that the child is an "
            "affected son.",
        ]
        used_verified = a_verified
    else:
        if b_value is None:
            out["reason"] = (
                "The second carrier risk is unknown, so no joint risk can be "
                f"computed ({b_label}). For autosomal recessive inheritance both "
                "parents' risks are needed and neither can be assumed."
            )
            return out
        point = a_value * b_value * 0.25
        ceiling = 0.25
        out["assumptions"] = [
            "Both parents' carrier risks are independent, which holds unless "
            "they are related. Consanguinity raises the joint risk and this "
            "calculation does not model it.",
            "Two carriers of an autosomal recessive condition have a one in four "
            "chance per pregnancy of an affected child.",
            "The figure is a chance per pregnancy and does not change with the "
            "outcome of previous pregnancies.",
        ]
        used_verified = a_verified and b_verified

    out["point"] = point
    out["point_as_fraction"] = as_fraction(point)

    if used_verified:
        out["low"] = point
        out["high"] = point
        out["range_basis"] = (
            "Both inputs are recorded as verified in this build, so the point "
            "estimate is reported without an uncertainty band. It still rests on "
            "population averages rather than on this couple."
        )
    else:
        low = point / UNCERTAINTY_FACTOR
        high = min(ceiling, point * UNCERTAINTY_FACTOR)
        out["low"] = low
        out["high"] = high
        out["range_basis"] = (
            f"At least one input rests on a carrier frequency or detection rate "
            f"that this build records as unverified, so the answer is a range "
            f"spanning a factor of {UNCERTAINTY_FACTOR:g} either side rather "
            f"than a single number. A single number here would be false "
            f"precision."
        )
    out["range_as_fraction"] = (
        f"{as_fraction(out['high'])} to {as_fraction(out['low'])}"
        if out["low"] != out["high"] else as_fraction(point)
    )
    out["reason"] = (
        f"Chance per pregnancy of an affected child: {out['range_as_fraction']}."
        if out["low"] != out["high"] else
        f"Chance per pregnancy of an affected child: {out['point_as_fraction']}."
    )
    out["caveats"] = [
        "This is arithmetic on population averages, not a clinical risk "
        "assessment. A genetic counsellor works from family history, ancestry "
        "and clinical-grade testing, all of which move this number.",
        "A residual risk that came from a consumer array is a lower bound, so "
        "any joint risk built on it is also a lower bound.",
    ]
    return out


# ---------------------------------------------------------------------------
# ACMG secondary findings
#
# ACMG SF v3.2, the 2023 revision of the secondary findings list. The list below
# was HAND-ENCODED for this build and was NOT reconciled item by item against
# the published table, which is why ACMG_SF_LIST_VERIFIED is False. A gene being
# absent from this copy must not be read as evidence that ACMG does not list it.
#
# The published v3.2 table carries 81 genes. This copy carries whatever
# len(ACMG_SF_GENES) reports, and the coverage report states both numbers so the
# discrepancy is visible rather than buried.
#
# The important fact about all of them is the same: a consumer SNP array reads
# essentially none of the clinically relevant positions in these genes. The
# coverage report says so in plain words for every single gene rather than
# rendering a reassuring zero.
# ---------------------------------------------------------------------------

ACMG_SF_VERSION = "v3.2 (2023), hand-encoded"
ACMG_SF_LIST_VERIFIED = False

ACMG_SF_GENES: tuple[str, ...] = (
    "ACTA2", "ACTC1", "ACTN2", "ACVRL1", "APC", "APOB", "ATP7B", "BAG3",
    "BMPR1A", "BRCA1", "BRCA2", "BTD", "CACNA1S", "CALM1", "CALM2", "CALM3",
    "CASQ2", "COL3A1", "DES", "DSC2", "DSG2", "DSP", "ENG", "FBN1", "FLNC",
    "GAA", "GLA", "HNF1A", "KCNH2", "KCNQ1", "LDLR", "LMNA", "MAX", "MEN1",
    "MLH1", "MSH2", "MSH6", "MUTYH", "MYBPC3", "MYH7", "MYH11", "MYL2", "MYL3",
    "NF2", "OTC", "PALB2", "PCSK9", "PKP2", "PLN", "PMS2", "PRKAG2", "PTEN",
    "RB1", "RBM20", "RET", "RPE65", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB",
    "SDHC", "SDHD", "SMAD3", "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TMEM127",
    "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TP53", "TPM1", "TRDN", "TSC1", "TSC2",
    "TTN", "TTR", "VHL", "WT1",
)

ACMG_SF_PUBLISHED_COUNT = 81

# The rsIDs a consumer array plausibly carries for an ACMG SF gene. This is a
# very short list on purpose: it is short because the truth is short. Every
# entry is unverified in this build.
ACMG_ARRAY_PROBES: dict[str, tuple[str, ...]] = {
    "BRCA1": ("rs80357713", "rs80357906"),
    "BRCA2": ("rs80359550",),
    "TTR": ("rs76992529",),
    "ATP7B": ("rs76151636",),
}

ACMG_PROBE_NOTE = (
    "These are founder or population-specific variants that some consumer "
    "arrays carry. The three BRCA1 and BRCA2 entries are the Ashkenazi founder "
    "variants. Carrying a probe for three BRCA variants is not BRCA testing: "
    "thousands of pathogenic BRCA variants exist and a negative array result "
    "excludes three of them. All rsIDs here are recorded as unverified in this "
    "build."
)

# Approximate counts of clinically relevant positions, used only as a
# denominator to make the fraction concrete. Every one is unverified and most
# genes have no entry at all, in which case the fraction is reported as None
# rather than as a made-up ratio.
ACMG_KNOWN_RELEVANT_POSITIONS: dict[str, dict] = {
    "BRCA1": _fig(4000, False, "Order of magnitude for BRCA1 variants "
                               "classified pathogenic or likely pathogenic in "
                               "ClinVar. Recalled, not re-verified."),
    "BRCA2": _fig(5000, False, "Order of magnitude for BRCA2 variants "
                               "classified pathogenic or likely pathogenic in "
                               "ClinVar. Recalled, not re-verified."),
    "LDLR": _fig(3000, False, "Order of magnitude for LDLR pathogenic variants. "
                              "Recalled, not re-verified."),
    "ATP7B": _fig(800, False, "Order of magnitude for ATP7B pathogenic variants. "
                              "Recalled, not re-verified."),
}


def acmg_coverage_report(genotypes: dict) -> dict:
    """Report how much of each ACMG secondary-findings gene this file actually reads.

    For nearly every gene the answer is zero, and this function says zero in
    words rather than rendering an empty row that a reader could mistake for a
    clean result. That distinction is the whole point: an ACMG SF gene absent
    from a consumer array has not been screened, and "no secondary findings" from
    array data is a statement about the array, not about the person.
    """
    g = _index(genotypes)
    rows: list[dict] = []
    genes_with_any = 0

    for gene in ACMG_SF_GENES:
        probes = ACMG_ARRAY_PROBES.get(gene, ())
        read = [rsid for rsid in probes if _read(g, rsid) is not None]
        if read:
            genes_with_any += 1
        denominator = ACMG_KNOWN_RELEVANT_POSITIONS.get(gene)
        known = denominator["value"] if denominator else None
        fraction = (len(read) / known) if (known and known > 0) else None

        if not probes:
            assessment = (
                f"Zero. This file reads no position in {gene} that this module "
                f"knows to be clinically relevant, so {gene} was not screened at "
                f"all."
            )
        elif not read:
            assessment = (
                f"Zero. This module knows {len(probes)} {gene} position(s) that "
                f"some arrays carry, and this file reads none of them, so {gene} "
                f"was not screened."
            )
        elif known:
            assessment = (
                f"Effectively zero. This file reads {len(read)} {gene} "
                f"position(s) out of roughly {known} clinically relevant ones, "
                f"which is {fraction:.2%} of the gene's known pathogenic "
                f"variants. A negative result covers only those {len(read)}."
            )
        else:
            assessment = (
                f"Effectively zero. This file reads {len(read)} {gene} "
                f"position(s) out of an unknown but large number of clinically "
                f"relevant ones. A negative result covers only those {len(read)}."
            )

        rows.append({
            "gene": gene,
            "positions_read": len(read),
            "positions_read_rsids": read,
            "probes_known": len(probes),
            "clinically_relevant_positions": known,
            "clinically_relevant_positions_verified": (
                bool(denominator["verified"]) if denominator else False),
            "fraction": round(fraction, 6) if fraction is not None else None,
            "assessment": assessment,
        })

    total = len(ACMG_SF_GENES)
    return {
        "version": ACMG_SF_VERSION,
        "list_verified": ACMG_SF_LIST_VERIFIED,
        "genes_encoded": total,
        "genes_published": ACMG_SF_PUBLISHED_COUNT,
        "list_discrepancy_note": (
            f"The published ACMG SF v3.2 table lists {ACMG_SF_PUBLISHED_COUNT} "
            f"genes. This hand-encoded copy holds {total} and was not reconciled "
            f"item by item, so a gene missing from this list is not evidence "
            f"that ACMG does not list it."
        ),
        "genes": rows,
        "genes_with_any_coverage": genes_with_any,
        "genes_with_zero_coverage": total - genes_with_any,
        "fraction_of_genes_with_any_coverage": round(genes_with_any / total, 6),
        "probe_note": ACMG_PROBE_NOTE,
        "summary": (
            f"This file reads at least one clinically relevant position in "
            f"{genes_with_any} of {total} ACMG secondary-findings genes. For the "
            f"other {total - genes_with_any} the coverage is zero: not low, not "
            f"partial, zero. A consumer array is not a screen for these "
            f"conditions, and finding nothing here means nothing was looked at."
        ),
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------
# Introspection for the documentation of known gaps
# ---------------------------------------------------------------------------

def unverified_figures() -> list[dict]:
    """Every figure and mapping in this module that was not verified at source.

    Generated from the tables so the documented gap list cannot drift away from
    the code, the same arrangement backend/diplotype.unverified_entries uses.
    """
    out: list[dict] = []
    for gene, entry in CARRIER_PANEL.items():
        for variant in entry["tested_variants"]:
            if not variant.get("verified"):
                out.append({"gene": gene, "kind": "variant",
                            "item": f"{variant['name']} ({variant['rsid']})",
                            "note": variant.get("note", "")})
        total = entry["total_known_pathogenic"]
        if not total["verified"]:
            out.append({"gene": gene, "kind": "total_known_pathogenic",
                        "item": str(total["value"]), "note": total["basis"]})
        for pop, fig in entry["carrier_frequency"].items():
            if not fig["verified"]:
                out.append({"gene": gene, "kind": "carrier_frequency",
                            "item": f"{pop} = {as_fraction(fig['value'])}",
                            "note": fig["basis"]})
        for pop, fig in entry["detection_rate"].items():
            if not fig["verified"]:
                out.append({"gene": gene, "kind": "detection_rate",
                            "item": f"{pop} = {fig['value']}",
                            "note": fig["basis"]})
    if not ACMG_SF_LIST_VERIFIED:
        out.append({"gene": "", "kind": "acmg_sf_gene_list",
                    "item": ACMG_SF_VERSION,
                    "note": "Hand-encoded and not reconciled item by item "
                            "against the published table."})
    for gene, probes in ACMG_ARRAY_PROBES.items():
        out.append({"gene": gene, "kind": "acmg_array_probe",
                    "item": ", ".join(probes), "note": ACMG_PROBE_NOTE})
    for gene, fig in ACMG_KNOWN_RELEVANT_POSITIONS.items():
        if not fig["verified"]:
            out.append({"gene": gene, "kind": "acmg_relevant_position_count",
                        "item": str(fig["value"]), "note": fig["basis"]})
    return out
