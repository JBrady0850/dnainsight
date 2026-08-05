"""
diplotype.py -- star-allele diplotype calling and CPIC phenotype translation.

WHY THIS MODULE EXISTS
----------------------
The existing PGx view groups per-SNP findings by drug. That is a real
improvement on nothing, but it is not what a pharmacist reads. The clinical
standard is a DIPLOTYPE: two star alleles, one per chromosome, translated into a
metabolizer phenotype through a published CPIC translation table. A per-SNP list
cannot say "CYP2C19 *2/*17"; only an allele-definition table can.

What the competition does, and why this is worth building anyway:

  SelfDecode  ships PGx but only for its own kit and refuses uploaded files.
  23andMe     holds the only DTC pharmacogenetics authorisation and covers
              three genes.
  Xcode       advertises roughly 500 drug-gene interactions and names no source
              database, which means the claim cannot be checked.

A CPIC-sourced offline report that states its own blind spots out loud is more
useful than any of those, and it is the only one a user can audit.

THE THREE RULES THIS MODULE IS BUILT AROUND
-------------------------------------------
1. REFUSE TO GUESS. ``backend/traits.py`` set the precedent: when the decisive
   position was not read, the answer is "unknown", never a confident call. Here
   the same rule has teeth, because the failure mode is worse. If the array did
   not read the position that defines *2, then *2 is UNTESTABLE, not absent, and
   the diplotype cannot honestly be reported as *1/*1. Every allele carries a
   tri-state: present, absent, untestable.

2. INDETERMINATE IS THE DEFAULT. When coverage is incomplete the phenotype is
   Indeterminate, never Normal. Defaulting to Normal is precisely how these
   tools tell somebody they metabolise a drug fine when nobody actually checked.

3. CYP2D6 IS PROVISIONAL BY CONSTRUCTION. Consumer arrays cannot resolve copy
   number, whole-gene deletions, duplications or CYP2D6-CYP2D7 hybrid alleles.
   Those are not rare curiosities; the *5 deletion and the xN duplications are
   common, and they move the activity score more than any SNP does. So every
   CYP2D6 result carries provisional=True, a plain-English reason, and a capped
   confidence. This is enforced in code, not left to the report writer.

WHAT ``verified`` MEANS IN THE TABLES BELOW
-------------------------------------------
Star-allele definitions are hand-encoded from CPIC allele definition tables,
which are CC0 and therefore safe to ship under this project's MIT licence. The
per-allele ``verified`` flag has one precise meaning and is not a vibe:

    verified True   The rsID-to-star-allele mapping is textbook-level standard
                    AND the variant base recorded here on the GRCh37 plus strand
                    is corroborated by ``data/evidence_overlay.py`` inside this
                    repository, which already carries a reviewed plus-strand risk
                    allele for that rsID.

    verified False  Recalled from the CPIC definition tables but NOT corroborated
                    in-tree at build time. The mapping may be right; it has not
                    been checked here. Every one of these carries a note saying
                    what specifically is unconfirmed.

An unverified allele is still used for calling, because refusing to call it
would silently reduce coverage without telling anybody. Instead it is reported
in the payload under ``unverified_alleles_used`` so the caller can show it. That
is the same trade the frequency module makes for ambiguous strands: use it, flag
it, never hide it.

STRAND CONVENTION
-----------------
Every variant base below is on the GRCh37 PLUS strand, which is what 23andMe and
AncestryDNA report, matching ``data/evidence_overlay.py``. Several of these
genes are transcribed from the minus strand (CYP2D6, TPMT, DPYD, VKORC1,
NUDT15), so the base here deliberately differs from the c. notation used in the
literature. Where that happens it is stated in the allele's note, because a
silently flipped base inverts a diplotype and is worse than no call at all.

OFFLINE CONTRACT
----------------
No network access on any path in this module, at import or at call time.
"""

from __future__ import annotations

from typing import Any, Iterable, Sequence

__all__ = [
    "GENES", "ALLELE_DEFINITIONS", "REFERENCE_ALLELE", "UNTESTABLE_ALLELES",
    "GENE_META", "PROVISIONAL_GENES", "PHENOTYPE_SCALES", "PHENOTYPES",
    "ACTIVITY_BANDS", "DRUG_GENE_PAIRS", "DRUG_ALIASES", "CATEGORY_TEXT",
    "BANNED_PRESCRIPTIVE_PHRASES", "DISCLAIMER", "PRESCRIBER_FRAMING",
    "CONFIDENCE_ORDER", "STATUS_PRESENT", "STATUS_ABSENT", "STATUS_UNTESTABLE",
    "call_diplotype", "translate_phenotype", "prescription_guard",
    "allele_table", "unverified_entries", "known_drugs",
    "contains_prescriptive_language", "audit_language",
]


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------

CONFIDENCE_ORDER: tuple[str, ...] = ("none", "low", "moderate", "high")

STATUS_PRESENT = "present"
STATUS_ABSENT = "absent"
STATUS_UNTESTABLE = "untestable"

# The metabolizer scale. Indeterminate is a first-class member, not an error
# state, because "we could not establish this" is a real and common answer.
PHENOTYPES: tuple[str, ...] = (
    "Poor", "Intermediate", "Normal", "Rapid", "Ultrarapid", "Indeterminate",
)

# Not every CPIC gene is scored on the metabolizer scale. SLCO1B1 is a
# transporter and CPIC reports function, not metabolism. VKORC1 is not an enzyme
# at all: it is the warfarin target, and its phenotype is a dose-sensitivity
# band. Forcing all three onto one vocabulary would be a lie of convenience, so
# each gene declares its scale and the caller is told which one it got.
PHENOTYPE_SCALES: dict[str, tuple[str, ...]] = {
    "metabolizer": PHENOTYPES,
    "transporter_function": (
        "Poor Function", "Decreased Function", "Normal Function",
        "Increased Function", "Indeterminate",
    ),
    "warfarin_sensitivity": (
        "High sensitivity", "Intermediate sensitivity", "Normal sensitivity",
        "Indeterminate",
    ),
}

_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "?", ".", "NN"}
NOCALL_KEY = "NN"


# ---------------------------------------------------------------------------
# Star-allele definitions.
#
# Shape, per allele:
#   variants   rsID -> variant base on the GRCh37 plus strand. Empty for the
#              reference allele.
#   function   CPIC function term, for display.
#   activity   CPIC activity value contributed by ONE copy, or None when the
#              function is uncertain. None propagates to an Indeterminate
#              phenotype rather than being treated as zero.
#   verified   see the module docstring. This is a claim about THIS FILE, not
#              about CPIC.
#   note       why an entry is unverified, or what strand trap applies.
# ---------------------------------------------------------------------------

ALLELE_DEFINITIONS: dict[str, dict[str, dict]] = {

    # -- CYP2C19. Transcribed from the chr10 plus strand, so the c. notation in
    #    the literature and the plus-strand base recorded here agree. That makes
    #    this the least strand-hazardous gene in the table.
    "CYP2C19": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True,
            "note": "Reference allele. Assigned only when the defining variants "
                    "of the other alleles were actually read and found absent.",
        },
        "*2": {
            "variants": {"rs4244285": "A"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.681G>A, the splice defect. The full CPIC definition also "
                    "lists rs12769205 on the same haplotype; rs4244285 alone is "
                    "what consumer arrays carry and is decisive for function.",
        },
        "*3": {
            "variants": {"rs4986893": "A"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.636G>A, W212X. Common in East Asian ancestry, rare "
                    "elsewhere.",
        },
        "*4": {
            "variants": {"rs28399504": "G"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.1A>G, loss of the initiation codon. *4B additionally "
                    "carries rs12248560; this file does not separate *4A from "
                    "*4B and reports both as *4.",
        },
        "*8": {
            "variants": {"rs41291556": "C"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.358T>C, W120R.",
        },
        "*17": {
            "variants": {"rs12248560": "T"}, "function": "increased function",
            "activity": 1.5, "verified": True,
            "note": "-806C>T promoter variant. The only increased-function "
                    "allele in this gene's table.",
        },
        "*6": {
            "variants": {"rs72552267": "A"}, "function": "no function",
            "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: rsID recalled from the CPIC allele definition "
                    "table, not corroborated in-tree. Rare and almost never "
                    "genotyped on a consumer array.",
        },
        "*9": {
            "variants": {"rs17884712": "A"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: rsID recalled, not corroborated in-tree.",
        },
        "*35": {
            "variants": {"rs12769205": "G"}, "function": "no function",
            "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: *35 is reported as rs12769205 without "
                    "rs4244285, which makes it the direct discriminator against "
                    "*2. Neither the rsID nor the plus-strand base was "
                    "corroborated in-tree, so a *35 call from this file should "
                    "be treated as a hypothesis.",
        },
    },

    # -- CYP2C9. chr10 plus strand, same low strand hazard as CYP2C19.
    "CYP2C9": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True, "note": "Reference allele.",
        },
        "*2": {
            "variants": {"rs1799853": "T"}, "function": "decreased function",
            "activity": 0.5, "verified": True,
            "note": "c.430C>T, R144C.",
        },
        "*3": {
            "variants": {"rs1057910": "C"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.1075A>C, I359L. The largest single-allele effect on "
                    "warfarin requirement in this gene.",
        },
        "*5": {
            "variants": {"rs28371686": "G"}, "function": "no function",
            "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: rsID recalled, not corroborated in-tree.",
        },
        "*8": {
            "variants": {"rs7900194": "A"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED, and there is an active conflict worth naming: "
                    "data/evidence_overlay.py files rs28371706 under CYP2C9, "
                    "while rs28371706 is elsewhere reported as the CYP2D6*17 "
                    "defining variant. rs7900194 is recorded here as CYP2C9*8 "
                    "on recall. One of the two attributions is wrong and this "
                    "build did not resolve which.",
        },
        "*11": {
            "variants": {"rs28371685": "T"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: rsID recalled, not corroborated in-tree.",
        },
    },

    # -- VKORC1. Not an enzyme and not a star-allele gene. CPIC treats warfarin
    #    dosing through the single -1639G>A promoter variant, so this table
    #    encodes haplotypes named for that variant rather than inventing star
    #    numbers that CPIC does not use.
    #
    #    VKORC1 is transcribed from the chr16 MINUS strand. The literature says
    #    -1639G>A; the plus-strand base an array reports for the low-dose
    #    haplotype is T, which is what data/evidence_overlay.py already records
    #    for rs9923231. Getting this backwards would invert every warfarin
    #    sensitivity call in the report.
    "VKORC1": {
        "-1639G": {
            "variants": {}, "function": "reference", "activity": 0.0,
            "verified": True,
            "note": "Reference promoter haplotype, normal VKORC1 expression.",
        },
        "-1639A": {
            "variants": {"rs9923231": "T"}, "function": "reduced expression",
            "activity": 1.0, "verified": True,
            "note": "Quoted in the literature as -1639G>A on the minus strand. "
                    "The plus-strand base recorded by an array is T. Activity "
                    "here counts SENSITIVITY copies, not enzyme activity.",
        },
    },

    # -- SLCO1B1. A hepatic uptake transporter, not a metaboliser. CPIC reports
    #    function, so this gene uses the transporter_function scale.
    "SLCO1B1": {
        "*1": {
            "variants": {}, "function": "normal function", "activity": 1.0,
            "verified": True, "note": "Reference allele.",
        },
        "*5": {
            "variants": {"rs4149056": "C"}, "function": "decreased function",
            "activity": 0.5, "verified": True,
            "note": "c.521T>C, V174A. The statin myopathy variant.",
        },
        "*1B": {
            "variants": {"rs2306283": "G"}, "function": "normal function",
            "activity": 1.0, "verified": False,
            "note": "UNVERIFIED: c.388A>G, N130D. rsID and plus-strand base "
                    "recalled, not corroborated in-tree.",
        },
        "*15": {
            "variants": {"rs2306283": "G", "rs4149056": "C"},
            "function": "decreased function", "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: the c.388A>G plus c.521T>C haplotype. Because "
                    "array data is unphased, *15 cannot be distinguished from "
                    "*1B plus *5 in trans when both sites are heterozygous.",
        },
        "*17": {
            "variants": {"rs2306283": "G", "rs4149056": "C", "rs4149015": "A"},
            "function": "decreased function", "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: adds the -11187G>A promoter variant to *15. "
                    "All three rsIDs recalled, none corroborated in-tree.",
        },
    },

    # -- TPMT. Transcribed from the chr6 MINUS strand, so every base below is
    #    the complement of the c. notation in the thiopurine literature.
    "TPMT": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True, "note": "Reference allele.",
        },
        "*3C": {
            "variants": {"rs1142345": "C"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.719A>G, Y240C. Plus-strand base is C because TPMT is on "
                    "the minus strand.",
        },
        "*3B": {
            "variants": {"rs1800460": "T"}, "function": "no function",
            "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: c.460G>A, A154T. The plus-strand base was not "
                    "corroborated in-tree and this gene is minus strand, which "
                    "is exactly where a flip goes unnoticed.",
        },
        "*3A": {
            "variants": {"rs1800460": "T", "rs1142345": "C"},
            "function": "no function", "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: the two-variant haplotype and the most common "
                    "non-functional TPMT allele in European ancestry. Inherits "
                    "the *3B uncertainty. Unphased array data cannot separate "
                    "*3A from *3B plus *3C in trans.",
        },
        "*2": {
            "variants": {"rs1800462": "C"}, "function": "no function",
            "activity": 0.0, "verified": False,
            "note": "UNVERIFIED: c.238G>C, A80P. Plus-strand base not "
                    "corroborated in-tree.",
        },
    },

    # -- NUDT15. chr13 minus strand.
    "NUDT15": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True, "note": "Reference allele.",
        },
        "*3": {
            "variants": {"rs116855232": "T"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.415C>T, R139C. The dominant NUDT15 risk allele in East "
                    "Asian and Hispanic ancestry.",
        },
        "*5": {
            "variants": {"rs186364861": "A"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: c.52G>A, V18I. rsID recalled, not corroborated "
                    "in-tree.",
        },
    },

    # -- DPYD. chr1 minus strand. CPIC names DPYD alleles by HGVS rather than by
    #    star number for most of the actionable set, so both namings appear here
    #    and the HGVS name is authoritative.
    "DPYD": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True, "note": "Reference allele.",
        },
        "*2A": {
            "variants": {"rs3918290": "A"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.1905+1G>A, the IVS14+1 splice variant. A homozygote can "
                    "suffer life-threatening fluoropyrimidine toxicity, which "
                    "is why DPYD carries the highest consequence-per-call ratio "
                    "in this table.",
        },
        "*13": {
            "variants": {"rs55886062": "C"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.1679T>G, I560S.",
        },
        "c.2846A>T": {
            "variants": {"rs67376798": "A"}, "function": "decreased function",
            "activity": 0.5, "verified": True,
            "note": "D949V. Named by HGVS because CPIC does not give it a star "
                    "number.",
        },
        "HapB3": {
            "variants": {"rs56038477": "T"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: HapB3 is defined by the deep intronic "
                    "c.1129-5923C>G, which arrays do not read. rs56038477 "
                    "(c.1236G>A) is used as a linked tag. Both the rsID and the "
                    "plus-strand base are uncorroborated in-tree, and a tag is "
                    "not the variant.",
        },
    },

    # -- UGT1A1. chr2 plus strand. The clinically important allele, *28, is a
    #    promoter TA dinucleotide repeat, not a substitution, so an array cannot
    #    read it directly. rs887829 is the standard linked tag.
    "UGT1A1": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True, "note": "Reference allele, (TA)6 promoter.",
        },
        "*28": {
            "variants": {"rs887829": "T"}, "function": "decreased function",
            "activity": 0.0, "verified": True, "tag_only": True,
            "note": "TAG, NOT THE VARIANT. *28 is the (TA)7 promoter repeat. "
                    "An array cannot count TA repeats, so rs887829 is used as a "
                    "linked tag. Linkage is strong in European and African "
                    "ancestry and weaker elsewhere, so a *28 call from this "
                    "file is an inference about a repeat nobody measured.",
        },
        "*6": {
            "variants": {"rs4148323": "A"}, "function": "decreased function",
            "activity": 0.0, "verified": True,
            "note": "c.211G>A, G71R. Common in East Asian ancestry, where it "
                    "matters more than *28.",
        },
    },

    # -- CYP2D6. chr22 MINUS strand, and the hardest gene here by a wide margin.
    #    Read PROVISIONAL_GENES and the UNTESTABLE_ALLELES entry below before
    #    trusting anything this gene produces.
    "CYP2D6": {
        "*1": {
            "variants": {}, "function": "normal", "activity": 1.0,
            "verified": True,
            "note": "Reference allele. On an array this really means 'none of "
                    "the SNP-defined alleles below were detected', which is a "
                    "much weaker statement than *1 normally implies.",
        },
        "*4": {
            "variants": {"rs3892097": "A"}, "function": "no function",
            "activity": 0.0, "verified": True,
            "note": "c.506-1G>A, the splice defect. The most common "
                    "non-functional CYP2D6 allele in European ancestry.",
        },
        "*10": {
            "variants": {"rs1065852": "A"}, "function": "decreased function",
            "activity": 0.25, "verified": True,
            "note": "c.100C>T, P34S. The full CPIC definition pairs it with "
                    "rs1135840; rs1065852 alone is what array callers use, so "
                    "*10 here may capture some *4 haplotypes that also carry "
                    "rs1065852.",
        },
        "*17": {
            "variants": {"rs28371706": "T"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED and conflicted: data/evidence_overlay.py files "
                    "rs28371706 under CYP2C9 with medicines warfarin and "
                    "phenytoin, while this entry uses it as the CYP2D6*17 "
                    "c.1023C>T defining variant. Both cannot be right. Until "
                    "that is resolved at source, treat any *17 call as a "
                    "hypothesis. *17 matters mainly in African ancestry.",
        },
        "*41": {
            "variants": {"rs28371725": "T"}, "function": "decreased function",
            "activity": 0.5, "verified": False,
            "note": "UNVERIFIED: c.985+39G>A splicing variant. rsID and "
                    "plus-strand base recalled, not corroborated in-tree.",
        },
        "*2": {
            "variants": {"rs16947": "A", "rs1135840": "C"},
            "function": "normal", "activity": 1.0, "verified": False,
            "note": "UNVERIFIED: neither rsID nor plus-strand base corroborated "
                    "in-tree. *2 matters mostly as the backbone of the *2xN "
                    "duplications, which an array cannot see at all.",
        },
    },
}

GENES: tuple[str, ...] = tuple(ALLELE_DEFINITIONS)

# The allele assigned to a chromosome when no defining variant was detected on
# it. Named per gene because VKORC1 does not use star numbers.
REFERENCE_ALLELE: dict[str, str] = {
    "CYP2C19": "*1", "CYP2C9": "*1", "VKORC1": "-1639G", "SLCO1B1": "*1",
    "TPMT": "*1", "NUDT15": "*1", "DPYD": "*1", "UGT1A1": "*1", "CYP2D6": "*1",
}

# Alleles that a consumer array CANNOT test, ever, regardless of which positions
# the chip carries. Structural variants, indels and repeats. These are listed
# separately from ALLELE_DEFINITIONS because they are not "we did not read the
# position" -- they are "this assay cannot see this class of variant". They are
# always reported as untestable and always downgrade confidence. Omitting them
# would make coverage look better than it is, which is the whole failure this
# module exists to prevent.
UNTESTABLE_ALLELES: dict[str, tuple[dict, ...]] = {
    "CYP2D6": (
        {"allele": "*5", "reason": "Whole-gene deletion. A SNP array reads bases, "
                                   "not gene copy number, so a deleted CYP2D6 is "
                                   "invisible and looks like a homozygous call."},
        {"allele": "*1xN / *2xN and other duplications",
         "reason": "Gene duplications raise the activity score and are the usual "
                   "cause of an ultrarapid metabolizer phenotype. Copy number is "
                   "not measurable from array genotypes."},
        {"allele": "*13, *68 and other CYP2D6-CYP2D7 hybrids",
         "reason": "Hybrid alleles formed by recombination with the CYP2D7 "
                   "pseudogene require long-range or long-read assays."},
        {"allele": "*3", "reason": "c.2549delA, a single-base deletion. Arrays "
                                   "report substitutions, not indels."},
        {"allele": "*6", "reason": "c.1707delT, a single-base deletion."},
        {"allele": "*9", "reason": "c.841_843delAAG, an in-frame deletion."},
    ),
    "CYP2C9": (
        {"allele": "*6", "reason": "c.818delA, a single-base deletion that an "
                                   "array cannot call."},
    ),
    "NUDT15": (
        {"allele": "*2", "reason": "Carries the c.36_37insGGAGTC insertion in "
                                   "addition to the *3 variant. The insertion is "
                                   "not array-readable, so *2 and *3 cannot be "
                                   "separated from this file."},
    ),
    "UGT1A1": (
        {"allele": "*36 (TA)5", "reason": "A promoter dinucleotide repeat length. "
                                          "Arrays cannot count repeats."},
        {"allele": "*37 (TA)8", "reason": "A promoter dinucleotide repeat length. "
                                          "Arrays cannot count repeats."},
    ),
    "DPYD": (
        {"allele": "HapB3 (c.1129-5923C>G)",
         "reason": "The causal variant is deep intronic and is not on any "
                   "consumer array. Only a linked tag is available, and a tag is "
                   "not the variant."},
    ),
}

# Which scale each gene is scored on, plus how a summed activity score maps to a
# phenotype. Bands are (lower_inclusive, upper_exclusive, phenotype); None means
# unbounded on that side.
GENE_META: dict[str, dict] = {
    "CYP2C19": {"scale": "metabolizer", "cpic_guideline": True,
                "chromosome": "10", "strand": "+",
                "summary": "Activates clopidogrel and clears several PPIs and SSRIs."},
    "CYP2C9": {"scale": "metabolizer", "cpic_guideline": True,
               "chromosome": "10", "strand": "+",
               "summary": "Clears warfarin, phenytoin and several NSAIDs."},
    "VKORC1": {"scale": "warfarin_sensitivity", "cpic_guideline": True,
               "chromosome": "16", "strand": "-",
               "summary": "The warfarin target. Not a metaboliser; its phenotype "
                          "is a dose-sensitivity band."},
    "SLCO1B1": {"scale": "transporter_function", "cpic_guideline": True,
                "chromosome": "12", "strand": "+",
                "summary": "Hepatic statin uptake transporter."},
    "TPMT": {"scale": "metabolizer", "cpic_guideline": True,
             "chromosome": "6", "strand": "-",
             "summary": "Thiopurine methyltransferase."},
    "NUDT15": {"scale": "metabolizer", "cpic_guideline": True,
               "chromosome": "13", "strand": "-",
               "summary": "Thiopurine nucleotide diphosphatase."},
    "DPYD": {"scale": "metabolizer", "cpic_guideline": True,
             "chromosome": "1", "strand": "-",
             "summary": "Clears fluoropyrimidines."},
    "UGT1A1": {"scale": "metabolizer", "cpic_guideline": True,
               "chromosome": "2", "strand": "+",
               "summary": "Glucuronidation of bilirubin and several drugs."},
    "CYP2D6": {"scale": "metabolizer", "cpic_guideline": True,
               "chromosome": "22", "strand": "-",
               "summary": "Metabolises roughly a quarter of prescribed drugs. "
                          "Array coverage of this gene is structurally incomplete."},
}

ACTIVITY_BANDS: dict[str, tuple[tuple[float | None, float | None, str], ...]] = {
    # CPIC CYP2C19 has no published activity score, so these bands are a faithful
    # re-expression of its published diplotype table: *17/*17 ultrarapid,
    # *1/*17 rapid, *1/*1 normal, one no-function allele intermediate, two poor.
    "CYP2C19": ((3.0, None, "Ultrarapid"), (2.5, 3.0, "Rapid"),
                (2.0, 2.5, "Normal"), (0.0001, 2.0, "Intermediate"),
                (None, 0.0001, "Poor")),
    # CPIC CYP2C9 activity score: 2 normal, 1.5 and 1 intermediate, 0.5 and 0 poor.
    "CYP2C9": ((2.0, None, "Normal"), (1.0, 2.0, "Intermediate"),
               (None, 1.0, "Poor")),
    # CPIC CYP2D6 consensus bands: 0 poor, >0 to <1.25 intermediate,
    # 1.25 to 2.25 normal, above 2.25 ultrarapid.
    "CYP2D6": ((2.2501, None, "Ultrarapid"), (1.25, 2.2501, "Normal"),
               (0.0001, 1.25, "Intermediate"), (None, 0.0001, "Poor")),
    "TPMT": ((2.0, None, "Normal"), (0.0001, 2.0, "Intermediate"),
             (None, 0.0001, "Poor")),
    "NUDT15": ((2.0, None, "Normal"), (0.0001, 2.0, "Intermediate"),
               (None, 0.0001, "Poor")),
    "DPYD": ((2.0, None, "Normal"), (1.0, 2.0, "Intermediate"),
             (None, 1.0, "Poor")),
    "UGT1A1": ((2.0, None, "Normal"), (0.0001, 2.0, "Intermediate"),
               (None, 0.0001, "Poor")),
    "SLCO1B1": ((2.0, None, "Normal Function"), (1.25, 2.0, "Decreased Function"),
                (None, 1.25, "Poor Function")),
    # VKORC1 counts sensitivity copies rather than activity: 2 copies of the
    # -1639A haplotype is the high-sensitivity, low-dose-requirement band.
    "VKORC1": ((2.0, None, "High sensitivity"),
               (1.0, 2.0, "Intermediate sensitivity"),
               (None, 1.0, "Normal sensitivity")),
}

# Genes whose array-derived diplotype is provisional by construction. This is a
# property of the ASSAY, not of the person's data, so it is a constant and not
# something computed per file.
PROVISIONAL_GENES: set[str] = {"CYP2D6"}

PROVISIONAL_REASON: dict[str, str] = {
    "CYP2D6": (
        "This CYP2D6 result is provisional. A consumer DNA array reads single "
        "bases. It cannot count how many copies of CYP2D6 you have, so it cannot "
        "see a deleted gene (*5), a duplicated gene (the xN alleles), or a "
        "CYP2D6-CYP2D7 hybrid. Those are common and they change the answer more "
        "than any single letter does: a duplication can turn a normal result "
        "into an ultrarapid one, and a deletion can turn a normal result into a "
        "poor one. Treat this as a starting point for a conversation, not as a "
        "CYP2D6 test result."
    ),
}


# ---------------------------------------------------------------------------
# Genotype reading. Same tri-state discipline as backend/traits.py: absent from
# the file and no-called are both "we did not read this", and neither is ever
# allowed to mean "reference".
# ---------------------------------------------------------------------------

def _index(genotypes: Any) -> dict:
    """Lower-case the rsID keys of a caller-supplied genotype mapping."""
    if not isinstance(genotypes, dict):
        return {}
    return {str(k).strip().lower(): v for k, v in genotypes.items()}


def _read(genotypes: dict, rsid: str) -> str | None:
    """Canonical sorted genotype key for one rsID, or None when unreadable.

    None covers both "the array does not carry this position" and "the probe
    failed". Both mean the same thing for allele calling: nothing was measured.
    """
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
    """Copies (0, 1 or 2) of one variant base, or None when nothing was read."""
    key = _read(genotypes, rsid)
    if key is None:
        return None
    wanted = str(allele or "").strip().upper()
    if not wanted:
        return None
    return sum(1 for base in key if base == wanted)


def _weaker(first: str, second: str) -> str:
    try:
        return CONFIDENCE_ORDER[min(CONFIDENCE_ORDER.index(first),
                                    CONFIDENCE_ORDER.index(second))]
    except ValueError:
        return "none"


def _downgrade(confidence: str, steps: int = 1) -> str:
    try:
        index = CONFIDENCE_ORDER.index(confidence)
    except ValueError:
        return "none"
    return CONFIDENCE_ORDER[max(0, index - steps)]


def _canonical_gene(gene: Any) -> str:
    """Match a gene name case-insensitively without inventing one."""
    raw = str(gene or "").strip().upper().replace(" ", "")
    for known in ALLELE_DEFINITIONS:
        if known.upper() == raw:
            return known
    return ""


def _allele_sort_key(gene: str, allele: str) -> tuple:
    """Order alleles the way CPIC writes a diplotype: reference first, then by
    star number, then lexically. Keeps *1/*17 rather than *17/*1."""
    ref = REFERENCE_ALLELE.get(gene, "*1")
    if allele == ref:
        return (0, 0, "")
    if allele.startswith("*"):
        digits = ""
        for char in allele[1:]:
            if char.isdigit():
                digits += char
            else:
                break
        if digits:
            return (1, int(digits), allele)
    return (2, 0, allele)


# ---------------------------------------------------------------------------
# Diplotype calling
# ---------------------------------------------------------------------------

def call_diplotype(gene: str, genotypes: dict) -> dict:
    """Call a star-allele diplotype for one gene from array genotypes.

    ``genotypes`` maps rsID to either a two-character genotype string or a
    two-item sequence of alleles, exactly as ``backend/traits.py`` accepts.

    The algorithm is the standard one, with one addition that most
    implementations skip:

      1. For every defined allele, read its defining positions and give the
         allele a TRI-STATE:

             present     every defining variant base is carried at least once
             absent      every defining position was read and at least one
                         defining base is carried zero times
             untestable  at least one defining position was not read at all

         Untestable is the point of this function. An array that does not carry
         rs4244285 has not shown that CYP2C19*2 is absent; it has shown nothing.
         Reporting that as *1/*1 is the single most dangerous thing a consumer
         PGx tool can do, and it is what happens whenever "missing" is quietly
         folded into "reference".

      2. Assign detected alleles to the two chromosome slots MOST SPECIFIC
         FIRST, so an allele defined by three variants is consumed before the
         one-variant allele that is a subset of it. Component variants are
         consumed as they are used, so *15 does not also produce a *5.

      3. Fill any remaining slot with the reference allele, and record
         ``reference_inferred`` when that happened while untestable alleles
         existed. That flag is what makes ``translate_phenotype`` return
         Indeterminate instead of Normal.

    Returns a dict. ``diplotype`` is None when nothing could be called at all.
    """
    canonical = _canonical_gene(gene)
    if not canonical:
        return {
            "gene": str(gene or ""), "known_gene": False, "diplotype": None,
            "alleles": [], "allele_status": {}, "untestable_alleles": [],
            "detected_alleles": [], "positions_read": 0, "positions_expected": 0,
            "coverage": 0.0, "confidence": "none", "provisional": False,
            "reference_inferred": False, "unverified_alleles_used": [],
            "caveats": [f"{gene!r} is not a gene this module has an allele "
                        f"definition table for."],
            "disclaimer": DISCLAIMER,
        }

    defs = ALLELE_DEFINITIONS[canonical]
    ref_name = REFERENCE_ALLELE.get(canonical, "*1")
    g = _index(genotypes)

    allele_status: dict[str, dict] = {}
    copies_at: dict[str, int] = {}
    positions_expected: set[str] = set()
    positions_read: set[str] = set()

    for name, spec in defs.items():
        variants: dict = spec.get("variants") or {}
        if not variants:
            continue                                   # the reference allele
        per_variant: list[int] = []
        missing: list[str] = []
        for rsid, base in variants.items():
            positions_expected.add(rsid)
            copies = _copies(g, rsid, base)
            if copies is None:
                missing.append(rsid)
            else:
                positions_read.add(rsid)
                copies_at[rsid] = copies
                per_variant.append(copies)
        if missing:
            allele_status[name] = {
                "status": STATUS_UNTESTABLE, "copies": None,
                "missing_positions": sorted(missing),
                "verified": bool(spec.get("verified")),
            }
            continue
        # A haplotype requires ALL of its defining variants, so the number of
        # copies it can account for is bounded by its scarcest component.
        copies = min(per_variant) if per_variant else 0
        allele_status[name] = {
            "status": STATUS_PRESENT if copies > 0 else STATUS_ABSENT,
            "copies": copies, "missing_positions": [],
            "verified": bool(spec.get("verified")),
        }

    caveats: list[str] = []

    # Structural alleles this assay can never see. Always reported, never
    # silently dropped.
    structural = [dict(entry) for entry in UNTESTABLE_ALLELES.get(canonical, ())]
    for entry in structural:
        allele_status.setdefault(entry["allele"], {
            "status": STATUS_UNTESTABLE, "copies": None,
            "missing_positions": [], "verified": True,
            "structural": True, "reason": entry["reason"],
        })

    untestable = sorted(n for n, s in allele_status.items()
                        if s["status"] == STATUS_UNTESTABLE)
    detected = [n for n, s in allele_status.items() if s["status"] == STATUS_PRESENT]

    # Most specific first, so a multi-variant haplotype consumes its components.
    order = sorted(detected,
                   key=lambda n: (-len(defs.get(n, {}).get("variants") or {}),
                                  _allele_sort_key(canonical, n)))
    consumed: dict[str, int] = {}
    slots: list[str] = []
    for name in order:
        variants = defs.get(name, {}).get("variants") or {}
        if not variants:
            continue
        available = min(
            (copies_at.get(rsid, 0) - consumed.get(rsid, 0)) for rsid in variants
        )
        take = max(0, min(available, 2 - len(slots)))
        for _ in range(take):
            slots.append(name)
            for rsid in variants:
                consumed[rsid] = consumed.get(rsid, 0) + 1
        if len(slots) >= 2:
            break

    if len(order) > len(slots) and len(slots) >= 2:
        caveats.append(
            "More distinct alleles were detected than two chromosomes can carry. "
            "Array data is unphased, so which variants sit together on the same "
            "chromosome is not knowable from this file."
        )

    reference_inferred = len(slots) < 2
    filled = list(slots)
    while len(filled) < 2:
        filled.append(ref_name)
    filled.sort(key=lambda a: _allele_sort_key(canonical, a))

    tested_any = any(s["status"] in (STATUS_PRESENT, STATUS_ABSENT)
                     for s in allele_status.values())
    if not tested_any:
        diplotype = None
    else:
        diplotype = "/".join(filled)

    # Confidence. Starts high and is only ever pushed down.
    confidence = "high"
    if not tested_any:
        confidence = "none"
    else:
        real_untestable = [n for n in untestable
                           if not allele_status[n].get("structural")]
        if real_untestable and reference_inferred:
            # The worst combination: we are about to call a reference allele on a
            # chromosome whose defining positions were never read.
            confidence = "low"
            caveats.append(
                "At least one chromosome is reported as the reference allele "
                "while " + ", ".join(real_untestable) + " could not be tested at "
                "all. Reference here means 'no listed variant was detected', not "
                "'no variant is present'."
            )
        elif real_untestable:
            confidence = _downgrade(confidence)
            caveats.append(
                "Both chromosomes carry a positively detected allele, but "
                + ", ".join(real_untestable) + " could not be tested, so a rarer "
                "allele cannot be excluded."
            )
        if structural:
            confidence = _downgrade(confidence)
            caveats.append(
                "This gene has alleles that no SNP array can read: "
                + "; ".join(f"{e['allele']} ({e['reason']})" for e in structural)
            )

    provisional = canonical in PROVISIONAL_GENES
    if provisional:
        # A cap, not a downgrade. However complete the SNP coverage looks, a
        # CYP2D6 call from array data cannot be better than low confidence,
        # because the variant class that matters most was never measured.
        confidence = _weaker(confidence, "low")
        caveats.append(PROVISIONAL_REASON[canonical])

    unverified_used = sorted({
        name for name in filled
        if name in defs and not defs[name].get("verified", False)
    })
    if unverified_used:
        confidence = _downgrade(confidence)
        caveats.append(
            "This call uses allele definitions that are not corroborated inside "
            "this repository: " + ", ".join(unverified_used) + ". See the note on "
            "each allele in backend/diplotype.py."
        )

    tag_only = sorted({
        name for name in filled
        if name in defs and defs[name].get("tag_only")
    })
    if tag_only:
        caveats.append(
            "Called through a linked tag SNP rather than the causal variant: "
            + ", ".join(tag_only) + ". Linkage varies by ancestry, so the tag can "
            "be right about the population and wrong about the person."
        )

    coverage = (len(positions_read) / len(positions_expected)
                if positions_expected else 0.0)

    return {
        "gene": canonical,
        "known_gene": True,
        "diplotype": diplotype,
        "alleles": filled if diplotype else [],
        "allele_status": allele_status,
        "detected_alleles": sorted(set(detected)),
        "untestable_alleles": untestable,
        "structural_untestable": structural,
        "reference_inferred": bool(reference_inferred and diplotype is not None),
        "positions_read": len(positions_read),
        "positions_expected": len(positions_expected),
        "coverage": round(coverage, 3),
        "confidence": confidence,
        "provisional": provisional,
        "provisional_reason": PROVISIONAL_REASON.get(canonical, ""),
        "unverified_alleles_used": unverified_used,
        "scale": GENE_META.get(canonical, {}).get("scale", "metabolizer"),
        "caveats": caveats,
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------
# Phenotype translation
# ---------------------------------------------------------------------------

def _band_phenotype(gene: str, score: float) -> str:
    for low, high, phenotype in ACTIVITY_BANDS.get(gene, ()):
        if (low is None or score >= low) and (high is None or score < high):
            return phenotype
    return "Indeterminate"


def _phenotype_label(gene: str, phenotype: str) -> str:
    scale = GENE_META.get(gene, {}).get("scale", "metabolizer")
    if phenotype == "Indeterminate":
        return "Indeterminate"
    if scale == "metabolizer":
        return f"{phenotype} Metabolizer"
    return phenotype


def translate_phenotype(gene: str, diplotype: Any) -> dict:
    """Translate a diplotype into a CPIC activity score and phenotype.

    ``diplotype`` may be the dict returned by :func:`call_diplotype`, which is
    the useful case because it carries coverage information, or a bare string
    such as ``"*1/*2"`` when a caller already knows the diplotype from elsewhere.

    INDETERMINATE IS THE DEFAULT, and that is the entire point of this function.
    It is returned whenever any of the following is true:

      * nothing could be called;
      * a chromosome was filled with the reference allele while some allele of
        that gene was untestable, so "reference" means "we did not look";
      * any allele in the diplotype has an unknown activity value;
      * the allele name is not in this module's table.

    Returning Normal in any of those situations would be a statement nobody
    made. The cost of Indeterminate is a user who has to ask a pharmacist. The
    cost of a wrong Normal is a user who does not.
    """
    canonical = _canonical_gene(gene)
    scale = GENE_META.get(canonical, {}).get("scale", "metabolizer")
    call: dict | None = None
    caveats: list[str] = []

    if isinstance(diplotype, dict):
        call = diplotype
        canonical = _canonical_gene(call.get("gene") or gene) or canonical
        scale = GENE_META.get(canonical, {}).get("scale", scale)
        text = call.get("diplotype")
        caveats = list(call.get("caveats") or [])
    else:
        text = diplotype
        caveats.append(
            "Translated from a diplotype string supplied by the caller. This "
            "function cannot see how much of the gene was actually read, so the "
            "coverage caveats that normally apply are absent, not satisfied."
        )

    base = {
        "gene": canonical or str(gene or ""),
        "diplotype": text if isinstance(text, str) else None,
        "scale": scale,
        "allowed_phenotypes": list(PHENOTYPE_SCALES.get(scale, PHENOTYPES)),
        "activity_score": None,
        "phenotype": "Indeterminate",
        "phenotype_label": "Indeterminate",
        "confidence": (call or {}).get("confidence", "none"),
        "provisional": canonical in PROVISIONAL_GENES,
        "provisional_reason": PROVISIONAL_REASON.get(canonical, ""),
        "reason": "",
        "caveats": caveats,
        "disclaimer": DISCLAIMER,
    }

    if not canonical:
        base["reason"] = (f"{gene!r} is not a gene this module has a translation "
                          f"table for.")
        return base

    if not isinstance(text, str) or not text.strip():
        base["reason"] = (
            "No diplotype could be called, so no phenotype is reported. Nothing "
            "was found and nothing was excluded."
        )
        return base

    parts = [p.strip() for p in text.split("/") if p.strip()]
    if len(parts) != 2:
        base["reason"] = (f"{text!r} is not a two-allele diplotype, so it cannot "
                          f"be translated.")
        return base

    defs = ALLELE_DEFINITIONS[canonical]
    unknown = [p for p in parts if p not in defs]
    if unknown:
        base["reason"] = (
            "Diplotype contains allele(s) with no definition in this module: "
            + ", ".join(unknown) + ". Guessing an activity value for an unknown "
            "allele would invent evidence, so the phenotype stays Indeterminate."
        )
        return base

    activities = [defs[p].get("activity") for p in parts]
    if any(a is None for a in activities):
        base["reason"] = (
            "At least one allele in this diplotype has no assigned activity "
            "value, so the activity score cannot be summed."
        )
        return base

    score = round(float(sum(activities)), 4)
    base["activity_score"] = score

    # The coverage veto. This is where Indeterminate stops being a default and
    # starts being a decision.
    #
    # It fires only when a chromosome was filled in with the REFERENCE allele,
    # because that is the only slot whose content was inferred rather than
    # observed. A diplotype where both chromosomes carry a positively detected
    # allele is not weakened by an untested rare allele elsewhere in the gene.
    #
    # Two things trigger it:
    #   * a defining position of some allele was never read, so "reference"
    #     means "we did not look";
    #   * the gene is provisional, which for CYP2D6 means an inferred reference
    #     chromosome might really be a deletion or a duplication and the array
    #     has no way to tell.
    # A rare indel allele that no array can ever read is NOT enough on its own.
    # If it were, Normal would be permanently unreachable for half these genes
    # and the word Indeterminate would stop carrying information.
    if call is not None and call.get("reference_inferred"):
        real_untestable = [
            name for name in (call.get("untestable_alleles") or [])
            if not (call.get("allele_status") or {}).get(name, {}).get("structural")
        ]
        if real_untestable or canonical in PROVISIONAL_GENES:
            why = (
                "some defining positions of this gene were never read from this "
                "file (" + ", ".join(real_untestable) + " could not be tested)"
                if real_untestable else
                "this gene's copy number, deletions and duplications cannot be "
                "measured from array data at all"
            )
            base["reason"] = (
                "A chromosome was reported as the reference allele, and " + why +
                ". An inferred reference chromosome is not evidence of normal "
                "function, so the phenotype is Indeterminate rather than Normal."
            )
            base["confidence"] = _weaker(base["confidence"], "low")
            return base

    phenotype = _band_phenotype(canonical, score)
    allowed = PHENOTYPE_SCALES.get(scale, PHENOTYPES)
    if phenotype not in allowed:
        phenotype = "Indeterminate"
    base["phenotype"] = phenotype
    base["phenotype_label"] = _phenotype_label(canonical, phenotype)
    base["reason"] = (
        f"Activity score {score} from {text} places this in the {phenotype} band "
        f"of the CPIC {scale.replace('_', ' ')} scale."
    )
    if base["provisional"]:
        base["confidence"] = _weaker(base["confidence"], "low")
    return base


# ---------------------------------------------------------------------------
# Prescription guard
#
# Everything below produces text that a user will read next to the name of a
# medicine they are taking. It must never read as an instruction. The
# categories are descriptions of what CPIC DOCUMENTS; what to do about it is a
# decision for a prescriber who can see the whole person.
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "DNAInsight is not a medical device and this page is not medical advice. "
    "It reports what published CPIC guidance documents for a gene and a "
    "medicine at a given phenotype, and nothing more. Consumer DNA arrays do "
    "not read every position that defines a star allele, so a result that looks "
    "normal here can simply be a position that was never read. Nothing on this "
    "page is an instruction about any medicine. Any decision to begin, continue "
    "or change a medicine belongs to you together with a licensed prescriber or "
    "pharmacist, working from clinically validated testing. Bring this page to "
    "that conversation."
)

PRESCRIBER_FRAMING = (
    "Discuss this with your prescriber or pharmacist. This is information for "
    "that conversation, not a decision about your medicines."
)

# Neutral, non-imperative descriptions. Each says what CPIC RECORDS. None of
# them says what anybody should do.
CATEGORY_TEXT: dict[str, str] = {
    "standard_prescribing": (
        "CPIC records no genotype-based change for this gene and medicine at "
        "this phenotype."
    ),
    "altered_exposure": (
        "CPIC documents altered drug exposure for this gene and medicine at "
        "this phenotype."
    ),
    "reduced_activation": (
        "CPIC documents reduced conversion of this medicine to its active form "
        "at this phenotype."
    ),
    "increased_activation": (
        "CPIC documents increased conversion of this medicine to its active "
        "form at this phenotype, and a correspondingly higher exposure to the "
        "active drug."
    ),
    "toxicity_risk_documented": (
        "CPIC documents an increased risk of drug toxicity for this gene and "
        "medicine at this phenotype."
    ),
    "efficacy_risk_documented": (
        "CPIC documents a risk of reduced therapeutic effect for this gene and "
        "medicine at this phenotype."
    ),
    "insufficient_evidence": (
        "The phenotype could not be established from this file, so CPIC's "
        "genotype-based guidance cannot be applied to this pair. This is an "
        "absence of information, not a reassuring result."
    ),
}

# Phrases that must never appear anywhere in a prescription_guard payload.
# Checked by a test, and by audit_language() here, because a guideline nobody
# enforces is a comment.
#
# Every entry is lower case and is matched as a substring. "avoid " carries a
# trailing space so it cannot fire on the word "avoidance" in a citation title.
BANNED_PRESCRIPTIVE_PHRASES: tuple[str, ...] = (
    "stop taking", "start taking", "do not take", "should take", "must take",
    "you should stop", "you should start", "discontinue", "switch to",
    "avoid ", "should avoid", "do not use", "use instead", "take instead",
    "reduce the dose", "reduce your dose", "increase the dose",
    "increase your dose", "lower the dose", "raise the dose",
    "adjust the dose", "adjust your dose", "dose reduction", "dose increase",
    "starting dose", "recommended dose", "standard dose", "half the dose",
    "double the dose", "mg/kg", "mg per", "twice daily", "once daily",
    "three times a day", "contraindicated", "we recommend", "i recommend",
    "you must", "you need to", "prescribe ", "alternative therapy",
    "alternative agent", "consider an alternative", "use an alternative",
)


def contains_prescriptive_language(text: Any) -> list[str]:
    """Return every banned phrase found in ``text``. Empty means clean."""
    lowered = str(text or "").lower()
    return [phrase for phrase in BANNED_PRESCRIPTIVE_PHRASES if phrase in lowered]


def audit_language(payload: Any) -> list[dict]:
    """Walk any nested payload and report strings carrying imperative language.

    Used by the test suite. Kept in the module rather than in the test so the
    banned list and the walker cannot drift apart.
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
            found = contains_prescriptive_language(node)
            if found:
                hits.append({"path": path, "phrases": found, "text": node})

    walk(payload, "")
    return hits


# Drug-gene pairs. ``effects`` maps a phenotype on that gene's scale to a
# category key from CATEGORY_TEXT. Phenotypes not listed fall through to
# "insufficient_evidence", which is the safe direction.
#
# ``cpic_level`` uses the same vocabulary as backend/scoring.CPIC_LEVELS. An
# empty string means CPIC publishes no guideline for that pair, which is stated
# rather than hidden behind a plausible-looking letter.
DRUG_GENE_PAIRS: tuple[dict, ...] = (
    {"drug": "clopidogrel", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC clopidogrel and CYP2C19",
     "mechanism": "Clopidogrel is a prodrug that CYP2C19 converts to its active "
                  "antiplatelet form.",
     "effects": {"Poor": "reduced_activation", "Intermediate": "reduced_activation",
                 "Normal": "standard_prescribing", "Rapid": "standard_prescribing",
                 "Ultrarapid": "standard_prescribing"}},
    {"drug": "omeprazole", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC proton pump inhibitors and CYP2C19",
     "mechanism": "CYP2C19 clears omeprazole.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing", "Rapid": "efficacy_risk_documented",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "pantoprazole", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC proton pump inhibitors and CYP2C19",
     "mechanism": "CYP2C19 clears pantoprazole.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing", "Rapid": "efficacy_risk_documented",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "escitalopram", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC SSRIs and CYP2C19",
     "mechanism": "CYP2C19 is the main clearance route for escitalopram.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing", "Rapid": "efficacy_risk_documented",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "citalopram", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC SSRIs and CYP2C19",
     "mechanism": "CYP2C19 is the main clearance route for citalopram.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing", "Rapid": "efficacy_risk_documented",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "voriconazole", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC voriconazole and CYP2C19",
     "mechanism": "CYP2C19 clears voriconazole, and its exposure window is narrow.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Rapid": "efficacy_risk_documented",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "warfarin", "gene": "CYP2C9", "cpic_level": "A",
     "guideline": "CPIC warfarin, CYP2C9 and VKORC1",
     "mechanism": "CYP2C9 clears the more active S-enantiomer of warfarin.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "warfarin", "gene": "VKORC1", "cpic_level": "A",
     "guideline": "CPIC warfarin, CYP2C9 and VKORC1",
     "mechanism": "VKORC1 is the enzyme warfarin inhibits. Promoter genotype "
                  "shifts how much warfarin is needed to reach the same effect.",
     "effects": {"High sensitivity": "altered_exposure",
                 "Intermediate sensitivity": "altered_exposure",
                 "Normal sensitivity": "standard_prescribing"}},
    {"drug": "phenytoin", "gene": "CYP2C9", "cpic_level": "A",
     "guideline": "CPIC phenytoin, CYP2C9 and HLA-B",
     "mechanism": "CYP2C9 clears phenytoin, which has a narrow therapeutic window.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "celecoxib", "gene": "CYP2C9", "cpic_level": "A",
     "guideline": "CPIC NSAIDs and CYP2C9",
     "mechanism": "CYP2C9 clears celecoxib.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "ibuprofen", "gene": "CYP2C9", "cpic_level": "A",
     "guideline": "CPIC NSAIDs and CYP2C9",
     "mechanism": "CYP2C9 clears ibuprofen.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "simvastatin", "gene": "SLCO1B1", "cpic_level": "A",
     "guideline": "CPIC statins and SLCO1B1, ABCG2 and CYP2C9",
     "mechanism": "SLCO1B1 carries simvastatin acid into the liver. Reduced "
                  "transport leaves more in the circulation, where muscle "
                  "toxicity happens.",
     "effects": {"Poor Function": "toxicity_risk_documented",
                 "Decreased Function": "toxicity_risk_documented",
                 "Normal Function": "standard_prescribing",
                 "Increased Function": "standard_prescribing"}},
    {"drug": "atorvastatin", "gene": "SLCO1B1", "cpic_level": "A",
     "guideline": "CPIC statins and SLCO1B1, ABCG2 and CYP2C9",
     "mechanism": "SLCO1B1 carries atorvastatin into the liver.",
     "effects": {"Poor Function": "toxicity_risk_documented",
                 "Decreased Function": "altered_exposure",
                 "Normal Function": "standard_prescribing",
                 "Increased Function": "standard_prescribing"}},
    {"drug": "rosuvastatin", "gene": "SLCO1B1", "cpic_level": "A",
     "guideline": "CPIC statins and SLCO1B1, ABCG2 and CYP2C9",
     "mechanism": "SLCO1B1 carries rosuvastatin into the liver.",
     "effects": {"Poor Function": "toxicity_risk_documented",
                 "Decreased Function": "altered_exposure",
                 "Normal Function": "standard_prescribing",
                 "Increased Function": "standard_prescribing"}},
    {"drug": "azathioprine", "gene": "TPMT", "cpic_level": "A",
     "guideline": "CPIC thiopurines, TPMT and NUDT15",
     "mechanism": "TPMT inactivates thiopurine metabolites. Less TPMT activity "
                  "means more active metabolite in the marrow.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "azathioprine", "gene": "NUDT15", "cpic_level": "A",
     "guideline": "CPIC thiopurines, TPMT and NUDT15",
     "mechanism": "NUDT15 removes active thiopurine nucleotides before they are "
                  "incorporated into DNA.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "mercaptopurine", "gene": "TPMT", "cpic_level": "A",
     "guideline": "CPIC thiopurines, TPMT and NUDT15",
     "mechanism": "TPMT inactivates mercaptopurine metabolites.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "mercaptopurine", "gene": "NUDT15", "cpic_level": "A",
     "guideline": "CPIC thiopurines, TPMT and NUDT15",
     "mechanism": "NUDT15 removes active thiopurine nucleotides.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "fluorouracil", "gene": "DPYD", "cpic_level": "A",
     "guideline": "CPIC fluoropyrimidines and DPYD",
     "mechanism": "DPYD is the rate-limiting step in clearing fluoropyrimidines. "
                  "Very low activity has caused fatal toxicity.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "capecitabine", "gene": "DPYD", "cpic_level": "A",
     "guideline": "CPIC fluoropyrimidines and DPYD",
     "mechanism": "Capecitabine is converted to fluorouracil, which DPYD clears.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "toxicity_risk_documented",
                 "Normal": "standard_prescribing"}},
    {"drug": "atazanavir", "gene": "UGT1A1", "cpic_level": "A",
     "guideline": "CPIC atazanavir and UGT1A1",
     "mechanism": "Atazanavir inhibits UGT1A1. Low baseline UGT1A1 activity is "
                  "associated with visible jaundice.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "irinotecan", "gene": "UGT1A1", "cpic_level": "",
     "guideline": "No CPIC guideline. The FDA label carries UGT1A1 information "
                  "and DPWG publishes separate guidance.",
     "mechanism": "UGT1A1 glucuronidates SN-38, the active irinotecan metabolite.",
     "effects": {"Poor": "toxicity_risk_documented",
                 "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing"}},
    {"drug": "codeine", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC codeine and tramadol with CYP2D6",
     "mechanism": "Codeine is a prodrug that CYP2D6 converts to morphine.",
     "effects": {"Poor": "reduced_activation", "Intermediate": "reduced_activation",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "increased_activation"}},
    {"drug": "tramadol", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC codeine and tramadol with CYP2D6",
     "mechanism": "Tramadol is a prodrug that CYP2D6 converts to its active "
                  "opioid metabolite.",
     "effects": {"Poor": "reduced_activation", "Intermediate": "reduced_activation",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "increased_activation"}},
    {"drug": "tamoxifen", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC tamoxifen and CYP2D6",
     "mechanism": "CYP2D6 converts tamoxifen to endoxifen, the active form.",
     "effects": {"Poor": "reduced_activation", "Intermediate": "reduced_activation",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "standard_prescribing"}},
    {"drug": "nortriptyline", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC tricyclic antidepressants and CYP2D6 and CYP2C19",
     "mechanism": "CYP2D6 clears nortriptyline.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "amitriptyline", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC tricyclic antidepressants and CYP2D6 and CYP2C19",
     "mechanism": "CYP2D6 clears amitriptyline and its active metabolite.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "amitriptyline", "gene": "CYP2C19", "cpic_level": "A",
     "guideline": "CPIC tricyclic antidepressants and CYP2D6 and CYP2C19",
     "mechanism": "CYP2C19 demethylates amitriptyline to nortriptyline.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Rapid": "altered_exposure", "Ultrarapid": "altered_exposure"}},
    {"drug": "paroxetine", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC SSRIs and CYP2D6",
     "mechanism": "CYP2D6 clears paroxetine.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "ondansetron", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC ondansetron and tropisetron with CYP2D6",
     "mechanism": "CYP2D6 clears ondansetron.",
     "effects": {"Poor": "standard_prescribing",
                 "Intermediate": "standard_prescribing",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "atomoxetine", "gene": "CYP2D6", "cpic_level": "A",
     "guideline": "CPIC atomoxetine and CYP2D6",
     "mechanism": "CYP2D6 clears atomoxetine.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
    {"drug": "metoprolol", "gene": "CYP2D6", "cpic_level": "",
     "guideline": "No CPIC guideline. DPWG publishes guidance for this pair.",
     "mechanism": "CYP2D6 clears metoprolol.",
     "effects": {"Poor": "altered_exposure", "Intermediate": "altered_exposure",
                 "Normal": "standard_prescribing",
                 "Ultrarapid": "efficacy_risk_documented"}},
)

# A deliberately small brand-name map. Brand names are regional and a long list
# would rot; these are the handful a user is most likely to type.
DRUG_ALIASES: dict[str, str] = {
    "plavix": "clopidogrel",
    "coumadin": "warfarin",
    "jantoven": "warfarin",
    "zocor": "simvastatin",
    "lipitor": "atorvastatin",
    "crestor": "rosuvastatin",
    "prilosec": "omeprazole",
    "protonix": "pantoprazole",
    "lexapro": "escitalopram",
    "celexa": "citalopram",
    "vfend": "voriconazole",
    "dilantin": "phenytoin",
    "celebrex": "celecoxib",
    "imuran": "azathioprine",
    "purinethol": "mercaptopurine",
    "6-mp": "mercaptopurine",
    "5-fu": "fluorouracil",
    "xeloda": "capecitabine",
    "reyataz": "atazanavir",
    "camptosar": "irinotecan",
    "ultram": "tramadol",
    "nolvadex": "tamoxifen",
    "pamelor": "nortriptyline",
    "elavil": "amitriptyline",
    "paxil": "paroxetine",
    "zofran": "ondansetron",
    "strattera": "atomoxetine",
    "lopressor": "metoprolol",
    "toprol": "metoprolol",
}


def known_drugs() -> list[str]:
    """Every medicine this module has a gene pairing for, sorted."""
    return sorted({pair["drug"] for pair in DRUG_GENE_PAIRS})


def _normalise_drug(name: Any) -> str:
    raw = str(name or "").strip().lower()
    if not raw:
        return ""
    raw = raw.split("(")[0].strip()
    # Strip a trailing strength or form so "simvastatin 20mg tablets" matches.
    tokens = [t for t in raw.replace(",", " ").split() if t]
    for token in tokens:
        cleaned = token.strip(".").strip()
        if cleaned in DRUG_ALIASES:
            return DRUG_ALIASES[cleaned]
        if cleaned in {p["drug"] for p in DRUG_GENE_PAIRS}:
            return cleaned
    if raw in DRUG_ALIASES:
        return DRUG_ALIASES[raw]
    return raw


def _phenotype_for(gene: str, value: Any) -> dict:
    """Resolve whatever the caller passed for a gene into a phenotype block."""
    if isinstance(value, dict) and "phenotype" in value and "activity_score" in value:
        return value                                   # already translated
    if isinstance(value, dict):
        return translate_phenotype(gene, value)        # a call_diplotype payload
    return translate_phenotype(gene, value)            # a diplotype string


def prescription_guard(medications: Iterable[Any],
                       diplotypes: dict[str, Any]) -> dict:
    """Return only the gene-drug pairs that apply to the medicines supplied.

    ``medications`` is what the user says they are currently taking: a list of
    strings, or dicts carrying a ``name`` key. ``diplotypes`` maps gene name to
    either a :func:`call_diplotype` payload, an already-translated
    :func:`translate_phenotype` payload, or a bare diplotype string.

    The output is deliberately narrow. A user taking two medicines gets the
    pairs for those two medicines, not a 500-row interaction table they have to
    search. Every entry describes what CPIC DOCUMENTS and frames the next step
    as a conversation with a prescriber. No string in the payload is an
    instruction to begin, continue or change a medicine, and
    :func:`audit_language` exists so that claim can be tested rather than
    believed.
    """
    requested: list[dict] = []
    for item in (medications or []):
        raw = item.get("name") if isinstance(item, dict) else item
        text = str(raw or "").strip()
        if not text:
            continue
        requested.append({"as_entered": text, "normalised": _normalise_drug(text)})

    resolved: dict[str, dict] = {}
    for gene, value in (diplotypes or {}).items():
        canonical = _canonical_gene(gene)
        if not canonical:
            continue
        resolved[canonical] = _phenotype_for(canonical, value)

    matches: list[dict] = []
    matched_names: set[str] = set()
    for entry in requested:
        for pair in DRUG_GENE_PAIRS:
            if pair["drug"] != entry["normalised"]:
                continue
            gene = pair["gene"]
            if gene not in resolved:
                # No diplotype for this gene, so there is nothing to say about
                # this pair. Silence beats a reassuring blank row.
                continue
            matched_names.add(entry["as_entered"])
            block = resolved[gene]
            phenotype = block.get("phenotype") or "Indeterminate"
            category = pair["effects"].get(phenotype, "insufficient_evidence")
            if phenotype == "Indeterminate":
                category = "insufficient_evidence"
            matches.append({
                "drug": pair["drug"],
                "drug_as_entered": entry["as_entered"],
                "gene": gene,
                "diplotype": block.get("diplotype"),
                "phenotype": phenotype,
                "phenotype_label": block.get("phenotype_label", phenotype),
                "activity_score": block.get("activity_score"),
                "scale": block.get("scale", "metabolizer"),
                "cpic_level": pair["cpic_level"],
                "cpic_level_note": (
                    "CPIC publishes no guideline for this pair."
                    if not pair["cpic_level"] else ""
                ),
                "guideline": pair["guideline"],
                "mechanism": pair["mechanism"],
                "category": category,
                "description": CATEGORY_TEXT[category],
                "next_step": PRESCRIBER_FRAMING,
                "confidence": block.get("confidence", "none"),
                "provisional": bool(block.get("provisional")),
                "provisional_reason": block.get("provisional_reason", ""),
                "disclaimer": DISCLAIMER,
            })

    # Sort by how much a prescriber would want to see it first: CPIC A pairs,
    # then documented risk categories, then everything else.
    risk_order = {
        "toxicity_risk_documented": 0, "reduced_activation": 1,
        "increased_activation": 1, "efficacy_risk_documented": 2,
        "altered_exposure": 3, "insufficient_evidence": 4,
        "standard_prescribing": 5,
    }
    matches.sort(key=lambda m: (m["cpic_level"] != "A",
                                risk_order.get(m["category"], 9),
                                m["drug"], m["gene"]))

    unmatched = [e["as_entered"] for e in requested
                 if e["as_entered"] not in matched_names]

    return {
        "medications_supplied": [e["as_entered"] for e in requested],
        "matches": matches,
        "count": len(matches),
        "unmatched_medications": unmatched,
        "unmatched_note": (
            "These medicines have no gene pairing in this module, or the gene "
            "they pair with was not called from this file. That is not a "
            "statement that they are unaffected by your genotype."
        ) if unmatched else "",
        "genes_considered": sorted(resolved),
        "framing": PRESCRIBER_FRAMING,
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------
# Introspection, for the report and for the documentation of known gaps
# ---------------------------------------------------------------------------

def allele_table(gene: str) -> list[dict]:
    """Flat, serialisable view of one gene's allele definitions."""
    canonical = _canonical_gene(gene)
    if not canonical:
        return []
    rows: list[dict] = []
    for name, spec in ALLELE_DEFINITIONS[canonical].items():
        rows.append({
            "gene": canonical, "allele": name,
            "variants": dict(spec.get("variants") or {}),
            "function": spec.get("function", ""),
            "activity": spec.get("activity"),
            "verified": bool(spec.get("verified")),
            "tag_only": bool(spec.get("tag_only")),
            "note": spec.get("note", ""),
            "array_testable": True,
        })
    for entry in UNTESTABLE_ALLELES.get(canonical, ()):
        rows.append({
            "gene": canonical, "allele": entry["allele"], "variants": {},
            "function": "", "activity": None, "verified": True,
            "tag_only": False, "note": entry["reason"], "array_testable": False,
        })
    rows.sort(key=lambda r: _allele_sort_key(canonical, r["allele"]))
    return rows


def unverified_entries() -> list[dict]:
    """Every allele definition this build could NOT corroborate in-tree.

    This is the list that belongs in the documentation as known gaps. Producing
    it from the table rather than from a hand-maintained prose list means it
    cannot silently fall out of date.
    """
    out: list[dict] = []
    for gene, defs in ALLELE_DEFINITIONS.items():
        for name, spec in defs.items():
            if spec.get("verified"):
                continue
            out.append({
                "gene": gene, "allele": name,
                "variants": dict(spec.get("variants") or {}),
                "note": spec.get("note", ""),
            })
    out.sort(key=lambda r: (r["gene"], _allele_sort_key(r["gene"], r["allele"])))
    return out
