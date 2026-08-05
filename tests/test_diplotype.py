"""Tests for backend.diplotype star-allele calling and CPIC phenotype translation.

The behaviours under test here are safety behaviours, not features. Every one of
them exists because the alternative is a consumer PGx tool telling somebody they
metabolise a drug normally when nobody checked.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.diplotype import (
    ACTIVITY_BANDS,
    ALLELE_DEFINITIONS,
    BANNED_PRESCRIPTIVE_PHRASES,
    CONFIDENCE_ORDER,
    DISCLAIMER,
    DRUG_ALIASES,
    DRUG_GENE_PAIRS,
    GENE_META,
    GENES,
    PHENOTYPES,
    PHENOTYPE_SCALES,
    PROVISIONAL_GENES,
    REFERENCE_ALLELE,
    STATUS_ABSENT,
    STATUS_PRESENT,
    STATUS_UNTESTABLE,
    UNTESTABLE_ALLELES,
    allele_table,
    audit_language,
    call_diplotype,
    contains_prescriptive_language,
    known_drugs,
    prescription_guard,
    translate_phenotype,
    unverified_entries,
)

# The em dash character itself, written as an escape so this file
# contains no literal one. House style forbids them in source.
EM_DASH = "\u2014"

MODULE_PATH = Path(__file__).parent.parent / "backend" / "diplotype.py"


# ---------------------------------------------------------------------------
# Synthetic genotypes.
#
# Each "clear" mapping reads every defining position of the gene and carries no
# variant base, which is the only state in which a reference call is honest.
# ---------------------------------------------------------------------------

CYP2C19_CLEAR = {
    "rs4244285": "GG", "rs4986893": "GG", "rs28399504": "AA",
    "rs41291556": "TT", "rs12248560": "CC", "rs72552267": "GG",
    "rs17884712": "GG", "rs12769205": "AA",
}
CYP2C9_CLEAR = {
    "rs1799853": "CC", "rs1057910": "AA", "rs28371686": "CC",
    "rs7900194": "GG", "rs28371685": "CC",
}
CYP2D6_CLEAR = {
    "rs3892097": "GG", "rs1065852": "GG", "rs28371706": "CC",
    "rs28371725": "CC", "rs16947": "GG", "rs1135840": "GG",
}
SLCO1B1_CLEAR = {"rs4149056": "TT", "rs2306283": "AA", "rs4149015": "GG"}
TPMT_CLEAR = {"rs1142345": "TT", "rs1800460": "CC", "rs1800462": "GG"}
NUDT15_CLEAR = {"rs116855232": "CC", "rs186364861": "GG"}
DPYD_CLEAR = {"rs3918290": "GG", "rs55886062": "AA", "rs67376798": "GG",
              "rs56038477": "CC"}
UGT1A1_CLEAR = {"rs887829": "CC", "rs4148323": "GG"}
VKORC1_CLEAR = {"rs9923231": "CC"}

CLEAR_BY_GENE = {
    "CYP2C19": CYP2C19_CLEAR, "CYP2C9": CYP2C9_CLEAR, "CYP2D6": CYP2D6_CLEAR,
    "SLCO1B1": SLCO1B1_CLEAR, "TPMT": TPMT_CLEAR, "NUDT15": NUDT15_CLEAR,
    "DPYD": DPYD_CLEAR, "UGT1A1": UGT1A1_CLEAR, "VKORC1": VKORC1_CLEAR,
}


def with_variant(clear, **overrides):
    """A clear background with specific positions carrying the variant base."""
    merged = dict(clear)
    merged.update(overrides)
    return merged


# ---------------------------------------------------------------------------
# Table integrity
# ---------------------------------------------------------------------------

def test_every_required_gene_is_present():
    for gene in ("CYP2C19", "CYP2C9", "VKORC1", "SLCO1B1", "TPMT", "NUDT15",
                 "DPYD", "UGT1A1", "CYP2D6"):
        assert gene in ALLELE_DEFINITIONS
        assert gene in GENES


def test_every_gene_has_a_reference_allele_that_exists_in_its_table():
    for gene in GENES:
        ref = REFERENCE_ALLELE[gene]
        assert ref in ALLELE_DEFINITIONS[gene]
        assert ALLELE_DEFINITIONS[gene][ref]["variants"] == {}


def test_every_allele_definition_carries_the_required_fields():
    for gene, defs in ALLELE_DEFINITIONS.items():
        for name, spec in defs.items():
            assert isinstance(spec.get("variants"), dict), (gene, name)
            assert isinstance(spec.get("verified"), bool), (gene, name)
            assert spec.get("note"), (gene, name)
            assert "activity" in spec, (gene, name)


def test_unverified_alleles_say_so_in_their_note():
    """An unverified entry that does not admit it is worse than no entry."""
    for gene, defs in ALLELE_DEFINITIONS.items():
        for name, spec in defs.items():
            if not spec["verified"]:
                assert "UNVERIFIED" in spec["note"], (gene, name)


def test_defining_variants_use_rsids_and_single_plus_strand_bases():
    for gene, defs in ALLELE_DEFINITIONS.items():
        for name, spec in defs.items():
            for rsid, base in spec["variants"].items():
                assert rsid.startswith("rs") and rsid[2:].isdigit(), (gene, name, rsid)
                assert base in ("A", "C", "G", "T"), (gene, name, rsid, base)


def test_every_gene_declares_a_scale_and_activity_bands():
    for gene in GENES:
        scale = GENE_META[gene]["scale"]
        assert scale in PHENOTYPE_SCALES
        assert gene in ACTIVITY_BANDS
        for _low, _high, phenotype in ACTIVITY_BANDS[gene]:
            assert phenotype in PHENOTYPE_SCALES[scale], (gene, phenotype)


def test_unverified_entries_are_reported_and_each_carries_a_note():
    rows = unverified_entries()
    assert rows, "the honest gap list should not be empty for this build"
    for row in rows:
        assert row["gene"] in ALLELE_DEFINITIONS
        assert "UNVERIFIED" in row["note"]


def test_allele_table_includes_structurally_untestable_alleles():
    rows = allele_table("CYP2D6")
    untestable = [r for r in rows if not r["array_testable"]]
    assert untestable, "CYP2D6 must expose its structural blind spots"
    assert any("*5" in r["allele"] for r in untestable)


def test_module_source_contains_no_em_dash():
    text = MODULE_PATH.read_text(encoding="utf-8")
    assert EM_DASH not in text


# ---------------------------------------------------------------------------
# Diplotype calling
# ---------------------------------------------------------------------------

def test_unknown_gene_is_refused_rather_than_guessed():
    call = call_diplotype("NOTAGENE", {})
    assert call["known_gene"] is False
    assert call["diplotype"] is None
    assert call["confidence"] == "none"


def test_heterozygous_star_two_gives_one_over_two():
    call = call_diplotype("CYP2C19", with_variant(CYP2C19_CLEAR, rs4244285="AG"))
    assert call["diplotype"] == "*1/*2"
    assert call["allele_status"]["*2"]["status"] == STATUS_PRESENT
    assert call["allele_status"]["*2"]["copies"] == 1


def test_homozygous_star_two_gives_two_over_two_and_is_not_reference_inferred():
    call = call_diplotype("CYP2C19", with_variant(CYP2C19_CLEAR, rs4244285="AA"))
    assert call["diplotype"] == "*2/*2"
    assert call["reference_inferred"] is False


def test_positions_read_and_found_absent_are_reported_absent():
    call = call_diplotype("CYP2C19", CYP2C19_CLEAR)
    assert call["allele_status"]["*2"]["status"] == STATUS_ABSENT
    assert call["untestable_alleles"] == []
    assert call["diplotype"] == "*1/*1"


def test_a_position_the_array_never_read_is_untestable_not_absent():
    """The tri-state. This is the whole point of the module."""
    partial = {k: v for k, v in CYP2C19_CLEAR.items() if k != "rs4244285"}
    call = call_diplotype("CYP2C19", partial)
    assert call["allele_status"]["*2"]["status"] == STATUS_UNTESTABLE
    assert call["allele_status"]["*2"]["copies"] is None
    assert "*2" in call["untestable_alleles"]
    assert call["allele_status"]["*2"]["missing_positions"] == ["rs4244285"]


def test_a_no_call_probe_is_untestable_not_absent():
    """A failed probe must never read as a reference allele."""
    for token in ("--", "NN", "", "00"):
        call = call_diplotype("CYP2C19",
                              with_variant(CYP2C19_CLEAR, rs4244285=token))
        assert call["allele_status"]["*2"]["status"] == STATUS_UNTESTABLE, token


def test_untestable_allele_plus_inferred_reference_downgrades_confidence():
    partial = {k: v for k, v in CYP2C19_CLEAR.items() if k != "rs4244285"}
    call = call_diplotype("CYP2C19", partial)
    assert call["reference_inferred"] is True
    assert CONFIDENCE_ORDER.index(call["confidence"]) <= CONFIDENCE_ORDER.index("low")


def test_full_coverage_reference_call_keeps_high_confidence():
    call = call_diplotype("CYP2C19", CYP2C19_CLEAR)
    assert call["confidence"] == "high"


def test_nothing_readable_produces_no_diplotype_at_all():
    call = call_diplotype("CYP2C19", {})
    assert call["diplotype"] is None
    assert call["confidence"] == "none"
    assert call["coverage"] == 0.0


def test_most_specific_allele_wins_and_consumes_its_component_variants():
    """SLCO1B1*15 is *1B plus *5. It must not also produce a separate *5."""
    call = call_diplotype("SLCO1B1",
                          with_variant(SLCO1B1_CLEAR, rs2306283="GG",
                                       rs4149056="CC"))
    assert call["diplotype"] == "*15/*15"
    assert call["alleles"].count("*15") == 2
    assert "*5" not in call["alleles"]


def test_single_copy_of_a_compound_allele_pairs_with_the_reference():
    call = call_diplotype("SLCO1B1",
                          with_variant(SLCO1B1_CLEAR, rs2306283="AG",
                                       rs4149056="CT"))
    assert call["diplotype"] == "*1/*15"
    assert call["reference_inferred"] is True


def test_compound_tpmt_allele_is_called_over_its_components():
    call = call_diplotype("TPMT",
                          with_variant(TPMT_CLEAR, rs1800460="TT",
                                       rs1142345="CC"))
    assert call["diplotype"] == "*3A/*3A"


def test_using_an_unverified_allele_definition_is_declared_and_costs_confidence():
    call = call_diplotype("TPMT",
                          with_variant(TPMT_CLEAR, rs1800460="TT",
                                       rs1142345="CC"))
    assert "*3A" in call["unverified_alleles_used"]
    assert CONFIDENCE_ORDER.index(call["confidence"]) < CONFIDENCE_ORDER.index("high")


def test_tag_only_alleles_are_flagged_in_the_caveats():
    call = call_diplotype("UGT1A1", with_variant(UGT1A1_CLEAR, rs887829="TT"))
    assert call["diplotype"] == "*28/*28"
    assert any("tag" in c.lower() for c in call["caveats"])


def test_every_call_carries_the_disclaimer():
    for gene, clear in CLEAR_BY_GENE.items():
        assert call_diplotype(gene, clear)["disclaimer"] == DISCLAIMER


# ---------------------------------------------------------------------------
# The CYP2D6 caveat
# ---------------------------------------------------------------------------

def test_provisional_genes_is_exactly_cyp2d6():
    assert PROVISIONAL_GENES == {"CYP2D6"}


@pytest.mark.parametrize("genotypes", [
    {},
    CYP2D6_CLEAR,
    with_variant(CYP2D6_CLEAR, rs3892097="AG"),
    with_variant(CYP2D6_CLEAR, rs3892097="AA"),
    with_variant(CYP2D6_CLEAR, rs1065852="AA"),
])
def test_every_cyp2d6_call_is_marked_provisional(genotypes):
    call = call_diplotype("CYP2D6", genotypes)
    assert call["provisional"] is True
    assert call["provisional_reason"]
    assert "copies" in call["provisional_reason"] or "duplicat" in call["provisional_reason"]


def test_cyp2d6_confidence_is_capped_even_with_complete_snp_coverage():
    call = call_diplotype("CYP2D6", with_variant(CYP2D6_CLEAR, rs3892097="AA"))
    assert CONFIDENCE_ORDER.index(call["confidence"]) <= CONFIDENCE_ORDER.index("low")


def test_cyp2d6_lists_the_alleles_no_array_can_ever_read():
    call = call_diplotype("CYP2D6", CYP2D6_CLEAR)
    names = [entry["allele"] for entry in call["structural_untestable"]]
    assert "*5" in names
    assert any("xN" in n for n in names)
    assert any("hybrid" in n.lower() for n in names)


def test_cyp2d6_untestable_alleles_include_the_structural_ones():
    call = call_diplotype("CYP2D6", CYP2D6_CLEAR)
    assert "*5" in call["untestable_alleles"]
    assert UNTESTABLE_ALLELES["CYP2D6"]


def test_cyp2d6_phenotype_is_provisional_in_the_translation_too():
    call = call_diplotype("CYP2D6", with_variant(CYP2D6_CLEAR, rs3892097="AA"))
    phenotype = translate_phenotype("CYP2D6", call)
    assert phenotype["provisional"] is True
    assert phenotype["provisional_reason"]


def test_cyp2d6_string_translation_is_still_marked_provisional():
    assert translate_phenotype("CYP2D6", "*1/*1")["provisional"] is True


def test_cyp2d6_reference_call_never_returns_normal():
    """A *1/*1 from array data could be a deletion or a duplication."""
    phenotype = translate_phenotype("CYP2D6", call_diplotype("CYP2D6", CYP2D6_CLEAR))
    assert phenotype["phenotype"] == "Indeterminate"


# ---------------------------------------------------------------------------
# Phenotype translation, including every activity-score boundary
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("diplotype,expected,score", [
    ("*2/*2", "Poor", 0.0),
    ("*2/*3", "Poor", 0.0),
    ("*1/*2", "Intermediate", 1.0),
    ("*2/*17", "Intermediate", 1.5),
    ("*1/*1", "Normal", 2.0),
    ("*1/*17", "Rapid", 2.5),
    ("*17/*17", "Ultrarapid", 3.0),
])
def test_cyp2c19_activity_score_boundaries(diplotype, expected, score):
    result = translate_phenotype("CYP2C19", diplotype)
    assert result["activity_score"] == score
    assert result["phenotype"] == expected
    assert result["phenotype_label"] == f"{expected} Metabolizer"


@pytest.mark.parametrize("diplotype,expected,score", [
    ("*3/*3", "Poor", 0.0),
    ("*2/*3", "Poor", 0.5),
    ("*1/*3", "Intermediate", 1.0),
    ("*2/*2", "Intermediate", 1.0),
    ("*1/*2", "Intermediate", 1.5),
    ("*1/*1", "Normal", 2.0),
])
def test_cyp2c9_activity_score_boundaries(diplotype, expected, score):
    result = translate_phenotype("CYP2C9", diplotype)
    assert result["activity_score"] == score
    assert result["phenotype"] == expected


@pytest.mark.parametrize("diplotype,expected,score", [
    ("*4/*4", "Poor", 0.0),
    ("*10/*10", "Intermediate", 0.5),
    ("*1/*4", "Intermediate", 1.0),
    ("*1/*10", "Normal", 1.25),
    ("*1/*1", "Normal", 2.0),
])
def test_cyp2d6_activity_score_boundaries(diplotype, expected, score):
    result = translate_phenotype("CYP2D6", diplotype)
    assert result["activity_score"] == score
    assert result["phenotype"] == expected


@pytest.mark.parametrize("diplotype,expected", [
    ("*1/*1", "Normal Function"),
    ("*1/*5", "Decreased Function"),
    ("*5/*5", "Poor Function"),
    ("*1/*15", "Decreased Function"),
])
def test_slco1b1_uses_the_transporter_function_scale(diplotype, expected):
    result = translate_phenotype("SLCO1B1", diplotype)
    assert result["scale"] == "transporter_function"
    assert result["phenotype"] == expected
    assert result["phenotype"] in PHENOTYPE_SCALES["transporter_function"]


@pytest.mark.parametrize("diplotype,expected", [
    ("-1639G/-1639G", "Normal sensitivity"),
    ("-1639G/-1639A", "Intermediate sensitivity"),
    ("-1639A/-1639A", "High sensitivity"),
])
def test_vkorc1_uses_the_warfarin_sensitivity_scale(diplotype, expected):
    result = translate_phenotype("VKORC1", diplotype)
    assert result["scale"] == "warfarin_sensitivity"
    assert result["phenotype"] == expected


@pytest.mark.parametrize("gene,diplotype,expected", [
    ("TPMT", "*1/*1", "Normal"),
    ("TPMT", "*1/*3C", "Intermediate"),
    ("TPMT", "*3C/*3C", "Poor"),
    ("NUDT15", "*1/*3", "Intermediate"),
    ("NUDT15", "*3/*3", "Poor"),
    ("DPYD", "*1/*1", "Normal"),
    ("DPYD", "*1/c.2846A>T", "Intermediate"),
    ("DPYD", "*1/*2A", "Intermediate"),
    ("DPYD", "*2A/*2A", "Poor"),
    ("UGT1A1", "*1/*28", "Intermediate"),
    ("UGT1A1", "*28/*28", "Poor"),
])
def test_remaining_gene_boundaries(gene, diplotype, expected):
    assert translate_phenotype(gene, diplotype)["phenotype"] == expected


def test_indeterminate_is_the_default_when_a_position_was_never_read():
    partial = {k: v for k, v in CYP2C19_CLEAR.items() if k != "rs4244285"}
    result = translate_phenotype("CYP2C19", call_diplotype("CYP2C19", partial))
    assert result["phenotype"] == "Indeterminate"
    assert "reference" in result["reason"].lower()


def test_normal_is_only_returned_when_the_reference_call_was_actually_earned():
    result = translate_phenotype("CYP2C19", call_diplotype("CYP2C19", CYP2C19_CLEAR))
    assert result["phenotype"] == "Normal"


def test_a_fully_detected_diplotype_survives_an_untested_rare_allele():
    """Both chromosomes observed, so an untested rare allele does not veto."""
    partial = {k: v for k, v in CYP2C19_CLEAR.items() if k != "rs72552267"}
    partial["rs4244285"] = "AA"
    result = translate_phenotype("CYP2C19", call_diplotype("CYP2C19", partial))
    assert result["phenotype"] == "Poor"


def test_empty_call_translates_to_indeterminate_not_normal():
    result = translate_phenotype("CYP2C19", call_diplotype("CYP2C19", {}))
    assert result["phenotype"] == "Indeterminate"
    assert result["activity_score"] is None


def test_unknown_allele_in_a_diplotype_string_is_indeterminate():
    result = translate_phenotype("CYP2C19", "*1/*9999")
    assert result["phenotype"] == "Indeterminate"
    assert "*9999" in result["reason"]


def test_malformed_diplotype_string_is_indeterminate():
    for text in ("*1", "", "*1/*2/*3", None, 17):
        assert translate_phenotype("CYP2C19", text)["phenotype"] == "Indeterminate"


def test_unknown_gene_translates_to_indeterminate():
    assert translate_phenotype("NOTAGENE", "*1/*1")["phenotype"] == "Indeterminate"


def test_every_phenotype_returned_is_in_that_genes_declared_scale():
    for gene, clear in CLEAR_BY_GENE.items():
        result = translate_phenotype(gene, call_diplotype(gene, clear))
        assert result["phenotype"] in PHENOTYPE_SCALES[result["scale"]]


def test_indeterminate_is_a_member_of_every_scale():
    for terms in PHENOTYPE_SCALES.values():
        assert "Indeterminate" in terms
    assert "Indeterminate" in PHENOTYPES


# ---------------------------------------------------------------------------
# Prescription guard
# ---------------------------------------------------------------------------

def build_diplotypes():
    return {
        "CYP2C19": call_diplotype("CYP2C19",
                                  with_variant(CYP2C19_CLEAR, rs4244285="AA")),
        "CYP2D6": call_diplotype("CYP2D6",
                                 with_variant(CYP2D6_CLEAR, rs3892097="AA")),
        "SLCO1B1": call_diplotype("SLCO1B1",
                                  with_variant(SLCO1B1_CLEAR, rs4149056="CC")),
    }


def test_guard_returns_only_pairs_for_the_medicines_supplied():
    result = prescription_guard(["clopidogrel"], build_diplotypes())
    assert result["count"] == 1
    assert result["matches"][0]["drug"] == "clopidogrel"
    assert result["matches"][0]["gene"] == "CYP2C19"


def test_guard_ignores_medicines_with_no_gene_pairing():
    result = prescription_guard(["aspirin", "paracetamol"], build_diplotypes())
    assert result["matches"] == []
    assert sorted(result["unmatched_medications"]) == ["aspirin", "paracetamol"]
    assert result["unmatched_note"]


def test_guard_skips_pairs_whose_gene_was_not_called():
    result = prescription_guard(["warfarin"], {"CYP2C19": "*1/*1"})
    assert result["matches"] == []


def test_guard_resolves_brand_names():
    result = prescription_guard(["Plavix"], build_diplotypes())
    assert result["count"] == 1
    assert result["matches"][0]["drug"] == "clopidogrel"
    assert result["matches"][0]["drug_as_entered"] == "Plavix"


def test_guard_tolerates_a_strength_suffix_on_the_medicine_name():
    result = prescription_guard(["simvastatin 20"], build_diplotypes())
    assert result["count"] == 1
    assert result["matches"][0]["drug"] == "simvastatin"


def test_guard_accepts_dicts_as_well_as_strings():
    result = prescription_guard([{"name": "clopidogrel"}], build_diplotypes())
    assert result["count"] == 1


def test_guard_reports_the_documented_category_not_an_instruction():
    result = prescription_guard(["clopidogrel"], build_diplotypes())
    match = result["matches"][0]
    assert match["category"] == "reduced_activation"
    assert match["next_step"].startswith("Discuss")


def test_guard_returns_insufficient_evidence_for_an_indeterminate_phenotype():
    partial = {k: v for k, v in CYP2C19_CLEAR.items() if k != "rs4244285"}
    result = prescription_guard(["clopidogrel"],
                                {"CYP2C19": call_diplotype("CYP2C19", partial)})
    match = result["matches"][0]
    assert match["phenotype"] == "Indeterminate"
    assert match["category"] == "insufficient_evidence"
    assert "not a reassuring result" in match["description"]


def test_guard_marks_cyp2d6_pairs_provisional():
    result = prescription_guard(["codeine"], build_diplotypes())
    assert result["matches"][0]["provisional"] is True
    assert result["matches"][0]["provisional_reason"]


def test_guard_names_pairs_with_no_cpic_guideline_rather_than_inventing_a_level():
    result = prescription_guard(["metoprolol"], build_diplotypes())
    match = result["matches"][0]
    assert match["cpic_level"] == ""
    assert "no guideline" in match["cpic_level_note"].lower()


def test_guard_attaches_the_disclaimer_to_the_payload_and_to_every_match():
    result = prescription_guard(["clopidogrel", "codeine"], build_diplotypes())
    assert result["disclaimer"] == DISCLAIMER
    assert result["matches"]
    for match in result["matches"]:
        assert match["disclaimer"] == DISCLAIMER


def test_guard_sorts_cpic_level_a_and_documented_risk_first():
    result = prescription_guard(["metoprolol", "clopidogrel"], build_diplotypes())
    assert result["matches"][0]["cpic_level"] == "A"


def test_guard_accepts_bare_diplotype_strings():
    result = prescription_guard(["clopidogrel"], {"CYP2C19": "*2/*2"})
    assert result["count"] == 1
    assert result["matches"][0]["phenotype"] == "Poor"


def test_guard_handles_an_empty_medication_list():
    result = prescription_guard([], build_diplotypes())
    assert result["matches"] == []
    assert result["disclaimer"] == DISCLAIMER


# ---------------------------------------------------------------------------
# The banned-imperative-language assertion.
#
# Every string in every prescription_guard payload is checked against the banned
# phrase list. This is the test that stops a well-meaning edit from turning an
# information page into a prescription.
# ---------------------------------------------------------------------------

ALL_MEDICATIONS = sorted(set(known_drugs()) | set(DRUG_ALIASES))


def all_gene_diplotypes():
    """One call per gene, at several phenotypes, to exercise every text path."""
    out = {}
    for gene, clear in CLEAR_BY_GENE.items():
        out[gene] = call_diplotype(gene, clear)
    return out


def test_no_output_string_contains_imperative_dosing_language():
    payload = prescription_guard(ALL_MEDICATIONS, all_gene_diplotypes())
    assert payload["count"] > 0
    assert audit_language(payload) == []


def test_no_imperative_language_at_any_phenotype_for_any_pair():
    for gene in CLEAR_BY_GENE:
        for diplotype in ALLELE_DEFINITIONS[gene]:
            call = {gene: f"{diplotype}/{REFERENCE_ALLELE[gene]}"}
            payload = prescription_guard(ALL_MEDICATIONS, call)
            assert audit_language(payload) == [], (gene, diplotype)


def test_no_imperative_language_in_the_indeterminate_path():
    payload = prescription_guard(ALL_MEDICATIONS,
                                 {"CYP2D6": call_diplotype("CYP2D6", {})})
    assert audit_language(payload) == []


def test_the_banned_phrase_list_is_not_empty_and_is_lower_case():
    assert BANNED_PRESCRIPTIVE_PHRASES
    for phrase in BANNED_PRESCRIPTIVE_PHRASES:
        assert phrase == phrase.lower()


def test_the_language_checker_actually_catches_imperative_text():
    bad = "You should stop taking simvastatin and reduce your dose."
    assert contains_prescriptive_language(bad)
    assert audit_language({"advice": bad})


def test_the_language_checker_passes_clean_descriptive_text():
    good = ("CPIC documents an increased risk of drug toxicity for this gene "
            "and medicine at this phenotype.")
    assert contains_prescriptive_language(good) == []


def test_every_drug_gene_pair_declares_a_gene_this_module_can_call():
    for pair in DRUG_GENE_PAIRS:
        assert pair["gene"] in ALLELE_DEFINITIONS
        assert pair["guideline"]
        assert pair["mechanism"]
