"""Tests for backend.carrier: carrier status wording, residual risk and ACMG coverage.

Two of these tests are the reason the module exists. The residual-risk
arithmetic is checked against values computed by hand from the Bayes expression
in the docstring, and every string every code path can produce is grepped for
the phrase "not a carrier" used without a scope. A carrier result that drops the
scope is a false reassurance delivered to somebody planning a pregnancy.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.carrier import (
    ACMG_ARRAY_PROBES,
    ACMG_SF_GENES,
    ACMG_SF_LIST_VERIFIED,
    ACMG_SF_PUBLISHED_COUNT,
    ACMG_SF_VERSION,
    CARRIER_PANEL,
    CARRIER_POPULATIONS,
    DISCLAIMER,
    FORBIDDEN_PHRASES,
    PANEL_GENES,
    STATUS_CARRIER,
    STATUS_NEGATIVE_FOR_TESTED,
    STATUS_UNTESTABLE,
    UNCERTAINTY_FACTOR,
    acmg_coverage_report,
    as_fraction,
    audit_wording,
    carrier_report,
    carrier_status,
    has_forbidden_phrasing,
    joint_reproductive_risk,
    negative_wording,
    normalise_population,
    residual_risk,
    residual_risk_value,
    unverified_figures,
)

# The em dash character itself, written as an escape so this file
# contains no literal one. House style forbids them in source.
EM_DASH = "\u2014"

MODULE_PATH = Path(__file__).parent.parent / "backend" / "carrier.py"


# ---------------------------------------------------------------------------
# Synthetic genotypes. "CLEAR" reads every array-testable position in the panel
# and carries none of the variant bases.
# ---------------------------------------------------------------------------

CFTR_CLEAR = {"rs113993959": "GG", "rs75527207": "GG",
              "rs77010898": "GG", "rs78655421": "GG"}
G6PD_CLEAR = {"rs1050828": "CC", "rs1050829": "TT", "rs5030868": "GG"}
HBB_CLEAR = {"rs334": "TT", "rs33930165": "CC"}
HEXA_CLEAR = {"rs147324677": "TT", "rs121907954": "GG"}
GJB2_CLEAR = {"rs72474224": "GG"}
PAH_CLEAR = {"rs5030858": "GG", "rs75193786": "GG"}
ATP7B_CLEAR = {"rs76151636": "CC"}
GALT_CLEAR = {"rs75391579": "AA", "rs2070074": "AA"}
ACADM_CLEAR = {"rs77931234": "AA"}
BTD_CLEAR = {"rs13078881": "TT"}

ALL_CLEAR = {}
for _block in (CFTR_CLEAR, G6PD_CLEAR, HBB_CLEAR, HEXA_CLEAR, GJB2_CLEAR,
               PAH_CLEAR, ATP7B_CLEAR, GALT_CLEAR, ACADM_CLEAR, BTD_CLEAR):
    ALL_CLEAR.update(_block)

CFTR_CARRIER = dict(CFTR_CLEAR, rs113993959="GT")
CFTR_HOMOZYGOUS = dict(CFTR_CLEAR, rs113993959="TT")
CFTR_PARTIAL = {"rs113993959": "GG"}          # one of four positions read


# ---------------------------------------------------------------------------
# Panel integrity
# ---------------------------------------------------------------------------

def test_every_required_gene_is_in_the_panel():
    for gene in ("CFTR", "HEXA", "SMN1", "HBB", "GJB2", "PAH", "ATP7B",
                 "GALT", "ACADM", "BTD", "G6PD"):
        assert gene in CARRIER_PANEL
        assert gene in PANEL_GENES


def test_every_panel_entry_carries_the_required_fields():
    for gene, entry in CARRIER_PANEL.items():
        assert entry["condition"], gene
        assert entry["inheritance"], gene
        assert isinstance(entry["tested_variants"], tuple), gene
        assert "total_known_pathogenic" in entry, gene
        assert isinstance(entry["carrier_frequency"], dict), gene
        assert isinstance(entry["detection_rate"], dict), gene


def test_every_figure_declares_whether_it_was_verified_and_states_its_basis():
    for gene, entry in CARRIER_PANEL.items():
        figures = [entry["total_known_pathogenic"]]
        figures += list(entry["carrier_frequency"].values())
        figures += list(entry["detection_rate"].values())
        for figure in figures:
            assert isinstance(figure["verified"], bool), gene
            assert figure["basis"], gene


def test_every_tested_variant_declares_array_testability_and_a_note():
    for gene, entry in CARRIER_PANEL.items():
        for variant in entry["tested_variants"]:
            assert variant["rsid"].startswith("rs"), gene
            assert isinstance(variant["array_testable"], bool), gene
            assert isinstance(variant["verified"], bool), gene
            assert variant["note"], (gene, variant["rsid"])


def test_unverified_figures_are_enumerable_for_the_documentation():
    rows = unverified_figures()
    assert rows
    for row in rows:
        assert row["kind"]
        assert row["note"]


def test_smn1_declares_that_it_cannot_be_tested_on_an_array():
    entry = CARRIER_PANEL["SMN1"]
    assert entry["tested_variants"] == ()
    assert entry["detection_rate"]["GENERAL"]["value"] == 0.0
    assert any("UNDETECTABLE" in note for note in entry["notes"])


def test_gjb2_common_european_allele_is_not_array_testable():
    entry = CARRIER_PANEL["GJB2"]
    deletion = [v for v in entry["tested_variants"] if "35delG" in v["name"]]
    assert deletion and deletion[0]["array_testable"] is False


def test_module_source_contains_no_em_dash():
    assert EM_DASH not in MODULE_PATH.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Wording
# ---------------------------------------------------------------------------

def test_negative_wording_is_always_scoped_to_a_count():
    assert negative_wording(4) == "not a carrier for the 4 variants tested"
    assert negative_wording(1) == "not a carrier for the 1 variant tested"
    assert negative_wording(0) == "not a carrier for the 0 variants tested"


def test_negative_wording_passes_the_forbidden_phrase_check():
    for n in range(0, 6):
        assert has_forbidden_phrasing(negative_wording(n)) == []


def test_the_bare_phrase_is_caught():
    assert "not a carrier" in has_forbidden_phrasing("You are not a carrier.")
    assert has_forbidden_phrasing("This is a non-carrier result")
    assert has_forbidden_phrasing("noncarrier")


def test_other_forbidden_reassurances_are_caught():
    assert has_forbidden_phrasing("This rules out cystic fibrosis")
    assert has_forbidden_phrasing("You have no risk")
    assert FORBIDDEN_PHRASES


def test_forbidden_phrasing_ignores_clean_text():
    clean = ("No variant was detected at the 4 positions this file could read, "
             "so this result is not a carrier for the 4 variants tested.")
    assert has_forbidden_phrasing(clean) == []


@pytest.mark.parametrize("gene,genotypes", [
    ("CFTR", CFTR_CLEAR), ("CFTR", CFTR_CARRIER), ("CFTR", CFTR_HOMOZYGOUS),
    ("CFTR", CFTR_PARTIAL), ("CFTR", {}),
    ("SMN1", ALL_CLEAR), ("SMN1", {}),
    ("HBB", HBB_CLEAR), ("HBB", dict(HBB_CLEAR, rs334="AT")),
    ("GJB2", GJB2_CLEAR), ("GJB2", {}),
    ("G6PD", G6PD_CLEAR), ("G6PD", dict(G6PD_CLEAR, rs1050828="CT")),
    ("HEXA", HEXA_CLEAR), ("PAH", PAH_CLEAR), ("ATP7B", ATP7B_CLEAR),
    ("GALT", GALT_CLEAR), ("ACADM", ACADM_CLEAR), ("BTD", BTD_CLEAR),
    ("NOTAGENE", ALL_CLEAR),
])
def test_no_output_path_ever_emits_the_forbidden_phrasing(gene, genotypes):
    assert audit_wording(carrier_status(gene, genotypes)) == []


def test_the_full_report_never_emits_the_forbidden_phrasing():
    for population in ("EUROPEAN", "AFRICAN", "ASHKENAZI", "", "nonsense"):
        report = carrier_report(ALL_CLEAR, population=population)
        assert audit_wording(report) == [], population


def test_residual_and_joint_payloads_never_emit_the_forbidden_phrasing():
    risk = residual_risk("CFTR", True, "EUROPEAN")
    assert audit_wording(risk) == []
    assert audit_wording(joint_reproductive_risk("CFTR", risk, risk,
                                                "autosomal recessive")) == []
    assert audit_wording(joint_reproductive_risk("G6PD", 0.1, None,
                                                 "X-linked")) == []


# ---------------------------------------------------------------------------
# Carrier status
# ---------------------------------------------------------------------------

def test_a_detected_variant_is_a_carrier_result():
    result = carrier_status("CFTR", CFTR_CARRIER)
    assert result["status"] == STATUS_CARRIER
    assert result["copies"] == 1
    assert result["detected"][0]["name"].startswith("G542X")


def test_two_copies_is_reported_as_needing_confirmation_not_as_carrier_wording():
    result = carrier_status("CFTR", CFTR_HOMOZYGOUS)
    assert result["status"] == STATUS_CARRIER
    assert result["copies"] == 2
    assert "Two copies" in result["statement"]
    assert "confirmation" in result["statement"]


def test_everything_read_and_nothing_found_is_scoped_to_the_count_read():
    result = carrier_status("CFTR", CFTR_CLEAR)
    assert result["status"] == STATUS_NEGATIVE_FOR_TESTED
    assert result["tested_count"] == 4
    assert "not a carrier for the 4 variants tested" in result["statement"]
    assert "2100" in result["statement"]


def test_nothing_readable_is_untestable_and_says_nothing_was_excluded():
    result = carrier_status("CFTR", {})
    assert result["status"] == STATUS_UNTESTABLE
    assert result["tested_count"] == 0
    assert "nothing was excluded" in result["statement"]


def test_untestable_and_negative_are_different_statuses():
    negative = carrier_status("CFTR", CFTR_CLEAR)
    untestable = carrier_status("CFTR", {})
    assert negative["status"] != untestable["status"]
    assert negative["tested_count"] > untestable["tested_count"]


def test_a_partial_read_scopes_the_negative_to_what_was_actually_read():
    result = carrier_status("CFTR", CFTR_PARTIAL)
    assert result["status"] == STATUS_NEGATIVE_FOR_TESTED
    assert result["tested_count"] == 1
    assert "not a carrier for the 1 variant tested" in result["statement"]
    assert len(result["not_read"]) == 3


def test_variants_no_array_can_read_are_listed_separately_from_unread_ones():
    result = carrier_status("CFTR", CFTR_CLEAR)
    names = [v["name"] for v in result["not_array_testable"]]
    assert any("F508del" in n for n in names)
    assert result["not_read"] == []


def test_smn1_is_always_untestable_however_complete_the_file_is():
    result = carrier_status("SMN1", ALL_CLEAR)
    assert result["status"] == STATUS_UNTESTABLE
    assert result["panel_size"] == 0


def test_unknown_gene_is_refused_rather_than_answered():
    result = carrier_status("NOTAGENE", ALL_CLEAR)
    assert result["known_gene"] is False
    assert result["status"] == STATUS_UNTESTABLE


def test_unverified_variant_mappings_are_declared_in_the_caveats():
    result = carrier_status("CFTR", CFTR_CLEAR)
    assert result["unverified_variants_used"]
    assert any("not corroborated" in c for c in result["caveats"])


def test_every_status_payload_carries_the_disclaimer():
    for gene in PANEL_GENES:
        assert carrier_status(gene, ALL_CLEAR)["disclaimer"] == DISCLAIMER


def test_the_fraction_of_known_variants_tested_is_reported_and_is_tiny():
    result = carrier_status("CFTR", CFTR_CLEAR)
    assert result["fraction_of_known_variants_tested"] == pytest.approx(
        4 / 2100, abs=1e-6)
    assert result["fraction_of_known_variants_tested"] < 0.01


def test_carrier_report_covers_every_panel_gene():
    report = carrier_report(ALL_CLEAR, population="EUROPEAN")
    assert len(report["results"]) == len(PANEL_GENES)
    assert sum(report["counts"].values()) == len(PANEL_GENES)
    assert report["summary"]


# ---------------------------------------------------------------------------
# Residual risk, against values computed by hand
# ---------------------------------------------------------------------------

def test_residual_risk_matches_the_hand_computed_cftr_european_value():
    # f = 1/25 = 0.04, DR = 0.88.
    # numerator   = 0.04 * (1 - 0.88) = 0.0048
    # denominator = 1 - (0.04 * 0.88) = 0.9648
    # residual    = 0.0048 / 0.9648   = 0.0049751243781...
    result = residual_risk("CFTR", True, "EUROPEAN")
    assert result["residual_risk"] == pytest.approx(0.0048 / 0.9648, rel=1e-12)
    assert result["residual_risk"] == pytest.approx(0.004975124, abs=1e-9)
    assert result["residual_risk_as_fraction"] == "1 in 201"


def test_residual_risk_matches_a_second_hand_computed_case():
    # f = 0.5, DR = 0.5 gives 0.25 / 0.75 = one third exactly.
    result = residual_risk("CFTR", True, "EUROPEAN",
                           detection_rate=0.5, carrier_frequency=0.5)
    assert result["residual_risk"] == pytest.approx(1 / 3, rel=1e-12)


def test_residual_risk_matches_a_third_hand_computed_case():
    # f = 0.1, DR = 0.9 gives 0.01 / 0.91.
    result = residual_risk("CFTR", True, "EUROPEAN",
                           detection_rate=0.9, carrier_frequency=0.1)
    assert result["residual_risk"] == pytest.approx(0.01 / 0.91, rel=1e-12)


def test_residual_risk_is_always_lower_than_the_prior():
    result = residual_risk("CFTR", True, "EUROPEAN")
    assert result["residual_risk"] < result["prior_risk"]


def test_a_perfect_detection_rate_drives_residual_risk_to_zero():
    result = residual_risk("CFTR", True, "EUROPEAN", detection_rate=1.0)
    assert result["residual_risk"] == pytest.approx(0.0)


def test_residual_risk_is_declared_a_lower_bound():
    result = residual_risk("CFTR", True, "EUROPEAN")
    assert result["is_lower_bound"] is True
    assert "LOWER BOUND" in result["reason"]


def test_residual_risk_returns_none_with_a_reason_when_the_frequency_is_unknown():
    result = residual_risk("CFTR", True, "MIDDLE_EASTERN")
    assert result["residual_risk"] is None
    assert "carrier frequency" in result["reason"]
    assert "MIDDLE_EASTERN" in result["reason"]


def test_residual_risk_returns_none_with_a_reason_when_the_detection_rate_is_unknown():
    result = residual_risk("PAH", True, "EUROPEAN")
    assert result["residual_risk"] is None
    assert "detection rate" in result["reason"]
    assert result["prior_risk"] == pytest.approx(1 / 50)


def test_residual_risk_returns_none_when_the_population_cannot_be_mapped():
    result = residual_risk("CFTR", True, "Klingon")
    assert result["residual_risk"] is None
    assert "population" in result["reason"]


def test_residual_risk_returns_none_for_an_unknown_gene():
    result = residual_risk("NOTAGENE", True, "EUROPEAN")
    assert result["residual_risk"] is None
    assert "not a gene" in result["reason"]


def test_residual_risk_refuses_when_there_was_no_negative_result():
    result = residual_risk("CFTR", False, "EUROPEAN")
    assert result["residual_risk"] is None
    assert "only after a negative result" in result["reason"]


def test_zero_variants_tested_leaves_the_risk_at_the_population_baseline():
    """Nothing tested means nothing subtracted, whatever the panel achieves."""
    result = residual_risk("CFTR", 0, "EUROPEAN")
    assert result["residual_risk"] == pytest.approx(1 / 25)
    assert result["detection_rate"] == 0.0
    assert "nothing was tested" in result["reason"]


def test_smn1_residual_risk_equals_the_carrier_frequency_because_dr_is_zero():
    result = residual_risk("SMN1", True, "EUROPEAN")
    assert result["residual_risk"] == pytest.approx(1 / 47)
    assert result["detection_rate"] == 0.0


def test_residual_risk_accepts_a_list_of_negative_variants():
    result = residual_risk("CFTR", ["rs113993959", "rs75527207"], "EUROPEAN")
    assert result["residual_risk"] == pytest.approx(0.0048 / 0.9648, rel=1e-12)


def test_residual_risk_flags_unverified_inputs():
    result = residual_risk("CFTR", True, "EUROPEAN")
    assert result["verified"] is False
    assert "unverified" in result["reason"]


def test_residual_risk_value_returns_a_bare_float_or_none():
    assert residual_risk_value("CFTR", True, "EUROPEAN") == pytest.approx(
        0.0048 / 0.9648, rel=1e-12)
    assert residual_risk_value("PAH", True, "EUROPEAN") is None


def test_residual_risk_records_the_basis_of_both_inputs():
    result = residual_risk("CFTR", True, "EUROPEAN")
    assert result["basis"]["carrier_frequency"]
    assert result["basis"]["detection_rate"]
    assert "ACMG" in result["basis"]["panel"] or result["basis"]["panel"]


def test_population_aliases_map_onto_the_carrier_groups():
    assert normalise_population("CEU") == "EUROPEAN"
    assert normalise_population("YRI") == "AFRICAN"
    assert normalise_population("Ashkenazi Jewish") == "ASHKENAZI"
    assert normalise_population("CHB") == "EAST_ASIAN"
    assert normalise_population("Klingon") == ""
    for code in CARRIER_POPULATIONS:
        assert normalise_population(code) == code


def test_as_fraction_renders_risks_the_way_people_read_them():
    assert as_fraction(0.04) == "1 in 25"
    assert as_fraction(0.004975124) == "1 in 201"
    assert as_fraction(None) == ""
    assert as_fraction(0) == "0"


# ---------------------------------------------------------------------------
# Joint reproductive risk
# ---------------------------------------------------------------------------

def test_autosomal_recessive_joint_risk_is_the_product_times_a_quarter():
    result = joint_reproductive_risk("CFTR", 0.04, 0.04, "autosomal recessive")
    assert result["point"] == pytest.approx(0.04 * 0.04 * 0.25)
    assert result["point"] == pytest.approx(0.0004)


def test_two_obligate_carriers_give_exactly_one_in_four():
    result = joint_reproductive_risk("CFTR", 1.0, 1.0, "autosomal recessive")
    assert result["point"] == pytest.approx(0.25)
    assert result["high"] == pytest.approx(0.25)


def test_joint_risk_accepts_residual_risk_payloads_directly():
    risk = residual_risk("CFTR", True, "EUROPEAN")
    result = joint_reproductive_risk("CFTR", risk, risk, "autosomal recessive")
    assert result["point"] == pytest.approx(risk["residual_risk"] ** 2 * 0.25)


def test_joint_risk_returns_a_range_when_an_input_is_unverified():
    result = joint_reproductive_risk("CFTR", 0.04, 0.04, "autosomal recessive")
    assert result["low"] == pytest.approx(result["point"] / UNCERTAINTY_FACTOR)
    assert result["high"] == pytest.approx(result["point"] * UNCERTAINTY_FACTOR)
    assert result["low"] < result["high"]
    assert "false precision" in result["range_basis"]
    assert " to " in result["range_as_fraction"]


def test_x_linked_uses_a_different_formula_and_ignores_the_second_parent():
    result = joint_reproductive_risk("G6PD", 0.11, 0.99, "X-linked")
    assert result["point"] == pytest.approx(0.11 * 0.5 * 0.5)
    assert any("not used" in a for a in result["assumptions"])


def test_x_linked_risk_does_not_change_when_the_second_input_changes():
    a = joint_reproductive_risk("G6PD", 0.11, 0.0, "X-linked")
    b = joint_reproductive_risk("G6PD", 0.11, 1.0, "X-linked")
    assert a["point"] == b["point"]


def test_x_linked_works_with_no_second_input_at_all():
    result = joint_reproductive_risk("G6PD", 0.11, None, "X-linked")
    assert result["point"] == pytest.approx(0.0275)


def test_autosomal_recessive_refuses_when_the_second_input_is_missing():
    result = joint_reproductive_risk("CFTR", 0.04, None, "autosomal recessive")
    assert result["point"] is None
    assert "second carrier risk is unknown" in result["reason"]


def test_joint_risk_refuses_when_the_first_input_is_missing():
    result = joint_reproductive_risk("CFTR", None, 0.04, "autosomal recessive")
    assert result["point"] is None
    assert "first carrier risk is unknown" in result["reason"]


def test_joint_risk_refuses_an_inheritance_pattern_it_does_not_model():
    result = joint_reproductive_risk("CFTR", 0.04, 0.04, "autosomal dominant")
    assert result["point"] is None
    assert "implemented for" in result["reason"]


def test_joint_risk_defaults_to_the_panels_inheritance_when_none_is_given():
    result = joint_reproductive_risk("CFTR", 0.04, 0.04)
    assert result["inheritance"] == "autosomal recessive"
    assert result["point"] == pytest.approx(0.0004)


def test_joint_risk_high_bound_is_clamped_at_the_arithmetic_maximum():
    result = joint_reproductive_risk("CFTR", 0.9, 0.9, "autosomal recessive")
    assert result["high"] <= 0.25


def test_joint_risk_states_its_assumptions_and_caveats():
    result = joint_reproductive_risk("CFTR", 0.04, 0.04, "autosomal recessive")
    assert result["assumptions"]
    assert result["caveats"]
    assert result["disclaimer"] == DISCLAIMER


# ---------------------------------------------------------------------------
# ACMG secondary findings
# ---------------------------------------------------------------------------

def test_the_acmg_list_states_its_version_and_admits_it_is_unverified():
    assert "v3.2" in ACMG_SF_VERSION
    assert ACMG_SF_LIST_VERIFIED is False
    assert len(ACMG_SF_GENES) > 50


def test_the_acmg_list_has_no_duplicates_and_is_upper_case():
    assert len(set(ACMG_SF_GENES)) == len(ACMG_SF_GENES)
    for gene in ACMG_SF_GENES:
        assert gene == gene.upper()


def test_the_acmg_report_states_the_discrepancy_with_the_published_count():
    report = acmg_coverage_report({})
    assert report["genes_published"] == ACMG_SF_PUBLISHED_COUNT
    assert str(len(ACMG_SF_GENES)) in report["list_discrepancy_note"]
    assert report["list_verified"] is False


def test_acmg_coverage_is_zero_for_a_file_with_no_relevant_positions():
    report = acmg_coverage_report({})
    assert report["genes_with_any_coverage"] == 0
    assert report["genes_with_zero_coverage"] == len(ACMG_SF_GENES)
    assert report["fraction_of_genes_with_any_coverage"] == 0.0


def test_acmg_coverage_says_zero_in_words_for_every_uncovered_gene():
    report = acmg_coverage_report({})
    for row in report["genes"]:
        assert row["assessment"].startswith("Zero"), row["gene"]
        assert row["positions_read"] == 0


def test_acmg_coverage_is_still_effectively_zero_when_a_probe_is_present():
    report = acmg_coverage_report({"rs80357713": "GG"})
    assert report["genes_with_any_coverage"] == 1
    brca1 = next(r for r in report["genes"] if r["gene"] == "BRCA1")
    assert brca1["positions_read"] == 1
    assert "Effectively zero" in brca1["assessment"]
    assert brca1["fraction"] is not None and brca1["fraction"] < 0.001


def test_acmg_summary_says_plainly_that_nothing_was_looked_at():
    report = acmg_coverage_report({})
    assert "nothing was looked at" in report["summary"]
    assert "not low, not partial, zero" in report["summary"]


def test_acmg_probe_note_denies_that_three_variants_is_brca_testing():
    report = acmg_coverage_report({})
    assert "not BRCA testing" in report["probe_note"]
    assert set(ACMG_ARRAY_PROBES) <= set(ACMG_SF_GENES)


def test_acmg_report_never_emits_forbidden_phrasing():
    assert audit_wording(acmg_coverage_report({"rs80357713": "GG"})) == []


def test_a_no_call_probe_does_not_count_as_acmg_coverage():
    report = acmg_coverage_report({"rs80357713": "--"})
    assert report["genes_with_any_coverage"] == 0
