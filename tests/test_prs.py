"""Tests for backend.prs additive polygenic scoring."""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

from backend import prs
from backend.prs import (
    BAND_MAGNITUDE,
    BANDS,
    CAVEAT_ANCESTRY,
    CAVEAT_BOTH_DIRECTIONS,
    CAVEAT_NOT_DIAGNOSTIC,
    MANDATORY_CAVEATS,
    RELIABLE_COVERAGE,
    band_for_percentile,
    build_caveats,
    compute_all,
    compute_model,
    coverage_caveat,
    dosage,
    get_model,
    list_models,
    load_models,
    model_ids,
    normal_cdf,
    percentile_from_z,
    to_findings,
)
from data.build_reference import REFERENCE

MODELS_PATH = Path(__file__).parent.parent / "data" / "prs_models.json"


def bundled_models() -> dict:
    """Return the bundled models, skipping the test when they are not built."""
    if not MODELS_PATH.exists():
        pytest.skip("data/prs_models.json has not been built")
    models = load_models()
    if not models:
        pytest.skip("data/prs_models.json carries no models")
    return models


def synthetic_model() -> dict:
    """Return a four-variant model with hand-checkable weights.

    Weights and frequencies are round numbers so every expected score in the
    tests below can be verified by hand rather than by rerunning the code.
    """
    return {
        "id": "synthetic",
        "trait": "Synthetic trait",
        "license": "MIT",
        "citation": "DNAInsight test fixture",
        "variants": [
            {"rsid": "rs1", "effect_allele": "A", "other_allele": "G",
             "weight": 1.0, "effect_allele_frequency": 0.5},
            {"rsid": "rs2", "effect_allele": "C", "other_allele": "T",
             "weight": 2.0, "effect_allele_frequency": 0.25},
            {"rsid": "rs3", "effect_allele": "T", "other_allele": "A",
             "weight": -0.5, "effect_allele_frequency": 0.1},
            {"rsid": "rs4", "effect_allele": "G", "other_allele": "C",
             "weight": 0.25, "effect_allele_frequency": 0.4},
        ],
        "reference": {"population": "EUR", "mean": 0.0, "sd": 1.0},
    }


def model_without_frequencies() -> dict:
    """Return a two-variant model that supplies no effect allele frequencies."""
    return {
        "id": "no_freq",
        "trait": "No frequency trait",
        "variants": [
            {"rsid": "rs1", "effect_allele": "A", "other_allele": "G",
             "weight": 1.0, "effect_allele_frequency": None},
            {"rsid": "rs2", "effect_allele": "C", "other_allele": "T",
             "weight": 2.0},
        ],
        "reference": {"population": "EUR", "mean": 0.0, "sd": 1.0},
    }


def ten_variant_model() -> dict:
    """Return a ten-variant model, so coverage lands on exact tenths."""
    return {
        "id": "ten",
        "trait": "Ten variant trait",
        "variants": [
            {"rsid": f"rs{n}", "effect_allele": "A", "other_allele": "G",
             "weight": 0.1, "effect_allele_frequency": 0.5}
            for n in range(1, 11)
        ],
        "reference": {"population": "EUR", "mean": 1.0, "sd": 0.5},
    }


def recompute_moments(variants: list) -> tuple:
    """Recompute a model's analytic Hardy-Weinberg mean and sd from scratch.

    Deliberately a second, independent implementation of the formula in
    data/build_prs.py so that the two have to agree.
    """
    mean = 0.0
    variance = 0.0
    for variant in variants:
        f = float(variant["effect_allele_frequency"])
        w = float(variant["weight"])
        mean += 2.0 * f * w
        variance += 2.0 * f * (1.0 - f) * w * w
    return mean, math.sqrt(variance)


class TestDosage:
    def test_homozygous_effect_allele_is_two(self):
        assert dosage("A", "A", "A") == 2

    def test_heterozygous_is_one(self):
        assert dosage("A", "G", "A") == 1

    def test_heterozygous_either_order_is_one(self):
        assert dosage("G", "A", "A") == 1

    def test_homozygous_other_allele_is_zero(self):
        assert dosage("G", "G", "A") == 0

    def test_zero_is_not_none(self):
        # A real 0 and an absent call are different facts and must not collapse.
        assert dosage("G", "G", "A") is not None

    def test_dash_no_call_is_none(self):
        assert dosage("-", "-", "A") is None

    def test_single_sided_no_call_is_none(self):
        assert dosage("A", "N", "A") is None
        assert dosage("N", "A", "A") is None

    def test_every_no_call_token_is_none(self):
        for token in ("", "N", "-", "--", "0", "D", "I"):
            assert dosage(token, "A", "A") is None

    def test_missing_effect_allele_is_none(self):
        assert dosage("A", "A", "") is None
        assert dosage("A", "A", "N") is None

    def test_comparison_is_case_insensitive(self):
        assert dosage("a", "a", "A") == 2
        assert dosage("A", "G", "a") == 1

    def test_whitespace_is_tolerated(self):
        assert dosage(" A ", "A", "A") == 2


class TestNormalCdf:
    def test_cdf_at_zero_is_one_half(self):
        assert abs(normal_cdf(0.0) - 0.5) < 1e-12

    def test_cdf_is_symmetric(self):
        for z in (0.25, 1.0, 1.96, 3.0):
            assert abs(normal_cdf(z) + normal_cdf(-z) - 1.0) < 1e-12

    def test_cdf_is_monotonic(self):
        values = [normal_cdf(z) for z in (-3.0, -1.0, 0.0, 1.0, 3.0)]
        assert values == sorted(values)

    def test_cdf_stays_inside_the_unit_interval(self):
        for z in (-40.0, -5.0, 0.0, 5.0, 40.0):
            assert 0.0 <= normal_cdf(z) <= 1.0

    def test_cdf_matches_the_one_sigma_value(self):
        assert abs(normal_cdf(1.0) - 0.8413447) < 1e-6


class TestPercentileFromZ:
    def test_z_zero_is_about_fifty(self):
        assert abs(percentile_from_z(0.0) - 50.0) < 0.05

    def test_z_one_point_six_four_five_is_about_ninety_five(self):
        assert abs(percentile_from_z(1.645) - 95.0) < 0.05

    def test_z_minus_one_point_two_eight_two_is_about_ten(self):
        assert abs(percentile_from_z(-1.282) - 10.0) < 0.05

    def test_none_z_gives_none_percentile(self):
        assert percentile_from_z(None) is None

    def test_percentile_is_rounded_to_one_decimal(self):
        value = percentile_from_z(0.5)
        assert value == round(value, 1)

    def test_percentile_never_leaves_zero_to_one_hundred(self):
        for z in (-12.0, -2.0, 0.0, 2.0, 12.0):
            assert 0.0 <= percentile_from_z(z) <= 100.0

    def test_percentile_uses_the_modules_own_cdf(self):
        # The percentile must be nothing more than the erf-based CDF scaled,
        # so no separate approximation can drift away from it.
        for z in (-1.5, -0.3, 0.0, 0.8, 2.2):
            assert percentile_from_z(z) == round(normal_cdf(z) * 100.0, 1)


class TestBandBoundaries:
    def test_zero_percentile_is_low(self):
        assert band_for_percentile(0.0) == "low"

    def test_just_under_ten_is_low(self):
        assert band_for_percentile(9.9) == "low"

    def test_exactly_ten_is_below_average(self):
        assert band_for_percentile(10.0) == "below_average"

    def test_just_under_thirty_is_below_average(self):
        assert band_for_percentile(29.9) == "below_average"

    def test_exactly_thirty_is_average(self):
        assert band_for_percentile(30.0) == "average"

    def test_fifty_is_average(self):
        assert band_for_percentile(50.0) == "average"

    def test_exactly_seventy_is_average(self):
        assert band_for_percentile(70.0) == "average"

    def test_just_over_seventy_is_above_average(self):
        assert band_for_percentile(70.1) == "above_average"

    def test_exactly_ninety_is_above_average(self):
        assert band_for_percentile(90.0) == "above_average"

    def test_just_over_ninety_is_high(self):
        assert band_for_percentile(90.1) == "high"

    def test_one_hundred_is_high(self):
        assert band_for_percentile(100.0) == "high"

    def test_none_is_unknown(self):
        assert band_for_percentile(None) == "unknown"

    def test_unparsable_percentile_is_unknown(self):
        assert band_for_percentile("not a number") == "unknown"

    def test_every_band_returned_is_a_declared_band(self):
        for percentile in (0.0, 9.9, 10.0, 29.9, 30.0, 70.0, 90.0, 100.0, None):
            assert band_for_percentile(percentile) in BANDS

    def test_every_band_has_a_magnitude(self):
        for band in BANDS:
            assert band in BAND_MAGNITUDE


class TestComputeModelArithmetic:
    def test_full_coverage_score_is_the_hand_computed_sum(self):
        # rs1 AA = 2 x 1.0, rs2 CT = 1 x 2.0, rs3 AA = 0, rs4 GC = 1 x 0.25
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["raw_score"] == 4.25

    def test_full_coverage_counts_every_variant(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["variants_used"] == 4
        assert result["variants_total"] == 4
        assert result["coverage"] == 1.0

    def test_three_of_four_variants_gives_coverage_zero_point_seven_five(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["variants_used"] == 3
        assert result["coverage"] == 0.75

    def test_no_genotypes_gives_zero_coverage(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        assert result["variants_used"] == 0
        assert result["coverage"] == 0.0

    def test_missing_rsids_are_reported(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["missing_rsids"] == ["rs4"]

    def test_no_call_genotype_counts_as_missing(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "--"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["variants_used"] == 3
        assert "rs4" in result["missing_rsids"]

    def test_rsid_lookup_is_case_insensitive(self):
        genotypes = {"RS1": "AA", "RS2": "CT", "RS3": "AA", "RS4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["coverage"] == 1.0

    def test_tuple_genotypes_are_accepted(self):
        genotypes = {"rs1": ("A", "A"), "rs2": ("C", "T"),
                     "rs3": ("A", "A"), "rs4": ("G", "C")}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["raw_score"] == 4.25

    def test_malformed_genotype_is_ignored_not_crashed_on(self):
        genotypes = {"rs1": "AAA", "rs2": None, "rs3": 7, "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["variants_used"] == 1

    def test_negative_weight_lowers_the_score(self):
        # rs3 carries weight -0.5, so two copies must subtract 1.0.
        base = compute_model("synthetic", {"rs3": "AA"},
                             model=synthetic_model())["raw_score"]
        carrier = compute_model("synthetic", {"rs3": "TT"},
                                model=synthetic_model())["raw_score"]
        assert carrier < base
        assert abs((carrier - base) + 1.0) < 1e-9

    def test_z_score_uses_the_stored_reference(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["z"] == 4.25

    def test_missing_reference_gives_no_percentile_and_unknown_band(self):
        model = synthetic_model()
        model["reference"] = {"population": "EUR"}
        result = compute_model("synthetic", {"rs1": "AA"}, model=model)
        assert result["percentile"] is None
        assert result["z"] is None
        assert result["band"] == "unknown"

    def test_zero_standard_deviation_gives_no_percentile(self):
        model = synthetic_model()
        model["reference"] = {"population": "EUR", "mean": 1.0, "sd": 0.0}
        result = compute_model("synthetic", {"rs1": "AA"}, model=model)
        assert result["percentile"] is None

    def test_result_carries_every_documented_key(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        for key in ("id", "trait", "raw_score", "variants_used",
                    "variants_total", "coverage", "missing_rsids",
                    "mean_imputed", "percentile", "band", "z",
                    "population_reference", "license", "citation",
                    "caveats", "reliable"):
            assert key in result

    def test_empty_variant_list_does_not_divide_by_zero(self):
        model = synthetic_model()
        model["variants"] = []
        result = compute_model("synthetic", {"rs1": "AA"}, model=model)
        assert result["coverage"] == 0.0
        assert result["variants_total"] == 0


class TestMeanImputation:
    def test_missing_variant_with_a_frequency_sets_mean_imputed(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["mean_imputed"] is True

    def test_full_coverage_does_not_set_mean_imputed(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["mean_imputed"] is False

    def test_imputation_does_not_inflate_coverage(self):
        # rs4 is imputed at 2 x 0.4 x 0.25 = 0.2, but coverage counts only the
        # three variants that were actually genotyped.
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["coverage"] == 0.75
        assert result["variants_used"] == 3
        assert result["raw_score"] == 4.2

    def test_imputed_variant_is_still_listed_as_missing(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert "rs4" in result["missing_rsids"]

    def test_a_model_without_frequencies_does_not_impute(self):
        result = compute_model("no_freq", {}, model=model_without_frequencies())
        assert result["mean_imputed"] is False
        assert result["raw_score"] == 0.0
        assert result["coverage"] == 0.0

    def test_fully_imputed_score_equals_the_analytic_mean(self):
        # Imputing every variant at 2f must reproduce sum(2 f w), which is the
        # same quantity data/build_prs.py stores as the reference mean.
        model = synthetic_model()
        result = compute_model("synthetic", {}, model=model)
        expected_mean, _ = recompute_moments(model["variants"])
        assert abs(result["raw_score"] - expected_mean) < 1e-6

    def test_imputation_pulls_the_score_toward_the_middle(self):
        model = synthetic_model()
        model["reference"] = {"population": "EUR", "mean": 2.1, "sd": 1.0}
        result = compute_model("synthetic", {}, model=model)
        assert abs(result["percentile"] - 50.0) < 0.5


class TestReliability:
    def test_reliable_threshold_constant_is_ninety_percent(self):
        assert RELIABLE_COVERAGE == 0.90

    def test_full_coverage_is_reliable(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["reliable"] is True

    def test_seventy_five_percent_coverage_is_not_reliable(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert result["reliable"] is False

    def test_exactly_ninety_percent_coverage_is_reliable(self):
        genotypes = {f"rs{n}": "AA" for n in range(1, 10)}
        result = compute_model("ten", genotypes, model=ten_variant_model())
        assert result["coverage"] == 0.9
        assert result["reliable"] is True

    def test_eighty_percent_coverage_is_not_reliable(self):
        genotypes = {f"rs{n}": "AA" for n in range(1, 9)}
        result = compute_model("ten", genotypes, model=ten_variant_model())
        assert result["coverage"] == 0.8
        assert result["reliable"] is False

    def test_zero_coverage_is_not_reliable(self):
        result = compute_model("ten", {}, model=ten_variant_model())
        assert result["reliable"] is False


class TestCaveats:
    def test_build_caveats_returns_four_strings(self):
        caveats = build_caveats(0.5)
        assert len(caveats) == 4
        assert all(isinstance(c, str) and c for c in caveats)

    def test_the_three_mandatory_caveats_are_present(self):
        caveats = build_caveats(0.5)
        for mandatory in MANDATORY_CAVEATS:
            assert mandatory in caveats

    def test_mandatory_caveats_name_the_three_structural_limits(self):
        assert CAVEAT_NOT_DIAGNOSTIC in MANDATORY_CAVEATS
        assert CAVEAT_BOTH_DIRECTIONS in MANDATORY_CAVEATS
        assert CAVEAT_ANCESTRY in MANDATORY_CAVEATS

    def test_coverage_caveat_carries_the_real_percentage(self):
        assert "75.0 percent" in coverage_caveat(0.75)
        assert "100.0 percent" in coverage_caveat(1.0)
        assert "0.0 percent" in coverage_caveat(0.0)

    def test_result_caveats_quote_the_results_own_coverage(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        assert any("75.0 percent" in c for c in result["caveats"])

    def test_every_result_carries_all_mandatory_caveats(self):
        for genotypes in ({}, {"rs1": "AA"},
                          {"rs1": "AA", "rs2": "CT", "rs3": "AA", "rs4": "GC"}):
            result = compute_model("synthetic", genotypes,
                                   model=synthetic_model())
            for mandatory in MANDATORY_CAVEATS:
                assert mandatory in result["caveats"]
            assert len(result["caveats"]) == 4

    def test_ancestry_caveat_is_stated_not_implied(self):
        assert "ancestry" in CAVEAT_ANCESTRY.lower()
        assert "European" in CAVEAT_ANCESTRY

    def test_caveats_survive_a_model_with_no_reference(self):
        model = synthetic_model()
        model.pop("reference")
        result = compute_model("synthetic", {"rs1": "AA"}, model=model)
        assert len(result["caveats"]) == 4


class TestToFindings:
    def test_entity_type_is_prs_for_every_finding(self):
        results = [compute_model("synthetic", {}, model=synthetic_model()),
                   compute_model("ten", {}, model=ten_variant_model())]
        findings = to_findings(results)
        assert findings
        assert all(f["entity_type"] == "prs" for f in findings)

    def test_repute_is_empty_for_every_finding(self):
        results = [compute_model("synthetic", {}, model=synthetic_model()),
                   compute_model("ten", {}, model=ten_variant_model())]
        findings = to_findings(results)
        assert all(f["repute"] == "" for f in findings)

    def test_rsid_field_carries_the_model_id(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        finding = to_findings([result])[0]
        assert finding["rsid"] == "synthetic"

    def test_silo_and_category_are_informational_prs(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        finding = to_findings([result])[0]
        assert finding["silo"] == "informational"
        assert finding["category"] == "PRS"

    def test_magnitude_follows_the_band(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        finding = to_findings([result])[0]
        assert finding["magnitude"] == BAND_MAGNITUDE[result["band"]]

    def test_full_result_is_carried_through_under_prs(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        finding = to_findings([result])[0]
        assert finding["prs"] is result
        assert finding["prs"]["caveats"] == result["caveats"]

    def test_interpretation_mentions_coverage(self):
        genotypes = {"rs1": "AA", "rs2": "CT", "rs3": "AA"}
        result = compute_model("synthetic", genotypes, model=synthetic_model())
        finding = to_findings([result])[0]
        assert "coverage" in finding["interpretation"]

    def test_low_coverage_interpretation_says_indicative_only(self):
        result = compute_model("synthetic", {}, model=synthetic_model())
        finding = to_findings([result])[0]
        assert "indicative only" in finding["interpretation"]

    def test_empty_and_none_input_give_an_empty_list(self):
        assert to_findings([]) == []
        assert to_findings(None) == []

    def test_non_dict_entries_are_skipped(self):
        findings = to_findings(["nonsense", 7, None])
        assert findings == []


def effect_homozygous_genotypes(models: dict) -> dict:
    """Return a genotype map that is effect-homozygous at every model variant.

    One rsID can appear in two models with opposite effect alleles, rs1801282
    being the real case, so the last model wins for that position. That is
    fine here: the tests using this map assert on coverage and on direction
    against the reference mean, not on an exact score.
    """
    genotypes: dict = {}
    for model in models.values():
        for variant in model.get("variants") or []:
            genotypes[variant["rsid"]] = variant["effect_allele"] * 2
    return genotypes


class TestUnknownModel:
    def test_unknown_model_id_raises_key_error(self):
        with pytest.raises(KeyError):
            compute_model("no_such_model", {"rs1": "AA"})

    def test_error_message_names_the_unknown_id(self):
        with pytest.raises(KeyError) as excinfo:
            compute_model("no_such_model", {})
        assert "no_such_model" in str(excinfo.value)

    def test_error_message_lists_what_is_available(self):
        with pytest.raises(KeyError) as excinfo:
            compute_model("no_such_model", {})
        assert "Available ids" in str(excinfo.value)

    def test_get_model_returns_none_rather_than_raising(self):
        assert get_model("no_such_model") is None

    def test_compute_all_never_raises_on_odd_input(self):
        for genotypes in ({}, {"rs1": None}, {"": ""}, {"rs1": "AAA"}):
            assert isinstance(compute_all(genotypes), list)


class TestModelLoading:
    def test_load_models_returns_a_dict(self):
        assert isinstance(load_models(), dict)

    def test_model_ids_match_the_loaded_keys(self):
        assert model_ids() == list(load_models().keys())

    def test_list_models_omits_the_weights(self):
        bundled_models()
        for summary in list_models():
            assert "variants" not in summary

    def test_list_models_reports_the_variant_count(self):
        models = bundled_models()
        for summary in list_models():
            stored = models[summary["id"]]
            assert summary["variant_count"] == len(stored["variants"])

    def test_a_missing_models_file_returns_an_empty_dict(self):
        absent = Path(__file__).parent / "no_such_prs_models.json"
        try:
            assert load_models(absent) == {}
        finally:
            # Restore the module cache for the tests that follow.
            load_models(prs.MODELS_FILE)

    def test_metadata_declares_a_version(self):
        bundled_models()
        assert prs.get_metadata().get("version")


class TestBundledModels:
    def test_at_least_seven_models_are_bundled(self):
        assert len(bundled_models()) >= 7

    def test_the_seven_expected_model_ids_are_present(self):
        models = bundled_models()
        for model_id in ("t2d", "cad", "bmi", "vte", "ldl",
                         "homocysteine", "inflammation"):
            assert model_id in models

    def test_stored_reference_mean_matches_a_recomputation(self):
        for model_id, model in bundled_models().items():
            expected_mean, _ = recompute_moments(model["variants"])
            stored = float(model["reference"]["mean"])
            assert abs(stored - expected_mean) < 1e-6, model_id

    def test_stored_reference_sd_matches_a_recomputation(self):
        for model_id, model in bundled_models().items():
            _, expected_sd = recompute_moments(model["variants"])
            stored = float(model["reference"]["sd"])
            assert abs(stored - expected_sd) < 1e-6, model_id

    def test_every_reference_sd_is_positive(self):
        for model_id, model in bundled_models().items():
            assert float(model["reference"]["sd"]) > 0.0, model_id

    def test_every_reference_declares_its_population(self):
        for model_id, model in bundled_models().items():
            assert model["reference"]["population"], model_id

    def test_every_variant_rsid_is_in_the_bundled_reference(self):
        known = {str(row[0]).strip().lower() for row in REFERENCE}
        for model_id, model in bundled_models().items():
            for variant in model["variants"]:
                assert variant["rsid"].lower() in known, \
                    f"{model_id} uses unbundled {variant['rsid']}"

    def test_no_model_license_is_noncommercial(self):
        for model_id, model in bundled_models().items():
            assert "NonCommercial" not in str(model.get("license", "")), model_id

    def test_no_model_license_is_no_derivatives(self):
        for model_id, model in bundled_models().items():
            assert "NoDerivatives" not in str(model.get("license", "")), model_id

    def test_variant_count_matches_the_variant_list(self):
        for model_id, model in bundled_models().items():
            assert model["variant_count"] == len(model["variants"]), model_id

    def test_every_variant_carries_the_required_keys(self):
        for model_id, model in bundled_models().items():
            for variant in model["variants"]:
                for key in ("rsid", "effect_allele", "other_allele", "weight",
                            "effect_allele_frequency", "af_source"):
                    assert key in variant, f"{model_id} {variant.get('rsid')} {key}"

    def test_effect_and_other_alleles_are_distinct_real_bases(self):
        for model_id, model in bundled_models().items():
            for variant in model["variants"]:
                effect, other = variant["effect_allele"], variant["other_allele"]
                assert effect in ("A", "C", "G", "T"), model_id
                assert other in ("A", "C", "G", "T"), model_id
                assert effect != other, f"{model_id} {variant['rsid']}"

    def test_effect_allele_frequencies_are_proper_fractions(self):
        for model_id, model in bundled_models().items():
            for variant in model["variants"]:
                freq = variant["effect_allele_frequency"]
                assert isinstance(freq, float), model_id
                assert 0.0 < freq < 1.0, f"{model_id} {variant['rsid']} {freq}"

    def test_every_weight_is_the_log_of_its_stored_odds_ratio(self):
        for model_id, model in bundled_models().items():
            for variant in model["variants"]:
                if "odds_ratio" not in variant:
                    continue
                expected = math.log(float(variant["odds_ratio"]))
                assert abs(float(variant["weight"]) - expected) < 1e-5, \
                    f"{model_id} {variant['rsid']}"

    def test_every_model_declares_the_grch37_build(self):
        for model_id, model in bundled_models().items():
            assert model["build"] == "GRCh37", model_id

    def test_every_model_names_a_citation_and_a_source(self):
        for model_id, model in bundled_models().items():
            assert model["citation"], model_id
            assert model["source"], model_id

    def test_every_model_carries_an_efo_string(self):
        # An empty string is allowed and means the id was not confirmed. A
        # missing key is not, because the UI reads it unconditionally.
        for model_id, model in bundled_models().items():
            assert isinstance(model["efo"], str), model_id

    def test_every_model_has_a_description(self):
        for model_id, model in bundled_models().items():
            assert len(model["description"]) > 40, model_id

    def test_compute_all_returns_one_result_per_bundled_model(self):
        models = bundled_models()
        results = compute_all(effect_homozygous_genotypes(models))
        assert len(results) == len(models)

    def test_every_bundled_result_is_fully_covered_by_a_complete_file(self):
        models = bundled_models()
        for result in compute_all(effect_homozygous_genotypes(models)):
            assert result["coverage"] == 1.0, result["id"]
            assert result["reliable"] is True, result["id"]
            assert result["mean_imputed"] is False, result["id"]

    def test_every_bundled_result_carries_the_mandatory_caveats(self):
        models = bundled_models()
        for result in compute_all(effect_homozygous_genotypes(models)):
            for mandatory in MANDATORY_CAVEATS:
                assert mandatory in result["caveats"], result["id"]
            assert len(result["caveats"]) == 4, result["id"]

    def test_effect_homozygous_scores_sit_above_the_reference_mean(self):
        models = bundled_models()
        for result in compute_all(effect_homozygous_genotypes(models)):
            stored_mean = float(models[result["id"]]["reference"]["mean"])
            assert result["raw_score"] > stored_mean, result["id"]

    def test_to_findings_marks_every_bundled_model_as_prs(self):
        models = bundled_models()
        findings = to_findings(compute_all(effect_homozygous_genotypes(models)))
        assert len(findings) == len(models)
        for finding in findings:
            assert finding["entity_type"] == "prs"
            assert finding["repute"] == ""

    def test_coverage_report_summarises_the_whole_bundle(self):
        models = bundled_models()
        report = prs.coverage_report(effect_homozygous_genotypes(models))
        assert report["models"] == len(models)
        assert report["reliable"] == len(models)
        assert report["mean_coverage"] == 1.0
        assert report["version"]
