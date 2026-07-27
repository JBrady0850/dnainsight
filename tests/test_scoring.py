"""Tests for backend.scoring magnitude, repute and confidence.

Covers docs/API_V2.md sections 2.2, 2.3 and 6. Two rules matter more than the
arithmetic: a no-call must score exactly 0.0, and a trait or a polygenic score
must never carry a repute.
"""

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.scoring import (
    BASE_SCORES,
    CLINVAR_SIG_CODES,
    CPIC_LEVELS,
    FDA_LABEL_TIERS,
    MAGNITUDE_MAX,
    MAGNITUDE_MIN,
    NEUTRAL_TERMS,
    PROTECTIVE_TERMS,
    REVIEW_STATUS_STARS,
    RISK_TERMS,
    UNSCORED_SORT_VALUE,
    base_magnitude,
    clinvar_sig_code,
    compute_confidence,
    compute_magnitude,
    compute_repute,
    evidence_label,
    normalize_cpic_level,
    review_stars,
    score_all,
    score_finding,
    sort_key_magnitude,
)

# The published documentation and the shipped data disagree about one review
# status string. Both spellings have to resolve, because a scan reads the data
# and a reader reads the docs.
DOCS_SPELLING = "no classification for the individual variant"
DATA_SPELLING = "no classification for the single variant"


class TestReviewStars:
    def test_every_documented_status_maps_to_its_star_count(self):
        for status, stars in REVIEW_STATUS_STARS.items():
            assert review_stars(status) == stars

    def test_review_status_stars_is_not_empty(self):
        assert len(REVIEW_STATUS_STARS) >= 16

    def test_practice_guideline_is_four_stars(self):
        assert review_stars("practice guideline") == 4

    def test_reviewed_by_expert_panel_is_three_stars(self):
        assert review_stars("reviewed by expert panel") == 3

    def test_multiple_submitters_no_conflicts_is_two_stars(self):
        assert review_stars(
            "criteria provided, multiple submitters, no conflicts") == 2

    def test_somatic_multiple_submitters_without_suffix_is_two_stars(self):
        assert review_stars("criteria provided, multiple submitters") == 2

    def test_conflicting_classifications_is_one_star(self):
        assert review_stars("criteria provided, conflicting classifications") == 1

    def test_conflicting_interpretations_is_one_star(self):
        assert review_stars("criteria provided, conflicting interpretations") == 1

    def test_single_submitter_is_one_star(self):
        assert review_stars("criteria provided, single submitter") == 1

    def test_no_assertion_criteria_provided_is_zero(self):
        assert review_stars("no assertion criteria provided") == 0

    def test_no_classification_provided_is_zero(self):
        assert review_stars("no classification provided") == 0

    def test_no_assertion_provided_is_zero(self):
        assert review_stars("no assertion provided") == 0

    def test_no_classifications_from_unflagged_records_is_zero(self):
        assert review_stars("no classifications from unflagged records") == 0

    def test_flagged_submission_is_zero(self):
        assert review_stars("flagged submission") == 0

    def test_bare_dash_is_zero(self):
        assert review_stars("-") == 0

    def test_empty_string_is_zero(self):
        assert review_stars("") == 0

    def test_data_spelling_of_the_discrepancy_is_zero(self):
        assert review_stars(DATA_SPELLING) == 0

    def test_docs_spelling_of_the_discrepancy_is_zero(self):
        assert review_stars(DOCS_SPELLING) == 0

    def test_both_spellings_of_the_discrepancy_are_known_keys(self):
        assert DATA_SPELLING in REVIEW_STATUS_STARS
        assert DOCS_SPELLING in REVIEW_STATUS_STARS

    def test_both_spellings_agree_with_each_other(self):
        assert review_stars(DATA_SPELLING) == review_stars(DOCS_SPELLING)

    def test_unknown_status_returns_zero_rather_than_raising(self):
        assert review_stars("reviewed by a passing stranger") == 0

    def test_future_clinvar_status_degrades_to_zero(self):
        assert review_stars("criteria provided, quantum submitters") == 0

    def test_none_returns_zero(self):
        assert review_stars(None) == 0

    def test_non_string_input_returns_zero(self):
        assert review_stars(17) == 0

    def test_case_is_ignored(self):
        assert review_stars("PRACTICE GUIDELINE") == 4

    def test_surrounding_whitespace_is_ignored(self):
        assert review_stars("  reviewed by expert panel  ") == 3

    def test_every_value_is_within_zero_to_four(self):
        for stars in REVIEW_STATUS_STARS.values():
            assert 0 <= stars <= 4


class TestClinvarSigCode:
    def test_every_documented_string_maps_to_its_code(self):
        for text, code in CLINVAR_SIG_CODES.items():
            assert clinvar_sig_code(text) == code

    def test_pathogenic_is_five(self):
        assert clinvar_sig_code("pathogenic") == 5

    def test_likely_pathogenic_is_four(self):
        assert clinvar_sig_code("likely pathogenic") == 4

    def test_pathogenic_likely_pathogenic_is_five(self):
        assert clinvar_sig_code("pathogenic/likely pathogenic") == 5

    def test_low_penetrance_pathogenic_stays_five(self):
        assert clinvar_sig_code("pathogenic, low penetrance") == 5

    def test_low_penetrance_likely_pathogenic_stays_four(self):
        assert clinvar_sig_code("likely pathogenic, low penetrance") == 4

    def test_likely_benign_is_three(self):
        assert clinvar_sig_code("likely benign") == 3

    def test_benign_is_two(self):
        assert clinvar_sig_code("benign") == 2

    def test_benign_likely_benign_is_two(self):
        assert clinvar_sig_code("benign/likely benign") == 2

    def test_uncertain_significance_is_one(self):
        assert clinvar_sig_code("uncertain significance") == 1

    def test_uncertain_risk_allele_is_one(self):
        assert clinvar_sig_code("uncertain risk allele") == 1

    def test_drug_response_is_six(self):
        assert clinvar_sig_code("drug response") == 6

    def test_histocompatibility_is_seven(self):
        assert clinvar_sig_code("histocompatibility") == 7

    def test_risk_factor_is_two_five_five(self):
        assert clinvar_sig_code("risk factor") == 255

    def test_protective_is_two_five_five(self):
        assert clinvar_sig_code("protective") == 255

    def test_not_provided_is_two_five_five(self):
        assert clinvar_sig_code("not provided") == 255

    def test_compound_string_takes_the_strongest_component(self):
        assert clinvar_sig_code("pathogenic/likely pathogenic, risk factor") == 5

    def test_compound_likely_pathogenic_with_risk_factor_is_four(self):
        assert clinvar_sig_code("likely pathogenic, risk factor") == 4

    def test_compound_pathogenic_with_drug_response_is_five(self):
        assert clinvar_sig_code("pathogenic, drug response") == 5

    def test_empty_input_is_none(self):
        assert clinvar_sig_code("") is None

    def test_none_input_is_none(self):
        assert clinvar_sig_code(None) is None

    def test_whitespace_only_input_is_none(self):
        assert clinvar_sig_code("   ") is None

    def test_unrecognised_string_is_two_five_five(self):
        assert clinvar_sig_code("association not found") == 255

    def test_junk_string_is_two_five_five(self):
        assert clinvar_sig_code("banana") == 255

    def test_case_and_whitespace_are_ignored(self):
        assert clinvar_sig_code("  PATHOGENIC ") == 5

    def test_every_code_is_documented_in_the_contract(self):
        documented = {1, 2, 3, 4, 5, 6, 7, 255}
        for code in CLINVAR_SIG_CODES.values():
            assert code in documented

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   clinvar_sig_code matches 'pathogenic' as a substring of 'pathogenicity',
    #   so a conflicting ClinVar record is coded 5 and is swept into the default
    #   clinvar_only={5,4} filter. A conflicting classification is none of the
    #   documented tiers, so it should code 255 (other).
    def test_conflicting_classifications_is_not_coded_pathogenic(self):
        assert clinvar_sig_code(
            "conflicting classifications of pathogenicity") == 255

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   Same substring fault as the classifications spelling: 'conflicting
    #   interpretations of pathogenicity' contains 'pathogenic' and is coded 5
    #   instead of 255.
    def test_conflicting_interpretations_is_not_coded_pathogenic(self):
        assert clinvar_sig_code(
            "conflicting interpretations of pathogenicity") == 255


class TestNormalizeCpicLevel:
    def test_cpic_levels_has_exactly_eight_entries(self):
        assert len(CPIC_LEVELS) == 8

    def test_cpic_levels_contains_the_three_splits(self):
        assert "A/B" in CPIC_LEVELS
        assert "B/C" in CPIC_LEVELS
        assert "C/D" in CPIC_LEVELS

    def test_cpic_levels_contains_the_non_letter_retired(self):
        assert "Retired" in CPIC_LEVELS

    def test_every_documented_level_round_trips(self):
        for level in CPIC_LEVELS:
            assert normalize_cpic_level(level) == level

    def test_plain_letters_are_accepted(self):
        assert normalize_cpic_level("A") == "A"
        assert normalize_cpic_level("D") == "D"

    def test_split_a_b_is_accepted(self):
        assert normalize_cpic_level("A/B") == "A/B"

    def test_split_b_c_is_accepted(self):
        assert normalize_cpic_level("B/C") == "B/C"

    def test_split_c_d_is_accepted(self):
        assert normalize_cpic_level("C/D") == "C/D"

    def test_retired_is_accepted(self):
        assert normalize_cpic_level("Retired") == "Retired"

    def test_lowercase_input_is_normalised(self):
        assert normalize_cpic_level("a/b") == "A/B"
        assert normalize_cpic_level("retired") == "Retired"

    def test_uppercase_retired_is_normalised(self):
        assert normalize_cpic_level("RETIRED") == "Retired"

    def test_surrounding_whitespace_is_stripped(self):
        assert normalize_cpic_level("  B/C  ") == "B/C"

    def test_internal_space_in_a_split_is_tolerated(self):
        assert normalize_cpic_level("A / B") == "A/B"

    def test_junk_returns_empty_string(self):
        assert normalize_cpic_level("Z") == ""

    def test_letters_without_the_slash_are_not_a_level(self):
        assert normalize_cpic_level("AB") == ""

    def test_none_and_empty_return_empty_string(self):
        assert normalize_cpic_level(None) == ""
        assert normalize_cpic_level("") == ""
        assert normalize_cpic_level("   ") == ""

    def test_the_splits_are_not_collapsed_onto_a_single_letter(self):
        # A four-value enum would silently fold these onto A, B or C and drop
        # every real split pair.
        for split in ("A/B", "B/C", "C/D"):
            assert normalize_cpic_level(split) == split
            assert len(normalize_cpic_level(split)) == 3


class TestFdaLabelTiers:
    def test_six_tiers_are_published_not_four(self):
        assert len(FDA_LABEL_TIERS) == 6

    def test_the_two_tiers_added_in_2024_are_present(self):
        assert "No Clinical PGx" in FDA_LABEL_TIERS
        assert "Criteria Not Met" in FDA_LABEL_TIERS

    def test_the_actionable_tiers_are_present(self):
        assert "Testing Required" in FDA_LABEL_TIERS
        assert "Testing Recommended" in FDA_LABEL_TIERS


# Minimal findings that reach exactly one base tier each.
BASE_FIXTURES = {
    "cpic_a": {"cpic_level": "A"},
    "clinvar_path_3star": {"clinical_sig": "pathogenic", "review_stars": 3},
    "cpic_b": {"cpic_level": "B"},
    "clinvar_path_2star": {"clinical_sig": "pathogenic", "review_stars": 2},
    "fda_testing": {"pgx_level": "Testing Required"},
    "clinvar_lp_2star": {"clinical_sig": "likely pathogenic", "review_stars": 2},
    "gwas_replicated": {"gwas_studies": 3},
    "clinvar_single": {"clinical_sig": "uncertain significance", "review_stars": 1},
    "default": {},
}


class TestBaseMagnitude:
    def test_every_base_scores_tier_is_reachable(self):
        for key, finding in BASE_FIXTURES.items():
            base, evidence_key = base_magnitude(dict(finding))
            assert evidence_key == key
            assert base == BASE_SCORES[key]

    def test_base_scores_has_the_nine_documented_tiers(self):
        assert set(BASE_SCORES) == set(BASE_FIXTURES)

    def test_cpic_a_scores_six(self):
        assert base_magnitude({"cpic_level": "A"}) == (6.0, "cpic_a")

    def test_clinvar_pathogenic_three_star_scores_six(self):
        finding = {"clinical_sig": "pathogenic", "review_stars": 3}
        assert base_magnitude(finding) == (6.0, "clinvar_path_3star")

    def test_clinvar_likely_pathogenic_three_star_also_scores_six(self):
        finding = {"clinical_sig": "likely pathogenic", "review_stars": 4}
        assert base_magnitude(finding) == (6.0, "clinvar_path_3star")

    def test_cpic_b_scores_four_and_a_half(self):
        assert base_magnitude({"cpic_level": "B"}) == (4.5, "cpic_b")

    def test_cpic_a_b_split_scores_as_b(self):
        assert base_magnitude({"cpic_level": "A/B"}) == (4.5, "cpic_b")

    def test_clinvar_pathogenic_two_star_scores_four_and_a_half(self):
        finding = {"clinical_sig": "pathogenic", "review_stars": 2}
        assert base_magnitude(finding) == (4.5, "clinvar_path_2star")

    def test_fda_testing_required_scores_four(self):
        finding = {"pgx_level": "Testing Required"}
        assert base_magnitude(finding) == (4.0, "fda_testing")

    def test_fda_testing_recommended_scores_four(self):
        finding = {"pgx_level": "Testing Recommended"}
        assert base_magnitude(finding) == (4.0, "fda_testing")

    def test_fda_informative_tier_does_not_reach_the_testing_base(self):
        finding = {"pgx_level": "Informative PGx"}
        assert base_magnitude(finding) == (1.0, "default")

    def test_clinvar_likely_pathogenic_two_star_scores_three_and_a_half(self):
        finding = {"clinical_sig": "likely pathogenic", "review_stars": 2}
        assert base_magnitude(finding) == (3.5, "clinvar_lp_2star")

    def test_replicated_gwas_scores_two_and_a_half(self):
        assert base_magnitude({"gwas_studies": 2}) == (2.5, "gwas_replicated")

    def test_a_gwas_source_string_counts_as_replicated(self):
        finding = {"sources": ["gwas_catalog"]}
        assert base_magnitude(finding) == (2.5, "gwas_replicated")

    def test_a_single_gwas_study_is_not_replicated(self):
        assert base_magnitude({"gwas_studies": 1}) == (1.0, "default")

    def test_clinvar_single_submitter_scores_one_and_a_half(self):
        finding = {"clinical_sig": "pathogenic", "review_stars": 1}
        assert base_magnitude(finding) == (1.5, "clinvar_single")

    def test_uncertain_significance_scores_one_and_a_half(self):
        finding = {"clinical_sig": "uncertain significance", "review_stars": 0}
        assert base_magnitude(finding) == (1.5, "clinvar_single")

    def test_no_evidence_at_all_scores_one(self):
        assert base_magnitude({}) == (1.0, "default")

    def test_review_status_string_is_used_when_stars_are_absent(self):
        finding = {"clinical_sig": "pathogenic",
                   "review_status": "reviewed by expert panel"}
        assert base_magnitude(finding) == (6.0, "clinvar_path_3star")

    def test_a_precomputed_sig_code_is_honoured(self):
        finding = {"clinvar_sig_code": 5, "review_stars": 3}
        assert base_magnitude(finding) == (6.0, "clinvar_path_3star")


class TestStrongestEvidenceWins:
    def test_cpic_a_beats_a_weak_clinvar_record(self):
        finding = {"cpic_level": "A", "clinical_sig": "uncertain significance",
                   "review_stars": 1}
        base, key = base_magnitude(finding)
        assert key == "cpic_a"
        assert base == 6.0

    def test_cpic_a_beats_an_fda_testing_label(self):
        finding = {"cpic_level": "A", "pgx_level": "Testing Required"}
        assert base_magnitude(finding)[1] == "cpic_a"

    def test_three_star_pathogenic_beats_cpic_b(self):
        finding = {"cpic_level": "B", "clinical_sig": "pathogenic",
                   "review_stars": 3}
        assert base_magnitude(finding)[1] == "clinvar_path_3star"

    def test_cpic_b_beats_a_two_star_likely_pathogenic(self):
        finding = {"cpic_level": "B", "clinical_sig": "likely pathogenic",
                   "review_stars": 2}
        assert base_magnitude(finding)[1] == "cpic_b"

    def test_fda_testing_beats_a_replicated_gwas(self):
        finding = {"pgx_level": "Testing Required", "gwas_studies": 9}
        assert base_magnitude(finding)[1] == "fda_testing"

    def test_gwas_beats_a_single_submitter_record(self):
        finding = {"gwas_studies": 4, "clinical_sig": "uncertain significance",
                   "review_stars": 1}
        assert base_magnitude(finding)[1] == "gwas_replicated"

    def test_a_cpic_a_variant_with_a_weak_record_scores_on_the_assignment(self):
        strong = {"cpic_level": "A", "clinical_sig": "uncertain significance",
                  "review_stars": 1, "variant_copies": 1}
        weak = {"clinical_sig": "uncertain significance", "review_stars": 1,
                "variant_copies": 1}
        assert compute_magnitude(strong)["magnitude"] == 6.0
        assert compute_magnitude(weak)["magnitude"] == 1.5


class TestComputeMagnitudeNoCall:
    def test_no_call_forces_magnitude_to_exactly_zero(self):
        result = compute_magnitude({"zygosity": "no_call"})
        assert result["magnitude"] == 0.0

    def test_no_call_sets_dubious_true(self):
        result = compute_magnitude({"zygosity": "no_call"})
        assert result["dubious"] is True

    def test_no_call_overrides_the_strongest_possible_evidence(self):
        finding = {"cpic_level": "A", "zygosity": "no_call",
                   "publications": 900, "variant_copies": 2,
                   "freq_band": "very_rare"}
        assert compute_magnitude(finding)["magnitude"] == 0.0

    def test_no_call_still_reports_the_base_it_would_have_had(self):
        result = compute_magnitude({"cpic_level": "A", "zygosity": "no_call"})
        assert result["base"] == 6.0
        assert result["evidence_key"] == "cpic_a"

    def test_no_call_records_the_reason_in_the_factors(self):
        result = compute_magnitude({"zygosity": "no_call"})
        assert any("no-call" in line for line in result["factors"])

    def test_no_call_is_case_insensitive(self):
        assert compute_magnitude({"zygosity": "NO_CALL"})["magnitude"] == 0.0


class TestComputeMagnitudeCarrier:
    def test_zero_copies_multiplies_by_a_quarter(self):
        finding = {"cpic_level": "A", "variant_copies": 0}
        assert compute_magnitude(finding)["magnitude"] == 1.5

    def test_zero_copies_names_the_multiplier(self):
        finding = {"cpic_level": "A", "variant_copies": 0}
        factors = compute_magnitude(finding)["factors"]
        assert any("x0.25" in line for line in factors)

    def test_two_copies_multiply_by_one_point_three(self):
        finding = {"cpic_level": "A", "variant_copies": 2}
        assert compute_magnitude(finding)["magnitude"] == 7.8

    def test_two_copies_name_the_multiplier(self):
        finding = {"cpic_level": "A", "variant_copies": 2}
        factors = compute_magnitude(finding)["factors"]
        assert any("x1.3" in line for line in factors)

    def test_one_copy_applies_no_carrier_multiplier(self):
        finding = {"cpic_level": "A", "variant_copies": 1}
        result = compute_magnitude(finding)
        assert result["magnitude"] == 6.0
        assert result["factors"] == ["base 6.00 from cpic_a"]

    def test_unknown_copy_number_applies_no_multiplier(self):
        finding = {"cpic_level": "A", "variant_copies": None}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_a_non_carrier_of_a_weak_variant_is_not_negative(self):
        finding = {"variant_copies": 0}
        assert compute_magnitude(finding)["magnitude"] == 0.25


class TestComputeMagnitudeRarity:
    def test_very_rare_adds_half_a_point(self):
        finding = {"cpic_level": "A", "freq_band": "very_rare"}
        assert compute_magnitude(finding)["magnitude"] == 6.5

    def test_rare_adds_a_quarter_point(self):
        finding = {"cpic_level": "A", "freq_band": "rare"}
        assert compute_magnitude(finding)["magnitude"] == 6.25

    def test_majority_subtracts_half_a_point(self):
        finding = {"cpic_level": "A", "freq_band": "majority"}
        assert compute_magnitude(finding)["magnitude"] == 5.5

    def test_uncommon_band_is_neutral(self):
        finding = {"cpic_level": "A", "freq_band": "uncommon"}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_common_band_is_neutral(self):
        finding = {"cpic_level": "A", "freq_band": "common"}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_unknown_band_is_neutral(self):
        finding = {"cpic_level": "A", "freq_band": "unknown"}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_missing_band_is_neutral(self):
        assert compute_magnitude({"cpic_level": "A"})["magnitude"] == 6.0

    def test_rarity_band_is_named_in_the_factors(self):
        finding = {"cpic_level": "A", "freq_band": "very_rare"}
        factors = compute_magnitude(finding)["factors"]
        assert any("very rare" in line for line in factors)


class TestComputeMagnitudePublications:
    def test_one_publication_adds_the_log_bump(self):
        finding = {"cpic_level": "A", "publications": 1}
        assert compute_magnitude(finding)["magnitude"] == 6.15

    def test_the_bump_matches_the_documented_formula(self):
        for count in (1, 2, 5, 12, 40, 99):
            expected = 6.0 + min(1.0, math.log10(1 + count) / 2.0)
            finding = {"cpic_level": "A", "publications": count}
            assert compute_magnitude(finding)["magnitude"] == round(expected, 2)

    def test_the_bump_saturates_at_one_point(self):
        finding = {"cpic_level": "A", "publications": 100}
        assert compute_magnitude(finding)["magnitude"] == 7.0

    def test_ten_thousand_publications_add_no_more_than_one(self):
        finding = {"cpic_level": "A", "publications": 10000}
        assert compute_magnitude(finding)["magnitude"] == 7.0

    def test_saturation_means_a_huge_count_ties_a_moderate_one(self):
        many = compute_magnitude({"cpic_level": "A", "publications": 10000})
        some = compute_magnitude({"cpic_level": "A", "publications": 100})
        assert many["magnitude"] == some["magnitude"]

    def test_zero_publications_add_nothing(self):
        finding = {"cpic_level": "A", "publications": 0}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_missing_publications_add_nothing(self):
        assert compute_magnitude({"cpic_level": "A"})["magnitude"] == 6.0

    def test_junk_publication_count_adds_nothing(self):
        finding = {"cpic_level": "A", "publications": "lots"}
        assert compute_magnitude(finding)["magnitude"] == 6.0

    def test_a_well_studied_common_variant_cannot_outrank_an_actionable_one(self):
        common = {"publications": 100000, "freq_band": "majority"}
        actionable = {"cpic_level": "A", "variant_copies": 1}
        assert (compute_magnitude(common)["magnitude"]
                < compute_magnitude(actionable)["magnitude"])

    def test_the_publication_count_is_named_in_the_factors(self):
        finding = {"cpic_level": "A", "publications": 42}
        factors = compute_magnitude(finding)["factors"]
        assert any("42 publications" in line for line in factors)


class TestComputeMagnitudeAmbiguity:
    def test_an_ambiguous_finding_is_capped_at_two(self):
        finding = {"cpic_level": "A", "ambiguous": True}
        assert compute_magnitude(finding)["magnitude"] == 2.0

    def test_a_freq_ambiguous_finding_is_capped_at_two(self):
        finding = {"cpic_level": "A", "freq_ambiguous": True}
        assert compute_magnitude(finding)["magnitude"] == 2.0

    def test_the_cap_applies_even_to_a_cpic_a_base(self):
        finding = {"cpic_level": "A", "variant_copies": 2,
                   "publications": 900, "freq_band": "very_rare",
                   "ambiguous": True}
        assert compute_magnitude(finding)["magnitude"] == 2.0

    def test_an_ambiguous_finding_sets_dubious(self):
        finding = {"cpic_level": "A", "ambiguous": True}
        assert compute_magnitude(finding)["dubious"] is True

    def test_a_freq_ambiguous_finding_sets_dubious(self):
        finding = {"freq_ambiguous": True}
        assert compute_magnitude(finding)["dubious"] is True

    def test_an_ambiguous_finding_below_the_cap_is_not_raised(self):
        finding = {"clinical_sig": "uncertain significance", "review_stars": 1,
                   "ambiguous": True}
        assert compute_magnitude(finding)["magnitude"] == 1.5

    def test_the_cap_is_recorded_in_the_factors(self):
        finding = {"cpic_level": "A", "ambiguous": True}
        factors = compute_magnitude(finding)["factors"]
        assert any("capped" in line for line in factors)

    def test_an_uncapped_palindrome_still_records_the_warning(self):
        finding = {"ambiguous": True}
        factors = compute_magnitude(finding)["factors"]
        assert any("strand not verifiable" in line for line in factors)

    def test_an_unverifiable_call_cannot_outrank_a_verifiable_one(self):
        unverifiable = {"cpic_level": "A", "ambiguous": True}
        verifiable = {"cpic_level": "A"}
        assert (compute_magnitude(unverifiable)["magnitude"]
                < compute_magnitude(verifiable)["magnitude"])


class TestComputeMagnitudeShape:
    def test_the_result_carries_the_five_documented_keys(self):
        result = compute_magnitude({})
        assert set(result) == {"magnitude", "base", "evidence_key",
                               "factors", "dubious"}

    def test_the_score_never_exceeds_the_documented_maximum(self):
        finding = {"cpic_level": "A", "variant_copies": 2,
                   "freq_band": "very_rare", "publications": 500000}
        assert compute_magnitude(finding)["magnitude"] <= MAGNITUDE_MAX

    def test_the_score_never_falls_below_the_documented_minimum(self):
        finding = {"variant_copies": 0, "freq_band": "majority"}
        assert compute_magnitude(finding)["magnitude"] >= MAGNITUDE_MIN

    def test_a_negative_intermediate_is_clamped_to_zero(self):
        finding = {"variant_copies": 0, "freq_band": "majority"}
        assert compute_magnitude(finding)["magnitude"] == 0.0

    def test_every_tier_stays_inside_the_documented_range(self):
        for finding in BASE_FIXTURES.values():
            for band in ("very_rare", "rare", "majority", "unknown"):
                for copies in (0, 1, 2):
                    probe = dict(finding)
                    probe.update({"freq_band": band, "variant_copies": copies,
                                  "publications": 5000})
                    magnitude = compute_magnitude(probe)["magnitude"]
                    assert MAGNITUDE_MIN <= magnitude <= MAGNITUDE_MAX

    def test_factors_are_never_empty(self):
        assert compute_magnitude({})["factors"]

    def test_factors_always_name_the_base(self):
        factors = compute_magnitude({})["factors"]
        assert factors[0].startswith("base 1.00 from default")

    def test_factors_name_every_step_that_fired(self):
        finding = {"cpic_level": "A", "variant_copies": 2,
                   "freq_band": "very_rare", "publications": 100}
        factors = compute_magnitude(finding)["factors"]
        assert len(factors) == 4
        assert "base 6.00 from cpic_a" in factors[0]
        assert "x1.3" in factors[1]
        assert "very rare" in factors[2]
        assert "100 publications" in factors[3]

    def test_the_combined_adjustments_reach_the_documented_total(self):
        finding = {"cpic_level": "A", "variant_copies": 2,
                   "freq_band": "very_rare", "publications": 100}
        assert compute_magnitude(finding)["magnitude"] == 9.3

    def test_a_pre_set_dubious_flag_survives(self):
        assert compute_magnitude({"dubious": True})["dubious"] is True

    def test_a_clean_finding_is_not_dubious(self):
        assert compute_magnitude({"cpic_level": "A"})["dubious"] is False


class TestComputeRepute:
    def test_a_trait_is_always_blank(self):
        finding = {"entity_type": "trait", "summary": "increased risk of harm",
                   "clinvar_sig_code": 5}
        assert compute_repute(finding) == ""

    def test_a_polygenic_score_is_always_blank(self):
        finding = {"entity_type": "prs", "summary": "increased risk of harm",
                   "clinvar_sig_code": 5}
        assert compute_repute(finding) == ""

    def test_a_trait_with_protective_text_is_still_blank(self):
        finding = {"entity_type": "trait", "summary": "protective effect"}
        assert compute_repute(finding) == ""

    def test_a_polygenic_score_with_protective_text_is_still_blank(self):
        finding = {"entity_type": "prs", "interpretation": "reduced risk"}
        assert compute_repute(finding) == ""

    def test_the_trait_rule_ignores_case(self):
        finding = {"entity_type": "TRAIT", "summary": "increased risk"}
        assert compute_repute(finding) == ""

    def test_a_no_call_is_blank(self):
        finding = {"zygosity": "no_call", "summary": "increased risk",
                   "clinvar_sig_code": 5}
        assert compute_repute(finding) == ""

    def test_a_confirmed_non_carrier_is_blank(self):
        finding = {"variant_copies": 0, "summary": "increased risk",
                   "clinvar_sig_code": 5}
        assert compute_repute(finding) == ""

    def test_a_carrier_of_one_copy_still_gets_a_repute(self):
        finding = {"variant_copies": 1, "summary": "increased risk"}
        assert compute_repute(finding) == "Bad"

    def test_a_homozygote_still_gets_a_repute(self):
        finding = {"variant_copies": 2, "summary": "protective effect"}
        assert compute_repute(finding) == "Good"

    def test_risk_vocabulary_is_bad(self):
        assert compute_repute({"summary": "increased risk of thrombosis"}) == "Bad"

    def test_poor_metabolizer_is_bad(self):
        assert compute_repute({"interpretation": "poor metabolizer"}) == "Bad"

    def test_contraindicated_is_bad(self):
        assert compute_repute({"summary": "contraindicated in this genotype"}) == "Bad"

    def test_loss_of_function_is_bad(self):
        assert compute_repute({"summary": "loss of function allele"}) == "Bad"

    def test_protective_vocabulary_is_good(self):
        assert compute_repute({"summary": "protective against infection"}) == "Good"

    def test_normal_metabolizer_is_good(self):
        assert compute_repute({"interpretation": "normal metabolizer"}) == "Good"

    def test_reduced_risk_is_good(self):
        assert compute_repute({"summary": "reduced risk of disease"}) == "Good"

    def test_a_neutral_term_forces_blank(self):
        assert compute_repute({"summary": "eye colour prediction"}) == ""

    def test_uncertain_significance_text_is_blank(self):
        assert compute_repute({"summary": "uncertain significance"}) == ""

    def test_conflicting_text_is_blank(self):
        assert compute_repute({"summary": "conflicting submissions"}) == ""

    def test_a_neutral_term_beats_risk_vocabulary_in_the_same_field(self):
        finding = {"summary": "eye colour, increased risk"}
        assert compute_repute(finding) == ""

    def test_risk_and_protective_terms_of_equal_length_tie_to_blank(self):
        finding = {"summary": "increased risk and decreased risk"}
        assert compute_repute(finding) == ""

    def test_the_longer_phrase_wins_a_mixed_field(self):
        finding = {"summary": "avoid, normal metabolizer"}
        assert compute_repute(finding) == "Good"

    def test_genotype_specific_text_beats_position_level_text(self):
        finding = {"summary": "protective effect",
                   "interpretation": "increased risk of harm"}
        assert compute_repute(finding) == "Good"

    def test_position_level_text_is_used_when_the_summary_is_silent(self):
        finding = {"summary": "", "interpretation": "increased risk of harm"}
        assert compute_repute(finding) == "Bad"

    def test_conditions_text_is_used_when_earlier_fields_are_silent(self):
        finding = {"conditions": "malignant hyperthermia"}
        assert compute_repute(finding) == "Bad"

    def test_a_pathogenic_code_with_no_text_is_bad(self):
        assert compute_repute({"clinvar_sig_code": 5}) == "Bad"

    def test_a_likely_pathogenic_code_with_no_text_is_bad(self):
        assert compute_repute({"clinvar_sig_code": 4}) == "Bad"

    def test_a_benign_code_with_no_text_is_good(self):
        assert compute_repute({"clinvar_sig_code": 2}) == "Good"

    def test_a_likely_benign_code_with_no_text_is_good(self):
        assert compute_repute({"clinvar_sig_code": 3}) == "Good"

    def test_an_uncertain_code_with_no_text_is_blank(self):
        assert compute_repute({"clinvar_sig_code": 1}) == ""

    def test_a_bare_finding_is_blank(self):
        assert compute_repute({}) == ""

    def test_the_code_is_derived_from_the_significance_string(self):
        assert compute_repute({"clinical_sig": "pathogenic"}) == "Bad"

    def test_the_three_vocabularies_are_populated(self):
        assert RISK_TERMS and PROTECTIVE_TERMS and NEUTRAL_TERMS

    def test_the_vocabularies_are_lowercase(self):
        for term in RISK_TERMS + PROTECTIVE_TERMS + NEUTRAL_TERMS:
            assert term == term.lower()


class TestComputeConfidence:
    def test_three_stars_is_high(self):
        assert compute_confidence({"review_stars": 3}) == "high"

    def test_four_stars_is_high(self):
        assert compute_confidence({"review_stars": 4}) == "high"

    def test_cpic_a_is_high(self):
        assert compute_confidence({"cpic_level": "A"}) == "high"

    def test_two_stars_is_moderate(self):
        assert compute_confidence({"review_stars": 2}) == "moderate"

    def test_cpic_b_is_moderate(self):
        assert compute_confidence({"cpic_level": "B"}) == "moderate"

    def test_cpic_a_b_split_is_moderate(self):
        assert compute_confidence({"cpic_level": "A/B"}) == "moderate"

    def test_one_star_is_low(self):
        assert compute_confidence({"review_stars": 1}) == "low"

    def test_a_replicated_gwas_is_low(self):
        assert compute_confidence({"gwas_studies": 2}) == "low"

    def test_a_gwas_source_is_low(self):
        assert compute_confidence({"sources": ["gwas_catalog"]}) == "low"

    def test_a_bare_significance_string_is_low(self):
        assert compute_confidence({"clinical_sig": "benign"}) == "low"

    def test_no_evidence_is_none(self):
        assert compute_confidence({}) == "none"

    def test_zero_stars_alone_is_none(self):
        assert compute_confidence({"review_stars": 0}) == "none"

    def test_a_no_call_is_none_however_strong_the_evidence(self):
        finding = {"zygosity": "no_call", "cpic_level": "A", "review_stars": 4}
        assert compute_confidence(finding) == "none"

    def test_an_ambiguous_strand_is_none(self):
        finding = {"ambiguous": True, "cpic_level": "A", "review_stars": 4}
        assert compute_confidence(finding) == "none"

    def test_a_freq_ambiguous_strand_is_none(self):
        finding = {"freq_ambiguous": True, "cpic_level": "A"}
        assert compute_confidence(finding) == "none"

    def test_cpic_a_outranks_a_two_star_record(self):
        finding = {"cpic_level": "A", "review_stars": 2}
        assert compute_confidence(finding) == "high"

    def test_the_review_status_string_is_used_when_stars_are_absent(self):
        finding = {"review_status": "practice guideline"}
        assert compute_confidence(finding) == "high"

    def test_every_answer_is_one_of_the_four_documented_words(self):
        allowed = {"high", "moderate", "low", "none"}
        for finding in BASE_FIXTURES.values():
            assert compute_confidence(dict(finding)) in allowed


class TestEvidenceLabel:
    def test_a_cpic_level_is_labelled(self):
        assert evidence_label({"cpic_level": "A"}) == "CPIC Level A"

    def test_a_split_cpic_level_is_labelled(self):
        assert evidence_label({"cpic_level": "B/C"}) == "CPIC Level B/C"

    def test_a_retired_level_is_not_used_as_the_label(self):
        assert evidence_label({"cpic_level": "Retired"}) == ""

    def test_a_clinvar_record_is_labelled_with_its_stars(self):
        finding = {"clinical_sig": "pathogenic", "review_stars": 3}
        assert evidence_label(finding) == "ClinVar pathogenic, 3 stars"

    def test_one_star_is_singular(self):
        finding = {"clinical_sig": "pathogenic", "review_stars": 1}
        assert evidence_label(finding) == "ClinVar pathogenic, 1 star"

    def test_an_fda_tier_is_labelled(self):
        finding = {"pgx_level": "Actionable PGx"}
        assert evidence_label(finding) == "FDA label: Actionable PGx"

    def test_a_gwas_study_count_is_labelled(self):
        assert evidence_label({"gwas_studies": 4}) == "GWAS, 4 studies"

    def test_no_evidence_gives_an_empty_label(self):
        assert evidence_label({}) == ""


DOCUMENTED_SCORE_KEYS = (
    "magnitude", "magnitude_source", "magnitude_base", "magnitude_factors",
    "repute", "confidence", "evidence", "review_stars", "clinvar_sig_code",
    "dubious",
)


class TestScoreFinding:
    def test_it_sets_every_documented_key(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        for key in DOCUMENTED_SCORE_KEYS:
            assert key in finding

    def test_entity_type_defaults_to_snp(self):
        finding = score_finding({"rsid": "rs1"})
        assert finding["entity_type"] == "snp"

    def test_magnitude_source_is_computed_by_default(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        assert finding["magnitude_source"] == "computed"

    def test_the_magnitude_matches_compute_magnitude(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "variant_copies": 2})
        assert finding["magnitude"] == 7.8

    def test_the_base_is_recorded(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        assert finding["magnitude_base"] == 6.0

    def test_the_audit_trail_is_recorded(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        assert finding["magnitude_factors"]

    def test_review_stars_are_derived_from_the_status_string(self):
        finding = score_finding({"rsid": "rs1",
                                 "review_status": "reviewed by expert panel"})
        assert finding["review_stars"] == 3

    def test_existing_review_stars_are_kept(self):
        finding = score_finding({"rsid": "rs1", "review_stars": 2,
                                 "review_status": "practice guideline"})
        assert finding["review_stars"] == 2

    def test_the_sig_code_is_derived_from_the_significance_string(self):
        finding = score_finding({"rsid": "rs1", "clinical_sig": "pathogenic"})
        assert finding["clinvar_sig_code"] == 5

    def test_an_existing_sig_code_is_kept(self):
        finding = score_finding({"rsid": "rs1", "clinical_sig": "pathogenic",
                                 "clinvar_sig_code": 6})
        assert finding["clinvar_sig_code"] == 6

    def test_confidence_is_set(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        assert finding["confidence"] == "high"

    def test_evidence_is_set(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A"})
        assert finding["evidence"] == "CPIC Level A"

    def test_an_existing_evidence_label_is_not_overwritten(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "evidence": "hand written"})
        assert finding["evidence"] == "hand written"

    def test_dubious_is_a_bool(self):
        finding = score_finding({"rsid": "rs1"})
        assert isinstance(finding["dubious"], bool)

    def test_a_no_call_scores_zero_and_is_dubious(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "zygosity": "no_call"})
        assert finding["magnitude"] == 0.0
        assert finding["dubious"] is True

    def test_the_finding_is_scored_in_place(self):
        finding = {"rsid": "rs1", "cpic_level": "A"}
        assert score_finding(finding) is finding

    def test_a_non_dict_is_returned_unchanged(self):
        assert score_finding("not a finding") == "not a finding"


class TestScoreFindingSnpediaOverlay:
    def test_a_cached_magnitude_wins(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "snpedia_magnitude": 3.7})
        assert finding["magnitude"] == 3.7

    def test_the_cached_source_is_recorded(self):
        finding = score_finding({"rsid": "rs1", "snpedia_magnitude": 3.7})
        assert finding["magnitude_source"] == "snpedia"

    def test_the_computed_value_is_kept_in_the_audit_trail(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "snpedia_magnitude": 3.7})
        assert any("local SNPedia cache" in line
                   for line in finding["magnitude_factors"])

    def test_a_cached_magnitude_is_clamped_to_the_documented_range(self):
        finding = score_finding({"rsid": "rs1", "snpedia_magnitude": 99.0})
        assert finding["magnitude"] == 10.0

    def test_a_negative_cached_magnitude_is_clamped(self):
        finding = score_finding({"rsid": "rs1", "snpedia_magnitude": -4.0})
        assert finding["magnitude"] == 0.0

    def test_a_cached_magnitude_of_zero_is_honoured(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "snpedia_magnitude": 0.0})
        assert finding["magnitude"] == 0.0

    def test_a_cached_repute_wins(self):
        finding = score_finding({"rsid": "rs1", "snpedia_repute": "Bad"})
        assert finding["repute"] == "Bad"

    def test_a_cached_repute_never_overrides_the_trait_rule(self):
        finding = score_finding({"rsid": "t1", "entity_type": "trait",
                                 "snpedia_repute": "Bad"})
        assert finding["repute"] == ""

    def test_a_cached_repute_never_overrides_the_prs_rule(self):
        finding = score_finding({"rsid": "p1", "entity_type": "prs",
                                 "snpedia_repute": "Good"})
        assert finding["repute"] == ""

    def test_a_junk_cached_repute_falls_back_to_the_computed_one(self):
        finding = score_finding({"rsid": "rs1", "snpedia_repute": "Maybe",
                                 "summary": "increased risk"})
        assert finding["repute"] == "Bad"

    def test_prefer_snpedia_false_ignores_the_cached_magnitude(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "snpedia_magnitude": 3.7},
                                prefer_snpedia=False)
        assert finding["magnitude"] == 6.0
        assert finding["magnitude_source"] == "computed"

    def test_prefer_snpedia_false_ignores_the_cached_repute(self):
        finding = score_finding({"rsid": "rs1", "snpedia_repute": "Bad"},
                                prefer_snpedia=False)
        assert finding["repute"] == ""

    def test_a_junk_cached_magnitude_is_ignored(self):
        finding = score_finding({"rsid": "rs1", "cpic_level": "A",
                                 "snpedia_magnitude": "very high"})
        assert finding["magnitude"] == 6.0
        assert finding["magnitude_source"] == "computed"


class TestScoreAll:
    def test_every_finding_is_scored(self):
        findings = [{"rsid": "rs1", "cpic_level": "A"}, {"rsid": "rs2"}]
        score_all(findings)
        assert all("magnitude" in f for f in findings)

    def test_the_same_list_is_returned(self):
        findings = [{"rsid": "rs1"}]
        assert score_all(findings) is findings

    def test_an_empty_list_is_tolerated(self):
        assert score_all([]) == []

    def test_none_is_tolerated(self):
        assert score_all(None) is None

    def test_prefer_snpedia_is_passed_through(self):
        findings = [{"rsid": "rs1", "cpic_level": "A", "snpedia_magnitude": 2.0}]
        score_all(findings, prefer_snpedia=False)
        assert findings[0]["magnitude_source"] == "computed"


class TestSortKeyMagnitude:
    def test_the_unscored_sort_value_is_one(self):
        assert UNSCORED_SORT_VALUE == 1.0

    def test_a_null_magnitude_sorts_as_one(self):
        assert sort_key_magnitude({"magnitude": None}) == UNSCORED_SORT_VALUE

    def test_a_null_magnitude_returns_exactly_one_point_zero(self):
        assert sort_key_magnitude({"magnitude": None}) == 1.0

    def test_a_missing_magnitude_sorts_as_one(self):
        assert sort_key_magnitude({}) == 1.0

    def test_a_null_magnitude_does_not_sort_as_zero(self):
        assert sort_key_magnitude({"magnitude": None}) != 0.0

    def test_a_zero_magnitude_sorts_as_zero(self):
        assert sort_key_magnitude({"magnitude": 0.0}) == 0.0

    def test_a_real_magnitude_is_returned_unchanged(self):
        assert sort_key_magnitude({"magnitude": 7.25}) == 7.25

    def test_a_string_magnitude_is_coerced(self):
        assert sort_key_magnitude({"magnitude": "3.5"}) == 3.5

    def test_junk_sorts_as_unscored(self):
        assert sort_key_magnitude({"magnitude": "high"}) == 1.0

    def test_an_unscored_finding_sits_between_two_and_zero(self):
        rows = [{"rsid": "zero", "magnitude": 0.0},
                {"rsid": "null", "magnitude": None},
                {"rsid": "two", "magnitude": 2.0}]
        ordered = sorted(rows, key=sort_key_magnitude, reverse=True)
        assert [r["rsid"] for r in ordered] == ["two", "null", "zero"]


class TestScoringConstants:
    def test_the_documented_magnitude_bounds(self):
        assert MAGNITUDE_MIN == 0.0
        assert MAGNITUDE_MAX == 10.0

    def test_the_base_scores_match_the_contract_table(self):
        assert BASE_SCORES["cpic_a"] == 6.0
        assert BASE_SCORES["clinvar_path_3star"] == 6.0
        assert BASE_SCORES["cpic_b"] == 4.5
        assert BASE_SCORES["clinvar_path_2star"] == 4.5
        assert BASE_SCORES["fda_testing"] == 4.0
        assert BASE_SCORES["clinvar_lp_2star"] == 3.5
        assert BASE_SCORES["gwas_replicated"] == 2.5
        assert BASE_SCORES["clinvar_single"] == 1.5
        assert BASE_SCORES["default"] == 1.0

    def test_every_base_score_is_inside_the_documented_range(self):
        for base in BASE_SCORES.values():
            assert MAGNITUDE_MIN <= base <= MAGNITUDE_MAX

    def test_the_clinvar_code_table_covers_the_nine_ui_codes(self):
        for text in ("pathogenic", "likely pathogenic", "likely benign",
                     "benign", "uncertain significance", "drug response",
                     "histocompatibility", "risk factor", "not provided"):
            assert text in CLINVAR_SIG_CODES
