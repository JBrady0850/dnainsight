"""Tests for backend.orientation strand handling."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.orientation import (
    AMBIGUOUS_PAIRS,
    COMPLEMENT,
    NOCALL_ALLELES,
    candidate_tokens,
    complement_allele,
    complement_genotype,
    flip_report,
    genotype_token,
    is_ambiguous_pair,
    is_no_call,
    normalize_orientation,
    orient_to_snpedia,
    sort_alleles,
)


class TestComplementAllele:
    def test_complement_a(self):
        assert complement_allele("A") == "T"

    def test_complement_t(self):
        assert complement_allele("T") == "A"

    def test_complement_c(self):
        assert complement_allele("C") == "G"

    def test_complement_g(self):
        assert complement_allele("G") == "C"

    def test_complement_is_case_insensitive(self):
        assert complement_allele("a") == "T"
        assert complement_allele("g") == "C"

    def test_complement_no_call_tokens_are_self(self):
        assert complement_allele("N") == "N"
        assert complement_allele("-") == "-"
        assert complement_allele("D") == "D"
        assert complement_allele("I") == "I"

    def test_complement_zero_and_empty_become_n(self):
        assert complement_allele("0") == "N"
        assert complement_allele("") == "N"

    def test_complement_unknown_character_becomes_n(self):
        assert complement_allele("X") == "N"
        assert complement_allele("?") == "N"

    def test_complement_is_involutive_for_real_bases(self):
        for base in ("A", "T", "C", "G"):
            assert complement_allele(complement_allele(base)) == base


class TestComplementGenotype:
    def test_complement_pair_preserves_order(self):
        assert complement_genotype("C", "T") == ("G", "A")

    def test_complement_pair_lowercase(self):
        assert complement_genotype("a", "g") == ("T", "C")

    def test_complement_pair_with_no_call(self):
        assert complement_genotype("A", "-") == ("T", "-")


class TestIsNoCall:
    def test_real_bases_are_calls(self):
        for base in ("A", "T", "C", "G"):
            assert is_no_call(base) is False

    def test_no_call_tokens(self):
        for token in ("", "N", "-", "0", "D", "I"):
            assert is_no_call(token) is True

    def test_no_call_is_case_insensitive(self):
        assert is_no_call("n") is True

    def test_nocall_alleles_membership(self):
        assert "N" in NOCALL_ALLELES
        assert "A" not in NOCALL_ALLELES


class TestIsAmbiguousPair:
    def test_at_heterozygote_is_ambiguous(self):
        assert is_ambiguous_pair("A", "T") is True
        assert is_ambiguous_pair("T", "A") is True

    def test_cg_heterozygote_is_ambiguous(self):
        assert is_ambiguous_pair("C", "G") is True
        assert is_ambiguous_pair("G", "C") is True

    def test_homozygous_aa_is_not_ambiguous(self):
        assert is_ambiguous_pair("A", "A") is False

    def test_homozygous_cc_is_not_ambiguous(self):
        assert is_ambiguous_pair("C", "C") is False

    def test_ag_heterozygote_is_not_ambiguous(self):
        assert is_ambiguous_pair("A", "G") is False

    def test_ct_heterozygote_is_not_ambiguous(self):
        assert is_ambiguous_pair("C", "T") is False

    def test_no_calls_are_not_ambiguous(self):
        assert is_ambiguous_pair("-", "-") is False
        assert is_ambiguous_pair("A", "N") is False

    def test_ambiguous_pairs_constant(self):
        assert frozenset({"A", "T"}) in AMBIGUOUS_PAIRS
        assert frozenset({"C", "G"}) in AMBIGUOUS_PAIRS
        assert frozenset({"A", "G"}) not in AMBIGUOUS_PAIRS


class TestNormalizeOrientation:
    def test_plus_words_and_symbols(self):
        assert normalize_orientation("plus") == "plus"
        assert normalize_orientation("PLUS") == "plus"
        assert normalize_orientation("+") == "plus"
        assert normalize_orientation("1") == "plus"

    def test_minus_words_and_symbols(self):
        assert normalize_orientation("minus") == "minus"
        assert normalize_orientation("MINUS") == "minus"
        assert normalize_orientation("-") == "minus"
        assert normalize_orientation("-1") == "minus"

    def test_none_and_empty_are_unknown(self):
        assert normalize_orientation(None) == ""
        assert normalize_orientation("") == ""
        assert normalize_orientation("   ") == ""

    def test_garbage_is_unknown(self):
        assert normalize_orientation("sideways") == ""
        assert normalize_orientation("2") == ""

    def test_surrounding_whitespace_is_tolerated(self):
        assert normalize_orientation("  Minus  ") == "minus"


class TestSortAlleles:
    def test_tc_sorts_to_ct(self):
        assert sort_alleles("T", "C") == ("C", "T")

    def test_already_sorted_pair_is_unchanged(self):
        assert sort_alleles("A", "G") == ("A", "G")

    def test_homozygote_is_unchanged(self):
        assert sort_alleles("G", "G") == ("G", "G")

    def test_no_call_is_normalised_and_pushed_last(self):
        assert sort_alleles("-", "A") == ("A", "N")
        assert sort_alleles("A", "0") == ("A", "N")

    def test_double_no_call_becomes_double_n(self):
        assert sort_alleles("-", "-") == ("N", "N")

    def test_lowercase_input_is_uppercased(self):
        assert sort_alleles("t", "c") == ("C", "T")


class TestGenotypeToken:
    def test_token_format_for_ct(self):
        assert genotype_token("C", "T") == "(C;T)"

    def test_token_is_sorted(self):
        assert genotype_token("T", "C") == "(C;T)"

    def test_token_for_homozygote(self):
        assert genotype_token("A", "A") == "(A;A)"

    def test_no_call_pair_token(self):
        assert genotype_token("-", "-") == "(-;-)"
        assert genotype_token("", "") == "(-;-)"

    def test_partial_no_call_token(self):
        assert genotype_token("A", "-") == "(A;N)"


class TestOrientToSnpedia:
    def test_returns_exactly_eight_keys(self):
        result = orient_to_snpedia("C", "T", "plus")
        assert set(result.keys()) == {
            "allele1",
            "allele2",
            "genotype",
            "token",
            "flipped",
            "ambiguous",
            "orientation",
            "source_genotype",
        }

    def test_rs1801133_minus_flips_ct_to_ag(self):
        result = orient_to_snpedia("C", "T", "minus")
        assert result["allele1"] == "A"
        assert result["allele2"] == "G"
        assert result["genotype"] == "AG"
        assert result["token"] == "(A;G)"
        assert result["flipped"] is True
        assert result["orientation"] == "minus"
        assert result["source_genotype"] == "CT"

    def test_rs1801133_plus_leaves_alleles_alone(self):
        result = orient_to_snpedia("C", "T", "plus")
        assert result["token"] == "(C;T)"
        assert result["flipped"] is False
        assert result["orientation"] == "plus"

    def test_unknown_orientation_does_not_flip(self):
        result = orient_to_snpedia("C", "T", None)
        assert result["token"] == "(C;T)"
        assert result["flipped"] is False
        assert result["orientation"] == ""

    def test_minus_homozygote_flips_both_alleles(self):
        result = orient_to_snpedia("T", "T", "minus")
        assert result["genotype"] == "AA"
        assert result["token"] == "(A;A)"
        assert result["flipped"] is True

    def test_result_is_sorted_after_flip(self):
        result = orient_to_snpedia("G", "A", "minus")
        assert (result["allele1"], result["allele2"]) == ("C", "T")
        assert result["token"] == "(C;T)"

    def test_falls_back_to_orientation_when_stabilized_absent(self):
        result = orient_to_snpedia("C", "T", None, "minus")
        assert result["orientation"] == "minus"
        assert result["token"] == "(A;G)"

    def test_stabilized_wins_over_orientation(self):
        result = orient_to_snpedia("C", "T", "plus", "minus")
        assert result["orientation"] == "plus"
        assert result["flipped"] is False
        assert result["token"] == "(C;T)"

    def test_empty_stabilized_falls_back_to_orientation(self):
        result = orient_to_snpedia("C", "T", "", "minus")
        assert result["orientation"] == "minus"
        assert result["flipped"] is True

    def test_ambiguous_flag_uses_input_pair(self):
        result = orient_to_snpedia("A", "T", "minus")
        assert result["ambiguous"] is True
        assert result["token"] == "(A;T)"
        assert result["source_genotype"] == "AT"

    def test_cg_heterozygote_is_flagged_ambiguous(self):
        result = orient_to_snpedia("G", "C", "minus")
        assert result["ambiguous"] is True
        assert result["token"] == "(C;G)"

    def test_unambiguous_pair_is_not_flagged(self):
        result = orient_to_snpedia("A", "G", "plus")
        assert result["ambiguous"] is False

    def test_no_call_never_becomes_a_real_base(self):
        result = orient_to_snpedia("-", "-", "minus")
        assert result["allele1"] == "N"
        assert result["allele2"] == "N"
        assert result["genotype"] == "NN"
        assert result["token"] == "(-;-)"
        assert result["flipped"] is False

    def test_partial_no_call_flips_only_the_real_base(self):
        result = orient_to_snpedia("A", "-", "minus")
        assert result["allele1"] == "T"
        assert result["allele2"] == "N"
        assert result["token"] == "(T;N)"
        assert result["flipped"] is True

    def test_lowercase_input_is_normalised(self):
        result = orient_to_snpedia("c", "t", "minus")
        assert result["token"] == "(A;G)"
        assert result["source_genotype"] == "CT"

    def test_symbol_orientation_values_are_accepted(self):
        assert orient_to_snpedia("C", "T", "-")["token"] == "(A;G)"
        assert orient_to_snpedia("C", "T", "+")["token"] == "(C;T)"


class TestCandidateTokens:
    def test_ag_returns_two_entries_unflipped_first(self):
        assert candidate_tokens("A", "G") == ["(A;G)", "(C;T)"]

    def test_unsorted_input_still_returns_two_entries(self):
        assert candidate_tokens("G", "A") == ["(A;G)", "(C;T)"]

    def test_palindromic_het_deduplicates_to_one(self):
        assert candidate_tokens("A", "T") == ["(A;T)"]
        assert candidate_tokens("C", "G") == ["(C;G)"]

    def test_homozygote_returns_two_entries(self):
        assert candidate_tokens("C", "C") == ["(C;C)", "(G;G)"]

    def test_no_call_pair_returns_empty_list(self):
        assert candidate_tokens("-", "-") == []
        assert candidate_tokens("", "") == []

    def test_partial_no_call_excludes_n_tokens(self):
        tokens = candidate_tokens("A", "-")
        assert tokens == []

    def test_never_returns_a_token_containing_n(self):
        for pair in (("A", "G"), ("C", "C"), ("A", "T"), ("A", "0"), ("D", "I")):
            for token in candidate_tokens(*pair):
                assert "N" not in token


class TestFlipReport:
    def test_empty_list_is_all_zeros(self):
        assert flip_report([]) == {
            "total": 0,
            "flipped": 0,
            "ambiguous": 0,
            "unknown_orientation": 0,
        }

    def test_aggregates_flipped_and_ambiguous(self):
        findings = [
            {"flipped": True, "ambiguous": False, "orientation": "minus"},
            {"flipped": False, "ambiguous": True, "orientation": "plus"},
            {"flipped": False, "ambiguous": False, "orientation": ""},
        ]
        assert flip_report(findings) == {
            "total": 3,
            "flipped": 1,
            "ambiguous": 1,
            "unknown_orientation": 1,
        }

    def test_missing_keys_count_as_unknown_orientation(self):
        report = flip_report([{}, {}])
        assert report["total"] == 2
        assert report["flipped"] == 0
        assert report["ambiguous"] == 0
        assert report["unknown_orientation"] == 2

    def test_counts_real_orient_results(self):
        findings = [
            orient_to_snpedia("C", "T", "minus"),
            orient_to_snpedia("A", "T", "plus"),
            orient_to_snpedia("A", "G", None),
        ]
        report = flip_report(findings)
        assert report["total"] == 3
        assert report["flipped"] == 1
        assert report["ambiguous"] == 1
        assert report["unknown_orientation"] == 1


class TestConstants:
    def test_complement_map_covers_documented_tokens(self):
        assert COMPLEMENT["A"] == "T"
        assert COMPLEMENT["T"] == "A"
        assert COMPLEMENT["C"] == "G"
        assert COMPLEMENT["G"] == "C"
        assert COMPLEMENT["0"] == "N"
        assert COMPLEMENT[""] == "N"

    def test_nocall_alleles_is_frozenset(self):
        assert isinstance(NOCALL_ALLELES, frozenset)
        assert len(NOCALL_ALLELES) == 6
