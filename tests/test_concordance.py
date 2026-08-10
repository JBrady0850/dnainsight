"""Tests for backend.concordance, cross-vendor agreement between the user's own kits.

WHY THIS FILE EXISTS, AND WHY IT WAS WRITTEN FIRST
--------------------------------------------------
Concordance publishes a number about a named company. "AncestryDNA and 23andMe
disagreed at 4,102 of your positions" is a factual claim about two businesses,
and if the arithmetic behind it is wrong it is a false accusation expressed as a
statistic. That is a different class of error from a wrong trait colour, so the
tests come before the implementation and they hold four lines specifically:

  1. A palindromic disagreement is never resolved. At an A/T or C/G heterozygote
     a strand flip and a real difference are indistinguishable, so the answer is
     "cannot tell", never "vendor error".
  2. An unparsed genotype can never become evidence against a company. The
     complement table in orientation.py deliberately tolerates no-call and indel
     tokens so that flipping a whole file does not destroy them. Inheriting that
     tolerance here would let NN become a disagreement.
  3. The three conflict buckets always sum to the conflict total. Folding
     indeterminate into genuine overstates vendor disagreement; folding it into
     artifact hides real disagreement.
  4. Every rate travels with the denominator it was computed over, and a pair
     with nothing shared reports None rather than 0.0 or 100.0.

Fixtures are built through merge.merge_sources rather than by hand wherever
possible, so the shapes these tests assert against are the shapes production
actually produces.
"""

import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.concordance import (
    AGREEMENT,
    CONFLICT_CLASSES,
    GENUINE,
    INDETERMINATE,
    NOT_COMPARABLE,
    ORIENTATION_ARTIFACT,
    analyse,
    classify_conflict,
    complement_genotype,
    coverage_by_provider,
    is_palindromic,
    pair_matrix,
    provider_by_label,
)
from backend.merge import merge_sources


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def snp(rsid, allele1, allele2, chromosome="1", position=1000):
    """One raw genotype row as parsers.py hands it over."""
    return {
        "rsid": rsid,
        "chromosome": chromosome,
        "position": position,
        "allele1": allele1,
        "allele2": allele2,
    }


def source(label, snps, provider="23andMe", role="self"):
    """One upload source."""
    return {"label": label, "role": role, "provider": provider, "snps": snps}


def two_kits():
    """Two pooled kits from different companies over four shared positions.

    rs1 agrees, rs2 is an orientation artifact (AG against CT), rs3 is a genuine
    disagreement (AA against GG) and rs4 is indeterminate (AT against AA).
    """
    a = source("kitA", [
        snp("rs1", "A", "G"),
        snp("rs2", "A", "G"),
        snp("rs3", "A", "A"),
        snp("rs4", "A", "T"),
    ], provider="23andMe")
    b = source("kitB", [
        snp("rs1", "A", "G"),
        snp("rs2", "C", "T"),
        snp("rs3", "G", "G"),
        snp("rs4", "A", "A"),
    ], provider="AncestryDNA")
    return merge_sources([a, b])


# ---------------------------------------------------------------------------
# complement_genotype
# ---------------------------------------------------------------------------

class TestComplementGenotype:
    def test_ag_complements_to_ct(self):
        assert complement_genotype("AG") == "CT"

    def test_the_complement_comes_back_in_canonical_order(self):
        # "GA" and "AG" are the same unordered call, so both must complement to
        # the same canonical string, otherwise a comparison against a sorted
        # merge key would miss.
        assert complement_genotype("GA") == "CT"

    def test_a_palindromic_pair_complements_to_itself(self):
        assert complement_genotype("AT") == "AT"
        assert complement_genotype("CG") == "CG"

    def test_a_no_call_is_refused_rather_than_complemented(self):
        assert complement_genotype("NN") is None

    def test_an_indel_token_is_refused_rather_than_inherited(self):
        # orientation.COMPLEMENT maps D to D and I to I on purpose, so that
        # flipping a whole file leaves indel rows intact. Validating against
        # that table instead of against ACGT would let "DD" through as a real
        # call and make it available as evidence of vendor disagreement.
        assert complement_genotype("DD") is None
        assert complement_genotype("II") is None

    def test_a_malformed_genotype_is_refused(self):
        for value in ("A", "", "AGG", None, "--", 17):
            assert complement_genotype(value) is None, value


# ---------------------------------------------------------------------------
# is_palindromic
# ---------------------------------------------------------------------------

class TestIsPalindromic:
    def test_at_is_palindromic(self):
        assert is_palindromic("AT") is True

    def test_allele_order_does_not_matter(self):
        assert is_palindromic("TA") is True

    def test_cg_is_palindromic(self):
        assert is_palindromic("CG") is True
        assert is_palindromic("GC") is True

    def test_a_homozygote_is_not_palindromic(self):
        # AA and TT are distinguishable from each other, so a homozygous call at
        # an A/T site is not the irreducible case. Only the heterozygote is.
        assert is_palindromic("AA") is False
        assert is_palindromic("TT") is False

    def test_a_non_palindromic_heterozygote_and_a_no_call_are_not_palindromic(self):
        assert is_palindromic("AG") is False
        assert is_palindromic("NN") is False


# ---------------------------------------------------------------------------
# classify_conflict
# ---------------------------------------------------------------------------

class TestClassifyConflict:
    def test_identical_genotypes_are_agreement(self):
        assert classify_conflict("AG", "AG")["classification"] == AGREEMENT

    def test_allele_order_does_not_create_a_conflict(self):
        assert classify_conflict("AG", "GA")["classification"] == AGREEMENT

    def test_ag_versus_ct_is_an_orientation_artifact(self):
        assert classify_conflict("AG", "CT")["classification"] == ORIENTATION_ARTIFACT

    def test_aa_versus_tt_is_an_orientation_artifact(self):
        # The classic whole-file strand flip. Neither call is an irreducible
        # heterozygote, so the complement explains it exactly.
        assert classify_conflict("AA", "TT")["classification"] == ORIENTATION_ARTIFACT

    def test_aa_versus_gg_is_genuine(self):
        assert classify_conflict("AA", "GG")["classification"] == GENUINE

    def test_ac_versus_ag_is_genuine(self):
        assert classify_conflict("AC", "AG")["classification"] == GENUINE

    def test_a_palindromic_call_against_a_homozygote_is_indeterminate(self):
        assert classify_conflict("AT", "AA")["classification"] == INDETERMINATE

    def test_a_palindromic_disagreement_is_never_resolved_either_way(self):
        # The load-bearing rule. Whatever it is compared against, a call whose
        # strand cannot be recovered may not be turned into a verdict.
        for partner in ("AA", "TT", "AC", "GG", "CT", "TG"):
            verdict = classify_conflict("AT", partner)["classification"]
            assert verdict != GENUINE, partner
            assert verdict != ORIENTATION_ARTIFACT, partner

    def test_an_unparsed_genotype_is_never_evidence_against_a_company(self):
        for pair in (("NN", "AA"), ("AA", "NN"), ("--", "AG"), ("A", "AG"), ("", "AA")):
            verdict = classify_conflict(*pair)["classification"]
            assert verdict == NOT_COMPARABLE, pair
            assert verdict not in CONFLICT_CLASSES, pair

    def test_an_indel_token_is_not_comparable(self):
        assert classify_conflict("DD", "AA")["classification"] == NOT_COMPARABLE
        assert classify_conflict("II", "AG")["classification"] == NOT_COMPARABLE

    def test_the_result_carries_both_genotypes_and_exactly_one_class(self):
        result = classify_conflict("AG", "CT")
        assert result["a"] == "AG"
        assert result["b"] == "CT"
        assert result["classification"] in (
            AGREEMENT, GENUINE, ORIENTATION_ARTIFACT, INDETERMINATE, NOT_COMPARABLE
        )
        assert result["reason"]


# ---------------------------------------------------------------------------
# provider_by_label
# ---------------------------------------------------------------------------

class TestProviderByLabel:
    def test_each_label_maps_to_its_declared_provider(self):
        mapping = provider_by_label(two_kits())
        assert mapping["kitA"] == "23andMe"
        assert mapping["kitB"] == "AncestryDNA"

    def test_an_undeclared_provider_is_empty_rather_than_guessed(self):
        merged = merge_sources([
            source("kitA", [snp("rs1", "A", "G")], provider=""),
            source("kitB", [snp("rs1", "A", "G")], provider="AncestryDNA"),
        ])
        mapping = provider_by_label(merged)
        assert mapping["kitA"] == ""
        assert mapping["kitB"] == "AncestryDNA"

    def test_comparison_sources_appear_in_the_map_too(self):
        merged = merge_sources([
            source("me", [snp("rs1", "A", "G")], provider="23andMe"),
            source("mum", [snp("rs1", "A", "G")], provider="MyHeritage", role="mother"),
        ])
        mapping = provider_by_label(merged)
        assert mapping["mum"] == "MyHeritage"

    def test_an_empty_merge_gives_an_empty_map(self):
        assert provider_by_label({}) == {}
        assert provider_by_label(None) == {}


# ---------------------------------------------------------------------------
# coverage_by_provider
# ---------------------------------------------------------------------------

class TestCoverageByProvider:
    def test_two_kits_from_one_provider_share_a_group(self):
        merged = merge_sources([
            source("kit2019", [snp("rs1", "A", "G")], provider="23andMe"),
            source("kit2024", [snp("rs2", "A", "G")], provider="23andMe"),
        ])
        groups = coverage_by_provider(merged)
        assert len(groups) == 1
        assert sorted(groups[0]["labels"]) == ["kit2019", "kit2024"]

    def test_two_providers_give_two_groups_with_their_own_counts(self):
        groups = {g["provider"]: g for g in coverage_by_provider(two_kits())}
        assert groups["23andMe"]["positions"] == 4
        assert groups["AncestryDNA"]["positions"] == 4

    def test_an_undeclared_provider_gets_its_own_group_per_kit(self):
        # Two kits with no declared provider are not evidence of one company.
        merged = merge_sources([
            source("kitA", [snp("rs1", "A", "G")], provider=""),
            source("kitB", [snp("rs1", "A", "G")], provider=""),
        ])
        groups = coverage_by_provider(merged)
        assert len(groups) == 2
        assert all(group["declared"] is False for group in groups)

    def test_the_groups_are_returned_in_a_deterministic_order(self):
        merged = two_kits()
        assert coverage_by_provider(merged) == coverage_by_provider(merged)
        assert [g["key"] for g in coverage_by_provider(merged)] == ["23andMe", "AncestryDNA"]


# ---------------------------------------------------------------------------
# pair_matrix
# ---------------------------------------------------------------------------

class TestPairMatrix:
    def test_two_kits_give_one_pair(self):
        pairs = pair_matrix(two_kits())
        assert len(pairs) == 1
        assert {pairs[0]["a"], pairs[0]["b"]} == {"kitA", "kitB"}

    def test_three_kits_give_three_pairs(self):
        merged = merge_sources([
            source("kitA", [snp("rs1", "A", "G")], provider="23andMe"),
            source("kitB", [snp("rs1", "A", "G")], provider="AncestryDNA"),
            source("kitC", [snp("rs1", "A", "G")], provider="MyHeritage"),
        ])
        assert len(pair_matrix(merged)) == 3

    def test_every_rate_travels_with_its_shared_denominator(self):
        pair = pair_matrix(two_kits())[0]
        assert pair["shared"] == 4
        for key in ("agreement_rate", "conflict_rate", "genuine_conflict_rate"):
            assert key in pair
        assert pair["agreement_rate"] == 0.25
        assert pair["comparable"] is True

    def test_a_pair_with_no_shared_positions_reports_none_not_zero(self):
        merged = merge_sources([
            source("kitA", [snp("rs1", "A", "G")], provider="23andMe"),
            source("kitB", [snp("rs2", "A", "G")], provider="AncestryDNA"),
        ])
        pair = pair_matrix(merged)[0]
        assert pair["shared"] == 0
        assert pair["comparable"] is False
        assert pair["agreement_rate"] is None
        assert pair["conflict_rate"] is None
        assert pair["genuine_conflict_rate"] is None

    def test_the_three_buckets_sum_to_the_conflict_total(self):
        merged = two_kits()
        # An unparsed call injected directly into the merge result, which is the
        # only way to get one this far: merge_sources drops no-calls before they
        # reach a call list. It must not become evidence of disagreement.
        merged["genotypes"]["rs5"] = {
            "rsid": "rs5", "chromosome": "1", "position": 5000,
            "allele1": "A", "allele2": "A", "genotype": "AA", "count": 2,
            "labels": ["kitA", "kitB"], "conflict": False,
            "calls": [
                {"label": "kitA", "allele1": "A", "allele2": "A", "genotype": "AA"},
                {"label": "kitB", "allele1": "N", "allele2": "N", "genotype": "NN"},
            ],
        }
        pair = pair_matrix(merged)[0]
        assert pair["genuine"] + pair["orientation_artifact"] + pair["indeterminate"] == pair["conflicts"]
        assert pair["genuine"] == 1
        assert pair["orientation_artifact"] == 1
        assert pair["indeterminate"] == 1
        assert pair["conflicts"] == 3
        assert pair["shared"] == 4          # rs5 excluded, not counted against kitB
        assert pair["not_comparable"] == 1

    def test_same_provider_pairs_are_compared_and_flagged_not_dropped(self):
        merged = merge_sources([
            source("kit2019", [snp("rs1", "A", "G")], provider="23andMe"),
            source("kit2024", [snp("rs1", "A", "A")], provider="23andMe"),
        ])
        pairs = pair_matrix(merged)
        assert len(pairs) == 1
        assert pairs[0]["same_provider"] is True
        assert pairs[0]["shared"] == 1

    def test_different_providers_are_not_flagged_same_provider(self):
        assert pair_matrix(two_kits())[0]["same_provider"] is False

    def test_only_the_users_own_pooled_kits_are_compared(self):
        merged = merge_sources([
            source("kitA", [snp("rs1", "A", "G")], provider="23andMe"),
            source("kitB", [snp("rs1", "A", "G")], provider="AncestryDNA"),
            source("mum", [snp("rs1", "A", "G")], provider="MyHeritage", role="mother"),
        ])
        pairs = pair_matrix(merged)
        assert len(pairs) == 1
        assert "mum" not in {pairs[0]["a"], pairs[0]["b"]}


# ---------------------------------------------------------------------------
# analyse
# ---------------------------------------------------------------------------

class TestAnalyse:
    def test_a_single_kit_is_an_absent_comparison_not_a_failed_one(self):
        merged = merge_sources([source("kitA", [snp("rs1", "A", "G")])])
        result = analyse(merged)
        assert result["available"] is False
        assert result["not_attempted"] is True
        assert result["pairs"] == []
        assert result["totals"] is None
        assert result["reason"]

    def test_two_kits_produce_an_available_result(self):
        result = analyse(two_kits())
        assert result["available"] is True
        assert result.get("not_attempted") is not True
        assert result["kits"] == 2
        assert len(result["pairs"]) == 1

    def test_findings_covered_is_none_when_no_findings_were_supplied(self):
        # None means nobody asked. Zero means someone asked and the answer was
        # none, which is a different statement.
        result = analyse(two_kits())
        assert result["findings_covered"] is None
        assert result["findings_total"] is None

    def test_findings_covered_is_zero_when_findings_were_supplied_and_none_are_covered(self):
        result = analyse(two_kits(), findings=[{"rsid": "rs999"}])
        assert result["findings_total"] == 1
        assert result["findings_covered"] == 0

    def test_analyse_never_mutates_the_merged_set(self):
        merged = two_kits()
        before = copy.deepcopy(merged)
        analyse(merged, findings=[{"rsid": "rs3"}])
        assert merged == before

    def test_the_totals_agree_with_the_sum_of_the_pairs(self):
        result = analyse(two_kits())
        for key in ("shared", "agree", "conflicts", "genuine", "orientation_artifact", "indeterminate"):
            assert result["totals"][key] == sum(pair[key] for pair in result["pairs"]), key
