"""Tests for backend.merge pooling, conflict retention, roles and trio checks."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.merge import (
    NOCALL_KEY,
    PRIMARY_ROLE,
    ROLES,
    TRIO_NOTE,
    MergeError,
    comparison_rows,
    dedupe_snps,
    genotype_key,
    is_pooled_role,
    is_real_call,
    mendelian_violation,
    merge_sources,
    normalize_role,
    transmission_probability,
    trio_annotate,
)


def snp(rsid, allele1, allele2, chromosome="1", position=1000):
    """One raw genotype row as parsers.py hands it over."""
    return {
        "rsid": rsid,
        "chromosome": chromosome,
        "position": position,
        "allele1": allele1,
        "allele2": allele2,
    }


def source(label, role, snps, provider="23andMe"):
    """One upload source."""
    return {"label": label, "role": role, "provider": provider, "snps": snps}


def _trio():
    """A merged trio: rs1 and rs3 are compatible, rs2 is impossible."""
    return merge_sources([
        source("child.txt", "self", [
            snp("rs1", "A", "G"),
            snp("rs2", "A", "A"),
            snp("rs3", "C", "T"),
            snp("rs4", "A", "A"),
        ]),
        source("mum.txt", "mother", [
            snp("rs1", "A", "A"),
            snp("rs2", "G", "G"),
            snp("rs3", "C", "C"),
        ]),
        source("dad.txt", "father", [
            snp("rs1", "G", "G"),
            snp("rs2", "G", "G"),
            snp("rs3", "T", "T"),
        ]),
    ])


class TestNormalizeRole:
    def test_known_roles_pass_through(self):
        for role in ROLES:
            assert normalize_role(role) == role

    def test_matching_is_case_insensitive(self):
        assert normalize_role("SELF") == "self"
        assert normalize_role("Mother") == "mother"

    def test_surrounding_whitespace_is_tolerated(self):
        assert normalize_role("  father  ") == "father"

    def test_unknown_role_becomes_other(self):
        assert normalize_role("cousin") == "other"

    def test_empty_role_becomes_other(self):
        assert normalize_role("") == "other"
        assert normalize_role("   ") == "other"

    def test_none_becomes_other(self):
        assert normalize_role(None) == "other"

    def test_non_string_becomes_other(self):
        assert normalize_role(42) == "other"

    def test_self_aliases(self):
        for alias in ("me", "myself", "primary", "proband"):
            assert normalize_role(alias) == "self"

    def test_parent_aliases(self):
        assert normalize_role("mom") == "mother"
        assert normalize_role("mum") == "mother"
        assert normalize_role("dad") == "father"

    def test_mate_and_sibling_aliases(self):
        assert normalize_role("spouse") == "mate"
        assert normalize_role("wife") == "mate"
        assert normalize_role("brother") == "sibling"
        assert normalize_role("sister") == "sibling"

    def test_ignore_aliases(self):
        for alias in ("ignore", "IGNORE", "skip", "exclude", "excluded"):
            assert normalize_role(alias) == "ignore"

    def test_spaces_and_hyphens_are_normalised(self):
        assert normalize_role("Half-Sibling") == "other"
        assert normalize_role("my self") == "other"

    def test_is_pooled_role_is_true_only_for_self(self):
        assert is_pooled_role("self") is True
        assert is_pooled_role("SELF") is True
        assert is_pooled_role("me") is True

    def test_is_pooled_role_is_false_for_everything_else(self):
        for role in ("mother", "father", "mate", "child", "sibling", "other", "ignore", ""):
            assert is_pooled_role(role) is False

    def test_primary_role_constant(self):
        assert PRIMARY_ROLE == "self"
        assert PRIMARY_ROLE in ROLES


class TestIsRealCall:
    def test_two_bases_are_a_real_call(self):
        assert is_real_call("A", "G") is True

    def test_lowercase_bases_are_a_real_call(self):
        assert is_real_call("a", "g") is True

    def test_whitespace_is_stripped(self):
        assert is_real_call(" A ", "G ") is True

    def test_no_call_tokens_on_the_left(self):
        for token in ("", "N", "-", "--", "0", "00", "D", "I", "?", "."):
            assert is_real_call(token, "G") is False

    def test_no_call_tokens_on_the_right(self):
        for token in ("", "N", "-", "--", "0", "00", "D", "I", "?", "."):
            assert is_real_call("A", token) is False

    def test_none_is_a_no_call(self):
        assert is_real_call(None, "A") is False
        assert is_real_call("A", None) is False

    def test_double_no_call(self):
        assert is_real_call("-", "-") is False


class TestGenotypeKey:
    def test_key_is_sorted(self):
        assert genotype_key("A", "G") == "AG"
        assert genotype_key("T", "C") == "CT"

    def test_allele_order_is_insignificant(self):
        assert genotype_key("A", "G") == genotype_key("G", "A")

    def test_homozygote_key(self):
        assert genotype_key("A", "A") == "AA"

    def test_lowercase_is_normalised(self):
        assert genotype_key("g", "a") == "AG"

    def test_no_call_becomes_the_no_call_key(self):
        assert genotype_key("-", "-") == NOCALL_KEY
        assert genotype_key("A", "-") == NOCALL_KEY
        assert genotype_key("", "") == NOCALL_KEY

    def test_indel_tokens_are_a_no_call(self):
        assert genotype_key("D", "I") == NOCALL_KEY
        assert genotype_key("D", "D") == NOCALL_KEY

    def test_no_call_key_constant(self):
        assert NOCALL_KEY == "NN"


class TestDedupeSnps:
    def test_repeated_rsid_keeps_the_first_real_call(self):
        result = dedupe_snps([snp("rs1", "A", "A"), snp("rs1", "G", "G")])
        assert len(result) == 1
        assert (result[0]["allele1"], result[0]["allele2"]) == ("A", "A")

    def test_a_later_real_call_replaces_an_earlier_no_call(self):
        result = dedupe_snps([snp("rs1", "-", "-"), snp("rs1", "A", "G")])
        assert len(result) == 1
        assert (result[0]["allele1"], result[0]["allele2"]) == ("A", "G")

    def test_two_no_calls_collapse_to_one_row(self):
        result = dedupe_snps([snp("rs1", "-", "-"), snp("rs1", "0", "0")])
        assert len(result) == 1

    def test_distinct_rsids_are_kept_in_order(self):
        result = dedupe_snps([snp("rs3", "A", "A"), snp("rs1", "G", "G")])
        assert [row["rsid"] for row in result] == ["rs3", "rs1"]

    def test_rsid_matching_is_case_insensitive(self):
        result = dedupe_snps([snp("RS1", "A", "A"), snp("rs1", "G", "G")])
        assert len(result) == 1

    def test_blank_rsids_are_dropped(self):
        result = dedupe_snps([snp("", "A", "A"), snp("rs1", "G", "G")])
        assert [row["rsid"] for row in result] == ["rs1"]

    def test_non_dict_entries_are_skipped(self):
        result = dedupe_snps(["nonsense", None, snp("rs1", "A", "A")])
        assert len(result) == 1

    def test_none_input_is_an_empty_list(self):
        assert dedupe_snps(None) == []
        assert dedupe_snps([]) == []

    def test_rows_are_copied_not_aliased(self):
        rows = [snp("rs1", "A", "A")]
        result = dedupe_snps(rows)
        assert result[0] is not rows[0]
        assert result[0] == rows[0]

    def test_a_duplicate_never_becomes_a_conflict(self):
        merged = merge_sources([
            source("one.txt", "self", [snp("rs1", "A", "A"), snp("rs1", "G", "G")]),
        ])
        assert merged["conflicts"] == []
        assert merged["genotypes"]["rs1"]["conflict"] is False
        assert merged["genotypes"]["rs1"]["count"] == 1


class TestSourceValidation:
    def test_none_raises(self):
        with pytest.raises(MergeError):
            merge_sources(None)

    def test_a_string_raises(self):
        with pytest.raises(MergeError):
            merge_sources("one.txt")

    def test_a_bare_dict_raises(self):
        with pytest.raises(MergeError):
            merge_sources({"label": "one.txt", "role": "self", "snps": []})

    def test_an_empty_list_raises(self):
        with pytest.raises(MergeError):
            merge_sources([])

    def test_a_non_dict_element_raises(self):
        with pytest.raises(MergeError):
            merge_sources(["one.txt"])

    def test_only_ignored_sources_raises(self):
        with pytest.raises(MergeError):
            merge_sources([source("one.txt", "ignore", [snp("rs1", "A", "A")])])

    def test_error_message_explains_the_empty_case(self):
        with pytest.raises(MergeError) as info:
            merge_sources([source("one.txt", "skip", [snp("rs1", "A", "A")])])
        assert "ignored" in str(info.value)


class TestRoleHandling:
    def test_an_ignored_source_is_dropped_entirely(self):
        merged = merge_sources([
            source("keep.txt", "self", [snp("rs1", "A", "A")]),
            source("drop.txt", "ignore", [snp("rs2", "G", "G")]),
        ])
        assert "rs2" not in merged["genotypes"]
        assert [s["label"] for s in merged["sources"]] == ["keep.txt"]
        assert merged["comparison"] == {}

    def test_an_ignored_source_never_reaches_comparison(self):
        merged = merge_sources([
            source("keep.txt", "self", [snp("rs1", "A", "A")]),
            source("drop.txt", "excluded", [snp("rs1", "G", "G")]),
        ])
        assert merged["comparison"] == {}
        assert merged["counts"]["comparison_sources"] == 0

    def test_non_self_roles_land_in_comparison_only(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [snp("rs2", "G", "G")]),
        ])
        assert set(merged["genotypes"]) == {"rs1"}
        assert set(merged["comparison"]) == {"mother"}
        assert set(merged["comparison"]["mother"]) == {"rs2"}

    def test_comparison_is_keyed_by_role(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("dad.txt", "dad", [snp("rs1", "A", "G")]),
            source("sis.txt", "sister", [snp("rs1", "G", "G")]),
        ])
        assert set(merged["comparison"]) == {"father", "sibling"}

    def test_two_sources_with_the_same_role_get_distinct_keys(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("kid1.txt", "child", [snp("rs1", "A", "G")]),
            source("kid2.txt", "child", [snp("rs1", "G", "G")]),
        ])
        assert len(merged["comparison"]) == 2
        assert "child" in merged["comparison"]

    def test_comparison_drops_no_calls(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [snp("rs1", "-", "-"), snp("rs2", "G", "G")]),
        ])
        assert set(merged["comparison"]["mother"]) == {"rs2"}

    def test_duplicate_labels_are_disambiguated(self):
        merged = merge_sources([
            source("raw.txt", "self", [snp("rs1", "A", "A")]),
            source("raw.txt", "self", [snp("rs2", "G", "G")]),
        ])
        labels = [s["label"] for s in merged["sources"]]
        assert labels == ["raw.txt", "raw.txt (2)"]
        assert merged["primary_labels"] == labels

    def test_missing_labels_are_generated_from_position(self):
        merged = merge_sources([
            {"role": "self", "snps": [snp("rs1", "A", "A")]},
            {"role": "self", "snps": [snp("rs2", "G", "G")]},
        ])
        assert [s["label"] for s in merged["sources"]] == ["file1", "file2"]

    def test_primary_labels_lists_only_pooled_sources(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [snp("rs1", "A", "G")]),
        ])
        assert merged["primary_labels"] == ["me.txt"]

    def test_source_counts(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A")]),
            source("b.txt", "self", [snp("rs2", "G", "G")]),
            source("mum.txt", "mother", [snp("rs1", "A", "A")]),
        ])
        assert merged["counts"]["pooled_sources"] == 2
        assert merged["counts"]["comparison_sources"] == 1


class TestPooling:
    def test_two_self_sources_union_their_positions(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs2", "C", "T")]),
            source("b.txt", "self", [snp("rs2", "C", "T"), snp("rs3", "G", "G")]),
        ])
        assert set(merged["genotypes"]) == {"rs1", "rs2", "rs3"}
        assert merged["counts"]["union"] == 3
        assert merged["counts"]["total_positions"] == 3

    def test_allele_order_is_not_a_conflict(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "G")]),
            source("b.txt", "self", [snp("rs1", "G", "A")]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert entry["genotype"] == "AG"
        assert entry["conflict"] is False
        assert entry["count"] == 2
        assert merged["conflicts"] == []

    def test_a_no_call_in_one_file_is_not_a_conflict(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "-", "-")]),
            source("b.txt", "self", [snp("rs1", "A", "G")]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert entry["conflict"] is False
        assert merged["conflicts"] == []

    def test_the_real_call_wins_over_a_no_call(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "0", "0")]),
            source("b.txt", "self", [snp("rs1", "A", "G")]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert (entry["allele1"], entry["allele2"]) == ("A", "G")
        assert entry["genotype"] == "AG"

    def test_a_no_call_in_both_files_stays_a_no_call(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "-", "-")]),
            source("b.txt", "self", [snp("rs1", "N", "N")]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert entry["genotype"] == NOCALL_KEY
        assert (entry["allele1"], entry["allele2"]) == ("N", "N")
        assert entry["count"] == 0

    def test_count_only_counts_real_calls(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "G")]),
            source("b.txt", "self", [snp("rs1", "-", "-")]),
        ])
        assert merged["genotypes"]["rs1"]["count"] == 1
        assert merged["genotypes"]["rs1"]["labels"] == ["a.txt"]

    def test_a_genuine_disagreement_sets_conflict_on_the_entry(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A")]),
            source("b.txt", "self", [snp("rs1", "G", "G")]),
        ])
        assert merged["genotypes"]["rs1"]["conflict"] is True

    def test_a_genuine_disagreement_is_recorded_in_conflicts(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A", "7", 123)]),
            source("b.txt", "self", [snp("rs1", "G", "G", "7", 123)]),
        ])
        assert len(merged["conflicts"]) == 1
        conflict = merged["conflicts"][0]
        assert conflict["rsid"] == "rs1"
        assert conflict["chromosome"] == "7"
        assert conflict["position"] == 123
        assert [c["genotype"] for c in conflict["calls"]] == ["AA", "GG"]
        assert [c["label"] for c in conflict["calls"]] == ["a.txt", "b.txt"]

    def test_both_disagreeing_calls_are_retained(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A")]),
            source("b.txt", "self", [snp("rs1", "G", "G")]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert entry["genotype"] == "AA"
        assert entry["count"] == 2
        assert [c["genotype"] for c in entry["calls"]] == ["AA", "GG"]

    def test_conflict_count_is_reported(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs2", "C", "C")]),
            source("b.txt", "self", [snp("rs1", "G", "G"), snp("rs2", "C", "C")]),
        ])
        assert merged["counts"]["conflicts"] == 1
        assert merged["counts"]["shared"] == 2

    def test_output_is_sorted_by_rsid(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs9", "A", "A"), snp("rs11", "C", "C")]),
            source("b.txt", "self", [snp("rs3", "G", "G")]),
        ])
        assert list(merged["genotypes"]) == sorted(merged["genotypes"])

    def test_comparison_rows_are_sorted_by_rsid(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [snp("rs9", "A", "A"), snp("rs11", "C", "C")]),
        ])
        rows = merged["comparison"]["mother"]
        assert list(rows) == sorted(rows)

    def test_repeated_runs_produce_identical_output(self):
        sources = [
            source("a.txt", "self", [snp("rs9", "A", "G"), snp("rs2", "C", "C")]),
            source("b.txt", "self", [snp("rs2", "C", "T"), snp("rs1", "G", "G")]),
            source("mum.txt", "mother", [snp("rs2", "C", "C")]),
        ]
        assert merge_sources(sources) == merge_sources(sources)

    def test_unique_positions_are_attributed_to_their_only_source(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs2", "C", "C")]),
            source("b.txt", "self", [snp("rs2", "C", "C"), snp("rs3", "G", "G")]),
        ])
        assert merged["counts"]["unique_by_label"] == {"a.txt": 1, "b.txt": 1}

    def test_chromosome_and_position_are_backfilled(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A", "", 0)]),
            source("b.txt", "self", [snp("rs1", "A", "A", "7", 123)]),
        ])
        entry = merged["genotypes"]["rs1"]
        assert entry["chromosome"] == "7"
        assert entry["position"] == 123

    def test_a_pool_of_one_source_still_works(self):
        merged = merge_sources([source("a.txt", "self", [snp("rs1", "A", "G")])])
        assert merged["genotypes"]["rs1"]["genotype"] == "AG"
        assert merged["counts"]["shared"] == 0
        assert merged["conflicts"] == []

    def test_no_self_source_leaves_the_primary_set_empty(self):
        merged = merge_sources([source("mum.txt", "mother", [snp("rs1", "A", "G")])])
        assert merged["genotypes"] == {}
        assert merged["primary_labels"] == []
        assert set(merged["comparison"]) == {"mother"}


class TestSourceStats:
    def test_stats_carry_the_documented_keys(self):
        merged = merge_sources([source("a.txt", "self", [snp("rs1", "A", "A")])])
        assert set(merged["sources"][0]) == {
            "label", "role", "provider", "snp_count",
            "contributed", "overlapped", "conflicting",
        }

    def test_role_and_provider_are_echoed(self):
        merged = merge_sources([
            source("a.txt", "me", [snp("rs1", "A", "A")], provider="AncestryDNA"),
        ])
        assert merged["sources"][0]["role"] == "self"
        assert merged["sources"][0]["provider"] == "AncestryDNA"

    def test_contributed_plus_overlapped_equals_the_snp_count(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs2", "C", "C")]),
            source("b.txt", "self", [snp("rs2", "C", "C"), snp("rs3", "G", "G")]),
        ])
        for stat in merged["sources"]:
            assert stat["contributed"] + stat["overlapped"] == stat["snp_count"]

    def test_contributed_and_overlapped_split_correctly(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs2", "C", "C")]),
            source("b.txt", "self", [snp("rs2", "C", "C"), snp("rs3", "G", "G")]),
        ])
        stats = {s["label"]: s for s in merged["sources"]}
        assert (stats["a.txt"]["contributed"], stats["a.txt"]["overlapped"]) == (2, 0)
        assert (stats["b.txt"]["contributed"], stats["b.txt"]["overlapped"]) == (1, 1)

    def test_snp_count_reflects_deduplication(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A"), snp("rs1", "G", "G")]),
        ])
        assert merged["sources"][0]["snp_count"] == 1

    def test_conflicting_is_counted_for_every_involved_label(self):
        merged = merge_sources([
            source("a.txt", "self", [snp("rs1", "A", "A")]),
            source("b.txt", "self", [snp("rs1", "G", "G")]),
        ])
        for stat in merged["sources"]:
            assert stat["conflicting"] == 1

    def test_comparison_stats_count_real_calls_and_overlap(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [
                snp("rs1", "A", "G"), snp("rs2", "C", "C"), snp("rs3", "-", "-"),
            ]),
        ])
        stats = {s["label"]: s for s in merged["sources"]}
        assert stats["mum.txt"]["contributed"] == 2
        assert stats["mum.txt"]["overlapped"] == 1
        assert stats["mum.txt"]["snp_count"] == 3


class TestTransmissionProbability:
    def test_aa_by_aa_gives_certainty_for_aa(self):
        assert transmission_probability("AA", "AA", "AA") == 1.0

    def test_aa_by_aa_gives_zero_for_ag(self):
        assert transmission_probability("AG", "AA", "AA") == 0.0

    def test_ag_by_ag_gives_a_quarter_for_aa(self):
        assert transmission_probability("AA", "AG", "AG") == 0.25

    def test_ag_by_ag_gives_a_half_for_ag(self):
        assert transmission_probability("AG", "AG", "AG") == 0.5

    def test_ag_by_ag_gives_a_quarter_for_gg(self):
        assert transmission_probability("GG", "AG", "AG") == 0.25

    def test_ag_by_ag_probabilities_sum_to_one(self):
        total = sum(
            transmission_probability(child, "AG", "AG")
            for child in ("AA", "AG", "GG")
        )
        assert total == 1.0

    def test_aa_by_gg_gives_certainty_for_ag(self):
        assert transmission_probability("AG", "AA", "GG") == 1.0

    def test_aa_by_gg_gives_zero_for_either_homozygote(self):
        assert transmission_probability("AA", "AA", "GG") == 0.0
        assert transmission_probability("GG", "AA", "GG") == 0.0

    def test_ag_by_gg_splits_evenly(self):
        assert transmission_probability("AG", "AG", "GG") == 0.5
        assert transmission_probability("GG", "AG", "GG") == 0.5
        assert transmission_probability("AA", "AG", "GG") == 0.0

    def test_a_no_call_child_returns_none(self):
        assert transmission_probability("--", "AA", "GG") is None
        assert transmission_probability("NN", "AA", "GG") is None

    def test_a_no_call_first_parent_returns_none(self):
        assert transmission_probability("AG", "--", "GG") is None

    def test_a_no_call_second_parent_returns_none(self):
        assert transmission_probability("AG", "AA", "") is None

    def test_none_inputs_return_none(self):
        assert transmission_probability(None, "AA", "GG") is None

    def test_allele_order_does_not_matter(self):
        assert transmission_probability("GA", "AA", "GG") == 1.0
        assert transmission_probability("AG", "GG", "AA") == 1.0

    def test_two_item_sequences_are_accepted(self):
        assert transmission_probability(("A", "G"), ("A", "A"), ("G", "G")) == 1.0

    def test_a_mendelian_impossible_child_gives_zero(self):
        assert transmission_probability("GG", "AA", "AG") == 0.0
        assert transmission_probability("CC", "AA", "AA") == 0.0

    def test_every_result_is_a_documented_value(self):
        allowed = {0.0, 0.25, 0.5, 0.75, 1.0}
        for child in ("AA", "AG", "GG"):
            for parent1 in ("AA", "AG", "GG"):
                for parent2 in ("AA", "AG", "GG"):
                    assert transmission_probability(child, parent1, parent2) in allowed

    def test_one_heterozygous_parent_splits_evenly(self):
        assert transmission_probability("AG", "AG", "AA") == 0.5
        assert transmission_probability("AA", "AG", "AA") == 0.5


class TestMendelianViolation:
    def test_an_impossible_trio_is_flagged(self):
        assert mendelian_violation("GG", "AA", "AG") is True

    def test_a_possible_trio_passes(self):
        assert mendelian_violation("AG", "AA", "GG") is False

    def test_a_certain_trio_passes(self):
        assert mendelian_violation("AA", "AA", "AA") is False

    def test_a_quarter_probability_trio_passes(self):
        assert mendelian_violation("AA", "AG", "AG") is False

    def test_a_no_call_is_not_a_violation(self):
        assert mendelian_violation("--", "AA", "GG") is False
        assert mendelian_violation("AG", "--", "GG") is False
        assert mendelian_violation("AG", "AA", "--") is False

    def test_a_wholly_foreign_allele_is_a_violation(self):
        assert mendelian_violation("CC", "AA", "GG") is True

    def test_allele_order_does_not_create_a_violation(self):
        assert mendelian_violation("GA", "AA", "GG") is False


class TestTrioAnnotate:
    def test_trio_is_available_with_both_parents(self):
        assert trio_annotate(_trio())["trio_available"] is True

    def test_summary_carries_the_documented_keys(self):
        summary = trio_annotate(_trio())
        assert set(summary) == {
            "trio_available", "compared", "violations",
            "violation_rsids", "violation_rate", "note",
        }

    def test_only_positions_with_both_parents_are_compared(self):
        assert trio_annotate(_trio())["compared"] == 3

    def test_a_compatible_position_is_annotated(self):
        merged = _trio()
        trio_annotate(merged)
        entry = merged["genotypes"]["rs1"]
        assert entry["parent1"] == "AA"
        assert entry["parent2"] == "GG"
        assert entry["probability"] == 1.0
        assert entry["mendelian_ok"] is True

    def test_parent1_is_the_mother_and_parent2_the_father(self):
        merged = _trio()
        trio_annotate(merged)
        entry = merged["genotypes"]["rs3"]
        assert entry["parent1"] == merged["comparison"]["mother"]["rs3"]["genotype"]
        assert entry["parent2"] == merged["comparison"]["father"]["rs3"]["genotype"]

    def test_an_impossible_position_is_annotated_as_not_ok(self):
        merged = _trio()
        trio_annotate(merged)
        entry = merged["genotypes"]["rs2"]
        assert entry["probability"] == 0.0
        assert entry["mendelian_ok"] is False

    def test_violations_are_counted(self):
        summary = trio_annotate(_trio())
        assert summary["violations"] == 1
        assert summary["violation_rsids"] == ["rs2"]

    def test_violation_rsids_are_sorted(self):
        summary = trio_annotate(_trio())
        assert summary["violation_rsids"] == sorted(summary["violation_rsids"])

    def test_violation_rate_is_rounded(self):
        assert trio_annotate(_trio())["violation_rate"] == 0.3333

    def test_positions_without_both_parents_are_left_alone(self):
        merged = _trio()
        trio_annotate(merged)
        entry = merged["genotypes"]["rs4"]
        assert "probability" not in entry
        assert "mendelian_ok" not in entry

    def test_annotation_happens_in_place(self):
        merged = _trio()
        entry = merged["genotypes"]["rs1"]
        trio_annotate(merged)
        assert entry["mendelian_ok"] is True

    def test_no_parents_means_no_trio(self):
        merged = merge_sources([source("me.txt", "self", [snp("rs1", "A", "G")])])
        summary = trio_annotate(merged)
        assert summary["trio_available"] is False
        assert summary["compared"] == 0
        assert summary["violations"] == 0
        assert summary["violation_rate"] == 0.0

    def test_one_parent_alone_means_no_trio(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "G")]),
            source("mum.txt", "mother", [snp("rs1", "A", "A")]),
        ])
        summary = trio_annotate(merged)
        assert summary["trio_available"] is False
        assert "probability" not in merged["genotypes"]["rs1"]

    def test_a_mate_is_not_treated_as_a_parent(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "G")]),
            source("wife.txt", "wife", [snp("rs1", "A", "A")]),
            source("mum.txt", "mother", [snp("rs1", "A", "A")]),
        ])
        assert trio_annotate(merged)["trio_available"] is False

    def test_an_empty_merge_result_is_handled(self):
        summary = trio_annotate({})
        assert summary["trio_available"] is False
        assert summary["compared"] == 0
        assert summary["note"] == TRIO_NOTE

    def test_none_is_handled(self):
        assert trio_annotate(None)["trio_available"] is False

    def test_the_note_is_always_included(self):
        assert trio_annotate(_trio())["note"] == TRIO_NOTE

    def test_the_note_blames_strand_and_chip_before_paternity(self):
        note = trio_annotate(_trio())["note"]
        assert "non-paternity" in note
        assert "strand" in note
        assert "chip" in note

    def test_a_clean_trio_reports_no_violations(self):
        merged = merge_sources([
            source("child.txt", "self", [snp("rs1", "A", "G")]),
            source("mum.txt", "mother", [snp("rs1", "A", "A")]),
            source("dad.txt", "father", [snp("rs1", "G", "G")]),
        ])
        summary = trio_annotate(merged)
        assert summary["compared"] == 1
        assert summary["violations"] == 0
        assert summary["violation_rate"] == 0.0
        assert summary["violation_rsids"] == []


class TestComparisonRows:
    def _merged(self):
        return merge_sources([
            source("me.txt", "self", [snp("rs1", "A", "G"), snp("rs2", "A", "A")]),
            source("mum.txt", "mother", [snp("rs1", "G", "A"), snp("rs2", "G", "G")]),
            source("dad.txt", "father", [snp("rs1", "A", "A")]),
        ])

    def test_a_matching_relative_is_marked_shared(self):
        rows = {row["role"]: row for row in comparison_rows(self._merged(), "rs1")}
        assert rows["mother"]["shared"] is True

    def test_a_differing_relative_is_not_marked_shared(self):
        rows = {row["role"]: row for row in comparison_rows(self._merged(), "rs1")}
        assert rows["father"]["shared"] is False

    def test_one_row_per_relative_with_a_call(self):
        assert len(comparison_rows(self._merged(), "rs1")) == 2
        assert len(comparison_rows(self._merged(), "rs2")) == 1

    def test_rows_carry_the_documented_keys(self):
        for row in comparison_rows(self._merged(), "rs1"):
            assert set(row) == {"label", "role", "genotype", "shared"}

    def test_rows_are_ordered_by_comparison_group(self):
        roles = [row["role"] for row in comparison_rows(self._merged(), "rs1")]
        assert roles == ["father", "mother"]

    def test_the_label_is_carried_through(self):
        rows = {row["role"]: row for row in comparison_rows(self._merged(), "rs1")}
        assert rows["mother"]["label"] == "mum.txt"

    def test_allele_order_does_not_break_sharing(self):
        rows = {row["role"]: row for row in comparison_rows(self._merged(), "rs1")}
        assert rows["mother"]["genotype"] == "AG"

    def test_rsid_lookup_is_case_insensitive(self):
        assert len(comparison_rows(self._merged(), "RS1")) == 2

    def test_an_unknown_rsid_gives_no_rows(self):
        assert comparison_rows(self._merged(), "rs999") == []

    def test_an_empty_rsid_gives_no_rows(self):
        assert comparison_rows(self._merged(), "") == []

    def test_an_empty_merge_result_gives_no_rows(self):
        assert comparison_rows({}, "rs1") == []
        assert comparison_rows(None, "rs1") == []

    def test_a_no_call_primary_is_never_shared(self):
        merged = merge_sources([
            source("me.txt", "self", [snp("rs1", "-", "-")]),
            source("mum.txt", "mother", [snp("rs1", "A", "G")]),
        ])
        rows = comparison_rows(merged, "rs1")
        assert len(rows) == 1
        assert rows[0]["shared"] is False
