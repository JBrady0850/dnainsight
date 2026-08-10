"""Tests for tools/audit_y_dbsnp.py, the Y backbone audit against NCBI dbSNP.

WHY THIS FILE EXISTS
--------------------
The audit produces claims about specific markers in a published tree, so a
wrong answer here is not a wrong number, it is a wrong correction applied to
genotyping logic.

One wrong answer has already been produced and discarded. A first draft split
the dbSNP `spdi` field on ":" and expected four parts. That field is a
COMMA-SEPARATED LIST when a site is multi-allelic, so the draft blanked the
allele pair for M45, M343, M269 and P312 and reported all four as conflicts
against the recorded table. Every one of them was fine: the recorded pair was a
subset of a larger observed set.

These tests pin the four behaviours that draft got wrong or nearly got wrong:

  1. Multi-allelic SPDI lists parse completely.
  2. A recorded pair inside a larger observed set is consistent, never CONFLICT.
  3. A non-substitution is CLASS, which is wrong in KIND, and never a conflict
     about allele values.
  4. The reference state is reported and never converted into an ancestral or
     derived assignment, because on the Y the reference carries the derived
     allele about half the time.

Every fixture below is a real dbSNP record retrieved on 2026-08-10, not an
invented shape.
"""

import copy
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import audit_y_dbsnp as A  # noqa: E402


# --- real dbSNP records, retrieved 2026-08-10 -------------------------------

RS2032595 = {  # M168, plain biallelic substitution
    "snp_class": "snv", "chr": "Y",
    "spdi": "NC_000024.10:12702061:T:C",
    "chrpos": "Y:12702062", "chrpos_prev_assm": "Y:14813991",
}
RS2032631 = {  # M45, tri-allelic
    "snp_class": "snv", "chr": "Y",
    "spdi": "NC_000024.10:19705900:A:G,NC_000024.10:19705900:A:T",
    "chrpos": "Y:19705901", "chrpos_prev_assm": "Y:21867787",
}
RS3908 = {  # M17, single base deletion in a G homopolymer
    "snp_class": "delins", "chr": "Y",
    "spdi": "NC_000024.10:19571278:GGGG:GGG",
    "chrpos": "Y:19571279", "chrpos_prev_assm": "Y:21733165",
}
RS3911 = {  # M20, genuinely disagrees with the recorded pair
    "snp_class": "snv", "chr": "Y",
    "spdi": "NC_000024.10:19571567:A:G",
    "chrpos": "Y:19571568", "chrpos_prev_assm": "Y:21733454",
}


def entry(marker, rsid, ancestral, derived):
    return {"marker": marker, "rsid": rsid,
            "ancestral": ancestral, "derived": derived}


class TestParseSpdi:
    def test_a_single_spdi_yields_its_reference_and_alternate(self):
        assert A.parse_spdi(RS2032595["spdi"]) == ("T", {"C"})

    def test_a_multi_allelic_list_yields_every_alternate(self):
        # The exact defect that produced the discarded draft.
        assert A.parse_spdi(RS2032631["spdi"]) == ("A", {"G", "T"})

    def test_an_indel_spdi_yields_whole_sequences(self):
        assert A.parse_spdi(RS3908["spdi"]) == ("GGGG", {"GGG"})

    def test_malformed_entries_are_skipped_not_guessed(self):
        assert A.parse_spdi("nonsense,NC_000024.10:1:A:G") == ("A", {"G"})

    def test_an_absent_field_yields_nothing(self):
        assert A.parse_spdi(None) == ("", set())
        assert A.parse_spdi("") == ("", set())

    def test_surrounding_whitespace_is_tolerated(self):
        assert A.parse_spdi(" NC_000024.10:1:A:G , NC_000024.10:1:A:T ") == ("A", {"G", "T"})


class TestReferenceState:
    def test_a_reference_matching_the_ancestral_state(self):
        assert A.reference_state("G", "G", "T", "snv") == "ancestral"

    def test_a_reference_matching_the_derived_state(self):
        # The case that would invert the tree if it were mapped onto ancestral.
        assert A.reference_state("T", "C", "T", "snv") == "derived"

    def test_an_opposite_strand_reference_is_named_as_such(self):
        assert A.reference_state("A", "T", "C", "snv") == "ancestral (opposite strand)"

    def test_a_non_substitution_has_no_determinable_state(self):
        assert A.reference_state("GGGG", "G", "A", "delins") == "undetermined"

    def test_an_absent_reference_has_no_determinable_state(self):
        assert A.reference_state("", "G", "A", "snv") == "undetermined"


class TestClassify:
    def test_an_exact_pair_is_consistent(self):
        row = A.classify(entry("M168", "rs2032595", "C", "T"), RS2032595)
        assert row["verdict"] == A.CONSISTENT
        assert row["ref_carries"] == "derived"

    def test_a_recorded_pair_inside_a_larger_set_is_consistent(self):
        row = A.classify(entry("M45", "rs2032631", "G", "A"), RS2032631)
        assert row["verdict"] == A.MULTI_ALLELIC

    def test_a_multi_allelic_record_is_never_reported_as_a_conflict(self):
        # The regression this file exists for. A false conflict here is a
        # published accusation about a marker in a real tree.
        for pair in (("G", "A"), ("A", "G"), ("A", "T"), ("T", "A")):
            row = A.classify(entry("M45", "rs2032631", *pair), RS2032631)
            assert row["verdict"] != A.CONFLICT, pair

    def test_a_non_substitution_is_a_class_error_not_an_allele_conflict(self):
        row = A.classify(entry("M17", "rs3908", "G", "A"), RS3908)
        assert row["verdict"] == A.CLASS_ERROR
        assert row["ref_carries"] == "undetermined"

    def test_a_genuine_allele_disagreement_is_a_conflict(self):
        row = A.classify(entry("M20", "rs3911", "A", "C"), RS3911)
        assert row["verdict"] == A.CONFLICT

    def test_a_complemented_pair_is_reported_as_strand_not_conflict(self):
        row = A.classify(entry("X", "rs1", "A", "G"), RS3911 | {"spdi": "NC_000024.10:1:T:C"})
        assert row["verdict"] == A.STRAND

    def test_an_empty_record_is_reported_as_absent_from_dbsnp(self):
        row = A.classify(entry("M31", "rs999999999", "G", "A"), {})
        assert row["verdict"] == A.NOT_FOUND

    def test_both_assembly_positions_are_carried_through(self):
        row = A.classify(entry("M168", "rs2032595", "C", "T"), RS2032595)
        assert row["grch38"] == "Y:12702062"
        assert row["grch37"] == "Y:14813991"

    def test_the_recorded_values_are_reported_alongside_the_observed_ones(self):
        row = A.classify(entry("M20", "rs3911", "A", "C"), RS3911)
        assert row["table_ancestral"] == "A"
        assert row["table_derived"] == "C"
        assert row["dbsnp_ref"] == "A"
        assert row["dbsnp_alts"] == ["G"]


class TestAuditOverTheBackbone:
    BACKBONE = {
        "root":  {"marker": None, "rsid": None},
        "CT":    entry("M168", "rs2032595", "C", "T"),
        "R1a1a": entry("M17", "rs3908", "G", "A"),
        "BT":    entry("M91", None, "A", "T"),
    }
    RECORDS = {"2032595": RS2032595, "3908": RS3908}

    def test_only_rows_with_an_rsid_are_audited(self):
        rows = A.audit(self.BACKBONE, self.RECORDS)
        assert {r["marker"] for r in rows} == {"M168", "M17"}

    def test_rows_without_an_rsid_are_reported_as_unresolved(self):
        assert A.unresolved(self.BACKBONE) == ["M91"]

    def test_the_marker_free_root_is_not_reported_as_unresolved(self):
        assert "None" not in A.unresolved(self.BACKBONE)

    def test_the_audit_never_modifies_the_backbone(self):
        before = copy.deepcopy(self.BACKBONE)
        A.audit(self.BACKBONE, self.RECORDS)
        A.unresolved(self.BACKBONE)
        assert self.BACKBONE == before

    def test_each_row_is_labelled_with_its_node(self):
        rows = A.audit(self.BACKBONE, self.RECORDS)
        assert {r["node"] for r in rows} == {"CT", "R1a1a"}
