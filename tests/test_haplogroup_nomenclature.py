"""Nomenclature equivalence and the reference-allele inversion guard.

WHY THIS FILE EXISTS
--------------------
Two distinct failures motivated it, both found on 2026-08-09 while checking a
public claim that consumer vendors "misclassify J1 as generic J-M267".

FAILURE ONE, THE READER'S: J-M267 *is* J1.
------------------------------------------
M267 is the SNP that defines haplogroup J1, so "J-M267" is not a vague or
degraded rendering of J1, it is the same node written in SNP-based nomenclature
instead of the older letter-number form. FamilyTreeDNA moved to SNP names
precisely because letter-number labels get renumbered whenever the tree is
revised, which is the same reason this module already stamps ``tree_name`` and
``tree_version`` on every payload.

Someone reading a DNAInsight result next to an FTDNA result had no way to know
the two names were the same node, so the obvious conclusion was that one of the
two tools was wrong. Reporting both names removes the ambiguity entirely.

Low resolution and a wrong call are different failures, and this project exists
to keep them apart. J1 and J2 are different branches under different SNPs
(M267 and M172). Calling J1 at low depth is not the same thing as calling it J2,
and the tests below hold that line: the equivalence sets for J1 and J2 must
never intersect.

FAILURE TWO, OURS: dbSNP ref/alt is NOT ancestral/derived.
----------------------------------------------------------
dbSNP reports rs2032595 (M168, the CT node) as chrY:12702062 T>C on the GRCh38
forward strand. ``Y_BACKBONE`` records M168 as ancestral C, derived T. That
reads as a flat contradiction and is not one: dbSNP reports REFERENCE over
ALTERNATE, and the GRCh38 reference Y chromosome comes from a non-African
lineage that carries the DERIVED allele at M168. Reference T is therefore the
derived state and alternate C the ancestral one.

So a builder that "verifies" ``ancestral``/``derived`` against dbSNP's
``ref``/``alt`` would invert roughly half this tree, and every existing test
would still pass, because the data would be internally consistent and
externally backwards. That is the same shape as the CPIC positive-strand
conflicts already recorded in docs/KNOWN_GAPS.md.

The structural fix is that an entry may not claim ``verified`` until it records
which assembly it was checked against and which state that assembly's reference
carries. The tests below make that unskippable rather than remembered.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from backend import haplogroups as h  # noqa: E402


MARKER_NODES = {name: entry for name, entry in h.Y_BACKBONE.items()
                if entry.get("marker")}


# ---------------------------------------------------------------------------
# SNP-based name construction
# ---------------------------------------------------------------------------

class TestSnpName:
    def test_j1_is_written_j_m267(self):
        # The exact claim that motivated this file.
        assert h.snp_name(h.Y_BACKBONE["J1"]) == "J-M267"

    def test_j2_is_written_j_m172(self):
        assert h.snp_name(h.Y_BACKBONE["J2"]) == "J-M172"

    def test_a_long_label_keeps_only_the_leading_letters(self):
        # R1b1a1b under M269 is written R-M269, not R1b1a1b-M269.
        assert h.snp_name(h.Y_BACKBONE["R-M269"]) == "R-M269"

    def test_a_multi_letter_internal_node_keeps_all_its_letters(self):
        # CT, DE, IJK and friends are genuinely multi-letter labels.
        assert h.snp_name(h.Y_BACKBONE["CT"]) == "CT-M168"

    def test_the_root_has_no_snp_name(self):
        # Every male is at the root by definition, so no marker defines it.
        assert h.snp_name(h.Y_BACKBONE["root"]) is None

    def test_every_marker_bearing_node_produces_a_name(self):
        missing = [name for name, entry in MARKER_NODES.items()
                   if not h.snp_name(entry)]
        assert not missing, f"no SNP name for: {missing}"

    def test_no_snp_name_contains_a_digit_before_the_hyphen(self):
        bad = []
        for name, entry in MARKER_NODES.items():
            prefix = h.snp_name(entry).split("-", 1)[0]
            if any(character.isdigit() for character in prefix):
                bad.append(f"{name} -> {h.snp_name(entry)}")
        assert not bad, "letter-number label leaked into a SNP name:\n  " + "\n  ".join(bad)


# ---------------------------------------------------------------------------
# Equivalence sets
# ---------------------------------------------------------------------------

class TestEquivalentNames:
    def test_j1_reports_both_of_its_names(self):
        assert set(h.equivalent_names(h.Y_BACKBONE["J1"])) == {"J1", "J-M267"}

    def test_the_label_is_always_present(self):
        for name, entry in h.Y_BACKBONE.items():
            assert entry["label"] in h.equivalent_names(entry), name

    def test_names_are_unique_within_an_entry(self):
        for name, entry in h.Y_BACKBONE.items():
            names = h.equivalent_names(entry)
            assert len(names) == len(set(names)), f"{name} repeats a name"

    def test_j1_and_j2_never_share_a_name(self):
        """The distinction the motivating claim collapsed.

        Under-resolution and misassignment are different failures. If these two
        sets ever intersect, a J1 result could render as a J2 name and the
        difference would stop being visible at all.
        """
        j1 = set(h.equivalent_names(h.Y_BACKBONE["J1"]))
        j2 = set(h.equivalent_names(h.Y_BACKBONE["J2"]))
        assert not (j1 & j2), f"J1 and J2 share names: {j1 & j2}"

    def test_no_name_is_claimed_by_two_different_nodes(self):
        seen: dict[str, str] = {}
        clashes = []
        for node, entry in h.Y_BACKBONE.items():
            for name in h.equivalent_names(entry):
                if name in seen and seen[name] != node:
                    clashes.append(f"{name!r} claimed by {seen[name]} and {node}")
                seen[name] = node
        assert not clashes, "ambiguous haplogroup name:\n  " + "\n  ".join(clashes)

    def test_the_node_key_is_included_when_it_differs_from_the_label(self):
        # "R-M269" is the dict key and "R1b1a1b" is the label. A caller holding
        # either string has to be able to find the other.
        names = set(h.equivalent_names(h.Y_BACKBONE["R-M269"]))
        assert {"R-M269", "R1b1a1b"} <= names


# ---------------------------------------------------------------------------
# The inversion guard
# ---------------------------------------------------------------------------

class TestReferenceOrientationIsRecorded:
    def test_every_marker_entry_carries_the_new_fields(self):
        missing = [name for name, entry in MARKER_NODES.items()
                   if "assembly" not in entry
                   or "ref_carries" not in entry
                   or "dbsnp_checked" not in entry]
        assert not missing, f"entry predates the inversion guard: {missing}"

    def test_ref_carries_is_ancestral_derived_or_unknown(self):
        allowed = {"ancestral", "derived", None}
        bad = [f"{name}={entry['ref_carries']!r}"
               for name, entry in MARKER_NODES.items()
               if entry["ref_carries"] not in allowed]
        assert not bad, f"unusable ref_carries value: {bad}"

    def test_nothing_may_be_verified_without_an_orientation(self):
        """The whole point. Verified plus unknown orientation is the bug.

        An entry that claims verification while nobody recorded which state the
        reference carries is exactly the state that lets a dbSNP-driven builder
        invert it silently.
        """
        offenders = [name for name, entry in MARKER_NODES.items()
                     if entry.get("verified") and entry.get("ref_carries") is None]
        assert not offenders, (
            "verified without a recorded reference orientation: " + str(offenders)
        )

    def test_nothing_may_be_verified_without_an_assembly(self):
        offenders = [name for name, entry in MARKER_NODES.items()
                     if entry.get("verified") and not entry.get("assembly")]
        assert not offenders, f"verified without an assembly: {offenders}"

    def test_nothing_may_be_verified_without_a_dbsnp_check(self):
        offenders = [name for name, entry in MARKER_NODES.items()
                     if entry.get("verified") and not entry.get("dbsnp_checked")]
        assert not offenders, f"verified without a dbSNP check: {offenders}"

    def test_dbsnp_checked_requires_an_rsid_to_check_against(self):
        offenders = [name for name, entry in MARKER_NODES.items()
                     if entry.get("dbsnp_checked") and not entry.get("rsid")]
        assert not offenders, f"claims a dbSNP check with no rsID: {offenders}"

    def test_the_ct_node_records_the_case_that_started_this(self):
        """rs2032595 is the worked example, so it is pinned as one.

        If someone later 'fixes' M168 to match dbSNP ref/alt, this fails and the
        docstring above explains why the fix is wrong.
        """
        ct = h.Y_BACKBONE["CT"]
        assert ct["rsid"] == "rs2032595"
        assert ct["ancestral"] == "C" and ct["derived"] == "T"


# ---------------------------------------------------------------------------
# The call payload carries the names through
# ---------------------------------------------------------------------------

class TestPayloadExposesBothNames:
    def _genotypes_reaching(self, node: str) -> dict:
        """Build a genotype map that walks the tree down to ``node``."""
        genotypes = {}
        for step in h.path_to(node):
            entry = h.Y_BACKBONE[step]
            if not entry.get("marker"):
                continue
            # marker_keys() lowercases, matching the project's genotype-map
            # convention, so a map keyed "M267" reads as not testable.
            key = (entry.get("rsid") or entry["marker"]).lower()
            genotypes[key] = entry["derived"]
        return genotypes

    def test_a_j1_call_reports_the_snp_name_too(self):
        call = h.call_y_backbone(self._genotypes_reaching("J1"), sex_hint="male")
        assert call["available"] is True
        assert call["haplogroup"] == "J1"
        assert "J-M267" in call["equivalent_names"]
        assert call["also_written_as"] == "J-M267"

    def test_an_unresolved_call_reports_no_names_rather_than_guessing(self):
        call = h.call_y_backbone({}, sex_hint="male")
        assert call["available"] is False
        assert call["equivalent_names"] == []
        assert call["also_written_as"] is None

    def test_the_mitochondrial_payload_has_the_same_keys(self):
        # Field-set parity between the two systems is a standing invariant.
        y = h.call_y_backbone({}, sex_hint="male")
        mt = h.call_mt_backbone({})
        assert set(y) - {"sex_hint"} == set(mt) - {"sex_hint"}

    @pytest.mark.parametrize("node", ["J1", "J2", "R-M269", "I1", "E-M2"])
    def test_named_calls_never_collide_across_branches(self, node):
        call = h.call_y_backbone(self._genotypes_reaching(node), sex_hint="male")
        assert call["haplogroup"] == node
        for other in ("J1", "J2", "R-M269", "I1", "E-M2"):
            if other == node:
                continue
            other_names = set(h.equivalent_names(h.Y_BACKBONE[other]))
            assert not (set(call["equivalent_names"]) & other_names)


class TestIndelMarkersAreWrongInKindNotInValue:
    """M17 and M91 are indels that were recorded as base substitutions.

    An indel recorded as a substitution is not a wrong value that a better
    allele would fix. A consumer array reports two base calls at a position and
    a deletion is not a base, so no genotype can ever satisfy the rule. The
    nodes those markers define are unreachable from array data, and that is a
    structural ceiling rather than missing coverage.

    M17 was confirmed against dbSNP rs3908 (delins, SPDI GGGG:GGG). M91 is
    stated as 9T to 8T in the Karafet et al. 2008 text and has no rsID, so it
    carries sequences from the publication and remains unverified.
    """

    def _indels(self):
        return {n: e for n, e in h.Y_BACKBONE.items()
                if e.get("variant_type") not in (None, "snv")}

    def test_the_two_known_indels_are_recorded_as_such(self):
        assert {e["marker"] for e in self._indels().values()} == {"M17", "M91"}

    def test_an_indel_clears_the_single_base_fields(self):
        # Widening these to hold "GGGG" would have been compared against an
        # array base call and silently never matched.
        for name, entry in self._indels().items():
            assert entry["ancestral"] is None, name
            assert entry["derived"] is None, name

    def test_an_indel_records_whole_sequences_instead(self):
        for name, entry in self._indels().items():
            assert entry["ancestral_seq"], name
            assert entry["derived_seq"], name
            assert len(entry["ancestral_seq"]) != len(entry["derived_seq"]), name

    def test_m17_carries_the_dbsnp_sequences(self):
        entry = h.Y_BACKBONE["R1a1a"]
        assert (entry["ancestral_seq"], entry["derived_seq"]) == ("GGGG", "GGG")

    def test_m91_carries_the_published_homopolymer_lengths(self):
        entry = h.Y_BACKBONE["BT"]
        assert len(entry["ancestral_seq"]) == 9
        assert len(entry["derived_seq"]) == 8

    def test_an_indel_may_never_be_marked_verified(self):
        for name, entry in self._indels().items():
            assert entry["verified"] is False, name

    def test_an_indel_has_no_reference_state(self):
        # An indel has no single reference base, so ref_carries cannot be
        # filled in and must not be guessed at.
        for name, entry in self._indels().items():
            assert entry["ref_carries"] is None, name

    def test_a_substitution_carries_no_sequences(self):
        for name, entry in h.Y_BACKBONE.items():
            if entry.get("variant_type") == "snv":
                assert entry["ancestral_seq"] is None, name
                assert entry["derived_seq"] is None, name

    def test_an_indel_can_never_be_called_from_an_array_genotype(self):
        # The load-bearing behaviour. Even handed a real base call at the site,
        # the marker must not resolve to derived or ancestral.
        state, _ = h.marker_state(h.Y_BACKBONE["R1a1a"], {"rs3908": ("G", "G")})
        assert state == h.UNUSABLE

    def test_untypeable_markers_reports_them_apart_from_unverified_ones(self):
        report = h.untypeable_markers()
        assert report["count"] == 2
        assert {row["marker"] for row in report["markers"]} == {"M17", "M91"}
        assert "structural ceiling" in report["reason"]


class TestTheDbsnpAuditIsRecordedWithoutOverclaiming:
    """v3.3.0 folded the dbSNP audit into the table. It settles class and
    position; it cannot settle ancestral against derived, and does not pretend to.
    """

    def _audited(self):
        return {n: e for n, e in h.Y_BACKBONE.items() if e.get("dbsnp_checked")}

    def test_every_audited_row_names_the_assembly_it_was_checked_against(self):
        for name, entry in self._audited().items():
            assert entry["assembly"] == "GRCh38", name

    def test_the_audit_marked_nothing_verified(self):
        # dbSNP answers reference over alternate. The table records ancestral
        # over derived. Recording one as the other is the inversion the
        # ref_carries guard exists to prevent.
        verified = [n for n, e in h.Y_BACKBONE.items()
                    if e.get("verified") and n != "root"]
        assert verified == []

    def test_the_reference_carries_the_derived_allele_more_often_than_not(self):
        # Measured, not asserted from memory: 10 of the 17 determinable nodes.
        carried = [e["ref_carries"] for e in h.Y_BACKBONE.values()
                   if e.get("ref_carries")]
        assert carried.count("derived") == 10
        assert carried.count("ancestral") == 7

    def test_the_unresolved_m20_conflict_is_recorded_rather_than_silently_fixed(self):
        note = h.Y_BACKBONE["L"]["note"]
        assert "CONFLICT" in note and "rs3911" in note

    def test_multi_allelic_sites_are_flagged_not_treated_as_conflicts(self):
        flagged = {e["marker"] for e in h.Y_BACKBONE.values()
                   if e.get("multi_allelic")}
        assert flagged == {"M45", "M343", "M269", "P312"}
