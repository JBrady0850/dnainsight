"""Pins what the Karafet 2008 supplement settled, and what it did not.

WHY THIS FILE EXISTS SEPARATELY
-------------------------------
``tests/test_haplogroups.py`` tests the calling machinery and
``tests/test_haplogroup_nomenclature.py`` tests naming and the dbSNP audit.
This file tests the third source layer, and specifically the three failure
modes a passing suite cannot otherwise catch:

1. An indel recorded as a base substitution. v3.3.0 found two. This source
   found two more, in markers nothing had flagged.
2. A published value quietly promoted to ``verified``. A publication gives
   ancestral over derived. It does not give the reference orientation, and
   without that the tree can be inverted while every test still passes.
3. A conflict resolved by picking a side. Two markers disagree with the
   supplement in a way no strand or direction operation reconciles. Those stay
   held, and this file fails the build if anyone writes a value into them.

Source: Karafet TM, Mendez FL, Meilerman MB, Underhill PA, Zegura SL, Hammer MF,
"New binary polymorphisms reshape and increase resolution of the human Y
chromosomal haplogroup tree", Genome Research 18:830-838, 2008,
Supplementary Table 1.
"""

from __future__ import annotations

import pytest

from backend import haplogroups as hg


# marker -> rsID, read off Supplementary Table 1 and spot-checked against the
# pages rendered at 170 dpi before any row was trusted.
KARAFET_RSIDS = {
    "M2": "rs3893",       "M3": "rs3894",        "M60": "rs2032623",
    "M69": "rs2032673",   "M91": "rs2032651",    "M130": "rs35284970",
    "M145": "rs3848982",  "M172": "rs2032604",   "M174": "rs2032602",
    "M175": "rs2032678",  "M184": "rs20320",     "M214": "rs2032674",
    "M217": "rs2032668",  "M231": "rs9341278",   "M267": "rs9341313",
    "M438": "rs17307294", "P143": "rs4141886",
}


def entry_for(marker):
    hits = [e for e in hg.Y_BACKBONE.values() if e.get("marker") == marker]
    assert len(hits) == 1, f"{marker} appears {len(hits)} times"
    return hits[0]


class TestResolvedRsids:

    @pytest.mark.parametrize("marker,rsid", sorted(KARAFET_RSIDS.items()))
    def test_the_marker_carries_the_rsid_the_supplement_gives_it(self, marker, rsid):
        assert entry_for(marker)["rsid"] == rsid

    def test_the_backbone_reaches_thirty_five_rsids(self):
        """17 before v3.4.0. The audit can only see rows that carry one."""
        assert sum(1 for e in hg.Y_BACKBONE.values() if e.get("rsid")) == 35

    @pytest.mark.parametrize("marker", sorted(KARAFET_RSIDS))
    def test_every_resolved_entry_cites_its_source(self, marker):
        assert "Karafet" in entry_for(marker)["note"]

    def test_a_citation_does_not_promote_an_entry_to_verified(self):
        for name, entry in hg.Y_BACKBONE.items():
            if name == "root":
                continue
            assert entry["verified"] is False, name


class TestMarkersTheSurveyLeftWithoutAnRsid:
    """Five markers were genotyped and assigned no RefSNP ID.

    That absence is a measured fact, not an oversight, and it is the reason
    those rows must never acquire a guessed rsID.
    """

    NO_RSID = ["M35", "P15", "P37.2", "M410", "M122"]

    @pytest.mark.parametrize("marker", NO_RSID)
    def test_the_rsid_stays_none_and_the_note_says_why(self, marker):
        entry = entry_for(marker)
        assert entry["rsid"] is None
        assert "no RefSNP ID" in entry["note"]


class TestTransposedPairs:
    """Three entries ran ancestral and derived the wrong way round."""

    @pytest.mark.parametrize("node,marker,pair", [
        ("DE", "M145", ("G", "A")),
        ("N-M178", "M178", ("T", "C")),
        ("J1", "M267", ("T", "G")),
    ])
    def test_the_pair_now_matches_the_source(self, node, marker, pair):
        entry = hg.Y_BACKBONE[node]
        assert entry["marker"] == marker
        assert (entry["ancestral"], entry["derived"]) == pair

    def test_m267_was_complemented_and_reversed_and_the_note_says_so(self):
        note = hg.Y_BACKBONE["J1"]["note"]
        assert "complemented AND reversed" in note


class TestTheTwoNewIndels:
    """M60 and M175 are the same class error v3.3.0 found in M17 and M91.

    Neither was flagged then. Both were stored as substitutions with allele
    pairs that cannot describe the variant, and nothing failed, because the
    data was internally consistent and externally wrong.
    """

    @pytest.mark.parametrize("node,marker,vtype,rsid", [
        ("B", "M60", "ins", "rs2032623"),
        ("O", "M175", "del", "rs2032678"),
    ])
    def test_the_marker_is_recorded_as_a_length_polymorphism(
            self, node, marker, vtype, rsid):
        entry = hg.Y_BACKBONE[node]
        assert entry["marker"] == marker
        assert entry["variant_type"] == vtype
        assert entry["rsid"] == rsid
        assert entry["ancestral"] is None
        assert entry["derived"] is None

    @pytest.mark.parametrize("node", ["B", "O"])
    def test_it_reads_as_untypeable_rather_than_ancestral(self, node):
        """The tri-state invariant, applied to a marker with no base call.

        Reporting an unreadable marker as ancestral turns "we could not look"
        into "we looked and found nothing", which is the single behaviour this
        module exists to prevent.
        """
        entry = hg.Y_BACKBONE[node]
        genotypes = {k: ("A", "A") for k in hg.marker_keys(entry)}
        state, _ = hg.marker_state(entry, genotypes)
        assert state in hg.NOT_TESTABLE_STATES
        assert state != hg.ANCESTRAL

    def test_both_cite_the_dbsnp_class_that_confirms_them(self):
        assert "snp_class ins" in hg.Y_BACKBONE["B"]["note"]
        assert "delins" in hg.Y_BACKBONE["O"]["note"]


class TestM91GainedAnRsid:

    def test_the_rsid_is_recorded_but_not_treated_as_checked(self):
        """NCBI esummary returns an empty record for rs2032651.

        It is a 2001-era accession and has most likely been merged forward.
        Recording the rsID is honest; recording a dbSNP check that never
        succeeded would not be.
        """
        entry = hg.Y_BACKBONE["BT"]
        assert entry["marker"] == "M91"
        assert entry["rsid"] == "rs2032651"
        assert entry["dbsnp_checked"] is False

    def test_the_clade_a_difference_is_recorded_and_not_acted_on(self):
        note = hg.Y_BACKBONE["BT"]["note"]
        assert "haplogroup A rather than BT" in note
        assert "NOT acted on" in note
        assert hg.Y_BACKBONE["BT"]["parent"] == "root"


class TestHeldConflicts:
    """M31 and M429 disagree with the supplement irreconcilably.

    Neither complementing nor transposing maps the stored pair onto the source
    pair. Exactly one side is wrong in each case and nothing available says
    which, so nothing is written. If a later session writes a value here
    without recording what arbitrated it, this fails.
    """

    @pytest.mark.parametrize("node,marker", [("A", "M31"), ("IJ", "M429")])
    def test_the_marker_is_held_with_no_rsid_and_a_stated_reason(self, node, marker):
        entry = hg.Y_BACKBONE[node]
        assert entry["marker"] == marker
        assert entry["rsid"] is None
        assert entry["note"].startswith("HELD.")
        assert entry["verified"] is False

    def test_m31_records_the_value_it_refused(self):
        assert "G->C" in hg.Y_BACKBONE["A"]["note"]

    def test_m429_records_the_rsid_it_refused_to_adopt(self):
        """The rsID travels with the alleles that conflict, so it is not
        written into the field either, only into the note."""
        assert "rs17306671" in hg.Y_BACKBONE["IJ"]["note"]
        assert hg.Y_BACKBONE["IJ"]["rsid"] is None


class TestMarkersTheSourceCannotSettle:
    """Six markers post-date the 2008 survey.

    Recorded so no future session spends effort re-reading a supplement that
    cannot contain the answer.
    """

    POST_2008 = {
        "GHIJK": "F1329", "HIJK": "F929", "IJK": "L15",
        "LT": "L298", "K2b": "P331", "R1a": "M420",
    }

    @pytest.mark.parametrize("node", sorted(POST_2008))
    def test_the_entry_says_the_source_can_never_resolve_it(self, node):
        entry = hg.Y_BACKBONE[node]
        assert entry["marker"] == self.POST_2008[node]
        assert entry["rsid"] is None
        assert "post-dates" in entry["note"]


class TestTheLayersStaySeparate:
    """Three sources of differing strength, kept apart on purpose.

    The literal table is literature recall, ``_apply_audit`` is what dbSNP
    measured, and ``_apply_karafet`` is what a publication states. Merging them
    would make it impossible to tell afterwards which rows rest on which
    source.
    """

    def test_the_karafet_pass_runs_after_the_dbsnp_audit(self):
        source = hg.__file__
        text = open(source, encoding="utf-8").read()
        assert text.index("_apply_audit(Y_BACKBONE)") < text.index("_apply_karafet(Y_BACKBONE)")

    def test_the_dbsnp_audit_still_owns_the_reference_orientation(self):
        """No Karafet row may set ref_carries. A publication cannot supply it."""
        for marker in list(KARAFET_RSIDS) + ["M31", "M429"]:
            entry = entry_for(marker)
            if entry.get("ref_carries") is not None:
                assert entry["dbsnp_checked"] is True, marker
