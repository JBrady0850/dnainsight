"""
Tests for backend.ledger: the reclassification ledger.

Every test runs against a throwaway SQLite file. The existing suite does this
by monkeypatching ``database.DB_PATH`` directly, because that module resolves
its path at IMPORT time, so setting the environment variable alone would land
every write in the developer's real database. Both are set here: the
environment variable so any subprocess or re-import also lands in the temp
directory, and the attribute so the already-imported module actually obeys.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import database as db
from backend import ledger


# ---------------------------------------------------------------------------
# Isolation
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolated_db(tmp_path, monkeypatch):
    target = tmp_path / "t.db"
    monkeypatch.setenv("DNAINSIGHT_DB_PATH", str(target))
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(db, "DB_PATH", target)
    ledger.init_ledger()
    yield target


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def snp(rsid="rs1801133", **over):
    finding = {
        "rsid": rsid,
        "entity_type": "snp",
        "gene": "MTHFR",
        "clinical_sig": "pathogenic",
        "clinvar_sig_code": 5,
        "review_stars": 2,
        "cpic_level": "",
        "magnitude": 5.0,
        "repute": "Bad",
        "variant_copies": 1,
        "zygosity": "heterozygous",
        "interpretation": "Some prose that changes constantly.",
        "summary": "More prose.",
        "conditions": "Homocystinuria",
    }
    finding.update(over)
    return finding


def vus(rsid="rs80357713", **over):
    return snp(rsid, clinical_sig="uncertain significance",
               clinvar_sig_code=1, repute="", magnitude=1.5, **over)


def genoset(name="gs100", coverage=1.0, **over):
    finding = {
        "rsid": name,
        "entity_type": "genoset",
        "gene": "",
        "coverage": coverage,
        "magnitude": 3.0,
        "repute": "",
        "clinical_sig": "",
        "review_stars": 0,
    }
    finding.update(over)
    return finding


def prs(name="PGS000018", reliable=True, **over):
    finding = {
        "rsid": name,
        "entity_type": "prs",
        "gene": "",
        "reliable": reliable,
        "coverage": 0.97 if reliable else 0.10,
        "magnitude": 2.0,
        "repute": "",
        "clinical_sig": "",
        "review_stars": 0,
    }
    finding.update(over)
    return finding


def kinds_of(diff):
    return {c["kind"] for c in diff["changes"]}


def change_of(diff, kind):
    for c in diff["changes"]:
        if c["kind"] == kind:
            return c
    raise AssertionError(f"no {kind} change in {sorted(kinds_of(diff))}")


# ---------------------------------------------------------------------------
# fingerprint
# ---------------------------------------------------------------------------

class TestFingerprint:

    def test_fingerprint_returns_a_digest(self):
        fp = ledger.fingerprint(snp())
        assert isinstance(fp["digest"], str) and len(fp["digest"]) == 64

    def test_fingerprint_key_combines_entity_type_and_rsid(self):
        assert ledger.fingerprint(snp("rs123"))["key"] == "snp:rs123"
        assert ledger.fingerprint(genoset("gs5"))["key"] == "genoset:gs5"

    def test_fingerprint_is_stable_across_interpretation_rewording(self):
        before = ledger.fingerprint(snp(interpretation="Old wording here."))
        after = ledger.fingerprint(snp(interpretation="Entirely new wording."))
        assert before["digest"] == after["digest"]

    def test_fingerprint_is_stable_across_summary_and_condition_churn(self):
        before = ledger.fingerprint(snp(summary="a", conditions="Condition A"))
        after = ledger.fingerprint(snp(summary="b", conditions="Condition A, B"))
        assert before["digest"] == after["digest"]

    def test_fingerprint_is_stable_across_gene_symbol_renaming(self):
        before = ledger.fingerprint(snp(gene="MTHFR"))
        after = ledger.fingerprint(snp(gene="MTHFR1"))
        assert before["digest"] == after["digest"]

    def test_fingerprint_is_stable_across_position_and_chromosome_edits(self):
        before = ledger.fingerprint(snp(chromosome="1", position=11856378))
        after = ledger.fingerprint(snp(chromosome="chr1", position=11796321))
        assert before["digest"] == after["digest"]

    def test_fingerprint_changes_when_significance_changes(self):
        before = ledger.fingerprint(snp(clinvar_sig_code=1))
        after = ledger.fingerprint(snp(clinvar_sig_code=5))
        assert before["digest"] != after["digest"]

    def test_fingerprint_rounds_magnitude_to_one_decimal_place(self):
        assert ledger.fingerprint(snp(magnitude=4.53))["magnitude"] == 4.5
        assert ledger.fingerprint(snp(magnitude=4.549))["magnitude"] == 4.5

    def test_magnitude_noise_below_one_decimal_place_does_not_move_the_digest(self):
        before = ledger.fingerprint(snp(magnitude=4.51))
        after = ledger.fingerprint(snp(magnitude=4.54))
        assert before["digest"] == after["digest"]

    def test_fingerprint_carries_no_prose_fields(self):
        fp = ledger.fingerprint(snp())
        for banned in ("interpretation", "summary", "conditions", "gene"):
            assert banned not in fp

    def test_fingerprint_computes_the_sig_code_when_absent(self):
        fp = ledger.fingerprint(snp(clinvar_sig_code=None,
                                    clinical_sig="likely pathogenic"))
        assert fp["clinvar_sig_code"] == 4

    def test_fingerprint_never_resolves_a_conflicting_record(self):
        fp = ledger.fingerprint(snp(
            clinvar_sig_code=None,
            clinical_sig="conflicting classifications of pathogenicity"))
        assert fp["clinvar_sig_code"] == 255

    def test_fingerprint_normalises_the_cpic_level(self):
        assert ledger.fingerprint(snp(cpic_level="a"))["cpic_level"] == "A"
        assert ledger.fingerprint(snp(cpic_level="junk"))["cpic_level"] == ""

    def test_fingerprint_records_carrier_copies(self):
        assert ledger.fingerprint(snp(variant_copies=0))["carrier"] == 0
        assert ledger.fingerprint(snp(variant_copies=2))["carrier"] == 2

    def test_fingerprint_carrier_is_none_when_the_risk_allele_is_unknown(self):
        assert ledger.fingerprint(snp(variant_copies=None))["carrier"] is None

    def test_fingerprint_rejects_a_non_repute_string(self):
        assert ledger.fingerprint(snp(repute="Maybe"))["repute"] == ""

    def test_a_no_call_is_not_evaluable(self):
        assert ledger.fingerprint(snp(zygosity="no_call"))["evaluable"] is False

    def test_a_fully_covered_genoset_is_evaluable(self):
        assert ledger.fingerprint(genoset(coverage=1.0))["evaluable"] is True

    def test_a_partially_covered_genoset_is_not_evaluable(self):
        assert ledger.fingerprint(genoset(coverage=0.75))["evaluable"] is False

    def test_a_reliable_prs_is_evaluable(self):
        assert ledger.fingerprint(prs(reliable=True))["evaluable"] is True

    def test_an_unreliable_prs_is_not_evaluable(self):
        assert ledger.fingerprint(prs(reliable=False))["evaluable"] is False

    def test_prs_evaluability_falls_back_to_the_nested_result_block(self):
        finding = {"rsid": "PGS1", "entity_type": "prs",
                   "prs": {"reliable": True, "coverage": 0.99}}
        assert ledger.fingerprint(finding)["evaluable"] is True

    def test_fingerprint_tolerates_a_non_dict(self):
        fp = ledger.fingerprint(None)
        assert fp["rsid"] == "" and isinstance(fp["digest"], str)


# ---------------------------------------------------------------------------
# Change vocabulary
# ---------------------------------------------------------------------------

class TestChangeVocabulary:

    REQUIRED = (
        "sig_upgraded", "sig_downgraded", "vus_resolved_pathogenic",
        "vus_resolved_benign", "stars_gained", "stars_lost",
        "cpic_level_changed", "magnitude_changed", "repute_changed",
        "newly_evaluable", "no_longer_evaluable", "finding_added",
        "finding_removed",
    )

    def test_every_required_kind_is_present(self):
        for kind in self.REQUIRED:
            assert kind in ledger.CHANGE_KINDS

    def test_every_kind_has_a_label_severity_and_direction(self):
        for kind, meta in ledger.CHANGE_KINDS.items():
            assert meta["label"], kind
            assert isinstance(meta["severity"], int), kind
            assert meta["direction"] in ("up", "down", "lateral"), kind

    def test_vus_resolved_pathogenic_is_the_highest_severity_in_the_system(self):
        top = max(ledger.CHANGE_KINDS.values(), key=lambda m: m["severity"])
        assert top is ledger.CHANGE_KINDS["vus_resolved_pathogenic"]

    def test_vus_resolved_pathogenic_severity_is_unique(self):
        peak = ledger.CHANGE_KINDS["vus_resolved_pathogenic"]["severity"]
        assert sum(1 for m in ledger.CHANGE_KINDS.values()
                   if m["severity"] == peak) == 1

    def test_severities_are_unique_so_ordering_is_deterministic(self):
        severities = [m["severity"] for m in ledger.CHANGE_KINDS.values()]
        assert len(severities) == len(set(severities))

    def test_material_kinds_is_a_frozenset(self):
        assert isinstance(ledger.MATERIAL_KINDS, frozenset)

    def test_magnitude_changed_is_not_material(self):
        assert "magnitude_changed" not in ledger.MATERIAL_KINDS

    def test_every_other_kind_is_material(self):
        assert ledger.MATERIAL_KINDS == frozenset(
            k for k in ledger.CHANGE_KINDS if k != "magnitude_changed")


# ---------------------------------------------------------------------------
# Snapshots
# ---------------------------------------------------------------------------

class TestSnapshots:

    def test_init_ledger_is_idempotent(self):
        ledger.init_ledger()
        ledger.init_ledger()
        assert ledger.list_snapshots(1) == []

    def test_init_ledger_does_not_disturb_existing_snapshots(self):
        first = ledger.snapshot(1, [snp()])
        ledger.init_ledger()
        assert ledger.get_snapshot(first) is not None
        assert len(ledger.list_snapshots(1)) == 1

    def test_snapshot_returns_an_id(self):
        assert isinstance(ledger.snapshot(1, [snp()]), int)

    def test_snapshot_records_the_finding_count(self):
        sid = ledger.snapshot(1, [snp("rs1"), snp("rs2")])
        assert ledger.get_snapshot(sid)["finding_count"] == 2

    def test_snapshot_stores_the_database_versions_in_force(self):
        sid = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        assert ledger.get_snapshot(sid)["db_versions"]["clinvar"] == "2026-07"

    def test_snapshot_fills_db_versions_from_provenance_when_omitted(self):
        sid = ledger.snapshot(1, [snp()])
        versions = ledger.get_snapshot(sid)["db_versions"]
        assert "clinvar" in versions and "dnainsight" in versions

    def test_snapshot_stores_the_label(self):
        sid = ledger.snapshot(1, [snp()], label="after August refresh")
        assert ledger.get_snapshot(sid)["label"] == "after August refresh"

    def test_snapshot_skips_a_finding_with_no_identity(self):
        sid = ledger.snapshot(1, [snp(), {"entity_type": "snp", "rsid": ""}])
        assert ledger.get_snapshot(sid)["finding_count"] == 1

    def test_snapshot_tolerates_a_non_dict_in_the_list(self):
        sid = ledger.snapshot(1, [snp(), None, "junk"])
        assert ledger.get_snapshot(sid)["finding_count"] == 1

    def test_snapshot_of_an_empty_finding_list_is_legal(self):
        sid = ledger.snapshot(1, [])
        assert ledger.get_snapshot(sid)["finding_count"] == 0

    def test_get_snapshot_returns_none_for_an_unknown_id(self):
        assert ledger.get_snapshot(99999) is None

    def test_get_snapshot_exposes_the_entries_by_key(self):
        sid = ledger.snapshot(1, [snp("rs1")])
        assert "snp:rs1" in ledger.get_snapshot(sid)["entries"]

    def test_get_snapshot_carries_gene_alongside_the_fingerprint(self):
        sid = ledger.snapshot(1, [snp("rs1", gene="MTHFR")])
        assert ledger.get_snapshot(sid)["genes"]["snp:rs1"] == "MTHFR"

    def test_latest_snapshot_is_none_before_any_scan(self):
        assert ledger.latest_snapshot(1) is None

    def test_latest_snapshot_returns_the_newest(self):
        ledger.snapshot(1, [snp()])
        second = ledger.snapshot(1, [snp()])
        assert ledger.latest_snapshot(1)["id"] == second

    def test_snapshots_are_isolated_per_profile(self):
        ledger.snapshot(1, [snp()])
        ledger.snapshot(2, [snp()])
        assert len(ledger.list_snapshots(1)) == 1
        assert len(ledger.list_snapshots(2)) == 1

    def test_list_snapshots_is_newest_first(self):
        ids = [ledger.snapshot(1, [snp()]) for _ in range(4)]
        assert [h["id"] for h in ledger.list_snapshots(1)] == list(reversed(ids))

    def test_list_snapshots_returns_headers_without_entries(self):
        ledger.snapshot(1, [snp()])
        assert "entries" not in ledger.list_snapshots(1)[0]

    def test_list_snapshots_is_empty_for_an_unknown_profile(self):
        assert ledger.list_snapshots(4242) == []

    def test_a_snapshot_is_never_rewritten_by_a_later_one(self):
        first = ledger.snapshot(1, [vus("rs9")], db_versions={"clinvar": "2026-07"})
        before = json.dumps(ledger.get_snapshot(first), sort_keys=True)
        ledger.snapshot(1, [snp("rs9")], db_versions={"clinvar": "2026-08"})
        after = json.dumps(ledger.get_snapshot(first), sort_keys=True)
        assert before == after

    def test_the_module_exposes_no_snapshot_mutator(self):
        # Additive by construction. If a writer for prior snapshots ever
        # appears, this test is the thing that has to be argued with first.
        assert not [name for name in dir(ledger)
                    if name.startswith(("update_", "edit_", "rewrite_",
                                        "amend_", "supersede_"))]


# ---------------------------------------------------------------------------
# Baseline
# ---------------------------------------------------------------------------

class TestBaseline:

    def test_diff_against_none_is_a_baseline_not_an_error(self):
        sid = ledger.snapshot(1, [snp()])
        diff = ledger.diff_snapshots(None, sid)
        assert diff["ok"] is True and diff["baseline"] is True

    def test_baseline_payload_has_no_changes(self):
        sid = ledger.snapshot(1, [snp()])
        diff = ledger.diff_snapshots(None, sid)
        assert diff["changes"] == [] and diff["material_count"] == 0

    def test_baseline_payload_counts_the_new_side(self):
        sid = ledger.snapshot(1, [snp("rs1"), snp("rs2")])
        diff = ledger.diff_snapshots(None, sid)
        assert diff["total_new"] == 2 and diff["total_old"] == 0

    def test_baseline_payload_says_baseline_established(self):
        sid = ledger.snapshot(1, [snp()])
        assert "aseline" in ledger.diff_snapshots(None, sid)["version_statement"]

    def test_addendum_on_a_single_snapshot_is_a_baseline(self):
        ledger.snapshot(1, [snp()])
        payload = ledger.addendum(1)
        assert payload["ok"] is True and payload["kind"] == "baseline"

    def test_changes_for_on_a_single_snapshot_is_a_baseline(self):
        ledger.snapshot(1, [snp()])
        assert ledger.changes_for(1)["baseline"] is True

    def test_changes_for_with_no_snapshots_at_all_is_well_formed(self):
        payload = ledger.changes_for(1)
        assert payload["ok"] is True and payload["snapshots"] == 0
        assert payload["changes"] == []

    def test_addendum_with_no_snapshots_reports_it_without_raising(self):
        payload = ledger.addendum(1)
        assert payload["ok"] is False and payload["error"] == "no_snapshots"

    def test_diff_with_an_unknown_new_id_names_the_field(self):
        diff = ledger.diff_snapshots(None, 99999)
        assert diff["ok"] is False and diff["field"] == "new_id"

    def test_diff_with_an_unknown_old_id_names_the_field(self):
        sid = ledger.snapshot(1, [snp()])
        diff = ledger.diff_snapshots(99999, sid)
        assert diff["ok"] is False and diff["field"] == "old_id"


# ---------------------------------------------------------------------------
# Change detection
# ---------------------------------------------------------------------------

class TestChangeDetection:

    def _diff(self, before, after):
        old = ledger.snapshot(1, before, db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, after, db_versions={"clinvar": "2026-08"})
        return ledger.diff_snapshots(old, new)

    def test_vus_to_pathogenic_is_detected(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=5)])
        assert "vus_resolved_pathogenic" in kinds_of(diff)

    def test_vus_to_pathogenic_direction_is_up(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=5)])
        assert change_of(diff, "vus_resolved_pathogenic")["direction"] == "up"

    def test_vus_to_likely_pathogenic_is_also_a_resolution(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=4)])
        assert "vus_resolved_pathogenic" in kinds_of(diff)

    def test_vus_to_benign_is_detected(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=2)])
        assert "vus_resolved_benign" in kinds_of(diff)

    def test_vus_to_benign_direction_is_down(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=2)])
        assert change_of(diff, "vus_resolved_benign")["direction"] == "down"

    def test_likely_pathogenic_to_pathogenic_is_an_upgrade(self):
        diff = self._diff([snp("rs9", clinvar_sig_code=4)],
                          [snp("rs9", clinvar_sig_code=5)])
        assert change_of(diff, "sig_upgraded")["direction"] == "up"

    def test_pathogenic_to_likely_benign_is_a_downgrade(self):
        diff = self._diff([snp("rs9", clinvar_sig_code=5)],
                          [snp("rs9", clinvar_sig_code=3)])
        assert change_of(diff, "sig_downgraded")["direction"] == "down"

    def test_pathogenic_to_conflicting_is_reclassified_not_a_downgrade(self):
        diff = self._diff([snp("rs9", clinvar_sig_code=5)],
                          [snp("rs9", clinvar_sig_code=255)])
        assert "sig_reclassified" in kinds_of(diff)
        assert "sig_downgraded" not in kinds_of(diff)

    def test_a_significance_change_names_both_positions(self):
        diff = self._diff([vus("rs9")], [snp("rs9", clinvar_sig_code=5)])
        change = change_of(diff, "vus_resolved_pathogenic")
        assert change["old_display"] == "uncertain significance"
        assert change["new_display"] == "pathogenic"

    def test_stars_gained_is_detected(self):
        diff = self._diff([snp("rs9", review_stars=1)],
                          [snp("rs9", review_stars=3)])
        change = change_of(diff, "stars_gained")
        assert change["old"] == 1 and change["new"] == 3
        assert change["direction"] == "up"

    def test_stars_lost_is_detected(self):
        diff = self._diff([snp("rs9", review_stars=3)],
                          [snp("rs9", review_stars=1)])
        assert change_of(diff, "stars_lost")["direction"] == "down"

    def test_cpic_level_change_is_detected_with_a_direction(self):
        diff = self._diff([snp("rs9", cpic_level="B")],
                          [snp("rs9", cpic_level="A")])
        change = change_of(diff, "cpic_level_changed")
        assert change["old"] == "B" and change["new"] == "A"
        assert change["direction"] == "up"

    def test_losing_a_cpic_level_is_a_downward_move(self):
        diff = self._diff([snp("rs9", cpic_level="A")],
                          [snp("rs9", cpic_level="")])
        assert change_of(diff, "cpic_level_changed")["direction"] == "down"

    def test_magnitude_change_is_detected(self):
        diff = self._diff([snp("rs9", magnitude=2.0)],
                          [snp("rs9", magnitude=6.0)])
        change = change_of(diff, "magnitude_changed")
        assert change["direction"] == "up" and change["material"] is False

    def test_repute_change_to_bad_is_an_upward_move(self):
        diff = self._diff([snp("rs9", repute="")], [snp("rs9", repute="Bad")])
        assert change_of(diff, "repute_changed")["direction"] == "up"

    def test_repute_cleared_from_bad_is_a_downward_move(self):
        diff = self._diff([snp("rs9", repute="Bad")], [snp("rs9", repute="")])
        assert change_of(diff, "repute_changed")["direction"] == "down"

    def test_newly_evaluable_is_detected(self):
        diff = self._diff([genoset("gs1", coverage=0.5)],
                          [genoset("gs1", coverage=1.0)])
        assert change_of(diff, "newly_evaluable")["direction"] == "up"

    def test_no_longer_evaluable_is_detected(self):
        diff = self._diff([snp("rs9", zygosity="heterozygous")],
                          [snp("rs9", zygosity="no_call")])
        assert change_of(diff, "no_longer_evaluable")["direction"] == "down"

    def test_finding_added_is_detected(self):
        diff = self._diff([snp("rs1")], [snp("rs1"), snp("rs2")])
        change = change_of(diff, "finding_added")
        assert change["rsid"] == "rs2" and change["direction"] == "up"

    def test_a_removed_finding_is_reported_rather_than_silently_vanishing(self):
        diff = self._diff([snp("rs1"), snp("rs2")], [snp("rs1")])
        change = change_of(diff, "finding_removed")
        assert change["rsid"] == "rs2"
        assert change["direction"] == "down"
        assert change["old_display"] == "present"
        assert change["new_display"] == "absent"

    def test_a_removed_finding_is_material(self):
        diff = self._diff([snp("rs1"), snp("rs2")], [snp("rs1")])
        assert change_of(diff, "finding_removed")["material"] is True

    def test_carrier_status_change_is_detected(self):
        diff = self._diff([snp("rs9", variant_copies=0)],
                          [snp("rs9", variant_copies=1)])
        change = change_of(diff, "carrier_status_changed")
        assert change["direction"] == "up"
        assert change["new_display"] == "1 copy"

    def test_every_required_change_kind_is_reachable_from_a_real_diff(self):
        before = [
            vus("rs_up"), vus("rs_down"),
            snp("rs_sigup", clinvar_sig_code=4),
            snp("rs_sigdown", clinvar_sig_code=5),
            snp("rs_stars_up", review_stars=1),
            snp("rs_stars_down", review_stars=3),
            snp("rs_cpic", cpic_level="B"),
            snp("rs_mag", magnitude=2.0),
            snp("rs_repute", repute=""),
            genoset("gs_eval", coverage=0.5),
            snp("rs_uneval", zygosity="heterozygous"),
            snp("rs_gone"),
        ]
        after = [
            snp("rs_up", clinvar_sig_code=5), snp("rs_down", clinvar_sig_code=2),
            snp("rs_sigup", clinvar_sig_code=5),
            snp("rs_sigdown", clinvar_sig_code=4),
            snp("rs_stars_up", review_stars=3),
            snp("rs_stars_down", review_stars=1),
            snp("rs_cpic", cpic_level="A"),
            snp("rs_mag", magnitude=7.0),
            snp("rs_repute", repute="Bad"),
            genoset("gs_eval", coverage=1.0),
            snp("rs_uneval", zygosity="no_call"),
            snp("rs_new"),
        ]
        found = kinds_of(self._diff(before, after))
        for kind in TestChangeVocabulary.REQUIRED:
            assert kind in found, kind

    def test_no_change_produces_no_change_records(self):
        diff = self._diff([snp("rs1")], [snp("rs1")])
        assert diff["changes"] == [] and diff["unchanged"] == 1

    def test_prose_churn_alone_produces_no_change_records(self):
        diff = self._diff([snp("rs1", interpretation="old", summary="old")],
                          [snp("rs1", interpretation="new", summary="new")])
        assert diff["changes"] == []

    def test_no_change_record_ever_says_only_changed(self):
        before = [vus("rs1"), snp("rs2", cpic_level="B"), snp("rs3", repute="")]
        after = [snp("rs1", clinvar_sig_code=5), snp("rs2", cpic_level="A"),
                 snp("rs3", repute="Bad")]
        for change in self._diff(before, after)["changes"]:
            assert change["direction"] in ("up", "down", "lateral")
            assert change["old_display"] and change["new_display"]

    def test_changes_are_sorted_by_severity_descending(self):
        before = [snp("rs_mag", magnitude=1.0), vus("rs_vus")]
        after = [snp("rs_mag", magnitude=9.0), snp("rs_vus", clinvar_sig_code=5)]
        severities = [c["severity"] for c in self._diff(before, after)["changes"]]
        assert severities == sorted(severities, reverse=True)

    def test_the_headline_leads_with_the_most_severe_material_change(self):
        before = [snp("rs_mag", magnitude=1.0), vus("rs_vus")]
        after = [snp("rs_mag", magnitude=9.0), snp("rs_vus", clinvar_sig_code=5)]
        headline = self._diff(before, after)["headline"]
        assert "rs_vus" in headline and "pathogenic" in headline.lower()

    def test_counts_tally_by_kind(self):
        diff = self._diff([vus("rs1"), vus("rs2")],
                          [snp("rs1", clinvar_sig_code=5),
                           snp("rs2", clinvar_sig_code=5)])
        assert diff["counts"]["vus_resolved_pathogenic"] == 2

    def test_material_excludes_the_magnitude_only_change(self):
        diff = self._diff([snp("rs1", magnitude=1.0)],
                          [snp("rs1", magnitude=9.0)])
        assert diff["counts"]["magnitude_changed"] == 1
        assert diff["material"] == []

    def test_diff_is_deterministic_across_repeated_calls(self):
        old = ledger.snapshot(1, [vus("rs1"), snp("rs2")])
        new = ledger.snapshot(1, [snp("rs1", clinvar_sig_code=5), snp("rs3")])
        first = ledger.diff_snapshots(old, new)["changes"]
        second = ledger.diff_snapshots(old, new)["changes"]
        assert first == second


# ---------------------------------------------------------------------------
# Database versions in the diff
# ---------------------------------------------------------------------------

class TestDatabaseVersionsInDiff:

    def test_a_version_bump_names_both_versions(self):
        old = ledger.snapshot(1, [vus("rs9")], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp("rs9")], db_versions={"clinvar": "2026-08"})
        statement = ledger.diff_snapshots(old, new)["version_statement"]
        assert "2026-07" in statement and "2026-08" in statement

    def test_the_version_statement_names_the_source(self):
        old = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-08"})
        assert "ClinVar" in ledger.diff_snapshots(old, new)["version_statement"]

    def test_the_changed_version_map_records_old_and_new(self):
        old = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-08"})
        changed = ledger.diff_snapshots(old, new)["db_versions"]["changed"]
        assert changed["clinvar"] == {"old": "2026-07", "new": "2026-08",
                                      "name": "ClinVar"}

    def test_an_unchanged_version_is_not_reported_as_changed(self):
        old = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        diff = ledger.diff_snapshots(old, new)
        assert diff["db_versions"]["changed"] == {}
        assert "No database version change" in diff["version_statement"]

    def test_multiple_source_bumps_are_all_named(self):
        old = ledger.snapshot(1, [snp()],
                              db_versions={"clinvar": "2026-07", "cpic": "1"})
        new = ledger.snapshot(1, [snp()],
                              db_versions={"clinvar": "2026-08", "cpic": "2"})
        statement = ledger.diff_snapshots(old, new)["version_statement"]
        assert "ClinVar" in statement and "CPIC" in statement

    def test_both_version_maps_survive_the_round_trip(self):
        old = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp()], db_versions={"clinvar": "2026-08"})
        diff = ledger.diff_snapshots(old, new)
        assert diff["db_versions"]["old"] == {"clinvar": "2026-07"}
        assert diff["db_versions"]["new"] == {"clinvar": "2026-08"}


# ---------------------------------------------------------------------------
# changes_for
# ---------------------------------------------------------------------------

class TestChangesFor:

    def _history(self):
        a = ledger.snapshot(1, [vus("rs1")], db_versions={"clinvar": "2026-06"})
        b = ledger.snapshot(1, [snp("rs1", clinvar_sig_code=5)],
                            db_versions={"clinvar": "2026-07"})
        c = ledger.snapshot(1, [snp("rs1", clinvar_sig_code=5, review_stars=4)],
                            db_versions={"clinvar": "2026-08"})
        return a, b, c

    def test_changes_for_walks_every_consecutive_pair(self):
        self._history()
        payload = ledger.changes_for(1)
        assert len(payload["comparisons"]) == 2

    def test_changes_for_finds_changes_from_both_pairs(self):
        self._history()
        found = {c["kind"] for c in ledger.changes_for(1)["changes"]}
        assert "vus_resolved_pathogenic" in found and "stars_gained" in found

    def test_each_change_records_the_snapshot_pair_it_came_from(self):
        a, b, _ = self._history()
        change = next(c for c in ledger.changes_for(1)["changes"]
                      if c["kind"] == "vus_resolved_pathogenic")
        assert change["old_snapshot_id"] == a and change["new_snapshot_id"] == b

    def test_each_change_carries_its_own_version_statement(self):
        self._history()
        change = next(c for c in ledger.changes_for(1)["changes"]
                      if c["kind"] == "vus_resolved_pathogenic")
        assert "2026-06" in change["version_statement"]
        assert "2026-07" in change["version_statement"]

    def test_since_filters_out_earlier_comparisons(self):
        _, b, _ = self._history()
        cutoff = ledger.get_snapshot(b)["created_at"]
        payload = ledger.changes_for(1, since=cutoff)
        assert len(payload["comparisons"]) == 1

    def test_since_in_the_future_returns_nothing(self):
        self._history()
        assert ledger.changes_for(1, since="2099-01-01")["changes"] == []

    def test_since_accepts_a_bare_date(self):
        self._history()
        assert ledger.changes_for(1, since="2000-01-01")["comparisons"]

    def test_limit_caps_the_returned_records(self):
        self._history()
        payload = ledger.changes_for(1, limit=1)
        assert len(payload["changes"]) == 1 and payload["truncated"] is True

    def test_limit_keeps_the_most_severe_record(self):
        self._history()
        payload = ledger.changes_for(1, limit=1)
        assert payload["changes"][0]["kind"] == "vus_resolved_pathogenic"

    def test_a_limit_above_the_count_does_not_truncate(self):
        self._history()
        assert ledger.changes_for(1, limit=500)["truncated"] is False

    def test_counts_are_computed_before_truncation(self):
        self._history()
        payload = ledger.changes_for(1, limit=1)
        assert sum(payload["counts"].values()) > len(payload["changes"])


# ---------------------------------------------------------------------------
# Addendum
# ---------------------------------------------------------------------------

class TestAddendum:

    def _two(self):
        old = ledger.snapshot(1, [vus("rs9")], db_versions={"clinvar": "2026-07"})
        new = ledger.snapshot(1, [snp("rs9", clinvar_sig_code=5)],
                              db_versions={"clinvar": "2026-08"})
        return old, new

    def test_addendum_defaults_to_the_two_newest_snapshots(self):
        old, new = self._two()
        payload = ledger.addendum(1)
        assert payload["old"]["id"] == old and payload["new"]["id"] == new

    def test_addendum_is_dated(self):
        self._two()
        assert ledger.addendum(1)["generated_at"].endswith("Z")

    def test_addendum_declares_itself_additive(self):
        self._two()
        assert ledger.addendum(1)["additive"] is True

    def test_addendum_supersedes_nothing(self):
        self._two()
        assert ledger.addendum(1)["supersedes"] is None

    def test_addendum_statement_says_it_does_not_replace_the_original(self):
        self._two()
        statement = ledger.addendum(1)["statement"].lower()
        assert "does not replace" in statement and "supersede" in statement

    def test_addendum_does_not_mutate_the_prior_snapshot(self):
        old, _ = self._two()
        before = json.dumps(ledger.get_snapshot(old), sort_keys=True)
        ledger.addendum(1)
        ledger.addendum(1)
        after = json.dumps(ledger.get_snapshot(old), sort_keys=True)
        assert before == after

    def test_addendum_does_not_add_or_remove_snapshots(self):
        self._two()
        before = len(ledger.list_snapshots(1))
        ledger.addendum(1)
        assert len(ledger.list_snapshots(1)) == before

    def test_addendum_carries_the_material_changes(self):
        self._two()
        payload = ledger.addendum(1)
        assert payload["material_count"] >= 1
        assert payload["material"][0]["kind"] == "vus_resolved_pathogenic"

    def test_addendum_accepts_explicit_ids(self):
        old, new = self._two()
        payload = ledger.addendum(1, old_id=old, new_id=new)
        assert payload["kind"] == "addendum"

    def test_addendum_for_the_oldest_snapshot_is_a_baseline(self):
        old, _ = self._two()
        assert ledger.addendum(1, new_id=old)["kind"] == "baseline"

    def test_addendum_names_both_database_versions(self):
        self._two()
        statement = ledger.addendum(1)["version_statement"]
        assert "2026-07" in statement and "2026-08" in statement


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

class TestAddendumHtml:

    def _html(self):
        ledger.snapshot(1, [vus("rs9"), snp("rs_mag", magnitude=1.0)],
                        db_versions={"clinvar": "2026-07"})
        ledger.snapshot(1, [snp("rs9", clinvar_sig_code=5),
                            snp("rs_mag", magnitude=9.0)],
                        db_versions={"clinvar": "2026-08"})
        return ledger.render_addendum_html(ledger.addendum(1))

    def test_html_contains_no_http_reference_at_all(self):
        assert "http" not in self._html().lower()

    def test_html_contains_no_external_resource_tags(self):
        html = self._html().lower()
        for tag in ("<link", "<script", "<img", "@import", "src=", "srcset="):
            assert tag not in html

    def test_html_is_a_fragment_not_a_document(self):
        html = self._html()
        assert html.startswith("<section")
        assert "<!doctype" not in html.lower()

    def test_html_carries_its_own_inline_style(self):
        assert "<style>" in self._html()

    def test_html_states_the_additive_rule(self):
        assert "does not replace" in self._html()

    def test_html_names_both_database_versions(self):
        html = self._html()
        assert "2026-07" in html and "2026-08" in html

    def test_html_shows_the_material_change(self):
        assert "rs9" in self._html()

    def test_html_notes_the_suppressed_magnitude_only_changes(self):
        assert "derived score" in self._html()

    def test_baseline_html_renders_without_a_table(self):
        ledger.snapshot(1, [snp()])
        html = ledger.render_addendum_html(ledger.addendum(1))
        assert "Baseline established" in html and "<table>" not in html

    def test_html_of_an_unchanged_diff_says_nothing_changed(self):
        ledger.snapshot(1, [snp("rs1")])
        ledger.snapshot(1, [snp("rs1")])
        html = ledger.render_addendum_html(ledger.addendum(1))
        assert "No finding changed" in html

    def test_html_escapes_markup_from_a_gene_symbol(self):
        ledger.snapshot(1, [vus("rs9", gene="<script>x</script>")])
        ledger.snapshot(1, [snp("rs9", clinvar_sig_code=5,
                                gene="<script>x</script>")])
        html = ledger.render_addendum_html(ledger.addendum(1))
        assert "<script>" not in html and "&lt;script&gt;" in html

    def test_html_tolerates_an_empty_payload(self):
        assert ledger.render_addendum_html({}).startswith("<section")

    def test_html_tolerates_a_non_dict_payload(self):
        assert ledger.render_addendum_html(None).startswith("<section")


# ---------------------------------------------------------------------------
# Pruning
# ---------------------------------------------------------------------------

class TestPruning:

    def test_pruning_keeps_the_newest_n(self):
        ids = [ledger.snapshot(1, [snp()]) for _ in range(10)]
        ledger.prune_snapshots(1, keep=3)
        assert [h["id"] for h in ledger.list_snapshots(1)] == ids[-3:][::-1]

    def test_pruning_returns_the_number_deleted(self):
        for _ in range(10):
            ledger.snapshot(1, [snp()])
        assert ledger.prune_snapshots(1, keep=4) == 6

    def test_pruning_below_the_threshold_deletes_nothing(self):
        for _ in range(3):
            ledger.snapshot(1, [snp()])
        assert ledger.prune_snapshots(1, keep=12) == 0

    def test_pruning_defaults_to_keeping_twelve(self):
        for _ in range(15):
            ledger.snapshot(1, [snp()])
        ledger.prune_snapshots(1)
        assert len(ledger.list_snapshots(1)) == 12

    def test_pruning_never_deletes_everything_even_when_asked_to(self):
        for _ in range(5):
            ledger.snapshot(1, [snp()])
        ledger.prune_snapshots(1, keep=0)
        assert len(ledger.list_snapshots(1)) == 1

    def test_pruning_with_a_negative_keep_still_leaves_one(self):
        for _ in range(5):
            ledger.snapshot(1, [snp()])
        ledger.prune_snapshots(1, keep=-100)
        assert len(ledger.list_snapshots(1)) == 1

    def test_pruning_one_profile_leaves_another_alone(self):
        for _ in range(5):
            ledger.snapshot(1, [snp()])
            ledger.snapshot(2, [snp()])
        ledger.prune_snapshots(1, keep=1)
        assert len(ledger.list_snapshots(2)) == 5

    def test_pruning_removes_the_entry_rows_too(self):
        ids = [ledger.snapshot(1, [snp()]) for _ in range(3)]
        ledger.prune_snapshots(1, keep=1)
        conn = db.get_connection()
        try:
            orphans = conn.execute(
                "SELECT COUNT(*) FROM reclass_entries WHERE snapshot_id IN (?,?)",
                (ids[0], ids[1]),
            ).fetchone()[0]
        finally:
            conn.close()
        assert orphans == 0

    def test_pruning_an_empty_profile_is_a_no_op(self):
        assert ledger.prune_snapshots(999, keep=3) == 0

    def test_the_surviving_snapshot_is_still_readable_after_pruning(self):
        for _ in range(5):
            ledger.snapshot(1, [snp("rs1")])
        ledger.prune_snapshots(1, keep=1)
        survivor = ledger.latest_snapshot(1)
        assert survivor is not None and "snp:rs1" in survivor["entries"]
