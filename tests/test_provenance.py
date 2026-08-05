"""
Tests for backend.provenance: the provenance graph, the licence contract and
the signed report manifest.

Isolation follows the same rule as the rest of the suite. ``database.DB_PATH``
is resolved at IMPORT time, so the environment variable alone would let a write
land in the developer's real database; the attribute is patched as well.
``DNAINSIGHT_HOME`` is redirected so no test can read, write or overwrite the
real HMAC key, which would be unrecoverable for any manifest already signed on
that machine.
"""

import hashlib
import json
import os
import stat
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import database as db
from backend import provenance as prov


HELLO_WORLD_SHA256 = "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
EMPTY_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


# ---------------------------------------------------------------------------
# Isolation
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    target = tmp_path / "t.db"
    monkeypatch.setenv("DNAINSIGHT_DB_PATH", str(target))
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(db, "DB_PATH", target)
    yield tmp_path


def finding(**over):
    base = {
        "rsid": "rs1801133",
        "entity_type": "snp",
        "gene": "MTHFR",
        "clinical_sig": "pathogenic",
        "clinvar_sig_code": 5,
        "review_stars": 2,
        "magnitude": 5.0,
        "silo": "actionable",
        "sources": ["bundled_reference"],
    }
    base.update(over)
    return base


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

class TestHashing:

    def test_content_hash_matches_a_known_value(self, tmp_path):
        path = tmp_path / "known.txt"
        path.write_bytes(b"hello world")
        assert prov.content_hash(path) == HELLO_WORLD_SHA256

    def test_content_hash_of_an_empty_file(self, tmp_path):
        path = tmp_path / "empty.txt"
        path.write_bytes(b"")
        assert prov.content_hash(path) == EMPTY_SHA256

    def test_content_hash_of_a_missing_file_returns_none(self, tmp_path):
        assert prov.content_hash(tmp_path / "nope.txt") is None

    def test_content_hash_of_a_directory_returns_none(self, tmp_path):
        assert prov.content_hash(tmp_path) is None

    def test_content_hash_of_none_returns_none(self):
        assert prov.content_hash(None) is None

    def test_content_hash_accepts_a_string_path(self, tmp_path):
        path = tmp_path / "known.txt"
        path.write_bytes(b"hello world")
        assert prov.content_hash(str(path)) == HELLO_WORLD_SHA256

    def test_content_hash_streams_a_file_larger_than_one_chunk(self, tmp_path):
        payload = b"x" * (3 * 1024 * 1024 + 7)
        path = tmp_path / "big.bin"
        path.write_bytes(payload)
        assert prov.content_hash(path) == hashlib.sha256(payload).hexdigest()

    def test_text_hash_matches_a_known_value(self):
        assert prov.text_hash("hello world") == HELLO_WORLD_SHA256

    def test_text_hash_of_an_empty_string(self):
        assert prov.text_hash("") == EMPTY_SHA256

    def test_text_hash_of_none_is_the_empty_hash(self):
        assert prov.text_hash(None) == EMPTY_SHA256

    def test_text_hash_is_utf8_not_platform_dependent(self):
        assert prov.text_hash("café") == hashlib.sha256(
            "café".encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# SOURCES
# ---------------------------------------------------------------------------

class TestSources:

    REQUIRED_KEYS = ("name", "licence", "spdx", "url", "role",
                     "redistributable", "commercial_ok", "verified")

    EXPECTED = ("cpic", "clinvar", "gnomad", "onekg_ensembl", "gwas_catalog",
                "pgs_catalog", "myvariant", "snpedia", "pharmgkb")

    def test_every_documented_source_is_present(self):
        for source_id in self.EXPECTED:
            assert source_id in prov.SOURCES, source_id

    def test_every_source_carries_every_required_key(self):
        for source_id, entry in prov.SOURCES.items():
            for key in self.REQUIRED_KEYS:
                assert key in entry, f"{source_id} missing {key}"

    def test_every_source_carries_a_verification_date(self):
        for source_id, entry in prov.SOURCES.items():
            assert len(entry["verified"]) == 10, source_id

    def test_cpic_is_cc0(self):
        assert prov.SOURCES["cpic"]["spdx"] == "CC0-1.0"

    def test_gnomad_is_cc0(self):
        assert prov.SOURCES["gnomad"]["spdx"] == "CC0-1.0"

    def test_clinvar_is_us_public_domain_with_no_spdx_identifier(self):
        entry = prov.SOURCES["clinvar"]
        assert entry["spdx"] == ""
        assert "17 USC 105" in entry["licence"]

    def test_onekg_is_open_with_citation_requested(self):
        assert "citation" in prov.SOURCES["onekg_ensembl"]["licence"].lower()

    def test_gwas_catalog_is_governed_by_embl_ebi_terms(self):
        assert "EMBL-EBI" in prov.SOURCES["gwas_catalog"]["licence"]

    def test_pgs_catalog_records_per_score_overrides(self):
        entry = prov.SOURCES["pgs_catalog"]
        assert entry["per_record"] is True
        assert "EMBL-EBI" in entry["licence"]

    def test_pgs_catalog_note_names_the_no_derivatives_problem(self):
        assert "CC BY-NC-ND" in prov.SOURCES["pgs_catalog"]["note"]

    def test_myvariant_code_is_apache_and_data_is_per_field(self):
        entry = prov.SOURCES["myvariant"]
        assert entry["spdx"] == "Apache-2.0"
        assert entry["per_record"] is True

    def test_myvariant_is_never_persisted(self):
        entry = prov.SOURCES["myvariant"]
        assert entry["never_bundle"] is True
        assert entry["redistributable"] is False

    def test_snpedia_is_cc_by_nc_sa(self):
        assert prov.SOURCES["snpedia"]["spdx"] == "CC-BY-NC-SA-3.0-US"

    def test_snpedia_is_never_bundled_and_not_commercial(self):
        entry = prov.SOURCES["snpedia"]
        assert entry["never_bundle"] is True
        assert entry["commercial_ok"] is False

    def test_pharmgkb_is_present_and_marked_not_used(self):
        entry = prov.SOURCES["pharmgkb"]
        assert entry["used"] is False
        assert entry["role"] == "not_used"

    def test_pharmgkb_records_the_reason_for_the_exclusion(self):
        note = prov.SOURCES["pharmgkb"]["note"].lower()
        licence = prov.SOURCES["pharmgkb"]["licence"].lower()
        assert "prohibiting sale" in licence
        assert "clinpgxlevel" in note

    def test_pharmgkb_publishes_no_usable_spdx_identifier(self):
        # It claims CC-BY-SA-4.0, but the data use agreement adds a term the
        # identifier cannot carry, so recording the identifier would be a lie.
        assert prov.SOURCES["pharmgkb"]["spdx"] == ""

    def test_pharmgkb_is_not_redistributable_or_commercial(self):
        entry = prov.SOURCES["pharmgkb"]
        assert entry["redistributable"] is False
        assert entry["commercial_ok"] is False


# ---------------------------------------------------------------------------
# source_record
# ---------------------------------------------------------------------------

class TestSourceRecord:

    def test_a_record_carries_the_licence_and_spdx(self):
        record = prov.source_record("cpic")
        assert record["licence"] and record["spdx"] == "CC0-1.0"

    def test_a_record_carries_the_name_url_and_role(self):
        record = prov.source_record("clinvar")
        assert record["name"] == "ClinVar"
        assert record["url"] and record["role"]

    def test_a_record_is_marked_known(self):
        assert prov.source_record("cpic")["known"] is True

    def test_an_explicit_version_wins_over_the_live_one(self):
        assert prov.source_record("clinvar", version="2026-08-24")["version"] == "2026-08-24"

    def test_an_omitted_version_falls_back_to_the_live_one(self):
        live = prov.database_versions()["clinvar"]
        assert prov.source_record("clinvar")["version"] == live

    def test_an_omitted_retrieval_date_defaults_to_today(self):
        assert len(prov.source_record("cpic")["retrieved"]) == 10

    def test_an_explicit_retrieval_date_is_kept(self):
        assert prov.source_record("cpic", retrieved="2026-01-01")["retrieved"] == "2026-01-01"

    def test_an_unknown_source_fails_closed(self):
        record = prov.source_record("some_new_database")
        assert record["known"] is False
        assert record["redistributable"] is False
        assert record["commercial_ok"] is False
        assert record["never_bundle"] is True

    def test_an_unknown_source_does_not_raise(self):
        assert prov.source_record("")["known"] is False

    def test_every_known_source_produces_a_record(self):
        for source_id in prov.SOURCES:
            assert prov.source_record(source_id)["known"] is True


# ---------------------------------------------------------------------------
# database_versions
# ---------------------------------------------------------------------------

class TestDatabaseVersions:

    def test_every_source_id_appears(self):
        versions = prov.database_versions()
        for source_id in prov.SOURCES:
            assert source_id in versions

    def test_every_value_is_a_string(self):
        assert all(isinstance(v, str) for v in prov.database_versions().values())

    def test_the_application_version_is_recorded(self):
        from backend import APP_VERSION
        assert prov.database_versions()["dnainsight"] == APP_VERSION

    def test_myvariant_is_recorded_as_live_not_as_a_date(self):
        assert prov.database_versions()["myvariant"] == "live"

    def test_pharmgkb_carries_no_version_because_it_is_not_used(self):
        assert prov.database_versions()["pharmgkb"] == ""

    def test_the_mapping_round_trips_through_json_unchanged(self):
        versions = prov.database_versions()
        assert json.loads(json.dumps(versions)) == versions

    def test_repeated_calls_agree(self):
        assert prov.database_versions() == prov.database_versions()


# ---------------------------------------------------------------------------
# attach
# ---------------------------------------------------------------------------

class TestAttach:

    def test_attach_adds_a_provenance_block(self):
        f = prov.attach(finding(), ["clinvar", "cpic"])
        assert len(f["provenance"]["sources"]) == 2

    def test_attach_returns_the_same_dict(self):
        original = finding()
        assert prov.attach(original, ["clinvar"]) is original

    def test_attach_records_the_source_ids(self):
        f = prov.attach(finding(), ["clinvar", "cpic"])
        assert f["source_ids"] == ["clinvar", "cpic"]

    def test_attach_accepts_a_single_string(self):
        assert prov.attach(finding(), "clinvar")["source_ids"] == ["clinvar"]

    def test_attach_deduplicates_source_ids(self):
        f = prov.attach(finding(), ["clinvar", "clinvar", "cpic"])
        assert f["source_ids"] == ["clinvar", "cpic"]

    def test_attach_does_not_clobber_the_existing_sources_list(self):
        f = prov.attach(finding(sources=["bundled_reference"]), ["clinvar"])
        assert f["sources"] == ["bundled_reference"]

    def test_attach_collects_the_licences(self):
        f = prov.attach(finding(), ["clinvar", "cpic"])
        assert len(f["provenance"]["licences"]) == 2

    def test_attach_marks_a_clean_finding_redistributable(self):
        assert prov.attach(finding(), ["clinvar", "cpic"])["provenance"]["redistributable"] is True

    def test_one_contaminating_source_marks_the_whole_finding(self):
        block = prov.attach(finding(), ["clinvar", "snpedia"])["provenance"]
        assert block["redistributable"] is False
        assert block["commercial_ok"] is False

    def test_a_never_bundle_source_makes_the_finding_unpersistable(self):
        block = prov.attach(finding(), ["clinvar", "myvariant"])["provenance"]
        assert block["persistable"] is False

    def test_attach_lists_unknown_sources(self):
        block = prov.attach(finding(), ["clinvar", "mystery_db"])["provenance"]
        assert block["unknown_sources"] == ["mystery_db"]

    def test_attach_honours_a_version_override(self):
        f = prov.attach(finding(), ["clinvar"], versions={"clinvar": "2026-08"})
        assert f["provenance"]["sources"][0]["version"] == "2026-08"

    def test_attach_records_a_timestamp(self):
        assert prov.attach(finding(), ["clinvar"])["provenance"]["attached_at"].endswith("Z")

    def test_attach_with_no_sources_is_legal(self):
        block = prov.attach(finding(), [])["provenance"]
        assert block["sources"] == [] and block["redistributable"] is True

    def test_attach_tolerates_a_non_dict(self):
        assert prov.attach("not a finding", ["clinvar"]) == "not a finding"


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------

class TestConflictDetection:

    def test_clinvar_benign_versus_a_replicated_gwas_hit_is_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2,
            gwas_studies=4, gwas_traits="Type 2 diabetes"))
        assert any(c["type"] == "clinvar_benign_vs_gwas_association"
                   for c in conflicts)

    def test_that_conflict_carries_both_positions(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, gwas_studies=4))
        positions = conflicts[0]["positions"]
        assert {p["source_id"] for p in positions} == {"clinvar", "gwas_catalog"}

    def test_that_conflict_carries_no_verdict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, gwas_studies=4))
        assert conflicts[0]["verdict"] is None
        assert conflicts[0]["resolved"] is False

    def test_every_conflict_record_is_display_only(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, gwas_studies=4,
            cpic_level="A"))
        assert conflicts and all(c["display_only"] for c in conflicts)

    def test_a_gwas_source_tag_counts_as_replication(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2,
            sources=["gwas_catalog"]))
        assert any(c["type"] == "clinvar_benign_vs_gwas_association"
                   for c in conflicts)

    def test_a_single_gwas_study_is_not_replication(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, gwas_studies=1,
            sources=["bundled_reference"]))
        assert conflicts == []

    def test_cpic_actionable_versus_clinvar_benign_is_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, cpic_level="A"))
        types = {c["type"] for c in conflicts}
        assert "cpic_actionable_vs_clinvar_benign" in types

    def test_that_conflict_names_cpic_and_clinvar(self):
        conflicts = [c for c in prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, cpic_level="A"))
            if c["type"] == "cpic_actionable_vs_clinvar_benign"]
        ids = {p["source_id"] for p in conflicts[0]["positions"]}
        assert ids == {"cpic", "clinvar"}

    def test_cpic_retired_versus_clinvar_pathogenic_is_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="pathogenic", clinvar_sig_code=5,
            cpic_level="Retired"))
        assert any(c["type"] == "cpic_no_action_vs_clinvar_pathogenic"
                   for c in conflicts)

    def test_cpic_and_clinvar_agreeing_is_not_a_conflict(self):
        assert prov.detect_conflicts(finding(
            clinical_sig="pathogenic", clinvar_sig_code=5,
            cpic_level="A")) == []

    def test_a_conflicting_clinvar_record_is_surfaced(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="conflicting classifications of pathogenicity",
            clinvar_sig_code=None))
        assert any(c["type"] == "clinvar_submitters_conflicting"
                   for c in conflicts)

    def test_a_conflicting_clinvar_record_is_not_resolved_to_pathogenic(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="conflicting classifications of pathogenicity",
            clinvar_sig_code=None))
        record = next(c for c in conflicts
                      if c["type"] == "clinvar_submitters_conflicting")
        assert record["verdict"] is None

    def test_snpedia_good_versus_clinvar_pathogenic_is_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="pathogenic", clinvar_sig_code=5,
            snpedia_repute="Good"))
        assert any(c["type"] == "snpedia_repute_vs_clinvar" for c in conflicts)

    def test_snpedia_bad_versus_clinvar_benign_is_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            clinical_sig="benign", clinvar_sig_code=2, snpedia_repute="Bad"))
        assert any(c["type"] == "snpedia_repute_vs_clinvar" for c in conflicts)

    def test_agreeing_snpedia_repute_is_not_a_conflict(self):
        assert prov.detect_conflicts(finding(
            clinical_sig="pathogenic", clinvar_sig_code=5,
            snpedia_repute="Bad")) == []

    def test_disagreeing_input_files_are_surfaced_as_a_conflict(self):
        conflicts = prov.detect_conflicts(finding(
            conflict=True, calls=["AG", "AA"]))
        record = next(c for c in conflicts if c["type"] == "input_files_disagree")
        assert len(record["positions"]) == 2

    def test_both_disagreeing_calls_are_kept(self):
        conflicts = prov.detect_conflicts(finding(
            conflict=True, calls=["AG", "AA"]))
        record = next(c for c in conflicts if c["type"] == "input_files_disagree")
        assert {p["claim"] for p in record["positions"]} == {"AG", "AA"}

    def test_the_input_conflict_states_there_is_no_automatic_winner(self):
        conflicts = prov.detect_conflicts(finding(conflict=True, calls=["AG"]))
        record = next(c for c in conflicts if c["type"] == "input_files_disagree")
        assert "no automatic winner" in record["note"]

    def test_no_conflict_on_a_plain_finding(self):
        assert prov.detect_conflicts(finding()) == []

    def test_conflicts_are_returned_in_a_stable_order(self):
        f = finding(clinical_sig="benign", clinvar_sig_code=2,
                    gwas_studies=3, cpic_level="A", snpedia_repute="Bad")
        assert prov.detect_conflicts(f) == prov.detect_conflicts(f)

    def test_multiple_conflicts_are_all_returned(self):
        f = finding(clinical_sig="benign", clinvar_sig_code=2,
                    gwas_studies=3, cpic_level="A", snpedia_repute="Bad")
        assert len(prov.detect_conflicts(f)) == 3

    def test_detect_conflicts_tolerates_a_non_dict(self):
        assert prov.detect_conflicts(None) == []

    def test_no_conflict_record_ever_carries_a_winner_key(self):
        f = finding(clinical_sig="benign", clinvar_sig_code=2,
                    gwas_studies=3, cpic_level="A")
        for record in prov.detect_conflicts(f):
            assert "winner" not in record and "preferred" not in record


# ---------------------------------------------------------------------------
# Licence audit
# ---------------------------------------------------------------------------

class TestLicenceAudit:

    def test_the_repository_as_declared_is_clean(self):
        audit = prov.licence_audit()
        assert audit["ok"] is True, audit["violations"]

    def test_the_audit_reports_how_much_it_checked(self):
        audit = prov.licence_audit()
        assert audit["artefacts"] >= 1 and audit["checked"] >= audit["artefacts"]

    def test_the_audit_lists_the_not_used_sources(self):
        assert "pharmgkb" in prov.licence_audit()["sources"]["not_used"]

    def test_the_audit_lists_the_never_bundle_sources(self):
        never = prov.licence_audit()["sources"]["never_bundle"]
        assert "snpedia" in never and "myvariant" in never

    def test_a_constructed_pharmgkb_violation_is_caught(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_pharmgkb.json", ["pharmgkb"])
        audit = prov.licence_audit()
        assert audit["ok"] is False
        rules = {v["rule"] for v in audit["violations"]}
        assert "not_used_source_bundled" in rules

    def test_the_pharmgkb_violation_names_the_artefact(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_pharmgkb.json", ["pharmgkb"])
        violation = prov.licence_audit()["violations"][0]
        assert violation["artefact"] == "data/oops_pharmgkb.json"
        assert violation["source_id"] == "pharmgkb"

    def test_a_constructed_snpedia_violation_is_caught(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_snpedia.json", ["snpedia"])
        audit = prov.licence_audit()
        assert audit["ok"] is False
        assert any(v["rule"] == "never_bundle" for v in audit["violations"])

    def test_a_constructed_myvariant_violation_is_caught(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_cache.json", ["myvariant"])
        assert prov.licence_audit()["ok"] is False

    def test_an_unassessed_source_in_a_bundled_artefact_is_a_violation(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/mystery.json", ["some_new_database"])
        audit = prov.licence_audit()
        assert audit["ok"] is False
        assert any(v["rule"] == "unknown_source" for v in audit["violations"])

    def test_a_violation_carries_the_licence_it_breached(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_snpedia.json", ["snpedia"])
        violation = next(v for v in prov.licence_audit()["violations"]
                         if v["source_id"] == "snpedia")
        assert "CC-BY-NC-SA" in violation["licence"]

    def test_a_violation_states_a_reason(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_snpedia.json", ["snpedia"])
        assert all(v["reason"] for v in prov.licence_audit()["violations"])

    def test_removing_the_constructed_violation_restores_a_clean_audit(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/oops_snpedia.json", ["snpedia"])
        assert prov.licence_audit()["ok"] is False
        monkeypatch.undo()
        assert prov.licence_audit()["ok"] is True

    def test_the_per_record_pgs_warning_is_a_warning_not_a_violation(self):
        audit = prov.licence_audit()
        assert any(w["source_id"] == "pgs_catalog" for w in audit["warnings"])
        assert audit["ok"] is True

    def test_a_clean_source_produces_no_violation(self, monkeypatch):
        monkeypatch.setitem(prov.BUNDLED_ARTEFACTS,
                            "data/fine.json", ["cpic", "gnomad"])
        assert prov.licence_audit()["ok"] is True


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

class TestManifest:

    def test_the_manifest_records_the_application_version(self):
        from backend import APP_VERSION
        manifest = prov.build_manifest(profile_id=1, findings=[])
        assert manifest["dnainsight_version"] == APP_VERSION

    def test_the_manifest_records_a_utc_timestamp(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        assert manifest["generated_at"].endswith("Z")

    def test_the_manifest_records_the_finding_count(self):
        manifest = prov.build_manifest(profile_id=1,
                                       findings=[finding(), finding(rsid="rs2")])
        assert manifest["finding_count"] == 2

    def test_the_manifest_records_every_database(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        recorded = {d["source_id"] for d in manifest["databases"]}
        assert recorded == set(prov.SOURCES)

    def test_every_database_entry_carries_name_version_retrieved_and_licence(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        for entry in manifest["databases"]:
            assert entry["name"] and entry["licence"]
            assert "version" in entry and entry["retrieved"]

    def test_the_manifest_records_pharmgkb_as_not_used(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        entry = next(d for d in manifest["databases"]
                     if d["source_id"] == "pharmgkb")
        assert entry["used"] is False

    def test_the_manifest_hashes_every_input_file(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_bytes(b"hello world")
        manifest = prov.build_manifest(profile_id=1, findings=[],
                                       input_files=[str(raw)])
        assert manifest["input_files"][0]["sha256"] == HELLO_WORLD_SHA256

    def test_an_input_file_record_carries_its_size(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_bytes(b"hello world")
        manifest = prov.build_manifest(profile_id=1, findings=[],
                                       input_files=[raw])
        assert manifest["input_files"][0]["bytes"] == 11

    def test_a_missing_input_file_is_recorded_not_dropped(self, tmp_path):
        manifest = prov.build_manifest(profile_id=1, findings=[],
                                       input_files=[tmp_path / "gone.txt"])
        record = manifest["input_files"][0]
        assert record["sha256"] is None and record["present"] is False

    def test_input_files_accept_dict_entries_with_extra_metadata(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_bytes(b"hello world")
        manifest = prov.build_manifest(
            profile_id=1, findings=[],
            input_files=[{"path": str(raw), "role": "self", "upload_id": 7}])
        record = manifest["input_files"][0]
        assert record["role"] == "self" and record["upload_id"] == 7

    def test_the_manifest_records_the_scan_parameters(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        assert manifest["scan_parameters"]["assembly"]

    def test_scan_parameters_can_be_overridden(self):
        manifest = prov.build_manifest(
            profile_id=1, findings=[],
            extra={"scan_parameters": {"network_enabled": True}})
        assert manifest["scan_parameters"]["network_enabled"] is True

    def test_extra_keys_pass_through(self):
        manifest = prov.build_manifest(profile_id=1, findings=[],
                                       extra={"operator": "clinic"})
        assert manifest["extra"]["operator"] == "clinic"

    def test_the_manifest_records_the_report_type(self):
        manifest = prov.build_manifest(profile_id=1, findings=[],
                                       report_type="doctor")
        assert manifest["report_type"] == "doctor"

    def test_the_manifest_counts_entity_types(self):
        manifest = prov.build_manifest(
            profile_id=1,
            findings=[finding(), finding(rsid="gs1", entity_type="genoset")])
        assert manifest["entity_counts"] == {"snp": 1, "genoset": 1}

    def test_the_findings_digest_is_order_independent(self):
        one = prov.build_manifest(profile_id=1,
                                  findings=[finding(), finding(rsid="rs2")])
        two = prov.build_manifest(profile_id=1,
                                  findings=[finding(rsid="rs2"), finding()])
        assert one["findings_digest"] == two["findings_digest"]

    def test_the_findings_digest_moves_when_a_classification_moves(self):
        one = prov.build_manifest(profile_id=1, findings=[finding()])
        two = prov.build_manifest(profile_id=1,
                                  findings=[finding(clinvar_sig_code=1)])
        assert one["findings_digest"] != two["findings_digest"]

    def test_the_manifest_is_json_serialisable(self):
        manifest = prov.build_manifest(profile_id=1, findings=[finding()])
        assert json.loads(json.dumps(manifest))["finding_count"] == 1


# ---------------------------------------------------------------------------
# Key handling
# ---------------------------------------------------------------------------

class TestManifestKey:

    def test_the_key_is_32_bytes(self):
        assert len(prov.manifest_key()) == 32

    def test_the_key_file_is_created_under_dnainsight_home(self, tmp_path):
        prov.manifest_key()
        assert (tmp_path / "home" / "manifest.key").exists()

    def test_the_key_is_created_once_and_reused(self):
        first = prov.manifest_key()
        second = prov.manifest_key()
        assert first == second

    def test_the_key_file_is_not_rewritten_on_reuse(self, tmp_path):
        prov.manifest_key()
        path = tmp_path / "home" / "manifest.key"
        before = path.read_bytes()
        prov.manifest_key()
        assert path.read_bytes() == before

    def test_an_existing_key_file_is_used_verbatim(self, tmp_path):
        home = tmp_path / "home"
        home.mkdir(parents=True, exist_ok=True)
        (home / "manifest.key").write_bytes(b"a" * 32)
        assert prov.manifest_key() == b"a" * 32

    def test_a_zero_length_key_file_is_replaced(self, tmp_path):
        home = tmp_path / "home"
        home.mkdir(parents=True, exist_ok=True)
        (home / "manifest.key").write_bytes(b"")
        assert len(prov.manifest_key()) == 32

    @pytest.mark.skipif(os.name == "nt",
                        reason="POSIX modes are not implemented on Windows")
    def test_the_key_file_is_owner_only(self, tmp_path):
        prov.manifest_key()
        mode = (tmp_path / "home" / "manifest.key").stat().st_mode
        assert not mode & (stat.S_IRWXG | stat.S_IRWXO)

    def test_two_homes_get_two_different_keys(self, tmp_path, monkeypatch):
        first = prov.manifest_key()
        monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "other"))
        assert prov.manifest_key() != first


# ---------------------------------------------------------------------------
# Signing and verification
# ---------------------------------------------------------------------------

class TestSigning:

    def _signed(self):
        return prov.sign_manifest(
            prov.build_manifest(profile_id=1, findings=[finding()],
                                report_type="doctor"))

    def test_a_round_trip_signs_and_verifies(self):
        assert prov.verify_manifest(self._signed())["ok"] is True

    def test_the_verdict_names_the_algorithm(self):
        assert prov.verify_manifest(self._signed())["algorithm"] == "HMAC-SHA256"

    def test_the_signature_block_states_its_scope(self):
        scope = self._signed()["signature"]["scope"].lower()
        assert "not a public-key attestation" in scope
        assert "does not prove authorship" in scope

    def test_the_signed_object_embeds_the_manifest_unchanged(self):
        manifest = prov.build_manifest(profile_id=1, findings=[finding()])
        assert prov.sign_manifest(manifest)["manifest"] == manifest

    def test_the_signature_is_hex(self):
        value = self._signed()["signature"]["value"]
        assert len(value) == 64 and int(value, 16) >= 0

    def test_the_signature_carries_a_signing_timestamp(self):
        assert self._signed()["signature"]["signed_at"].endswith("Z")

    def test_verification_survives_a_json_round_trip(self):
        signed = json.loads(json.dumps(self._signed()))
        assert prov.verify_manifest(signed)["ok"] is True

    def test_a_tampered_finding_count_fails_and_names_the_field(self):
        signed = self._signed()
        signed["manifest"]["finding_count"] = 999
        verdict = prov.verify_manifest(signed)
        assert verdict["ok"] is False
        assert verdict["field"] == "manifest.finding_count"

    def test_a_tampered_manifest_lists_the_changed_fields(self):
        signed = self._signed()
        signed["manifest"]["finding_count"] = 999
        assert prov.verify_manifest(signed)["changed_fields"] == ["finding_count"]

    def test_a_tampered_database_list_fails_and_names_the_field(self):
        signed = self._signed()
        signed["manifest"]["databases"] = []
        assert prov.verify_manifest(signed)["field"] == "manifest.databases"

    def test_a_tampered_input_hash_fails_and_names_the_field(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_bytes(b"hello world")
        signed = prov.sign_manifest(
            prov.build_manifest(profile_id=1, findings=[], input_files=[raw]))
        signed["manifest"]["input_files"][0]["sha256"] = "0" * 64
        assert prov.verify_manifest(signed)["field"] == "manifest.input_files"

    def test_two_tampered_fields_report_the_manifest_and_list_both(self):
        signed = self._signed()
        signed["manifest"]["finding_count"] = 999
        signed["manifest"]["report_type"] = "forged"
        verdict = prov.verify_manifest(signed)
        assert verdict["field"] == "manifest"
        assert verdict["changed_fields"] == ["finding_count", "report_type"]

    def test_an_added_manifest_field_is_detected(self):
        signed = self._signed()
        signed["manifest"]["injected"] = "value"
        verdict = prov.verify_manifest(signed)
        assert verdict["ok"] is False and "injected" in verdict["changed_fields"]

    def test_a_removed_manifest_field_is_detected(self):
        signed = self._signed()
        del signed["manifest"]["report_type"]
        verdict = prov.verify_manifest(signed)
        assert verdict["ok"] is False and "report_type" in verdict["changed_fields"]

    def test_a_tampered_signature_value_names_the_signature(self):
        signed = self._signed()
        signed["signature"]["value"] = "0" * 64
        verdict = prov.verify_manifest(signed)
        assert verdict["ok"] is False and verdict["field"] == "signature.value"

    def test_a_tampered_signing_timestamp_fails(self):
        signed = self._signed()
        signed["signature"]["signed_at"] = "1999-01-01T00:00:00Z"
        assert prov.verify_manifest(signed)["ok"] is False

    def test_a_downgraded_algorithm_is_rejected_by_name(self):
        signed = self._signed()
        signed["signature"]["algorithm"] = "none"
        verdict = prov.verify_manifest(signed)
        assert verdict["field"] == "signature.algorithm"
        assert verdict["reason"] == "unsupported_algorithm"

    def test_verification_fails_under_a_different_key(self, tmp_path, monkeypatch):
        signed = self._signed()
        monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "elsewhere"))
        verdict = prov.verify_manifest(signed)
        assert verdict["ok"] is False and verdict["field"] == "signature.value"

    def test_an_unsigned_payload_returns_a_structured_failure(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        verdict = prov.verify_manifest({"manifest": manifest})
        assert verdict["ok"] is False
        assert verdict["reason"] == "unsigned"
        assert verdict["field"] == "signature"

    def test_a_bare_manifest_does_not_raise(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        assert prov.verify_manifest(manifest)["ok"] is False

    def test_a_non_mapping_payload_returns_a_structured_failure(self):
        verdict = prov.verify_manifest("not a manifest")
        assert verdict["ok"] is False and verdict["field"] == "signed"

    def test_a_missing_manifest_names_the_manifest_field(self):
        verdict = prov.verify_manifest({"signature": {}})
        assert verdict["field"] == "manifest"

    def test_a_non_mapping_manifest_names_the_manifest_field(self):
        verdict = prov.verify_manifest({"manifest": [], "signature": {}})
        assert verdict["reason"] == "manifest_not_a_mapping"

    def test_a_non_mapping_signature_names_the_signature_field(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        verdict = prov.verify_manifest({"manifest": manifest, "signature": "x"})
        assert verdict["field"] == "signature"

    def test_a_signature_with_no_value_names_the_value_field(self):
        signed = self._signed()
        signed["signature"]["value"] = ""
        assert prov.verify_manifest(signed)["field"] == "signature.value"

    def test_the_verdict_is_never_a_bare_boolean(self):
        assert isinstance(prov.verify_manifest({}), dict)

    def test_every_verdict_carries_a_check_timestamp(self):
        assert prov.verify_manifest(self._signed())["checked_at"].endswith("Z")

    def test_every_verdict_carries_a_reason_and_a_detail(self):
        for payload in ({}, {"manifest": {}}, self._signed()):
            verdict = prov.verify_manifest(payload)
            assert verdict["reason"] and verdict["detail"]


# ---------------------------------------------------------------------------
# Text rendering
# ---------------------------------------------------------------------------

class TestManifestText:

    def _signed(self, **kw):
        return prov.sign_manifest(
            prov.build_manifest(profile_id=1, findings=[finding()],
                                report_type="doctor", **kw))

    def test_the_text_names_every_database(self):
        text = prov.render_manifest_text(self._signed())
        for entry in prov.SOURCES.values():
            assert entry["name"] in text

    def test_the_text_marks_pharmgkb_as_not_used(self):
        assert "[NOT USED]" in prov.render_manifest_text(self._signed())

    def test_the_text_carries_the_signature_value(self):
        signed = self._signed()
        assert signed["signature"]["value"] in prov.render_manifest_text(signed)

    def test_the_text_states_what_the_signature_does_not_prove(self):
        text = prov.render_manifest_text(self._signed()).lower()
        assert "does not prove authorship" in text

    def test_the_text_lists_the_input_file_hash(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_bytes(b"hello world")
        signed = prov.sign_manifest(
            prov.build_manifest(profile_id=1, findings=[], input_files=[raw]))
        assert HELLO_WORLD_SHA256 in prov.render_manifest_text(signed)

    def test_the_text_says_when_there_are_no_input_files(self):
        assert "none recorded" in prov.render_manifest_text(self._signed())

    def test_the_text_lists_the_scan_parameters(self):
        assert "assembly =" in prov.render_manifest_text(self._signed())

    def test_an_unsigned_payload_renders_as_unsigned(self):
        manifest = prov.build_manifest(profile_id=1, findings=[])
        assert "UNSIGNED" in prov.render_manifest_text({"manifest": manifest})

    def test_the_text_tolerates_an_empty_payload(self):
        assert "DNAINSIGHT REPORT MANIFEST" in prov.render_manifest_text({})

    def test_the_text_tolerates_a_non_dict_payload(self):
        assert "DNAINSIGHT REPORT MANIFEST" in prov.render_manifest_text(None)

    def test_no_line_of_the_signature_note_is_unreasonably_long(self):
        text = prov.render_manifest_text(self._signed())
        note_start = text.index("What this signature does")
        for line in text[note_start:].splitlines():
            assert len(line) <= 80


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

class TestSchema:

    def test_init_provenance_is_idempotent(self):
        prov.init_provenance()
        prov.init_provenance()
        conn = db.get_connection()
        try:
            names = {r["name"] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
        finally:
            conn.close()
        assert {"provenance_records", "report_manifests"} <= names

    def test_init_provenance_does_not_disturb_existing_rows(self):
        prov.init_provenance()
        conn = db.get_connection()
        try:
            conn.execute(
                "INSERT INTO provenance_records "
                "(profile_id, rsid, source_id) VALUES (1, 'rs1', 'clinvar')")
            conn.commit()
        finally:
            conn.close()
        prov.init_provenance()
        conn = db.get_connection()
        try:
            count = conn.execute(
                "SELECT COUNT(*) FROM provenance_records").fetchone()[0]
        finally:
            conn.close()
        assert count == 1

    def test_init_provenance_runs_on_a_database_that_already_has_tables(self):
        db.init_db()
        prov.init_provenance()
        prov.init_provenance()
        assert db.list_profiles() == []
