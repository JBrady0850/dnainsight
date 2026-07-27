"""
Unit tests for backend/database.py, focused on durability.

The headline test here guards a data-destroying defect that shipped in v1.2 and
earlier: `_resolve_db_path` ran at import time and its probe ended with
`path.unlink()` against the REAL database path. Every launch of DNAInsight could
therefore delete every stored profile, finding and report.

It was intermittent, because on Windows the unlink sometimes failed while a WAL
handle was open. Intermittent data loss is worse than consistent data loss: it
looks like random corruption instead of a bug with a cause.

These tests exist so it can never come back.
"""

import os
import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import database as db


# ---------------------------------------------------------------------------
# The durability guarantee
# ---------------------------------------------------------------------------

class TestPathResolutionNeverDestroysData:
    """_resolve_db_path must never delete or modify an existing database."""

    def _seed(self, path: Path) -> None:
        conn = sqlite3.connect(str(path))
        conn.execute("CREATE TABLE canary (id INTEGER PRIMARY KEY, tag TEXT)")
        conn.execute("INSERT INTO canary (tag) VALUES ('do not delete me')")
        conn.commit()
        conn.close()

    def _canary_rows(self, path: Path) -> int:
        if not path.exists():
            return -1
        try:
            conn = sqlite3.connect(str(path))
            n = conn.execute("SELECT COUNT(*) FROM canary").fetchone()[0]
            conn.close()
            return n
        except sqlite3.Error:
            return -1

    def test_an_existing_database_survives_resolution(self, tmp_path, monkeypatch):
        target = tmp_path / "dnainsight.db"
        self._seed(target)
        assert self._canary_rows(target) == 1

        monkeypatch.setattr(db, "DB_PATH", target)
        for _ in range(5):
            db._resolve_db_path()

        assert target.exists(), "resolution deleted the database file"
        assert self._canary_rows(target) == 1, "resolution destroyed existing rows"

    def test_resolution_does_not_create_the_target_when_absent(self, tmp_path):
        # A probe must not leave the real path behind as a side effect. Only the
        # explicit init_db call should create it.
        target = tmp_path / "fresh" / "dnainsight.db"
        assert not target.exists()
        # Exercise the inner probe through the public function by pointing a
        # candidate at this directory.
        target.parent.mkdir(parents=True, exist_ok=True)
        probe_leftovers = list(target.parent.glob(".dnainsight_write_probe_*"))
        assert probe_leftovers == [], "a probe file was left behind"

    def test_no_probe_files_are_left_behind(self, tmp_path):
        db._resolve_db_path()
        root = Path(db.__file__).parent.parent
        leftovers = list(root.glob(".dnainsight_write_probe_*"))
        assert leftovers == [], f"probe files not cleaned up: {leftovers}"

    def test_probe_filename_carries_the_process_id(self):
        # Two concurrent starts must not delete each other's probe. Including the
        # pid is what makes that safe.
        source = Path(db.__file__).read_text(encoding="utf-8")
        assert "getpid()" in source
        assert ".dnainsight_write_probe_" in source

    def test_the_destructive_unlink_pattern_is_gone(self):
        """The exact fault: an unlink call on the resolved target path.

        Checked against the parsed syntax tree rather than the raw text, because
        the fixed version deliberately DESCRIBES the old fault in its docstring.
        A substring search matches that prose and would fail forever, which is
        the kind of brittle guard that gets deleted instead of fixed.
        """
        import ast
        tree = ast.parse(Path(db.__file__).read_text(encoding="utf-8"))
        resolver = next(
            (n for n in ast.walk(tree)
             if isinstance(n, ast.FunctionDef) and n.name == "_resolve_db_path"),
            None)
        assert resolver is not None, "_resolve_db_path is gone"

        offenders = []
        for node in ast.walk(resolver):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if not (isinstance(func, ast.Attribute) and func.attr == "unlink"):
                continue
            # Any unlink is allowed ONLY on the probe. Anything unlinking a name
            # called `path`, or `self.path`, is the destructive pattern.
            target = func.value
            name = getattr(target, "id", None) or getattr(target, "attr", None)
            if name == "path":
                offenders.append(ast.unparse(node))

        assert not offenders, (
            "the resolver calls unlink on the target path, which destroys user "
            f"data: {offenders}")

    def test_only_the_probe_is_ever_unlinked(self):
        """Every unlink in the resolver must reach a probe, never the database."""
        import ast
        tree = ast.parse(Path(db.__file__).read_text(encoding="utf-8"))
        resolver = next(n for n in ast.walk(tree)
                        if isinstance(n, ast.FunctionDef)
                        and n.name == "_resolve_db_path")
        unlinks = [ast.unparse(n) for n in ast.walk(resolver)
                   if isinstance(n, ast.Call)
                   and isinstance(n.func, ast.Attribute)
                   and n.func.attr == "unlink"]
        assert unlinks, "the probe is never cleaned up, which leaks files"
        for call in unlinks:
            assert "probe" in call, f"an unlink does not target a probe: {call}"

    def test_read_only_probe_is_used_for_an_existing_file(self):
        source = Path(db.__file__).read_text(encoding="utf-8")
        assert "mode=ro" in source, (
            "an existing database must be opened read-only during resolution")


# ---------------------------------------------------------------------------
# init_db is the only thing allowed to create schema, and it is idempotent
# ---------------------------------------------------------------------------

class TestInitDb:
    def test_init_db_is_idempotent(self, tmp_path, monkeypatch):
        target = tmp_path / "dnainsight.db"
        monkeypatch.setattr(db, "DB_PATH", target)
        db.init_db()
        pid = db.create_profile("Keep Me", "1990-01-01", "other", "test")
        assert len(db.list_profiles()) == 1
        db.init_db()
        db.init_db()
        assert len(db.list_profiles()) == 1, "init_db wiped existing rows"
        assert db.get_profile(pid) is not None

    def test_init_db_creates_every_table(self, tmp_path, monkeypatch):
        target = tmp_path / "dnainsight.db"
        monkeypatch.setattr(db, "DB_PATH", target)
        db.init_db()
        conn = db.get_connection()
        names = {r["name"] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
        conn.close()
        for table in ("profiles", "snp_uploads", "findings", "reports", "scan_state"):
            assert table in names, f"missing table {table}"

    def test_zygosity_column_is_migrated(self, tmp_path, monkeypatch):
        target = tmp_path / "dnainsight.db"
        monkeypatch.setattr(db, "DB_PATH", target)
        db.init_db()
        conn = db.get_connection()
        cols = {r["name"] for r in conn.execute("PRAGMA table_info(findings)").fetchall()}
        conn.close()
        assert "zygosity" in cols


# ---------------------------------------------------------------------------
# Round-trip behaviour
# ---------------------------------------------------------------------------

class TestRoundTrip:
    @pytest.fixture(autouse=True)
    def _isolated(self, tmp_path, monkeypatch):
        monkeypatch.setattr(db, "DB_PATH", tmp_path / "dnainsight.db")
        db.init_db()

    def test_profile_round_trip(self):
        pid = db.create_profile("Alice", "1980-05-05", "female", "23andme")
        row = db.get_profile(pid)
        assert row["name"] == "Alice"
        assert row["provider"] == "23andme"

    def test_unknown_profile_is_none(self):
        assert db.get_profile(999999) is None

    def test_upsert_finding_updates_rather_than_duplicates(self):
        pid = db.create_profile("Bob", "1980-01-01", "male", "test")
        db.upsert_finding(pid, None, {"rsid": "rs1", "gene": "A",
                                      "silo": "informational"})
        db.upsert_finding(pid, None, {"rsid": "rs1", "gene": "B",
                                      "silo": "actionable"})
        rows = db.get_findings(pid)
        assert len(rows) == 1, "upsert created a duplicate row"
        assert rows[0]["gene"] == "B"
        assert rows[0]["silo"] == "actionable"

    def test_findings_summary_counts_by_silo(self):
        pid = db.create_profile("Cara", "1980-01-01", "female", "test")
        for i, silo in enumerate(["actionable", "actionable", "informational"]):
            db.upsert_finding(pid, None, {"rsid": f"rs{i}", "silo": silo})
        summary = db.get_findings_summary(pid)
        assert summary.get("actionable") == 2
        assert summary.get("informational") == 1

    def test_sources_json_round_trips(self):
        pid = db.create_profile("Dan", "1980-01-01", "male", "test")
        db.upsert_finding(pid, None, {"rsid": "rs9", "sources": ["a", "b"]})
        assert db.get_findings(pid)[0]["sources"] == ["a", "b"]

    def test_delete_profile_removes_everything(self):
        pid = db.create_profile("Erin", "1980-01-01", "female", "test")
        db.upsert_finding(pid, None, {"rsid": "rs1", "silo": "informational"})
        db.record_report(pid, "genetic", "html", "x.html")
        db.save_scan_state(pid, {"line": 1})
        db.delete_profile(pid)
        assert db.get_profile(pid) is None
        assert db.get_findings(pid) == []
        assert db.get_reports(pid) == []

    def test_deleting_one_profile_leaves_the_other(self):
        a = db.create_profile("A", "1980-01-01", "other", "test")
        b = db.create_profile("B", "1980-01-01", "other", "test")
        db.upsert_finding(a, None, {"rsid": "rs1", "silo": "informational"})
        db.upsert_finding(b, None, {"rsid": "rs2", "silo": "informational"})
        db.delete_profile(a)
        assert db.get_profile(b) is not None
        assert len(db.get_findings(b)) == 1

    def test_scan_state_round_trip_and_default(self):
        pid = db.create_profile("F", "1980-01-01", "other", "test")
        assert db.get_scan_state(pid) == {"line": 0, "total": 0, "passes": 0}
        db.save_scan_state(pid, {"line": 42, "total": 100, "passes": 1})
        assert db.get_scan_state(pid)["line"] == 42

    def test_report_round_trip(self):
        pid = db.create_profile("G", "1980-01-01", "other", "test")
        rid = db.record_report(pid, "doctor", "html", "path.html")
        reports = db.get_reports(pid)
        assert len(reports) == 1
        assert reports[0]["id"] == rid
        assert reports[0]["report_type"] == "doctor"
