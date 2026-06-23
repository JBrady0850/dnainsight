"""
database.py -- SQLite schema and data access layer for DNAInsight.
"""

import sqlite3
import json
import os
from pathlib import Path
from datetime import datetime


def _resolve_db_path() -> Path:
    """
    Resolve the SQLite database path.

    Priority:
      1. App directory (DNAInsight folder itself) -- standard install on Windows/Mac/Linux.
      2. User home directory ~/.dnainsight/ -- if the app folder is read-only.
      3. System temp directory -- last resort for FUSE/container environments.
    """
    import tempfile, sqlite3 as _sqlite3

    def _test_sqlite(path: Path) -> bool:
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            conn = _sqlite3.connect(str(path))
            # Test WAL mode -- this is what the app uses; FUSE mounts often fail here
            conn.execute("PRAGMA journal_mode=WAL")
            conn.executescript("""
                CREATE TABLE IF NOT EXISTS _write_test (x INTEGER);
                INSERT INTO _write_test VALUES (1);
                DROP TABLE _write_test;
            """)
            conn.close()
            # Clean up test file
            try:
                path.unlink(missing_ok=True)
            except Exception:
                pass
            return True
        except Exception:
            try:
                path.unlink(missing_ok=True)
            except Exception:
                pass
            return False

    candidates = [
        Path(__file__).parent.parent / "dnainsight.db",
        Path.home() / ".dnainsight" / "dnainsight.db",
        Path(tempfile.gettempdir()) / "dnainsight" / "dnainsight.db",
    ]
    for p in candidates:
        if _test_sqlite(p):
            return p
    # Absolute fallback -- should never reach here
    return candidates[-1]

DB_PATH = _resolve_db_path()


def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def init_db():
    """Create all tables if they don't exist."""
    conn = get_connection()
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS profiles (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            name        TEXT NOT NULL,
            dob         TEXT,
            sex         TEXT,
            provider    TEXT,
            created_at  TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS snp_uploads (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            profile_id  INTEGER NOT NULL REFERENCES profiles(id),
            filename    TEXT NOT NULL,
            snp_count   INTEGER,
            uploaded_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS findings (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            profile_id      INTEGER NOT NULL REFERENCES profiles(id),
            upload_id       INTEGER REFERENCES snp_uploads(id),
            rsid            TEXT NOT NULL,
            gene            TEXT,
            chromosome      TEXT,
            position        INTEGER,
            allele1         TEXT,
            allele2         TEXT,
            genotype        TEXT,
            clinical_sig    TEXT,
            conditions      TEXT,
            interpretation  TEXT,
            category        TEXT,
            silo            TEXT,
            sources         TEXT,
            discovered_at   TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_findings_profile ON findings(profile_id);
        CREATE INDEX IF NOT EXISTS idx_findings_rsid    ON findings(rsid);
        CREATE INDEX IF NOT EXISTS idx_findings_silo    ON findings(silo);

        CREATE TABLE IF NOT EXISTS reports (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            profile_id  INTEGER NOT NULL REFERENCES profiles(id),
            report_type TEXT NOT NULL,
            format      TEXT NOT NULL,
            filepath    TEXT NOT NULL,
            created_at  TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS scan_state (
            profile_id  INTEGER PRIMARY KEY REFERENCES profiles(id),
            checkpoint  TEXT NOT NULL DEFAULT '{}'
        );
    """)
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Profile operations
# ---------------------------------------------------------------------------

def create_profile(name: str, dob: str, sex: str, provider: str) -> int:
    conn = get_connection()
    cur = conn.execute(
        "INSERT INTO profiles (name, dob, sex, provider, created_at) VALUES (?, ?, ?, ?, ?)",
        (name, dob, sex, provider, datetime.utcnow().isoformat())
    )
    pid = cur.lastrowid
    conn.commit()
    conn.close()
    return pid


def get_profile(profile_id: int) -> dict | None:
    conn = get_connection()
    row = conn.execute("SELECT * FROM profiles WHERE id=?", (profile_id,)).fetchone()
    conn.close()
    return dict(row) if row else None


def list_profiles() -> list[dict]:
    conn = get_connection()
    rows = conn.execute("SELECT * FROM profiles ORDER BY created_at DESC").fetchall()
    conn.close()
    return [dict(r) for r in rows]


def delete_profile(profile_id: int):
    conn = get_connection()
    conn.execute("DELETE FROM findings WHERE profile_id=?", (profile_id,))
    conn.execute("DELETE FROM snp_uploads WHERE profile_id=?", (profile_id,))
    conn.execute("DELETE FROM scan_state WHERE profile_id=?", (profile_id,))
    conn.execute("DELETE FROM reports WHERE profile_id=?", (profile_id,))
    conn.execute("DELETE FROM profiles WHERE id=?", (profile_id,))
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Upload operations
# ---------------------------------------------------------------------------

def record_upload(profile_id: int, filename: str, snp_count: int) -> int:
    conn = get_connection()
    cur = conn.execute(
        "INSERT INTO snp_uploads (profile_id, filename, snp_count, uploaded_at) VALUES (?, ?, ?, ?)",
        (profile_id, filename, snp_count, datetime.utcnow().isoformat())
    )
    uid = cur.lastrowid
    conn.commit()
    conn.close()
    return uid


# ---------------------------------------------------------------------------
# Findings operations
# ---------------------------------------------------------------------------

def upsert_finding(profile_id: int, upload_id: int, finding: dict):
    """Insert or update a finding (keyed on profile_id + rsid)."""
    conn = get_connection()
    existing = conn.execute(
        "SELECT id FROM findings WHERE profile_id=? AND rsid=?",
        (profile_id, finding["rsid"])
    ).fetchone()

    now = datetime.utcnow().isoformat()
    fields = (
        finding.get("gene", ""),
        finding.get("chromosome", ""),
        finding.get("position", 0),
        finding.get("allele1", ""),
        finding.get("allele2", ""),
        finding.get("genotype", ""),
        finding.get("clinical_sig", ""),
        finding.get("conditions", ""),
        finding.get("interpretation", ""),
        finding.get("category", ""),
        finding.get("silo", "actionable"),
        json.dumps(finding.get("sources", [])),
    )

    if existing:
        conn.execute("""
            UPDATE findings SET
                gene=?, chromosome=?, position=?, allele1=?, allele2=?,
                genotype=?, clinical_sig=?, conditions=?, interpretation=?,
                category=?, silo=?, sources=?
            WHERE profile_id=? AND rsid=?
        """, fields + (profile_id, finding["rsid"]))
    else:
        conn.execute("""
            INSERT INTO findings
                (profile_id, upload_id, rsid, gene, chromosome, position,
                 allele1, allele2, genotype, clinical_sig, conditions,
                 interpretation, category, silo, sources, discovered_at)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (profile_id, upload_id, finding["rsid"]) + fields + (now,))

    conn.commit()
    conn.close()


def get_findings(profile_id: int, silo: str = None) -> list[dict]:
    conn = get_connection()
    if silo:
        rows = conn.execute(
            "SELECT * FROM findings WHERE profile_id=? AND silo=? ORDER BY silo, gene",
            (profile_id, silo)
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT * FROM findings WHERE profile_id=? ORDER BY silo, gene",
            (profile_id,)
        ).fetchall()
    conn.close()
    result = []
    for r in rows:
        d = dict(r)
        try:
            d["sources"] = json.loads(d.get("sources") or "[]")
        except Exception:
            d["sources"] = []
        result.append(d)
    return result


def get_findings_summary(profile_id: int) -> dict:
    conn = get_connection()
    rows = conn.execute(
        "SELECT silo, COUNT(*) as cnt FROM findings WHERE profile_id=? GROUP BY silo",
        (profile_id,)
    ).fetchall()
    conn.close()
    return {r["silo"]: r["cnt"] for r in rows}


# ---------------------------------------------------------------------------
# Scan state
# ---------------------------------------------------------------------------

def get_scan_state(profile_id: int) -> dict:
    conn = get_connection()
    row = conn.execute(
        "SELECT checkpoint FROM scan_state WHERE profile_id=?", (profile_id,)
    ).fetchone()
    conn.close()
    if not row:
        return {"line": 0, "total": 0, "passes": 0}
    try:
        return json.loads(row["checkpoint"])
    except Exception:
        return {"line": 0, "total": 0, "passes": 0}


def save_scan_state(profile_id: int, state: dict):
    conn = get_connection()
    conn.execute(
        "INSERT OR REPLACE INTO scan_state (profile_id, checkpoint) VALUES (?, ?)",
        (profile_id, json.dumps(state))
    )
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

def record_report(profile_id: int, report_type: str, fmt: str, filepath: str) -> int:
    conn = get_connection()
    cur = conn.execute(
        "INSERT INTO reports (profile_id, report_type, format, filepath, created_at) VALUES (?,?,?,?,?)",
        (profile_id, report_type, fmt, filepath, datetime.utcnow().isoformat())
    )
    rid = cur.lastrowid
    conn.commit()
    conn.close()
    return rid


def get_reports(profile_id: int) -> list[dict]:
    conn = get_connection()
    rows = conn.execute(
        "SELECT * FROM reports WHERE profile_id=? ORDER BY created_at DESC",
        (profile_id,)
    ).fetchall()
    conn.close()
    return [dict(r) for r in rows]
