"""
ledger.py -- the per-user reclassification ledger.

WHY THIS MODULE EXISTS
----------------------
ClinVar reclassifies monthly. ``data/DATA_SOURCES.md`` already warns that
"clinical assertions change between monthly releases, so a stale copy can be
actively misleading", and the Tier 2 ``meta`` table already records
``build_date`` and ``clinvar_version``. So the project already knows how old its
evidence is. What it has never had is the PER-USER DIFF: this person's variant
rs28897696 was Uncertain Significance in the July release and is Pathogenic in
the August one, and nobody told them.

Sequencing.com Premium advertises continuous reanalysis and shows no diff.
No consumer product gives a user a reclassification audit trail. A VUS that
became Pathogenic is the single highest-value event in personal genomics and
nobody surfaces it as an event. That is what this module does: it records a
fingerprint of the clinically comparable state of every finding at scan time,
and it can tell you exactly what changed between any two scans, in which
direction, and between which database releases.

WHY A FINGERPRINT RATHER THAN THE WHOLE FINDING
-----------------------------------------------
Interpretation prose churns constantly. A submitter rewords a condition name, a
summary gains a comma, our own template changes, and a naive whole-record diff
screams about a hundred findings that are clinically identical. Users stop
reading a report that cries wolf, and then the one real reclassification is
buried. ``fingerprint`` therefore captures only comparable clinical state:
significance code, review stars, CPIC level, magnitude to one decimal place,
repute, carrier status, entity type, and whether the entity was evaluable at
all. Prose is deliberately excluded.

THE ADDENDUM IS DATED AND ADDITIVE
----------------------------------
An addendum NEVER rewrites, replaces or supersedes the original report. This is
not a style preference. A clinician who acted on the January report has to be
able to see exactly what the January report said, because their decision is only
defensible against the evidence that existed at the time. Silently mutating the
old document would destroy that defence and would also destroy the user's
ability to notice that anything moved.

Enforced structurally: there is no function in this module that mutates a prior
snapshot or a prior report. Snapshots are insert-only. ``addendum()`` reads two
snapshots and returns a new dated payload whose ``supersedes`` field is always
None. The only destructive operation is ``prune_snapshots``, which drops the
OLDEST rows on explicit request and can never delete the newest one.

DURABILITY
----------
``backend/database.py`` carries a long comment about a defect that deleted every
user's database on launch. The lesson generalises: a helper that "just tidies
up" must never be able to reach live data. Everything here follows it.
``init_ledger`` is CREATE TABLE IF NOT EXISTS plus the ``_add_column_if_missing``
migration pattern only. No DROP, no rewrite, no table rebuild, ever, on any
path. ``prune_snapshots`` clamps ``keep`` to at least 1 so a caller passing 0
cannot wipe the ledger.

OFFLINE
-------
Nothing here touches the network. Standard library plus ``backend.database``.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from typing import Any

from .database import get_connection, _add_column_if_missing

__all__ = [
    "CHANGE_KINDS", "MATERIAL_KINDS",
    "init_ledger", "fingerprint", "snapshot", "list_snapshots",
    "latest_snapshot", "get_snapshot", "diff_snapshots", "changes_for",
    "addendum", "render_addendum_html", "prune_snapshots",
]


# ---------------------------------------------------------------------------
# Change vocabulary
#
# ``severity`` orders the addendum so it leads with what matters. It is a rank,
# not a probability and not a clinical grade. ``direction`` on the KIND is the
# static default; every emitted change record carries its own ``direction``
# computed from the actual old and new values, because "changed" on its own is
# useless to a reader who has to decide whether to worry.
#
# vus_resolved_pathogenic is the highest severity in the system by design. It is
# the event this module exists to surface.
#
# Two kinds beyond the required set:
#   sig_reclassified       for a move that cannot be placed on the benign to
#                          pathogenic axis at all (pathogenic to conflicting,
#                          uncertain to risk factor). Calling that an upgrade or
#                          a downgrade would be a lie, and dropping it would
#                          hide a real change. It still names both positions.
#   carrier_status_changed for when a reference update flips the recorded risk
#                          allele and a person who was told they carry nothing
#                          is now told they carry something. That is exactly the
#                          class of event this ledger exists for.
# ---------------------------------------------------------------------------

CHANGE_KINDS: dict[str, dict[str, Any]] = {
    "vus_resolved_pathogenic": {
        "label": "Uncertain significance resolved to pathogenic",
        "severity": 100, "direction": "up",
    },
    "sig_upgraded": {
        "label": "Clinical significance moved toward pathogenic",
        "severity": 90, "direction": "up",
    },
    "vus_resolved_benign": {
        "label": "Uncertain significance resolved to benign",
        "severity": 80, "direction": "down",
    },
    "sig_downgraded": {
        "label": "Clinical significance moved toward benign",
        "severity": 70, "direction": "down",
    },
    "sig_reclassified": {
        "label": "Clinical significance moved off the benign to pathogenic axis",
        "severity": 65, "direction": "lateral",
    },
    "carrier_status_changed": {
        "label": "Carrier status changed",
        "severity": 62, "direction": "lateral",
    },
    "cpic_level_changed": {
        "label": "CPIC actionability level changed",
        "severity": 60, "direction": "lateral",
    },
    "finding_added": {
        "label": "New finding appeared",
        "severity": 55, "direction": "up",
    },
    "finding_removed": {
        "label": "Finding no longer reported",
        "severity": 50, "direction": "down",
    },
    "stars_gained": {
        "label": "ClinVar review confidence increased",
        "severity": 45, "direction": "up",
    },
    "stars_lost": {
        "label": "ClinVar review confidence decreased",
        "severity": 40, "direction": "down",
    },
    "repute_changed": {
        "label": "Direction of effect changed",
        "severity": 35, "direction": "lateral",
    },
    "newly_evaluable": {
        "label": "Now evaluable on this data",
        "severity": 30, "direction": "up",
    },
    "no_longer_evaluable": {
        "label": "No longer evaluable on this data",
        "severity": 25, "direction": "down",
    },
    "magnitude_changed": {
        "label": "DNAInsight magnitude changed",
        "severity": 10, "direction": "lateral",
    },
}

# The kinds that warrant telling the user. magnitude_changed is the only
# exclusion, and the reason is that magnitude is OUR derived number: it moves
# when the scoring code changes as well as when the evidence changes, so a
# magnitude move on its own is not proof that anything upstream happened. Every
# other kind is driven by a database value, so every other kind is material.
MATERIAL_KINDS: frozenset[str] = frozenset(
    k for k in CHANGE_KINDS if k != "magnitude_changed"
)

# Position on the benign to pathogenic axis. Codes that are not on that axis at
# all (drug response 6, histocompatibility 7, risk factor / conflicting / other
# 255, and absent) map to 0, meaning "off axis", never to a middle rank. A
# conflicting record in particular must never be resolved to one of the
# positions it conflicts over, which is the same invariant scoring.py enforces.
_AXIS_RANK: dict[Any, int] = {
    2: 1,    # benign
    3: 2,    # likely benign
    1: 3,    # uncertain significance
    4: 4,    # likely pathogenic
    5: 5,    # pathogenic
}

_SIG_DISPLAY: dict[Any, str] = {
    None: "no ClinVar classification",
    1: "uncertain significance",
    2: "benign",
    3: "likely benign",
    4: "likely pathogenic",
    5: "pathogenic",
    6: "drug response",
    7: "histocompatibility",
    255: "other, conflicting or unranked",
}

# CPIC actionability strength. A is strongest. Absent is -1 rather than 0 so
# that gaining any level at all reads as an increase.
_CPIC_RANK: dict[str, int] = {
    "": -1, "Retired": 0, "D": 1, "C/D": 2, "C": 3,
    "B/C": 4, "B": 5, "A/B": 6, "A": 7,
}

# Repute as a concern ordering, so a repute move always has a direction.
_REPUTE_RANK: dict[str, int] = {"Good": 0, "": 1, "Bad": 2}

# Fallback display names for database sources, used only when
# backend.provenance cannot be imported. Duplication is deliberate and small:
# the ledger must remain usable without the provenance module present.
_SOURCE_LABELS: dict[str, str] = {
    "clinvar": "ClinVar",
    "cpic": "CPIC",
    "gnomad": "gnomAD",
    "onekg_ensembl": "1000 Genomes via Ensembl",
    "gwas_catalog": "GWAS Catalog",
    "pgs_catalog": "PGS Catalog",
    "myvariant": "MyVariant.info",
    "snpedia": "SNPedia",
    "pharmgkb": "PharmGKB",
    "dnainsight": "DNAInsight",
    "snp_reference": "bundled SNP reference",
    "genosets": "genoset corpus",
    "prs_models": "PRS models",
    "frequencies": "frequency table",
    "reference_db": "Tier 2 reference database",
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _utc_now() -> str:
    """UTC timestamp with microseconds, so two snapshots taken in the same
    second still order deterministically."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _text(value: Any) -> str:
    return str(value or "").strip()


def _canonical(payload: Any) -> str:
    """Stable JSON for hashing. Sorted keys and no incidental whitespace, so
    the same state always produces the same digest across processes."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, default=str)


def _digest(payload: Any) -> str:
    return hashlib.sha256(_canonical(payload).encode("utf-8")).hexdigest()


def _as_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _source_label(source_id: str) -> str:
    """Human name for a data source id, for the version statement."""
    try:
        from . import provenance as _prov
        entry = _prov.SOURCES.get(source_id)
        if isinstance(entry, dict) and entry.get("name"):
            return str(entry["name"])
    except Exception:
        # A missing or broken provenance module must not stop a diff from being
        # produced. The version statement degrades to the raw id, which is still
        # more informative than saying nothing.
        pass
    return _SOURCE_LABELS.get(source_id, source_id)


def _live_db_versions() -> dict:
    """Live database versions, or an empty dict when provenance is unavailable.

    Soft import on purpose. The ledger is the more fundamental of the two Wave 1
    modules and must keep working if provenance is absent or fails to load its
    bundled JSON.
    """
    try:
        from . import provenance as _prov
        versions = _prov.database_versions()
        return versions if isinstance(versions, dict) else {}
    except Exception:
        return {}


def _table_exists(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,)
    ).fetchone()
    return row is not None


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_SCHEMA = """
    CREATE TABLE IF NOT EXISTS reclass_snapshots (
        id            INTEGER PRIMARY KEY AUTOINCREMENT,
        profile_id    INTEGER NOT NULL,
        created_at    TEXT NOT NULL,
        label         TEXT NOT NULL DEFAULT '',
        db_versions   TEXT NOT NULL DEFAULT '{}',
        finding_count INTEGER NOT NULL DEFAULT 0,
        digest        TEXT NOT NULL DEFAULT ''
    );

    CREATE INDEX IF NOT EXISTS idx_reclass_snapshots_profile
        ON reclass_snapshots(profile_id, created_at);

    CREATE TABLE IF NOT EXISTS reclass_entries (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        snapshot_id INTEGER NOT NULL
                    REFERENCES reclass_snapshots(id) ON DELETE CASCADE,
        key         TEXT NOT NULL,
        rsid        TEXT NOT NULL DEFAULT '',
        gene        TEXT NOT NULL DEFAULT '',
        entity_type TEXT NOT NULL DEFAULT 'snp',
        state       TEXT NOT NULL DEFAULT '{}',
        digest      TEXT NOT NULL DEFAULT ''
    );

    CREATE INDEX IF NOT EXISTS idx_reclass_entries_snapshot
        ON reclass_entries(snapshot_id);

    CREATE UNIQUE INDEX IF NOT EXISTS idx_reclass_entries_key
        ON reclass_entries(snapshot_id, key);
"""


def _create_schema(conn: sqlite3.Connection) -> None:
    """Create the ledger tables and run forward migrations.

    CREATE TABLE IF NOT EXISTS only. Never DROP, never rebuild. If a table
    already exists with an older column set, ``_add_column_if_missing`` adds
    what is missing and leaves every existing row untouched, which is the
    pattern database.init_db uses and the only pattern that is safe against a
    database holding a user's real history.
    """
    conn.executescript(_SCHEMA)

    # Forward migrations for a ledger written by an earlier build. No-ops on a
    # fresh database. SQLite requires a DEFAULT on a NOT NULL added column, so
    # every one of these carries one.
    _add_column_if_missing(conn, "reclass_snapshots", "label", "TEXT NOT NULL DEFAULT ''")
    _add_column_if_missing(conn, "reclass_snapshots", "db_versions", "TEXT NOT NULL DEFAULT '{}'")
    _add_column_if_missing(conn, "reclass_snapshots", "finding_count", "INTEGER NOT NULL DEFAULT 0")
    _add_column_if_missing(conn, "reclass_snapshots", "digest", "TEXT NOT NULL DEFAULT ''")
    _add_column_if_missing(conn, "reclass_entries", "rsid", "TEXT NOT NULL DEFAULT ''")
    _add_column_if_missing(conn, "reclass_entries", "gene", "TEXT NOT NULL DEFAULT ''")
    _add_column_if_missing(conn, "reclass_entries", "entity_type", "TEXT NOT NULL DEFAULT 'snp'")
    _add_column_if_missing(conn, "reclass_entries", "digest", "TEXT NOT NULL DEFAULT ''")


def init_ledger() -> None:
    """Create the ledger tables if they do not exist.

    Idempotent and safe to call on an existing database, including one that
    already holds snapshots. Call it at startup next to ``database.init_db``.
    ``snapshot()`` also calls it, so a caller who forgets cannot lose data.
    """
    conn = get_connection()
    try:
        _create_schema(conn)
        conn.commit()
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Fingerprinting
# ---------------------------------------------------------------------------

def _evaluable(finding: dict) -> bool:
    """Was this entity actually evaluable against the user's data?

    A genoset needs every one of its required positions present, so anything
    short of full coverage means the rule could not be tested. A polygenic score
    carries its own ``reliable`` flag from prs.py, which is coverage against the
    documented threshold. Everything else is evaluable unless the probe failed.

    This matters for the ledger because "the array finally covers this genoset"
    and "ClinVar changed its mind" are different events and must not be reported
    as the same thing.
    """
    entity = _text(finding.get("entity_type")).lower() or "snp"

    if entity == "genoset":
        coverage = _as_float(finding.get("coverage"))
        return coverage is not None and coverage >= 1.0

    if entity == "prs":
        if isinstance(finding.get("reliable"), bool):
            return bool(finding["reliable"])
        block = finding.get("prs")
        if isinstance(block, dict) and isinstance(block.get("reliable"), bool):
            return bool(block["reliable"])
        coverage = _as_float(finding.get("coverage"))
        if coverage is None and isinstance(block, dict):
            coverage = _as_float(block.get("coverage"))
        return coverage is not None and coverage >= 0.90

    return _text(finding.get("zygosity")).lower() != "no_call"


def _sig_code_of(finding: dict) -> int | None:
    """The ClinVar significance code, computed if the finding lacks one."""
    code = finding.get("clinvar_sig_code")
    if isinstance(code, int) and not isinstance(code, bool):
        return code
    if code is not None:
        coerced = _as_int(code)
        if coerced is not None:
            return coerced
    raw = _text(finding.get("clinical_sig"))
    if not raw:
        return None
    try:
        from .scoring import clinvar_sig_code
        return clinvar_sig_code(raw)
    except Exception:
        return None


def _stars_of(finding: dict) -> int:
    stars = finding.get("review_stars")
    if isinstance(stars, int) and not isinstance(stars, bool):
        return stars
    coerced = _as_int(stars)
    if coerced is not None:
        return coerced
    try:
        from .scoring import review_stars
        return review_stars(finding.get("review_status"))
    except Exception:
        return 0


def _cpic_of(finding: dict) -> str:
    try:
        from .scoring import normalize_cpic_level
        return normalize_cpic_level(finding.get("cpic_level"))
    except Exception:
        return _text(finding.get("cpic_level"))


def fingerprint(finding: dict) -> dict:
    """Reduce one finding to the state worth diffing.

    Captures rsid, entity type, ClinVar significance code, review stars, CPIC
    level, magnitude rounded to one decimal place, repute, carrier status
    (copies of the variant allele, None when the risk allele is unknown) and an
    evaluability flag.

    Deliberately EXCLUDES interpretation, summary, conditions, gene and every
    other prose field. Those churn between releases without meaning anything: a
    reworded condition name is not a reclassification, and treating it as one
    trains the user to ignore the report. Gene is carried alongside the
    fingerprint for display but never enters the digest.

    Magnitude is rounded to one decimal place for the same reason. The scorer
    emits two decimals; a move from 4.53 to 4.52 caused by a publication count
    ticking over is noise.

    Returns a dict including a ``digest``, which is the sha256 of everything
    else. Two fingerprints are clinically identical when their digests match.
    """
    if not isinstance(finding, dict):
        finding = {}

    entity = _text(finding.get("entity_type")).lower() or "snp"
    rsid = _text(finding.get("rsid"))
    magnitude = _as_float(finding.get("magnitude"))
    copies = finding.get("variant_copies")
    copies = copies if isinstance(copies, int) and not isinstance(copies, bool) else None

    repute = _text(finding.get("repute"))
    if repute not in ("Good", "Bad"):
        repute = ""

    state = {
        "key": f"{entity}:{rsid}",
        "rsid": rsid,
        "entity_type": entity,
        "clinvar_sig_code": _sig_code_of(finding),
        "review_stars": _stars_of(finding),
        "cpic_level": _cpic_of(finding),
        "magnitude": round(magnitude, 1) if magnitude is not None else None,
        "repute": repute,
        "carrier": copies,
        "evaluable": _evaluable(finding),
    }
    state["digest"] = _digest({k: v for k, v in state.items()})
    return state


# ---------------------------------------------------------------------------
# Snapshots
# ---------------------------------------------------------------------------

def snapshot(profile_id: int, findings: list[dict], *,
             db_versions: dict | None = None, label: str = "") -> int:
    """Record the current clinical state of every finding for one profile.

    Insert-only. A snapshot is never updated and never rewritten, because a
    snapshot that can change is not evidence of anything.

    ``db_versions`` is the set of database versions in force AT SCAN TIME. It is
    stored with the snapshot rather than looked up later, so a diff can state
    "between ClinVar 2026-07-27 and ClinVar 2026-08-24" instead of the useless
    "between two scans". When omitted it is filled from
    ``provenance.database_versions()``.

    Returns the new snapshot id.
    """
    versions = db_versions if isinstance(db_versions, dict) else _live_db_versions()
    versions = {str(k): ("" if v is None else str(v)) for k, v in versions.items()}

    prints: dict[str, dict] = {}
    genes: dict[str, str] = {}
    for finding in (findings or []):
        if not isinstance(finding, dict):
            continue
        fp = fingerprint(finding)
        if not fp["rsid"]:
            # An entity with no identity cannot be tracked across scans. Skipping
            # it is honest; inventing a key would silently pair unrelated rows in
            # a later diff and report reclassifications that never happened.
            continue
        # Last write wins on a duplicate key. The pipeline already dedupes by
        # rsid, so this only fires on malformed input.
        prints[fp["key"]] = fp
        genes[fp["key"]] = _text(finding.get("gene"))

    created_at = _utc_now()
    snapshot_digest = _digest(sorted(fp["digest"] for fp in prints.values()))

    conn = get_connection()
    try:
        _create_schema(conn)
        cur = conn.execute(
            "INSERT INTO reclass_snapshots "
            "(profile_id, created_at, label, db_versions, finding_count, digest) "
            "VALUES (?,?,?,?,?,?)",
            (int(profile_id), created_at, _text(label), _canonical(versions),
             len(prints), snapshot_digest),
        )
        snapshot_id = int(cur.lastrowid)
        conn.executemany(
            "INSERT INTO reclass_entries "
            "(snapshot_id, key, rsid, gene, entity_type, state, digest) "
            "VALUES (?,?,?,?,?,?,?)",
            [
                (snapshot_id, fp["key"], fp["rsid"], genes.get(fp["key"], ""),
                 fp["entity_type"], _canonical(fp), fp["digest"])
                for fp in prints.values()
            ],
        )
        conn.commit()
    finally:
        conn.close()
    return snapshot_id


def _header_from_row(row: sqlite3.Row) -> dict:
    try:
        versions = json.loads(row["db_versions"] or "{}")
    except (TypeError, ValueError):
        versions = {}
    return {
        "id": int(row["id"]),
        "profile_id": int(row["profile_id"]),
        "created_at": row["created_at"],
        "label": row["label"] or "",
        "db_versions": versions if isinstance(versions, dict) else {},
        "finding_count": int(row["finding_count"] or 0),
        "digest": row["digest"] or "",
    }


def list_snapshots(profile_id: int) -> list[dict]:
    """Every snapshot header for a profile, NEWEST FIRST.

    Newest first matches ``database.list_profiles`` and matches what a user
    wants to see. Ties on ``created_at`` break on id descending, so two
    snapshots taken in the same microsecond still order deterministically
    instead of coming back in whatever order SQLite happens to scan.

    Headers only. Entries are large and are fetched by ``get_snapshot``.
    """
    conn = get_connection()
    try:
        if not _table_exists(conn, "reclass_snapshots"):
            return []
        rows = conn.execute(
            "SELECT * FROM reclass_snapshots WHERE profile_id=? "
            "ORDER BY created_at DESC, id DESC",
            (int(profile_id),),
        ).fetchall()
    finally:
        conn.close()
    return [_header_from_row(r) for r in rows]


def get_snapshot(snapshot_id: int) -> dict | None:
    """One snapshot, header plus its fingerprint entries. None when absent."""
    conn = get_connection()
    try:
        if not _table_exists(conn, "reclass_snapshots"):
            return None
        row = conn.execute(
            "SELECT * FROM reclass_snapshots WHERE id=?", (int(snapshot_id),)
        ).fetchone()
        if row is None:
            return None
        header = _header_from_row(row)
        entries = conn.execute(
            "SELECT * FROM reclass_entries WHERE snapshot_id=? ORDER BY key",
            (int(snapshot_id),),
        ).fetchall()
    finally:
        conn.close()

    prints: dict[str, dict] = {}
    genes: dict[str, str] = {}
    for entry in entries:
        try:
            state = json.loads(entry["state"] or "{}")
        except (TypeError, ValueError):
            state = {}
        if not isinstance(state, dict):
            state = {}
        key = entry["key"]
        prints[key] = state
        genes[key] = entry["gene"] or ""

    header["entries"] = prints
    header["genes"] = genes
    header["entry_count"] = len(prints)
    return header


def latest_snapshot(profile_id: int) -> dict | None:
    """The most recent snapshot for a profile, with entries. None when there
    has never been one."""
    headers = list_snapshots(profile_id)
    if not headers:
        return None
    return get_snapshot(headers[0]["id"])


def prune_snapshots(profile_id: int, keep: int = 12) -> int:
    """Delete all but the newest ``keep`` snapshots for a profile.

    Returns the number of snapshots deleted.

    ``keep`` is clamped to a minimum of 1. A caller passing 0 or a negative
    number almost certainly has a bug, and honouring it would erase the user's
    entire reclassification history, which is exactly the class of mistake
    database.py's path-probe defect was. Refusing costs nothing; obeying is
    unrecoverable.

    Entries are deleted explicitly rather than relying on ON DELETE CASCADE. A
    ledger created by an earlier build may not carry the cascade clause, and
    orphaned entry rows would grow forever without anything noticing.
    """
    keep = max(1, int(keep))
    conn = get_connection()
    try:
        if not _table_exists(conn, "reclass_snapshots"):
            return 0
        rows = conn.execute(
            "SELECT id FROM reclass_snapshots WHERE profile_id=? "
            "ORDER BY created_at DESC, id DESC",
            (int(profile_id),),
        ).fetchall()
        doomed = [int(r["id"]) for r in rows[keep:]]
        if not doomed:
            return 0
        marks = ",".join("?" for _ in doomed)
        conn.execute(f"DELETE FROM reclass_entries WHERE snapshot_id IN ({marks})", doomed)
        conn.execute(f"DELETE FROM reclass_snapshots WHERE id IN ({marks})", doomed)
        conn.commit()
    finally:
        conn.close()
    return len(doomed)


# ---------------------------------------------------------------------------
# Diffing
# ---------------------------------------------------------------------------

def _change(kind: str, fp: dict, gene: str, field: str, old: Any, new: Any,
            direction: str, detail: str,
            old_display: str = "", new_display: str = "") -> dict:
    meta = CHANGE_KINDS[kind]
    return {
        "kind": kind,
        "label": meta["label"],
        "severity": meta["severity"],
        "direction": direction,
        "material": kind in MATERIAL_KINDS,
        "key": fp.get("key", ""),
        "rsid": fp.get("rsid", ""),
        "gene": gene,
        "entity_type": fp.get("entity_type", "snp"),
        "field": field,
        "old": old,
        "new": new,
        "old_display": old_display if old_display else ("" if old is None else str(old)),
        "new_display": new_display if new_display else ("" if new is None else str(new)),
        "detail": detail,
    }


def _sig_display(code: Any) -> str:
    if code in _SIG_DISPLAY:
        return _SIG_DISPLAY[code]
    return f"ClinVar code {code}"


def _classify_significance(old_fp: dict, new_fp: dict, gene: str) -> list[dict]:
    """Classify a ClinVar significance move, always with a direction."""
    old_code = old_fp.get("clinvar_sig_code")
    new_code = new_fp.get("clinvar_sig_code")
    if old_code == new_code:
        return []

    old_rank = _AXIS_RANK.get(old_code, 0)
    new_rank = _AXIS_RANK.get(new_code, 0)
    old_show = _sig_display(old_code)
    new_show = _sig_display(new_code)

    def emit(kind: str, direction: str) -> list[dict]:
        return [_change(
            kind, new_fp, gene, "clinvar_sig_code", old_code, new_code,
            direction,
            f"ClinVar moved from {old_show} to {new_show}.",
            old_display=old_show, new_display=new_show,
        )]

    # The headline event. An uncertain call that has been resolved either way is
    # reported as a resolution, not as a generic up or down move, because
    # "resolved" is the word that makes a user act.
    if old_rank == 3:
        if new_rank >= 4:
            return emit("vus_resolved_pathogenic", "up")
        if 1 <= new_rank <= 2:
            return emit("vus_resolved_benign", "down")

    if old_rank and new_rank:
        return emit("sig_upgraded" if new_rank > old_rank else "sig_downgraded",
                    "up" if new_rank > old_rank else "down")

    if new_rank >= 4:
        # Arrived on the axis at the pathogenic end from somewhere off it.
        return emit("sig_upgraded", "up")
    if 1 <= new_rank <= 2:
        return emit("sig_downgraded", "down")

    # Either the record left the axis (pathogenic to conflicting), or it moved
    # between two off-axis states, or it landed on uncertain from off axis.
    # Calling any of those an upgrade or a downgrade would be a lie. The record
    # still names both positions, so the reader can judge.
    return emit("sig_reclassified", "lateral")


def _classify_pair(old_fp: dict, new_fp: dict, gene: str) -> list[dict]:
    """Every change between two fingerprints of the same entity."""
    changes: list[dict] = []
    changes.extend(_classify_significance(old_fp, new_fp, gene))

    old_stars = _as_int(old_fp.get("review_stars")) or 0
    new_stars = _as_int(new_fp.get("review_stars")) or 0
    if new_stars != old_stars:
        gained = new_stars > old_stars
        changes.append(_change(
            "stars_gained" if gained else "stars_lost", new_fp, gene,
            "review_stars", old_stars, new_stars, "up" if gained else "down",
            f"ClinVar review status moved from {old_stars} to {new_stars} stars.",
            old_display=f"{old_stars} stars", new_display=f"{new_stars} stars",
        ))

    old_cpic = _text(old_fp.get("cpic_level"))
    new_cpic = _text(new_fp.get("cpic_level"))
    if new_cpic != old_cpic:
        old_rank = _CPIC_RANK.get(old_cpic, -1)
        new_rank = _CPIC_RANK.get(new_cpic, -1)
        changes.append(_change(
            "cpic_level_changed", new_fp, gene, "cpic_level", old_cpic, new_cpic,
            "up" if new_rank > old_rank else "down",
            "CPIC actionability moved from "
            f"{old_cpic or 'no level'} to {new_cpic or 'no level'}.",
            old_display=old_cpic or "no level",
            new_display=new_cpic or "no level",
        ))

    old_carrier = old_fp.get("carrier")
    new_carrier = new_fp.get("carrier")
    if old_carrier != new_carrier:
        old_n = old_carrier if isinstance(old_carrier, int) else -1
        new_n = new_carrier if isinstance(new_carrier, int) else -1
        changes.append(_change(
            "carrier_status_changed", new_fp, gene, "carrier",
            old_carrier, new_carrier, "up" if new_n > old_n else "down",
            "Recorded copies of the variant allele moved from "
            f"{_carrier_words(old_carrier)} to {_carrier_words(new_carrier)}.",
            old_display=_carrier_words(old_carrier),
            new_display=_carrier_words(new_carrier),
        ))

    old_repute = _text(old_fp.get("repute"))
    new_repute = _text(new_fp.get("repute"))
    if new_repute != old_repute:
        old_rank = _REPUTE_RANK.get(old_repute, 1)
        new_rank = _REPUTE_RANK.get(new_repute, 1)
        changes.append(_change(
            "repute_changed", new_fp, gene, "repute", old_repute, new_repute,
            "up" if new_rank > old_rank else "down",
            "Direction of effect moved from "
            f"{old_repute or 'unset'} to {new_repute or 'unset'}.",
            old_display=old_repute or "unset",
            new_display=new_repute or "unset",
        ))

    old_eval = bool(old_fp.get("evaluable"))
    new_eval = bool(new_fp.get("evaluable"))
    if new_eval != old_eval:
        changes.append(_change(
            "newly_evaluable" if new_eval else "no_longer_evaluable", new_fp,
            gene, "evaluable", old_eval, new_eval,
            "up" if new_eval else "down",
            "This entity became evaluable against your data."
            if new_eval else
            "This entity is no longer evaluable against your data.",
            old_display="evaluable" if old_eval else "not evaluable",
            new_display="evaluable" if new_eval else "not evaluable",
        ))

    old_mag = _as_float(old_fp.get("magnitude"))
    new_mag = _as_float(new_fp.get("magnitude"))
    if old_mag != new_mag:
        old_n = old_mag if old_mag is not None else -1.0
        new_n = new_mag if new_mag is not None else -1.0
        changes.append(_change(
            "magnitude_changed", new_fp, gene, "magnitude", old_mag, new_mag,
            "up" if new_n > old_n else "down",
            "DNAInsight magnitude moved from "
            f"{'unscored' if old_mag is None else old_mag} to "
            f"{'unscored' if new_mag is None else new_mag}.",
            old_display="unscored" if old_mag is None else f"{old_mag:.1f}",
            new_display="unscored" if new_mag is None else f"{new_mag:.1f}",
        ))

    return changes


def _carrier_words(copies: Any) -> str:
    if copies == 0:
        return "no copies"
    if copies == 1:
        return "1 copy"
    if copies == 2:
        return "2 copies"
    return "unknown"


def _sort_changes(changes: list[dict]) -> list[dict]:
    """Highest severity first, then a stable alphabetical tiebreak so repeated
    runs produce identical output."""
    return sorted(changes, key=lambda c: (-c["severity"], c["key"], c["field"]))


def _version_diff(old_versions: dict, new_versions: dict) -> dict:
    changed: dict[str, dict] = {}
    for key in sorted(set(old_versions) | set(new_versions)):
        before = _text(old_versions.get(key))
        after = _text(new_versions.get(key))
        if before != after:
            changed[key] = {"old": before, "new": after,
                            "name": _source_label(key)}
    return changed


def _version_statement(changed: dict) -> str:
    """A sentence naming both database versions, not just "two scans".

    "which ClinVar release said this, and when" is the first question a real
    clinician asks, so the answer is written out rather than left implied.
    """
    if not changed:
        return ("No database version change was recorded between these two "
                "scans.")
    parts = []
    for key in sorted(changed):
        entry = changed[key]
        name = entry["name"]
        before = entry["old"] or "an unrecorded version"
        after = entry["new"] or "an unrecorded version"
        parts.append(f"{name} {before} and {name} {after}")
    return "Between " + "; ".join(parts) + "."


def _baseline_payload(new_snap: dict, profile_id: int) -> dict:
    """The first-ever snapshot has nothing to diff against.

    This is a normal state, not an error, and it has to be a well formed payload
    because the UI renders it the same way it renders a real diff. Returning an
    error here would mean the very first scan looked broken.
    """
    return {
        "ok": True,
        "baseline": True,
        "profile_id": profile_id,
        "old": None,
        "new": {k: v for k, v in new_snap.items()
                if k not in ("entries", "genes")},
        "db_versions": {
            "old": {},
            "new": dict(new_snap.get("db_versions") or {}),
            "changed": {},
        },
        "version_statement": (
            "Baseline established. There is no earlier scan to compare "
            "against yet."
        ),
        "changes": [],
        "material": [],
        "counts": {},
        "material_count": 0,
        "unchanged": 0,
        "total_old": 0,
        "total_new": int(new_snap.get("entry_count", 0)),
        "headline": "Baseline established.",
    }


def diff_snapshots(old_id: int, new_id: int) -> dict:
    """Compare two snapshots and return every classified change.

    Direction is always explicit. A change record never says only "changed": it
    names the field, both values in display form, and whether the move was up
    (more concerning or stronger evidence), down (less concerning or weaker
    evidence), or lateral (a move that genuinely does not sit on that axis, such
    as pathogenic becoming conflicting).

    ``old_id`` may be None, which produces the baseline payload rather than an
    error. Unknown ids return ``{"ok": False, "error": ...}`` rather than
    raising, because a stale id in a URL is a routine event and should not be a
    500.
    """
    new_snap = get_snapshot(new_id) if new_id is not None else None
    if new_snap is None:
        return {"ok": False, "baseline": False, "error": "snapshot_not_found",
                "field": "new_id", "changes": [], "counts": {},
                "material": [], "material_count": 0}

    if old_id is None:
        return _baseline_payload(new_snap, new_snap["profile_id"])

    old_snap = get_snapshot(old_id)
    if old_snap is None:
        return {"ok": False, "baseline": False, "error": "snapshot_not_found",
                "field": "old_id", "changes": [], "counts": {},
                "material": [], "material_count": 0}

    old_entries = old_snap.get("entries") or {}
    new_entries = new_snap.get("entries") or {}
    old_genes = old_snap.get("genes") or {}
    new_genes = new_snap.get("genes") or {}

    changes: list[dict] = []
    unchanged = 0

    for key in sorted(set(old_entries) | set(new_entries)):
        before = old_entries.get(key)
        after = new_entries.get(key)
        gene = new_genes.get(key) or old_genes.get(key) or ""

        if before is not None and after is None:
            # Reported, not silently dropped. A finding that vanishes between
            # scans is itself information: the variant may have been withdrawn
            # from ClinVar, or a pooled file may have been removed. Letting it
            # disappear without comment is how a user ends up unable to find
            # something they read last month and concludes the tool is lying.
            changes.append(_change(
                "finding_removed", before, gene, "presence", True, False, "down",
                "This finding was reported in the earlier scan and is not "
                "present in the later one.",
                old_display="present", new_display="absent",
            ))
            continue

        if before is None and after is not None:
            changes.append(_change(
                "finding_added", after, gene, "presence", False, True, "up",
                "This finding is present in the later scan and was not present "
                "in the earlier one.",
                old_display="absent", new_display="present",
            ))
            continue

        if before.get("digest") and before.get("digest") == after.get("digest"):
            unchanged += 1
            continue

        pair = _classify_pair(before, after, gene)
        if pair:
            changes.extend(pair)
        else:
            unchanged += 1

    changes = _sort_changes(changes)
    counts: dict[str, int] = {}
    for change in changes:
        counts[change["kind"]] = counts.get(change["kind"], 0) + 1

    material = [c for c in changes if c["material"]]
    changed_versions = _version_diff(old_snap.get("db_versions") or {},
                                     new_snap.get("db_versions") or {})

    return {
        "ok": True,
        "baseline": False,
        "profile_id": new_snap["profile_id"],
        "old": {k: v for k, v in old_snap.items() if k not in ("entries", "genes")},
        "new": {k: v for k, v in new_snap.items() if k not in ("entries", "genes")},
        "db_versions": {
            "old": dict(old_snap.get("db_versions") or {}),
            "new": dict(new_snap.get("db_versions") or {}),
            "changed": changed_versions,
        },
        "version_statement": _version_statement(changed_versions),
        "changes": changes,
        "material": material,
        "counts": counts,
        "material_count": len(material),
        "unchanged": unchanged,
        "total_old": len(old_entries),
        "total_new": len(new_entries),
        "headline": _headline(material, changes),
    }


def _headline(material: list[dict], changes: list[dict]) -> str:
    if material:
        top = material[0]
        gene = f" ({top['gene']})" if top["gene"] else ""
        return f"{top['label']}: {top['rsid']}{gene}."
    if changes:
        return "Only derived scores moved. No database classification changed."
    return "Nothing changed."


# ---------------------------------------------------------------------------
# Change history
# ---------------------------------------------------------------------------

def changes_for(profile_id: int, *, since: str | None = None,
                limit: int | None = None) -> dict:
    """Every change across a profile's whole snapshot history.

    Walks consecutive snapshot pairs in chronological order, so a variant that
    moved twice appears twice with the release that moved it each time.

    ``since`` filters on the LATER snapshot of each pair, compared as an ISO
    8601 string. A trailing Z is stripped from both sides first, so
    "2026-08-01" and "2026-08-01T00:00:00Z" both behave the way a caller
    expects instead of tripping over lexical ordering of the separator.

    ``limit`` caps the returned change records after sorting by severity, so a
    caller asking for 10 gets the 10 that matter, never the first 10 found.
    """
    headers = sorted(list_snapshots(profile_id), key=lambda h: (h["created_at"], h["id"]))
    if not headers:
        return {"ok": True, "profile_id": profile_id, "baseline": True,
                "snapshots": 0, "comparisons": [], "changes": [], "material": [],
                "counts": {}, "material_count": 0, "truncated": False,
                "since": since or "",
                "headline": "No scans have been recorded for this profile yet."}

    if len(headers) == 1:
        payload = _baseline_payload(get_snapshot(headers[0]["id"]), profile_id)
        payload.update({"snapshots": 1, "comparisons": [], "truncated": False,
                        "since": since or ""})
        return payload

    cutoff = _text(since).rstrip("Z") if since else ""

    all_changes: list[dict] = []
    comparisons: list[dict] = []
    for older, newer in zip(headers, headers[1:]):
        if cutoff and _text(newer["created_at"]).rstrip("Z") <= cutoff:
            continue
        diff = diff_snapshots(older["id"], newer["id"])
        if not diff.get("ok"):
            continue
        for change in diff["changes"]:
            enriched = dict(change)
            enriched["old_snapshot_id"] = older["id"]
            enriched["new_snapshot_id"] = newer["id"]
            enriched["observed_at"] = newer["created_at"]
            enriched["version_statement"] = diff["version_statement"]
            all_changes.append(enriched)
        comparisons.append({
            "old_id": older["id"],
            "new_id": newer["id"],
            "old_at": older["created_at"],
            "new_at": newer["created_at"],
            "counts": diff["counts"],
            "material_count": diff["material_count"],
            "version_statement": diff["version_statement"],
            "db_versions_changed": diff["db_versions"]["changed"],
        })

    all_changes = _sort_changes(all_changes)
    counts: dict[str, int] = {}
    for change in all_changes:
        counts[change["kind"]] = counts.get(change["kind"], 0) + 1
    material = [c for c in all_changes if c["material"]]

    truncated = False
    if isinstance(limit, int) and limit >= 0 and len(all_changes) > limit:
        all_changes = all_changes[:limit]
        truncated = True

    return {
        "ok": True,
        "baseline": False,
        "profile_id": profile_id,
        "since": since or "",
        "snapshots": len(headers),
        "comparisons": comparisons,
        "changes": all_changes,
        "material": material,
        "counts": counts,
        "material_count": len(material),
        "truncated": truncated,
        "headline": _headline(material, all_changes),
    }


# ---------------------------------------------------------------------------
# Addendum
# ---------------------------------------------------------------------------

ADDENDUM_STATEMENT = (
    "This addendum is additive. It does not replace, correct or supersede the "
    "original report, and the original report remains valid as a record of "
    "what the evidence said on the date it was produced. Read the two "
    "together."
)


def addendum(profile_id: int, *, old_id: int | None = None,
             new_id: int | None = None) -> dict:
    """Build the dated, additive addendum payload for a profile.

    Defaults compare the two most recent snapshots. With only one snapshot on
    record it returns the well formed baseline payload, because the first scan
    genuinely has nothing to compare against and that is not an error.

    ADDITIVE MEANS ADDITIVE. ``supersedes`` is always None and there is no code
    path anywhere in this module that writes to a prior snapshot or to a prior
    report. A clinician who acted on the original report must be able to see
    exactly what that report said, because their decision is only defensible
    against the evidence that existed at the time.
    """
    headers = list_snapshots(profile_id)

    if new_id is None:
        if not headers:
            return {
                "ok": False,
                "kind": "addendum",
                "error": "no_snapshots",
                "profile_id": profile_id,
                "generated_at": _utc_now(),
                "additive": True,
                "supersedes": None,
                "statement": ADDENDUM_STATEMENT,
                "changes": [], "material": [], "counts": {},
                "material_count": 0,
                "headline": "No scans have been recorded for this profile yet.",
            }
        new_id = headers[0]["id"]

    if old_id is None:
        ordered = [h["id"] for h in headers]
        try:
            position = ordered.index(new_id)
        except ValueError:
            position = -1
        # headers are newest first, so the predecessor sits one further along.
        if position >= 0 and position + 1 < len(ordered):
            old_id = ordered[position + 1]

    diff = diff_snapshots(old_id, new_id)
    payload = dict(diff)
    payload.update({
        "kind": "baseline" if diff.get("baseline") else "addendum",
        "profile_id": profile_id,
        "generated_at": _utc_now(),
        "additive": True,
        "supersedes": None,
        "statement": ADDENDUM_STATEMENT,
    })
    payload.setdefault("material", [])
    payload.setdefault("material_count", 0)
    return payload


# ---------------------------------------------------------------------------
# Offline HTML rendering
#
# The project rule is that a generated report is fully self contained: no CDN,
# no external font, no remote image, no tracking pixel, nothing that turns
# opening a genetic report into a network event that some third party can see.
# The fragment below carries its own inline style and nothing else. There is
# not a single URL in it.
# ---------------------------------------------------------------------------

_HTML_ESCAPES = {
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
}


def _esc(value: Any) -> str:
    """HTML-escape any value. A raw DNA file is attacker controllable in
    principle, so an interpretation string is never trusted as markup."""
    if value is None:
        return ""
    return "".join(_HTML_ESCAPES.get(c, c) for c in str(value))


_DIRECTION_STYLE = {
    "up": ("#b3261e", "#fdecea", "increase"),
    "down": ("#1b5e20", "#e8f5e9", "decrease"),
    "lateral": ("#7a4f01", "#fff4e0", "sideways"),
}

_ADDENDUM_CSS = """
.dnai-addendum{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;
 color:#1a1a1a;line-height:1.5;max-width:60rem}
.dnai-addendum h2{font-size:1.35rem;margin:0 0 .25rem}
.dnai-addendum .dnai-meta{font-size:.85rem;color:#555;margin:0 0 .75rem}
.dnai-addendum .dnai-note{background:#f2f4f8;border-left:4px solid #4a5568;
 padding:.6rem .8rem;margin:0 0 1rem;font-size:.9rem}
.dnai-addendum .dnai-versions{font-size:.9rem;margin:0 0 1rem;font-weight:600}
.dnai-addendum table{border-collapse:collapse;width:100%;font-size:.9rem}
.dnai-addendum th,.dnai-addendum td{border-bottom:1px solid #ddd;
 padding:.45rem .5rem;text-align:left;vertical-align:top}
.dnai-addendum th{background:#eef1f5;font-weight:600}
.dnai-addendum .dnai-kind{display:inline-block;border-radius:3px;
 padding:.05rem .4rem;font-size:.8rem;font-weight:600;white-space:nowrap}
.dnai-addendum .dnai-move{font-family:ui-monospace,Consolas,monospace;
 font-size:.85rem}
.dnai-addendum .dnai-gene{color:#555;font-size:.85rem}
.dnai-addendum .dnai-empty{color:#555;font-style:italic}
.dnai-addendum .dnai-counts{font-size:.85rem;color:#555;margin-top:.75rem}
"""


def _render_row(change: dict) -> str:
    colour, background, arrow_word = _DIRECTION_STYLE.get(
        change.get("direction", "lateral"), _DIRECTION_STYLE["lateral"])
    gene = change.get("gene") or ""
    label = _esc(change.get("rsid", ""))
    if gene:
        label += f" <span class=\"dnai-gene\">{_esc(gene)}</span>"
    return (
        "<tr>"
        f"<td><span class=\"dnai-kind\" style=\"color:{colour};"
        f"background:{background}\">{_esc(change.get('label', ''))}</span></td>"
        f"<td>{label}</td>"
        f"<td class=\"dnai-move\">{_esc(change.get('old_display', ''))} "
        f"to {_esc(change.get('new_display', ''))}</td>"
        f"<td>{_esc(arrow_word)}</td>"
        f"<td>{_esc(change.get('detail', ''))}</td>"
        "</tr>"
    )


def render_addendum_html(addendum_payload: dict) -> str:
    """Render an addendum as a self contained HTML fragment.

    Fragment, not a document: the caller drops it into an existing report shell.
    It carries its own inline style block and makes no network request of any
    kind. There is no CDN reference, no external font, no remote image and no
    link out. That is the same offline rule the interactive report follows, and
    it is the reason a user can open a genetic report on a plane, or on a
    machine that is deliberately air gapped, and see exactly what they would see
    online.
    """
    payload = addendum_payload if isinstance(addendum_payload, dict) else {}
    generated = payload.get("generated_at") or _utc_now()
    baseline = bool(payload.get("baseline"))
    changes = payload.get("changes") or []
    material = payload.get("material") or []
    counts = payload.get("counts") or {}

    heading = "Baseline established" if baseline else "Reclassification addendum"
    parts: list[str] = [
        "<section class=\"dnai-addendum\">",
        f"<style>{_ADDENDUM_CSS}</style>",
        f"<h2>{_esc(heading)}</h2>",
        f"<p class=\"dnai-meta\">Generated {_esc(generated)} (UTC). "
        f"Profile {_esc(payload.get('profile_id', ''))}.</p>",
        f"<p class=\"dnai-note\">{_esc(payload.get('statement') or ADDENDUM_STATEMENT)}</p>",
        f"<p class=\"dnai-versions\">{_esc(payload.get('version_statement') or '')}</p>",
    ]

    if baseline:
        parts.append(
            "<p class=\"dnai-empty\">This is the first recorded scan for this "
            "profile, so there is nothing to compare it against yet. Future "
            "scans will be compared to this one.</p>"
        )
    elif not changes:
        parts.append(
            "<p class=\"dnai-empty\">No finding changed between these two "
            "scans.</p>"
        )
    else:
        shown = material if material else changes
        parts.append(
            "<table><thead><tr><th>What changed</th><th>Finding</th>"
            "<th>Old to new</th><th>Direction</th><th>Detail</th></tr></thead>"
            "<tbody>"
        )
        parts.extend(_render_row(c) for c in shown)
        parts.append("</tbody></table>")
        if material and len(changes) > len(material):
            parts.append(
                "<p class=\"dnai-counts\">"
                f"{len(changes) - len(material)} further change(s) affected only "
                "the DNAInsight magnitude, which is a derived score rather than "
                "a database classification, and are not listed here.</p>"
            )

    if counts:
        summary = ", ".join(f"{_esc(kind)}: {int(n)}"
                            for kind, n in sorted(counts.items()))
        parts.append(f"<p class=\"dnai-counts\">{summary}</p>")

    parts.append("</section>")
    return "".join(parts)
