"""
provenance.py -- the provenance graph, the licence contract and the signed
report manifest.

WHY THIS MODULE EXISTS
----------------------
``backend/scoring.py`` already gives every finding a ``magnitude_factors``
audit trail, so a reader can see exactly why a number came out the way it did.
That trail stops at the score. It says "base 6.00 from clinvar_path_3star" and
it does not say WHICH ClinVar release asserted that, on what date it was
downloaded, or under what licence it arrived.

"Which ClinVar release said this, and when" is the first question a real
clinician asks. Today the honest answer is that nobody can reconstruct it,
because a report is a snapshot of databases that all move independently and
none of them are named in the output. A report that cannot be traced back to
its evidence is an opinion.

So two things live here:

1. A PROVENANCE GRAPH. Every finding can be annotated with the sources behind
   it, each carrying name, licence, SPDX identifier, version and retrieval
   date. ``detect_conflicts`` then surfaces disagreement between those sources
   rather than resolving it.

2. A SIGNED MANIFEST. Every report can carry a record of the DNAInsight
   version, every database name plus version plus retrieval date plus licence,
   the sha256 of every input DNA file, the finding count, the scan parameters
   and a UTC timestamp, HMAC-signed so drift is detectable. A clinician holding
   the report can reproduce it exactly or prove it has changed.

DISAGREEMENT IS DISPLAYED, NOT RESOLVED
---------------------------------------
``detect_conflicts`` emits both positions and NO verdict. When ClinVar says
benign and the GWAS Catalog carries a replicated risk association, the honest
output is "these two sources disagree, here is each one", not a silently
computed winner. When a CPIC actionability level and a ClinVar significance
point opposite ways, same answer.

This is the same invariant ``backend/merge.py`` enforces when two pooled files
disagree at a position: BOTH calls are retained and surfaced, with no voting,
no confidence weighting and no automatic winner. Disagreement is information,
not noise to be hidden. It is a deliberate design invariant of this project,
and picking a winner here would quietly manufacture a certainty that the
evidence does not contain.

THE LICENCE CONTRACT
--------------------
``SOURCES`` is the machine-readable form of ``data/DATA_SOURCES.md``. That
document is currently a promise kept by whoever last read it. ``licence_audit``
turns it into a runtime contract: any source marked non-redistributable, or
marked not used at all, that appears in a declared bundled artefact is a
violation with a named artefact and a named reason.

PharmGKB is in ``SOURCES`` precisely because it is NOT used. An exclusion that
exists only as prose is an exclusion nobody can test. Recording it as data means
``licence_audit`` can prove the exclusion still holds, and means anyone who adds
PharmGKB to a bundled artefact later trips a failing check instead of quietly
relicensing the repository away from MIT.

WHAT THE SIGNATURE DOES AND DOES NOT PROVE
------------------------------------------
The manifest is signed with HMAC-SHA256 using a secret generated on this
machine and stored at ``~/.dnainsight/manifest.key``. That proves ONE thing:
the report was not altered after generation on this machine, by anyone who does
not hold that key file.

It is NOT a public-key attestation. It does not prove authorship to a third
party. Anyone holding the key file can produce a valid signature over any
content, so a recipient who does not trust the sender learns nothing from it.
Overclaiming here would be worse than not signing at all, because a clinician
who believes a report is cryptographically attested to an author will stop
asking questions they should still be asking.

OFFLINE
-------
Nothing here touches the network, at import time or on any call path. Standard
library only, plus ``backend.database`` and ``backend.scoring``.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from . import APP_VERSION
from .database import get_connection, _add_column_if_missing

__all__ = [
    "SOURCES", "BUNDLED_ARTEFACTS", "SIGNATURE_SCOPE", "MANIFEST_VERSION",
    "init_provenance", "content_hash", "text_hash", "database_versions",
    "source_record", "attach", "detect_conflicts", "build_manifest",
    "manifest_key", "sign_manifest", "verify_manifest",
    "render_manifest_text", "licence_audit",
]


MANIFEST_VERSION = "1"

SIGNATURE_SCOPE = (
    "HMAC-SHA256 over the canonical manifest, keyed by a secret generated on "
    "and stored only on this machine. This proves the report was not altered "
    "after generation on this machine. It is NOT a public-key attestation, it "
    "does not prove authorship to a third party, and anyone holding the key "
    "file can sign anything."
)


# ---------------------------------------------------------------------------
# SOURCES: data/DATA_SOURCES.md, made machine readable.
#
# Keys per entry:
#   name             what the source calls itself
#   licence          the human statement of the terms, including the extra
#                    clauses an SPDX identifier cannot express
#   spdx             the SPDX identifier, or "" when the project publishes none.
#                    "" is itself information: it means the terms had to be read
#                    by hand, exactly as backend/external.py records
#   url              where to verify it
#   role             bundled / bundled_and_fetched / fetched_local / live_api /
#                    opt_in_local / not_used
#   redistributable  may this project ship it inside the repository
#   commercial_ok    may a downstream user sell something containing it
#   verified         the date a human last read the terms at source
#   used             is it used at all. False is a recorded decision, not an
#                    oversight
#   never_bundle     stronger than "not redistributable": bundling would
#                    relicense the repository
#   per_record       the licence varies per record or per field, so the file
#                    level statement is not sufficient on its own
#   note             the caveat that actually bites
#
# A licence that nobody re-reads drifts, so every entry carries its own
# verification date rather than relying on the document's header date.
# ---------------------------------------------------------------------------

def _source(*, name: str, licence: str, spdx: str, url: str, role: str,
            redistributable: bool, commercial_ok: bool, verified: str,
            used: bool, never_bundle: bool, per_record: bool,
            note: str) -> dict[str, Any]:
    """Build one SOURCES entry with every licence field stated.

    Keyword-only, and deliberately without a single default value. A default
    here would be the exact failure this table exists to prevent: a source added
    next year would inherit "redistributable, commercial use fine" from whoever
    wrote this function rather than from somebody who read the terms. Leaving a
    field out is a TypeError at import, which is a far better outcome than a
    licence claim nobody ever made.
    """
    return {
        "name": name,
        "licence": licence,
        "spdx": spdx,
        "url": url,
        "role": role,
        "redistributable": redistributable,
        "commercial_ok": commercial_ok,
        "verified": verified,
        "used": used,
        "never_bundle": never_bundle,
        "per_record": per_record,
        "note": note,
    }


SOURCES: dict[str, dict[str, Any]] = {
    "cpic": _source(
        name="CPIC",
        licence="CC0-1.0 public domain dedication. No conditions.",
        spdx="CC0-1.0",
        url="https://cpicpgx.org/",
        role="bundled", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=False,
        note=(
            "The CPIC dump carries clinpgxlevel and pgxtesting columns that are "
            "NOT CC0. They originate from PharmGKB and inherit its "
            "no-commercial-sale clause. Only genesymbol, drugname and cpiclevel "
            "are ever read, and the builder excludes the other two in the HTTP "
            "request itself so they are never transferred."
        ),
    ),
    "clinvar": _source(
        name="ClinVar",
        licence=(
            "US Government work, public domain in the United States under "
            "17 USC 105. Citation requested, not required."
        ),
        spdx="",
        url="https://www.ncbi.nlm.nih.gov/clinvar/",
        role="bundled_and_fetched", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=False,
        note=(
            "Clinical assertions change between monthly releases, so a stale "
            "copy can be actively misleading. This is the reason "
            "backend/ledger.py exists. Individual submitter records may also "
            "carry that submitter's own consent terms for onward use."
        ),
    ),
    "gnomad": _source(
        name="gnomAD",
        licence=(
            "CC0-1.0 for aggregate frequency data. Citation requested as a "
            "courtesy."
        ),
        spdx="CC0-1.0",
        url="https://gnomad.broadinstitute.org/",
        role="bundled", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=False,
        note=(
            "Only AGGREGATE frequencies are CC0 and only aggregates are used. "
            "Individual-level gnomAD data is not public and is not touched."
        ),
    ),
    "onekg_ensembl": _source(
        name="1000 Genomes via Ensembl",
        licence=(
            "Open with no restriction on use, per the project data reuse "
            "statement. Ensembl adds no usage restriction and asks only for "
            "citation."
        ),
        spdx="",
        url="https://www.internationalgenome.org/",
        role="bundled", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=False,
        note="Per-population frequency breakdowns, complementing gnomAD.",
    ),
    "gwas_catalog": _source(
        name="GWAS Catalog",
        licence=(
            "EMBL-EBI terms of use. Freely available, usable for research and "
            "commercial purposes, redistributable. Citation requested."
        ),
        spdx="",
        url="https://www.ebi.ac.uk/gwas/",
        role="fetched_local", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=False,
        note=(
            "Not committed for size reasons, not licence reasons. EMBL-EBI "
            "disclaims warranty and asks that the catalogue not be "
            "misrepresented as EBI-endorsed analysis."
        ),
    ),
    "pgs_catalog": _source(
        name="PGS Catalog",
        licence=(
            "EMBL-EBI terms of use for the catalogue as a whole, with PER-SCORE "
            "licence overrides set by submitting authors. The overrides win."
        ),
        spdx="",
        url="https://www.pgscatalog.org/",
        role="fetched_local", redistributable=True, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=False, per_record=True,
        note=(
            "Roughly 31 scores are CC BY-NC-ND, which forbids both commercial "
            "use and derivative works. A no-derivatives term is incompatible "
            "with reweighting or subsetting a score, which is exactly what a "
            "scoring engine does. Those scores must be filtered out before any "
            "PGS data is committed. Check each score's own license field; do "
            "not reason from the catalogue default."
        ),
    ),
    "myvariant": _source(
        name="MyVariant.info",
        licence=(
            "Apache-2.0 for the service code. The DATA is not under a single "
            "licence: each field keeps the licence of the upstream source it "
            "came from."
        ),
        spdx="Apache-2.0",
        url="https://myvariant.info/",
        role="live_api", redistributable=False, commercial_ok=True,
        verified="2026-07-27", used=True, never_bundle=True, per_record=True,
        note=(
            "NEVER PERSISTED into a bundled artefact. Because the licence is "
            "per field, a single response can blend US public domain ClinVar "
            "fields, CC0 gnomAD fields and fields with more restrictive terms. "
            "Use it live, show it to the user, do not bake it in. If a field is "
            "wanted permanently, go to that field's upstream source directly."
        ),
    ),
    "snpedia": _source(
        name="SNPedia",
        licence="CC-BY-NC-SA-3.0-US. Attribution, non-commercial and share-alike all bite.",
        spdx="CC-BY-NC-SA-3.0-US",
        url="https://www.snpedia.com/",
        role="opt_in_local", redistributable=False, commercial_ok=False,
        verified="2026-07-27", used=True, never_bundle=True, per_record=False,
        note=(
            "NEVER BUNDLED. Not in data/, not in a fixture, not in a test file, "
            "not in a docstring example. A harvested copy inside this repository "
            "would trigger share-alike and relicense the project away from MIT, "
            "and the non-commercial term would strip the right to sell anything "
            "built on it. Offered instead as an opt-in local harvest that writes "
            "to ~/.dnainsight/, outside the repository tree."
        ),
    ),
    "pharmgkb": _source(
        name="PharmGKB / ClinPGx",
        licence=(
            "Claims CC-BY-SA-4.0, but the accompanying data use agreement adds "
            "a term prohibiting sale of the data or of products containing it. "
            "That restriction is not part of CC-BY-SA-4.0 and cannot be added "
            "to it, so the effective licence is a bespoke non-commercial one "
            "and is NOT the standard identifier it appears to be."
        ),
        spdx="",
        url="https://www.pharmgkb.org/",
        role="not_used", redistributable=False, commercial_ok=False,
        verified="2026-07-27", used=False, never_bundle=True, per_record=False,
        note=(
            "USED FOR NOTHING. Recorded here so the exclusion is a decision "
            "rather than an oversight, and so licence_audit can prove it still "
            "holds. Bundling it would forbid downstream commercial use, which "
            "conflicts with the MIT grant this project makes. It leaks in from "
            "the CPIC dump via the clinpgxlevel and pgxtesting columns, which "
            "look safe by association and are not. Replaced by CPIC for "
            "actionability, ClinVar for significance and FDA label tiers for "
            "pharmacogenomic labelling."
        ),
    ),
}


# ---------------------------------------------------------------------------
# BUNDLED_ARTEFACTS: what is actually committed to this repository, and which
# sources each artefact derives from.
#
# This is the input to licence_audit. Adding a source to an artefact here and
# forgetting to check its licence is what the audit catches. Tests construct a
# violation by inserting an entry, which is the point: the contract has to be
# falsifiable or it is decoration.
# ---------------------------------------------------------------------------

BUNDLED_ARTEFACTS: dict[str, list[str]] = {
    "data/snp_reference.json": ["cpic", "clinvar", "gwas_catalog"],
    "data/evidence_overlay.py": ["cpic", "clinvar"],
    "data/frequencies.json": ["gnomad", "onekg_ensembl"],
    "data/prs_models.json": ["pgs_catalog"],
    "data/genosets.json": ["clinvar", "gwas_catalog"],
    "data/snp_reference_sources.csv": ["cpic", "clinvar", "gwas_catalog"],
}


# ---------------------------------------------------------------------------
# Time and hashing
# ---------------------------------------------------------------------------

def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _utc_date() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def _canonical(payload: Any) -> str:
    """Stable JSON for hashing and signing.

    Sorted keys, no incidental whitespace, ensure_ascii off so the bytes are
    the same on any platform. Signing a non-canonical serialisation is the
    classic way to build a signature that verifies on the machine that made it
    and nowhere else.
    """
    return json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, default=str)


def text_hash(text: str) -> str:
    """sha256 of a string, encoded as UTF-8."""
    return hashlib.sha256(str(text if text is not None else "").encode("utf-8")).hexdigest()


def content_hash(path) -> str | None:
    """sha256 of a file, streamed. None when the file is absent or unreadable.

    Streamed in 1 MiB chunks because a raw DNA export is routinely 20 to 200 MB
    and reading one into memory to hash it is how a hashing helper turns into an
    out-of-memory crash on a small machine.

    Returns None rather than raising for a missing, directory or permission
    denied path. A manifest that records "this input file could not be hashed"
    is more useful than a report generation that dies, and the None is visible
    in the manifest rather than silently omitted.
    """
    if path is None:
        return None
    try:
        target = Path(path)
    except TypeError:
        return None
    try:
        if not target.is_file():
            return None
        digest = hashlib.sha256()
        with open(target, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except (OSError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Schema
#
# The provenance table records which sources produced which finding for which
# profile, so a report can be traced after the fact even if the finding row has
# since been rescored. Same durability rules as everywhere else in this project:
# CREATE TABLE IF NOT EXISTS plus _add_column_if_missing, never a DROP, never a
# rebuild. See the long comment in backend/database.py about the defect that
# deleted every user's database.
# ---------------------------------------------------------------------------

_SCHEMA = """
    CREATE TABLE IF NOT EXISTS provenance_records (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        profile_id  INTEGER NOT NULL,
        rsid        TEXT NOT NULL DEFAULT '',
        source_id   TEXT NOT NULL,
        version     TEXT NOT NULL DEFAULT '',
        retrieved   TEXT NOT NULL DEFAULT '',
        licence     TEXT NOT NULL DEFAULT '',
        spdx        TEXT NOT NULL DEFAULT '',
        recorded_at TEXT NOT NULL DEFAULT ''
    );

    CREATE INDEX IF NOT EXISTS idx_provenance_profile
        ON provenance_records(profile_id, rsid);

    CREATE TABLE IF NOT EXISTS report_manifests (
        id          INTEGER PRIMARY KEY AUTOINCREMENT,
        profile_id  INTEGER NOT NULL,
        report_type TEXT NOT NULL DEFAULT '',
        created_at  TEXT NOT NULL,
        manifest    TEXT NOT NULL DEFAULT '{}',
        signature   TEXT NOT NULL DEFAULT '{}'
    );

    CREATE INDEX IF NOT EXISTS idx_report_manifests_profile
        ON report_manifests(profile_id, created_at);
"""


def init_provenance() -> None:
    """Create the provenance tables if they do not exist.

    Idempotent and safe to call on an existing database that already holds
    rows. Call it at startup next to ``database.init_db``.
    """
    conn = get_connection()
    try:
        conn.executescript(_SCHEMA)
        # Forward migrations for a database written by an earlier build. No-ops
        # on a fresh one. SQLite requires a DEFAULT on a NOT NULL added column.
        _add_column_if_missing(conn, "provenance_records", "version", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "provenance_records", "retrieved", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "provenance_records", "licence", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "provenance_records", "spdx", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "provenance_records", "recorded_at", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "report_manifests", "report_type", "TEXT NOT NULL DEFAULT ''")
        _add_column_if_missing(conn, "report_manifests", "signature", "TEXT NOT NULL DEFAULT '{}'")
        conn.commit()
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Live database versions
# ---------------------------------------------------------------------------

_BASE = Path(__file__).resolve().parent.parent
_DATA = _BASE / "data"


def _json_meta(filename: str) -> dict:
    """The _meta block of a bundled JSON artefact, or {} when unavailable.

    Reads the file directly rather than importing the module that owns it,
    because those modules cache aggressively and a cached corpus loaded before a
    rebuild would report a stale version, which is precisely the failure this
    function exists to prevent.
    """
    path = _DATA / filename
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
    except (OSError, ValueError):
        return {}
    if isinstance(raw, dict) and isinstance(raw.get("_meta"), dict):
        return raw["_meta"]
    return {}


def _reference_db_meta() -> dict:
    """The meta table of the Tier 2 reference database, or {} when absent.

    Opened read-only through a file: URI. This database belongs to the user's
    build, not to us, and a version probe must never be able to create, modify
    or lock it.
    """
    path = _DATA / "reference.db"
    if not path.exists():
        return {}
    try:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            rows = conn.execute("SELECT key, value FROM meta").fetchall()
        finally:
            conn.close()
    except sqlite3.Error:
        return {}
    return {str(k): str(v) for k, v in rows}


def _meta_version(meta: dict) -> str:
    """The identity version of a bundled ARTEFACT: its own semver if it has
    one, otherwise the date it was built."""
    for key in ("version", "built_at", "updated_at", "build_date"):
        value = str(meta.get(key) or "").strip()
        if value:
            return value
    return ""


def _meta_date(meta: dict) -> str:
    """The retrieval date behind a bundled artefact.

    Used for the upstream SOURCES rather than the artefact semver, because
    ClinVar, CPIC and the GWAS Catalog publish rolling files with no version
    string in the payload. The fetch date is the only honest version marker
    available, which is the same call data/build_full_reference.py made. An
    artefact semver of "2.0.0" says nothing at all about which ClinVar release
    is inside it, and reporting it as the ClinVar version would be worse than
    reporting nothing.
    """
    for key in ("build_date", "built_at", "updated_at"):
        value = str(meta.get(key) or "").strip()
        if value:
            return value
    return ""


def database_versions() -> dict:
    """The version of every data source currently in force, as flat strings.

    Returns ``{key: version}`` where every value is a string, so the mapping
    round-trips through JSON unchanged and can be compared directly by
    ``backend.ledger`` when it needs to say "between ClinVar 2026-07-27 and
    ClinVar 2026-08-24".

    Keys are every id in ``SOURCES`` plus the derived local artefacts and the
    application version. An empty string means "no version could be
    established", which is honest and distinguishable from a missing key.

    Version resolution order per source is Tier 2 reference database first
    (it is the freshest and it records the fetch date per source), then the
    bundled artefact that source feeds. ClinVar and the GWAS Catalog publish
    rolling files with no version string in the payload, so the fetch date is
    the only honest version marker available, which is the same decision
    data/build_full_reference.py made.
    """
    ref_meta = _reference_db_meta()
    snp_meta = _json_meta("snp_reference.json")
    freq_meta = _json_meta("frequencies.json")
    prs_meta = _json_meta("prs_models.json")
    geno_meta = _json_meta("genosets.json")

    snp_date = _meta_date(snp_meta)
    freq_date = _meta_date(freq_meta)

    versions: dict[str, str] = {
        "dnainsight": str(APP_VERSION),
        "clinvar": ref_meta.get("clinvar_version", "") or snp_date,
        "cpic": ref_meta.get("cpic_version", "") or snp_date,
        "gwas_catalog": ref_meta.get("gwas_version", "") or "",
        "gnomad": freq_date,
        "onekg_ensembl": freq_date,
        "pgs_catalog": _meta_date(prs_meta),
        # An API has no version, only a moment. Saying "live" is honest; saying
        # a date would imply a pinned snapshot that does not exist.
        "myvariant": "live",
        # Never bundled, and the local opt-in cache is the user's own copy, so
        # this project records no version for it.
        "snpedia": "",
        # Not used. Deliberately empty, and licence_audit proves it stays that
        # way.
        "pharmgkb": "",
        "snp_reference": _meta_version(snp_meta),
        "frequencies": _meta_version(freq_meta),
        "genosets": _meta_version(geno_meta),
        "prs_models": _meta_version(prs_meta),
        "reference_db": ref_meta.get("build_date", "") or "",
    }
    return {k: str(v or "") for k, v in versions.items()}


# ---------------------------------------------------------------------------
# Source records
# ---------------------------------------------------------------------------

def source_record(source_id: str, *, version: str = "", retrieved: str = "") -> dict:
    """A full provenance record for one source.

    Carries the licence and the SPDX identifier on every record, because a
    provenance entry that names a database without naming its terms is exactly
    the gap DATA_SOURCES.md exists to close.

    An unknown source id FAILS CLOSED: ``known`` is False, ``redistributable``
    and ``commercial_ok`` are both False, and the licence reads "unknown". A
    source nobody has assessed must never default to permissive, because the
    default is what a careless caller will ship.
    """
    key = str(source_id or "").strip()
    entry = SOURCES.get(key)

    if entry is None:
        return {
            "source_id": key,
            "name": key or "unknown source",
            "licence": "unknown, not assessed",
            "spdx": "",
            "url": "",
            "role": "unknown",
            "redistributable": False,
            "commercial_ok": False,
            "verified": "",
            "used": False,
            "never_bundle": True,
            "per_record": False,
            "note": (
                "This source id is not recorded in data/DATA_SOURCES.md. "
                "Treated as non-redistributable until it is assessed."
            ),
            "known": False,
            "version": str(version or ""),
            "retrieved": str(retrieved or _utc_date()),
        }

    live = database_versions()
    record = {
        "source_id": key,
        "name": entry["name"],
        "licence": entry["licence"],
        "spdx": entry["spdx"],
        "url": entry["url"],
        "role": entry["role"],
        "redistributable": bool(entry["redistributable"]),
        "commercial_ok": bool(entry["commercial_ok"]),
        "verified": entry["verified"],
        "used": bool(entry.get("used", True)),
        "never_bundle": bool(entry.get("never_bundle", False)),
        "per_record": bool(entry.get("per_record", False)),
        "note": entry.get("note", ""),
        "known": True,
        "version": str(version or live.get(key, "")),
        "retrieved": str(retrieved or _utc_date()),
    }
    return record


def _as_source_ids(source_ids: Any) -> list[str]:
    if source_ids is None:
        return []
    if isinstance(source_ids, str):
        candidates: Iterable[Any] = [source_ids]
    elif isinstance(source_ids, dict):
        candidates = list(source_ids.keys())
    elif isinstance(source_ids, Iterable):
        candidates = source_ids
    else:
        candidates = [source_ids]

    seen: list[str] = []
    for candidate in candidates:
        value = str(candidate or "").strip()
        if value and value not in seen:
            seen.append(value)
    return seen


def attach(finding: dict, source_ids, *, versions: dict | None = None) -> dict:
    """Attach a provenance block to a finding, in place, and return it.

    In place and returning the same dict, matching ``scoring.score_finding``.
    One convention across the codebase beats two.

    Sets ``finding["provenance"]`` and ``finding["source_ids"]``. It does NOT
    touch the existing ``sources`` list: that field is written by the scanner
    with its own vocabulary ("bundled_reference", "genoset", "prs") and is read
    by scoring, filters and both report generators. Overwriting it here to mean
    something slightly different would break all of them for no gain.

    ``versions`` overrides the live version per source id, which is what a
    replay or a regression fixture needs so the provenance block is
    deterministic rather than dependent on whatever is on disk today.
    """
    if not isinstance(finding, dict):
        return finding

    overrides = versions if isinstance(versions, dict) else {}
    ids = _as_source_ids(source_ids)
    records = [source_record(sid, version=str(overrides.get(sid, "")))
               for sid in ids]

    licences = sorted({r["licence"] for r in records if r["licence"]})
    finding["source_ids"] = ids
    finding["provenance"] = {
        "sources": records,
        "licences": licences,
        "spdx": sorted({r["spdx"] for r in records if r["spdx"]}),
        # All of these are ANDs. One non-redistributable source contaminates the
        # whole finding, which is the "contamination travels with the column"
        # rule from DATA_SOURCES.md applied at the record level.
        "redistributable": all(r["redistributable"] for r in records) if records else True,
        "commercial_ok": all(r["commercial_ok"] for r in records) if records else True,
        "persistable": not any(r["never_bundle"] for r in records),
        "unknown_sources": [r["source_id"] for r in records if not r["known"]],
        "attached_at": _utc_now(),
    }
    return finding


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------

def _position(source_id: str, claim: str, detail: str = "",
              evidence: str = "") -> dict:
    entry = SOURCES.get(source_id)
    return {
        "source_id": source_id,
        "source": entry["name"] if entry else source_id,
        "claim": claim,
        "detail": detail,
        "evidence": evidence,
    }


def _conflict(kind: str, positions: list[dict], note: str) -> dict:
    """A conflict record. Always two or more positions and never a verdict."""
    return {
        "type": kind,
        "positions": positions,
        # Explicitly None and explicitly False rather than absent. A consumer
        # that goes looking for a winner must find a stated refusal, not a
        # missing key it can interpret as "no conflict after all".
        "verdict": None,
        "resolved": False,
        "note": note,
        "display_only": True,
    }


def _gwas_replicated(finding: dict) -> bool:
    """Mirrors scoring's replicated-GWAS test. Two independent studies, or an
    explicit gwas source tag."""
    try:
        studies = int(finding.get("gwas_studies") or 0)
    except (TypeError, ValueError):
        studies = 0
    if studies >= 2:
        return True
    sources = finding.get("sources") or []
    if isinstance(sources, str):
        sources = [sources]
    return any("gwas" in str(s).lower() for s in sources)


def detect_conflicts(finding: dict) -> list[dict]:
    """Surface every disagreement between the sources behind one finding.

    DISPLAYS disagreement, never resolves it. Every returned record carries
    both positions, ``verdict`` None and ``resolved`` False. There is no code
    path in this function that picks a winner, and there should never be one.

    This mirrors ``backend/merge.py``: when two pooled files disagree at a
    position, BOTH calls are retained and surfaced, with no voting, no
    confidence weighting and no automatic winner. Disagreement is information,
    not noise to be hidden. It is a deliberate design invariant of this project.
    Silently choosing would manufacture a certainty the evidence does not
    contain, and the user would never know it had happened.

    Returns a list ordered by conflict type, so repeated calls produce
    identical output.
    """
    if not isinstance(finding, dict):
        return []

    from .scoring import clinvar_sig_code, normalize_cpic_level

    conflicts: list[dict] = []

    raw_sig = str(finding.get("clinical_sig") or "").strip()
    code = finding.get("clinvar_sig_code")
    if not isinstance(code, int) or isinstance(code, bool):
        code = clinvar_sig_code(raw_sig) if raw_sig else None
    cpic = normalize_cpic_level(finding.get("cpic_level"))
    stars = finding.get("review_stars")
    stars = stars if isinstance(stars, int) and not isinstance(stars, bool) else 0

    benign = code in (2, 3)
    pathogenic = code in (4, 5)

    # 1. ClinVar says benign, the GWAS Catalog carries a replicated association.
    if benign and _gwas_replicated(finding):
        traits = finding.get("gwas_traits") or ""
        if isinstance(traits, (list, tuple)):
            traits = "; ".join(str(t) for t in traits)
        conflicts.append(_conflict(
            "clinvar_benign_vs_gwas_association",
            [
                _position("clinvar", raw_sig or "benign",
                          "Clinically benign for the asserted condition.",
                          f"{stars} star review status"),
                _position("gwas_catalog", "replicated risk association",
                          str(traits) or "Replicated trait association.",
                          f"{finding.get('gwas_studies') or 2} independent studies"),
            ],
            "ClinVar assesses Mendelian pathogenicity; the GWAS Catalog "
            "records population-level association. Both can be correct at "
            "once. DNAInsight shows both and does not choose.",
        ))

    # 2. CPIC actionability and ClinVar significance pointing opposite ways.
    if cpic in ("A", "A/B", "B") and benign:
        conflicts.append(_conflict(
            "cpic_actionable_vs_clinvar_benign",
            [
                _position("cpic", f"CPIC Level {cpic}",
                          "Actionable for prescribing.",
                          "CPIC assignment"),
                _position("clinvar", raw_sig or "benign",
                          "Benign for disease causation.",
                          f"{stars} star review status"),
            ],
            "A variant can be benign for disease and still matter for a drug "
            "dose. DNAInsight shows both positions rather than collapsing them.",
        ))
    if cpic in ("D", "Retired") and pathogenic:
        conflicts.append(_conflict(
            "cpic_no_action_vs_clinvar_pathogenic",
            [
                _position("cpic", f"CPIC Level {cpic}",
                          "No prescribing action recommended.",
                          "CPIC assignment"),
                _position("clinvar", raw_sig or "pathogenic",
                          "Pathogenic for the asserted condition.",
                          f"{stars} star review status"),
            ],
            "CPIC grades prescribing actionability, not pathogenicity. Both "
            "positions are shown.",
        ))

    # 3. ClinVar disagreeing with itself. scoring.py already refuses to resolve
    #    a conflicting record to one of the positions it conflicts over; this
    #    surfaces the same fact to the reader.
    if "conflicting" in raw_sig.lower():
        conflicts.append(_conflict(
            "clinvar_submitters_conflicting",
            [
                _position("clinvar", raw_sig,
                          "Submitters do not agree with each other.",
                          f"{stars} star review status"),
                _position("clinvar", "no aggregate call",
                          "ClinVar publishes no consensus for this variant.",
                          "aggregate record"),
            ],
            "ClinVar itself records disagreement here. It is not resolved to "
            "either position, by ClinVar or by DNAInsight.",
        ))

    # 4. A curated SNPedia repute contradicting the ClinVar classification.
    snp_repute = str(finding.get("snpedia_repute") or "").strip()
    if snp_repute == "Good" and pathogenic:
        conflicts.append(_conflict(
            "snpedia_repute_vs_clinvar",
            [
                _position("snpedia", "Good",
                          "Hand-curated wiki judgement.",
                          "SNPedia editor"),
                _position("clinvar", raw_sig or "pathogenic",
                          "Pathogenic for the asserted condition.",
                          f"{stars} star review status"),
            ],
            "A curated human judgement and a clinical assertion disagree. "
            "Both are shown; neither is suppressed.",
        ))
    elif snp_repute == "Bad" and benign:
        conflicts.append(_conflict(
            "snpedia_repute_vs_clinvar",
            [
                _position("snpedia", "Bad",
                          "Hand-curated wiki judgement.",
                          "SNPedia editor"),
                _position("clinvar", raw_sig or "benign",
                          "Benign for the asserted condition.",
                          f"{stars} star review status"),
            ],
            "A curated human judgement and a clinical assertion disagree. "
            "Both are shown; neither is suppressed.",
        ))

    # 5. Pooled input files disagreeing at this position. merge.py has already
    #    kept both calls; this lifts them into the same conflict vocabulary so
    #    the UI has one shape to render rather than two.
    if finding.get("conflict"):
        calls = finding.get("calls") or []
        positions = [
            _position("dnainsight_input", str(call),
                      "Genotype call from one uploaded file.", "raw input file")
            for call in calls
        ] or [_position("dnainsight_input", "disagreement recorded",
                        "Uploaded files disagree at this position.",
                        "raw input file")]
        conflicts.append(_conflict(
            "input_files_disagree",
            positions,
            "Two uploaded files report different genotypes here. Both calls "
            "are retained. There is no voting and no automatic winner.",
        ))

    conflicts.sort(key=lambda c: c["type"])
    return conflicts


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

_DEFAULT_SCAN_PARAMETERS: dict[str, Any] = {
    "assembly": "GRCh37 / hg19",
    "magnitude_scale": "0 to 10, DNAInsight computed, not SNPedia magnitude",
    "network_enabled": False,
    "snpedia_enrichment": False,
    "carrier_aware": True,
    "palindromic_sites_capped": True,
}


def _finding_digest(finding: dict) -> str:
    """A stable digest of the clinically comparable part of one finding.

    Reuses the ledger's fingerprint so the manifest and the reclassification
    ledger agree on what "the same finding" means. Two modules with two
    different definitions of sameness would eventually disagree in front of a
    clinician.
    """
    try:
        from .ledger import fingerprint
        return str(fingerprint(finding).get("digest") or "")
    except Exception:
        return text_hash(_canonical(finding))


def build_manifest(*, profile_id, findings, input_files=None,
                   report_type: str = "", extra: dict | None = None) -> dict:
    """Assemble the reproducibility manifest for one report.

    Records the DNAInsight version, every database name plus version plus
    retrieval date plus licence, the sha256 of every input DNA file, the finding
    count, the scan parameters and a UTC timestamp.

    The point is that a clinician holding the report can either reproduce it
    exactly or prove it has drifted. Every field here exists because its absence
    would break one of those two.

    Input files whose hash cannot be computed are recorded with
    ``"sha256": None`` and ``"present": false`` rather than omitted. A silently
    missing input is indistinguishable from a report generated from nothing.
    """
    findings = list(findings or [])
    extra = extra if isinstance(extra, dict) else {}

    scan_parameters = dict(_DEFAULT_SCAN_PARAMETERS)
    if isinstance(extra.get("scan_parameters"), dict):
        scan_parameters.update(extra["scan_parameters"])

    versions = database_versions()
    retrieved = str(extra.get("retrieved") or _utc_date())

    databases = []
    for source_id in sorted(SOURCES):
        record = source_record(source_id,
                               version=versions.get(source_id, ""),
                               retrieved=retrieved)
        databases.append({
            "source_id": record["source_id"],
            "name": record["name"],
            "version": record["version"],
            "retrieved": record["retrieved"],
            "licence": record["licence"],
            "spdx": record["spdx"],
            "role": record["role"],
            "used": record["used"],
            "redistributable": record["redistributable"],
            "commercial_ok": record["commercial_ok"],
        })

    file_records = []
    for item in (input_files or []):
        path = item.get("path") if isinstance(item, dict) else item
        digest = content_hash(path)
        name = ""
        size = None
        try:
            resolved = Path(path)
            name = resolved.name
            if resolved.is_file():
                size = resolved.stat().st_size
        except (TypeError, OSError, ValueError):
            # An unusable path is recorded as unusable, never dropped.
            name = ""
        record = {
            "path": str(path),
            "name": name,
            "bytes": size,
            "sha256": digest,
            "present": digest is not None,
        }
        if isinstance(item, dict):
            for key in ("role", "provider", "label", "upload_id"):
                if item.get(key) is not None:
                    record[key] = item[key]
        file_records.append(record)

    entity_counts: dict[str, int] = {}
    silo_counts: dict[str, int] = {}
    for finding in findings:
        if not isinstance(finding, dict):
            continue
        entity = str(finding.get("entity_type") or "snp")
        entity_counts[entity] = entity_counts.get(entity, 0) + 1
        silo = str(finding.get("silo") or "")
        if silo:
            silo_counts[silo] = silo_counts.get(silo, 0) + 1

    digests = sorted(_finding_digest(f) for f in findings if isinstance(f, dict))

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "dnainsight_version": str(APP_VERSION),
        "generated_at": _utc_now(),
        "profile_id": profile_id,
        "report_type": str(report_type or ""),
        "finding_count": len(findings),
        "entity_counts": entity_counts,
        "silo_counts": silo_counts,
        "databases": databases,
        "database_versions": versions,
        "input_files": file_records,
        "scan_parameters": scan_parameters,
        "findings_digest": text_hash(_canonical(digests)),
        "licence_position": (
            "This report was produced from CC0 and US public domain data plus "
            "any opt-in local enrichment the user installed themselves. See "
            "data/DATA_SOURCES.md."
        ),
    }
    passthrough = {k: v for k, v in extra.items()
                   if k not in ("scan_parameters", "retrieved")}
    if passthrough:
        manifest["extra"] = passthrough
    return manifest


# ---------------------------------------------------------------------------
# Signing
# ---------------------------------------------------------------------------

_KEY_BYTES = 32


def _home_root() -> Path:
    """~/.dnainsight, honouring DNAINSIGHT_HOME exactly as external.py does."""
    override = os.environ.get("DNAINSIGHT_HOME")
    if override:
        return Path(override)
    return Path.home() / ".dnainsight"


def _key_path() -> Path:
    return _home_root() / "manifest.key"


def manifest_key() -> bytes:
    """The local HMAC key, generated on first use and reused thereafter.

    32 bytes from ``secrets.token_bytes``, written with 0600 permissions where
    the platform supports them, under ``~/.dnainsight/`` and therefore outside
    the repository tree, which is the same rule the SNPedia cache and the
    external tool registry follow.

    Created with O_CREAT plus O_EXCL so two processes starting at once cannot
    each write a key and leave one of them signing with a secret that has since
    been overwritten. The loser of that race re-reads the winner's key.

    An existing non-empty key file is used VERBATIM and is never regenerated,
    whatever its length. Regenerating would invalidate every manifest ever
    signed on this machine, turning "your report is intact" into "your report
    cannot be verified" with no way back. Only a zero-length file, which is what
    an interrupted first write leaves behind, is treated as absent.
    """
    path = _key_path()
    try:
        if path.exists() and path.stat().st_size > 0:
            return path.read_bytes()
    except OSError:
        pass

    path.parent.mkdir(parents=True, exist_ok=True)
    key = secrets.token_bytes(_KEY_BYTES)
    try:
        handle = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        # Either a concurrent writer won, or a zero-length file is in the way.
        try:
            existing = path.read_bytes()
        except OSError:
            existing = b""
        if existing:
            return existing
        with open(path, "wb") as fallback:
            fallback.write(key)
        _chmod_600(path)
        return key

    # os.fdopen takes ownership of the descriptor, so the with block closes it
    # on both the success and the failure path. No manual os.close here.
    with os.fdopen(handle, "wb") as fresh:
        fresh.write(key)
    _chmod_600(path)
    return key


def _chmod_600(path: Path) -> None:
    """Best effort 0600. Windows and some network filesystems do not implement
    POSIX modes, and a failure there must not stop a report being signed."""
    try:
        os.chmod(str(path), 0o600)
    except (OSError, NotImplementedError):
        pass


def _signing_payload(manifest: dict, signed_at: str, algorithm: str) -> bytes:
    """The exact bytes that get signed.

    ``signed_at`` and ``algorithm`` are inside the signed payload, not beside
    it. A signature that covers only the manifest body would let anyone rewrite
    the signing timestamp or downgrade the stated algorithm while the signature
    still verified, which is a real and well known way to make a signature
    meaningless.
    """
    return _canonical({
        "algorithm": algorithm,
        "signed_at": signed_at,
        "manifest": manifest,
    }).encode("utf-8")


def _field_digests(manifest: dict) -> dict:
    """Per top-level-field digests, so verification can NAME what drifted.

    Without these, a failed verification can only say "something changed",
    which sends a clinician hunting through a 200 line manifest by eye.
    """
    return {str(key): text_hash(_canonical(value))
            for key, value in manifest.items()}


def sign_manifest(manifest: dict) -> dict:
    """Sign a manifest with the local HMAC key.

    Returns ``{"manifest": ..., "signature": {...}}``. The manifest is embedded
    unchanged, not rewritten, so the signed object still reads as the report it
    describes.

    Read ``signature["scope"]``: this proves the report was not altered after
    generation ON THIS MACHINE. It is not a public-key attestation and it does
    not prove authorship to a third party.
    """
    body = manifest if isinstance(manifest, dict) else {}
    signed_at = _utc_now()
    algorithm = "HMAC-SHA256"
    payload = _signing_payload(body, signed_at, algorithm)
    value = hmac.new(manifest_key(), payload, hashlib.sha256).hexdigest()

    return {
        "manifest": body,
        "signature": {
            "algorithm": algorithm,
            "signed_at": signed_at,
            "value": value,
            "payload_sha256": hashlib.sha256(payload).hexdigest(),
            "manifest_sha256": text_hash(_canonical(body)),
            "field_digests": _field_digests(body),
            "scope": SIGNATURE_SCOPE,
        },
    }


def _verdict(ok: bool, reason: str, field: str, detail: str, **rest) -> dict:
    verdict = {
        "ok": bool(ok),
        "reason": reason,
        "field": field,
        "detail": detail,
        "checked_at": _utc_now(),
        "changed_fields": [],
        "algorithm": "",
        "signed_at": "",
    }
    verdict.update(rest)
    return verdict


def verify_manifest(signed: dict) -> dict:
    """Verify a signed manifest and return a structured verdict.

    NEVER a bare False. The verdict always names the field that failed, because
    "verification failed" tells a clinician nothing they can act on, while
    "manifest.finding_count was altered" tells them exactly what to distrust.

    An unsigned payload is a structured failure, not an exception. Being handed
    a plain manifest instead of a signed one is a routine mistake, and raising
    would turn it into a 500 in the route layer.

    Failure fields, in the order they are checked:
      signed                 not a mapping at all
      manifest               missing, not a mapping, or altered after signing
      manifest.<name>        the specific top-level field that drifted
      signature              missing or not a mapping (the unsigned case)
      signature.algorithm    absent or not one this build supports
      signature.signed_at    absent
      signature.value        absent, malformed, or does not match
    """
    if not isinstance(signed, dict):
        return _verdict(False, "not_a_mapping", "signed",
                        "Expected a mapping of manifest plus signature.")

    manifest = signed.get("manifest")
    if manifest is None:
        return _verdict(False, "missing_manifest", "manifest",
                        "The payload carries no manifest.")
    if not isinstance(manifest, dict):
        return _verdict(False, "manifest_not_a_mapping", "manifest",
                        "The manifest is not a mapping.")

    signature = signed.get("signature")
    if signature is None:
        return _verdict(False, "unsigned", "signature",
                        "The payload carries no signature block. It was never "
                        "signed, or the signature was stripped.")
    if not isinstance(signature, dict):
        return _verdict(False, "signature_not_a_mapping", "signature",
                        "The signature block is not a mapping.")

    algorithm = str(signature.get("algorithm") or "")
    if not algorithm:
        return _verdict(False, "missing_algorithm", "signature.algorithm",
                        "The signature block names no algorithm.")
    if algorithm != "HMAC-SHA256":
        return _verdict(False, "unsupported_algorithm", "signature.algorithm",
                        f"This build verifies HMAC-SHA256 only, not {algorithm}.",
                        algorithm=algorithm)

    signed_at = str(signature.get("signed_at") or "")
    if not signed_at:
        return _verdict(False, "missing_signed_at", "signature.signed_at",
                        "The signature block carries no signing timestamp, so "
                        "the timestamp cannot have been covered by the "
                        "signature.")

    value = str(signature.get("value") or "")
    if not value:
        return _verdict(False, "missing_signature_value", "signature.value",
                        "The signature block carries no value.",
                        algorithm=algorithm, signed_at=signed_at)

    # Name the drifted field before checking the HMAC. The HMAC only ever says
    # "something is wrong"; the per-field digests say what.
    recorded = signature.get("field_digests")
    changed: list[str] = []
    if isinstance(recorded, dict) and recorded:
        current = _field_digests(manifest)
        for key in sorted(set(recorded) | set(current)):
            if recorded.get(key) != current.get(key):
                changed.append(key)
        if changed:
            field = f"manifest.{changed[0]}" if len(changed) == 1 else "manifest"
            return _verdict(
                False, "manifest_modified", field,
                "The manifest was altered after signing. Altered field(s): "
                + ", ".join(changed) + ".",
                changed_fields=changed, algorithm=algorithm, signed_at=signed_at)

    payload = _signing_payload(manifest, signed_at, algorithm)
    expected = hmac.new(manifest_key(), payload, hashlib.sha256).hexdigest()

    # Constant time, so a caller cannot learn the correct signature one byte at
    # a time by measuring how long a rejection takes.
    if not hmac.compare_digest(expected, value):
        return _verdict(
            False, "signature_mismatch", "signature.value",
            "The signature does not match this manifest under the local key. "
            "Either the payload was altered or it was signed on another "
            "machine.",
            algorithm=algorithm, signed_at=signed_at)

    recorded_payload = str(signature.get("payload_sha256") or "")
    if recorded_payload and recorded_payload != hashlib.sha256(payload).hexdigest():
        return _verdict(False, "payload_digest_mismatch",
                        "signature.payload_sha256",
                        "The recorded payload digest does not match the "
                        "reconstructed payload.",
                        algorithm=algorithm, signed_at=signed_at)

    return _verdict(True, "verified", "",
                    "The manifest is byte for byte what was signed on this "
                    "machine.",
                    algorithm=algorithm, signed_at=signed_at)


# ---------------------------------------------------------------------------
# Human readable manifest
# ---------------------------------------------------------------------------

def render_manifest_text(signed: dict) -> str:
    """Render a signed manifest as plain text for the end of a report.

    Plain text on purpose. This block has to survive being printed, faxed,
    pasted into an electronic health record and retyped, which is what actually
    happens to clinical documents, and none of those preserve markup.
    """
    payload = signed if isinstance(signed, dict) else {}
    manifest = payload.get("manifest") if isinstance(payload.get("manifest"), dict) else {}
    signature = payload.get("signature") if isinstance(payload.get("signature"), dict) else {}

    lines: list[str] = []
    lines.append("DNAINSIGHT REPORT MANIFEST")
    lines.append("=" * 60)
    lines.append(f"Manifest version : {manifest.get('manifest_version', '')}")
    lines.append(f"DNAInsight       : {manifest.get('dnainsight_version', '')}")
    lines.append(f"Generated (UTC)  : {manifest.get('generated_at', '')}")
    lines.append(f"Profile          : {manifest.get('profile_id', '')}")
    lines.append(f"Report type      : {manifest.get('report_type', '')}")
    lines.append(f"Findings         : {manifest.get('finding_count', 0)}")
    lines.append(f"Findings digest  : {manifest.get('findings_digest', '')}")
    lines.append("")

    lines.append("DATABASES IN FORCE")
    lines.append("-" * 60)
    for entry in manifest.get("databases", []) or []:
        if not isinstance(entry, dict):
            continue
        status = "used" if entry.get("used") else "NOT USED"
        lines.append(f"  {entry.get('name', '')} [{status}]")
        lines.append(f"      version   : {entry.get('version', '') or 'not recorded'}")
        lines.append(f"      retrieved : {entry.get('retrieved', '') or 'not recorded'}")
        lines.append(f"      spdx      : {entry.get('spdx', '') or 'none published'}")
        lines.append(f"      licence   : {entry.get('licence', '')}")
    lines.append("")

    lines.append("INPUT FILES")
    lines.append("-" * 60)
    files = manifest.get("input_files", []) or []
    if not files:
        lines.append("  none recorded")
    for entry in files:
        if not isinstance(entry, dict):
            continue
        digest = entry.get("sha256") or "NOT HASHED, file absent or unreadable"
        lines.append(f"  {entry.get('name') or entry.get('path', '')}")
        lines.append(f"      sha256 : {digest}")
        if entry.get("bytes") is not None:
            lines.append(f"      bytes  : {entry.get('bytes')}")
    lines.append("")

    lines.append("SCAN PARAMETERS")
    lines.append("-" * 60)
    params = manifest.get("scan_parameters", {}) or {}
    for key in sorted(params):
        lines.append(f"  {key} = {params[key]}")
    lines.append("")

    lines.append("SIGNATURE")
    lines.append("-" * 60)
    lines.append(f"  algorithm : {signature.get('algorithm', '') or 'UNSIGNED'}")
    lines.append(f"  signed at : {signature.get('signed_at', '') or 'never'}")
    lines.append(f"  value     : {signature.get('value', '') or 'none'}")
    lines.append("")
    lines.append("  What this signature does and does not prove:")
    for chunk in _wrap(SIGNATURE_SCOPE, 66):
        lines.append(f"  {chunk}")
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    """Wrap without importing textwrap for one call. Words longer than the
    width are left intact rather than broken, because a hash split across two
    lines cannot be retyped."""
    words = str(text or "").split()
    lines: list[str] = []
    current = ""
    for word in words:
        if current and len(current) + 1 + len(word) > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(current)
    return lines


# ---------------------------------------------------------------------------
# Licence audit
# ---------------------------------------------------------------------------

def licence_audit() -> dict:
    """Check every declared bundled artefact against the licence rules.

    This is what turns ``data/DATA_SOURCES.md`` from a static document into an
    enforced runtime contract. The rules, each one taken straight from that
    document:

      1. A source marked ``used: False`` must not appear in any bundled
         artefact. PharmGKB is the live case, and its exclusion is the whole
         reason it has an entry at all.
      2. A source marked ``never_bundle`` must not appear in any bundled
         artefact. SNPedia would relicense the repository away from MIT through
         share-alike; MyVariant.info carries per-field licences that cannot be
         asserted over a persisted copy.
      3. A source that is not redistributable must not appear in any bundled
         artefact.
      4. A source id that is not in ``SOURCES`` at all is a violation, not a
         shrug. An unassessed source is exactly how contamination enters.

    Per-record licensing produces a WARNING rather than a violation, because
    the PGS Catalog is legitimately bundled after filtering. The warning exists
    so nobody forgets that the filtering is load bearing.

    Returns ``{"ok": bool, "violations": [...], "warnings": [...], ...}``.
    ``ok`` is False whenever there is at least one violation.
    """
    violations: list[dict] = []
    warnings: list[dict] = []

    for artefact in sorted(BUNDLED_ARTEFACTS):
        for source_id in BUNDLED_ARTEFACTS[artefact]:
            entry = SOURCES.get(source_id)

            if entry is None:
                violations.append({
                    "artefact": artefact,
                    "source_id": source_id,
                    "name": source_id,
                    "licence": "unknown, not assessed",
                    "spdx": "",
                    "rule": "unknown_source",
                    "reason": (
                        f"'{source_id}' is bundled in {artefact} but is not "
                        "recorded in data/DATA_SOURCES.md. An unassessed source "
                        "must never be shipped."
                    ),
                })
                continue

            base = {
                "artefact": artefact,
                "source_id": source_id,
                "name": entry["name"],
                "licence": entry["licence"],
                "spdx": entry["spdx"],
            }

            if not entry.get("used", True):
                violations.append({**base, "rule": "not_used_source_bundled",
                                   "reason": (
                                       f"{entry['name']} is recorded as NOT USED and "
                                       f"appears in the bundled artefact {artefact}.")})
            if entry.get("never_bundle"):
                violations.append({**base, "rule": "never_bundle",
                                   "reason": (
                                       f"{entry['name']} must never be bundled and "
                                       f"appears in {artefact}. {entry.get('note', '')}")})
            elif not entry.get("redistributable", False):
                violations.append({**base, "rule": "not_redistributable",
                                   "reason": (
                                       f"{entry['name']} is not redistributable and "
                                       f"appears in the bundled artefact {artefact}.")})
            if entry.get("per_record"):
                warnings.append({**base, "rule": "per_record_licence",
                                 "reason": (
                                     f"{entry['name']} licences per record, so "
                                     f"{artefact} is only clean if the "
                                     "per-record filter still runs. "
                                     f"{entry.get('note', '')}")})

    not_used = sorted(k for k, v in SOURCES.items() if not v.get("used", True))
    never_bundle = sorted(k for k, v in SOURCES.items() if v.get("never_bundle"))

    return {
        "ok": not violations,
        "generated_at": _utc_now(),
        "artefacts": len(BUNDLED_ARTEFACTS),
        "checked": sum(len(v) for v in BUNDLED_ARTEFACTS.values()),
        "violations": violations,
        "warnings": warnings,
        "sources": {
            "total": len(SOURCES),
            "not_used": not_used,
            "never_bundle": never_bundle,
            "redistributable": sorted(k for k, v in SOURCES.items()
                                      if v.get("redistributable")),
        },
    }
