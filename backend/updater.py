"""
updater.py -- Database refresh module for DNAInsight.

Re-fetches annotations from MyVariant.info for all bundled SNPs and
updates snp_reference.json in-place. Preserves the static interpretation
and category fields from build_reference.py; only overwrites clinical_sig
and conditions with fresh ClinVar data when the API returns a valid result.

Thread-safe: designed to run in a daemon thread; progress is reported via
an optional callback.
"""

import json
import requests
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent.parent / "data"
REFERENCE_PATH = DATA_DIR / "snp_reference.json"

# ---------------------------------------------------------------------------
# API config (mirrors scanner.py so they stay in sync)
# ---------------------------------------------------------------------------

MYVARIANT_URL = "https://myvariant.info/v1/variant"
CHUNK_SIZE    = 200
MAX_WORKERS   = 4
REQUEST_TIMEOUT = 30

KEEP_SIG = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "drug response",
    "risk factor",
    "protective",
    "pathogenic, low penetrance",
    "pathogenic low penetrance",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_reference() -> dict:
    """Load snp_reference.json. Returns the full raw dict (may include _meta + snps)."""
    with open(REFERENCE_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def _get_snps_dict(ref: dict) -> dict:
    """Return the flat rsID->entry dict regardless of file format (old flat or new nested)."""
    if "snps" in ref and "_meta" in ref:
        return ref["snps"]
    return {k: v for k, v in ref.items() if k.startswith("rs")}


def _save_reference(snps: dict, meta: dict):
    """Write snp_reference.json in the current versioned format."""
    payload = {
        "_meta": meta,
        "snps":  snps,
    }
    with open(REFERENCE_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def _extract_clinical_sig(hit: dict) -> str | None:
    """Pull the best ClinVar clinical significance string from a MyVariant hit."""
    clinvar = hit.get("clinvar", {})
    if not clinvar:
        return None

    # RCV list
    rcv = clinvar.get("rcv", [])
    if isinstance(rcv, dict):
        rcv = [rcv]
    for r in rcv:
        sig = (r.get("clinical_significance") or "").strip().lower()
        if sig in KEEP_SIG:
            return sig

    # Top-level clinical_significance dict
    sig_obj = clinvar.get("clinical_significance", {})
    if isinstance(sig_obj, dict):
        desc = (sig_obj.get("description") or "").strip().lower()
        if desc in KEEP_SIG:
            return desc
    elif isinstance(sig_obj, str):
        desc = sig_obj.strip().lower()
        if desc in KEEP_SIG:
            return desc

    return None


def _extract_conditions(hit: dict) -> str | None:
    """Extract disease/condition names from a MyVariant hit."""
    clinvar = hit.get("clinvar", {})
    rcv = clinvar.get("rcv", [])
    if isinstance(rcv, dict):
        rcv = [rcv]

    conditions = []
    for r in rcv:
        cond = r.get("conditions", {})
        if isinstance(cond, dict):
            name = cond.get("name", "")
            if name:
                conditions.append(name)
        elif isinstance(cond, list):
            for c in cond:
                if isinstance(c, dict):
                    name = c.get("name", "")
                    if name:
                        conditions.append(name)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for c in conditions:
        if c not in seen:
            seen.add(c)
            unique.append(c)

    return "; ".join(unique) if unique else None


def _fetch_chunk(rsids: list) -> dict:
    """
    POST one chunk of rsIDs to MyVariant.info.
    Returns {rsid: hit_dict} for successful hits only.
    On network error returns an empty dict (caller increments error counter).
    """
    try:
        payload = {
            "ids":      ",".join(rsids),
            "fields":   "clinvar,pharmgkb,dbsnp.gene",
            "assembly": "hg19",
        }
        resp = requests.post(MYVARIANT_URL, data=payload, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        hits = resp.json()
        result = {}
        for hit in hits:
            if isinstance(hit, dict) and not hit.get("notfound"):
                rsid = hit.get("_id", "")
                if rsid:
                    result[rsid] = hit
        return result
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def update_bundled_reference(progress_cb=None) -> dict:
    """
    Re-fetch all bundled rsIDs from MyVariant.info and update snp_reference.json.

    Strategy:
      - Only rsIDs (keys starting with "rs") are sent to the API.
      - The _meta and any non-rsid keys are preserved unchanged.
      - For each rsID hit: clinical_sig and conditions are updated if the API
        returns a KEEP_SIG-level result; otherwise the existing values are kept.
      - interpretation, gene, category are NEVER overwritten (they come from
        the curated build_reference.py source).

    Args:
        progress_cb: optional callable(processed_int, total_int)

    Returns:
        {success, updated, skipped, errors, timestamp, snp_count}
    """
    raw   = _load_reference()
    ref   = _get_snps_dict(raw)
    rsids = list(ref.keys())
    total = len(rsids)

    if total == 0:
        return {
            "success": False,
            "error":   "No rsIDs found in reference file.",
            "updated": 0, "skipped": 0, "errors": 0,
            "snp_count": 0,
        }

    if progress_cb:
        progress_cb(0, total)

    # Build chunks
    chunks = [rsids[i:i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)]

    updated   = 0
    skipped   = 0
    errors    = 0
    processed = 0
    all_hits: dict = {}

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_chunk = {executor.submit(_fetch_chunk, chunk): chunk for chunk in chunks}
        for future in as_completed(future_to_chunk):
            chunk = future_to_chunk[future]
            hits  = future.result()
            if not hits and chunk:
                # Network failure for this chunk
                errors += len(chunk)
            all_hits.update(hits)
            processed += len(chunk)
            if progress_cb:
                progress_cb(min(processed, total), total)

    # Merge fresh annotations into reference
    for rsid in rsids:
        entry = dict(ref.get(rsid, {}))
        hit   = all_hits.get(rsid)

        if not hit:
            skipped += 1
            ref[rsid] = entry
            continue

        fresh_sig  = _extract_clinical_sig(hit)
        fresh_cond = _extract_conditions(hit)

        if fresh_sig:
            entry["clinical_sig"] = fresh_sig
            updated += 1
        else:
            skipped += 1

        if fresh_cond:
            entry["conditions"] = fresh_cond

        ref[rsid] = entry

    # Build updated _meta block, preserving any existing fields
    now  = datetime.now(timezone.utc).isoformat()
    meta = raw.get("_meta", {})
    meta.update({
        "updated_at": now,
        "snp_count":  len(rsids),
        "source":     "MyVariant.info (ClinVar + PharmGKB)",
        "updated":    updated,
        "skipped":    skipped,
        "errors":     errors,
    })

    _save_reference(ref, meta)

    return {
        "success":   True,
        "updated":   updated,
        "skipped":   skipped,
        "errors":    errors,
        "timestamp": now,
        "snp_count": len(rsids),
    }


def get_reference_metadata() -> dict:
    """
    Return the current _meta block from snp_reference.json.
    Falls back gracefully if the file has no _meta yet (pre-first-update).
    """
    try:
        raw   = _load_reference()
        meta  = raw.get("_meta", {})
        snps  = _get_snps_dict(raw)
        return {
            "updated_at": meta.get("updated_at"),
            "built_at":   meta.get("built_at"),
            "version":    meta.get("version", "1.1.0"),
            "snp_count":  meta.get("snp_count", len(snps)),
            "source":     meta.get("source", "Bundled (offline)"),
            "updated":    meta.get("updated", 0),
            "skipped":    meta.get("skipped", 0),
            "errors":     meta.get("errors", 0),
        }
    except FileNotFoundError:
        return {"error": "snp_reference.json not found. Run install.bat to rebuild."}
    except Exception as e:
        return {"error": str(e)}
