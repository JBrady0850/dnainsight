"""
scanner.py -- SNP annotation engine for DNAInsight.

Annotation priority:
  1. Bundled reference database (data/snp_reference.json) -- instant, offline
  2. MyVariant.info API (hg19/GRCh37) -- for SNPs not in bundled set
  3. Local ClinVar SQLite DB (db/clinvar_local.db) -- if user ran the DB builder

Clinical significance filter (KEEP_SIG):
  pathogenic, likely pathogenic, drug response, risk factor, protective

Returns structured finding dicts consumed by database.upsert_finding().
"""

import json
import os
import re
import sqlite3
import urllib.request
import urllib.parse
import concurrent.futures
import tempfile
from pathlib import Path
from typing import Generator


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
API_URL      = "https://myvariant.info/v1/variant"
API_FIELDS   = "clinvar,pharmgkb,dbsnp.gene,dbsnp.ref,dbsnp.alt"
API_ASSEMBLY = "hg19"
API_TIMEOUT  = 30
CHUNK_SIZE   = 200
MAX_WORKERS  = 4

KEEP_SIG = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
    "drug response",
    "risk factor",
    "protective",
    "pathogenic, low penetrance",
    "likely pathogenic, low penetrance",
}

_BASE = Path(__file__).parent.parent
BUNDLED_REF  = _BASE / "data" / "snp_reference.json"
LOCAL_DB     = _BASE / "db" / "clinvar_local.db"


# ---------------------------------------------------------------------------
# Risk text classifier (adapted from batch_snp_scan.py)
# ---------------------------------------------------------------------------
_RISK_FRAGS = [
    'no function', 'decreased function', 'reduced function', 'poor function',
    'null function', 'no activity', 'loss of function',
    'increased risk', 'higher risk', 'elevated risk',
    'increased toxicity', 'severe toxicity', 'toxicity risk', 'toxicity',
    'contraindicated', 'avoid', 'do not use',
    'poor metabolizer', 'ultrarapid', 'ultra-rapid',
    'malignant hyperthermia', 'adverse reaction', 'no response',
    'decreased response', 'reduced efficacy', 'treatment failure',
    'pathogenic', 'disease causing', 'susceptibility',
    'associated with increased',
]
_SAFE_FRAGS = [
    'normal function', 'wildtype', 'standard response', 'no increased',
    'not at increased', 'lower risk', 'reduced risk', 'decreased risk',
    'protective', 'no significant', 'no effect', 'standard dose',
    'normal metabolizer',
]

def _conditions_risk(conditions: str, genotype: str) -> str:
    """Return 'risk', 'safe', or 'unknown'."""
    lc = conditions.lower()
    safe_score = sum(1 for f in _SAFE_FRAGS if f in lc)
    risk_score = sum(1 for f in _RISK_FRAGS if f in lc)
    if risk_score > safe_score:
        return "risk"
    if safe_score > risk_score:
        return "safe"
    return "unknown"


# ---------------------------------------------------------------------------
# Silo classification
# ---------------------------------------------------------------------------
PRE_RX_TERMS = [
    'contraindicated', 'avoid', 'dose adjustment', 'dose reduction',
    'poor metabolizer', 'ultrarapid', 'slow metabolizer', 'reduced metabolism',
    'myopathy risk', 'warfarin', 'fluorouracil', 'thiopurine', 'azathioprine',
    'fentanyl', 'codeine', 'tramadol', 'ssri', 'antidepressant', 'antipsychotic',
    'opioid', 'tacrolimus', 'cyclosporine', 'methotrexate', 'statin',
    'beta-blocker', 'phenytoin', 'carbamazepine',
]

def classify_silo(interp: str) -> str:
    il = interp.lower()
    if any(t in il for t in PRE_RX_TERMS):
        return "pre_prescription"
    if any(t in il for t in ['risk', 'susceptibility', 'elevated', 'increased risk']):
        return "actionable"
    return "informational"


# ---------------------------------------------------------------------------
# Bundled reference lookup
# ---------------------------------------------------------------------------
_bundled_cache: dict = {}

_bundled_meta: dict = {}

def _load_bundled() -> dict:
    global _bundled_cache, _bundled_meta
    if _bundled_cache:
        return _bundled_cache
    if BUNDLED_REF.exists():
        with open(BUNDLED_REF, "r", encoding="utf-8") as f:
            raw = json.load(f)
        # Support both old flat format and new versioned format
        if "snps" in raw and "_meta" in raw:
            _bundled_meta = raw["_meta"]
            _bundled_cache = raw["snps"]
        else:
            _bundled_cache = raw
    return _bundled_cache


def get_reference_metadata() -> dict:
    """Return _meta block from snp_reference.json (version, snp_count, built_at)."""
    _load_bundled()
    return _bundled_meta


def lookup_bundled(rsid: str) -> dict | None:
    ref = _load_bundled()
    return ref.get(rsid)


# ---------------------------------------------------------------------------
# Local ClinVar DB lookup
# ---------------------------------------------------------------------------

def lookup_local_db(rsid: str) -> list[dict]:
    if not LOCAL_DB.exists():
        return []
    try:
        conn = sqlite3.connect(str(LOCAL_DB))
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            "SELECT * FROM clinvar WHERE rsid=? LIMIT 10", (rsid,)
        ).fetchall()
        conn.close()
        return [dict(r) for r in rows]
    except Exception:
        return []


# ---------------------------------------------------------------------------
# MyVariant.info API
# ---------------------------------------------------------------------------

def _build_rsid_list(rsids: list[str]) -> str:
    return ",".join(rsids)


def _api_fetch_chunk(rsids: list[str]) -> dict:
    """POST a batch of rsIDs to MyVariant.info. Returns raw API response dict."""
    body = urllib.parse.urlencode({
        "ids":    _build_rsid_list(rsids),
        "fields": API_FIELDS,
        "assembly": API_ASSEMBLY,
    }).encode("utf-8")
    req = urllib.request.Request(
        API_URL,
        data=body,
        method="POST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    try:
        with urllib.request.urlopen(req, timeout=API_TIMEOUT) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return {}


def _extract_gene(hit: dict) -> str:
    try:
        g = hit.get("dbsnp", {}).get("gene", {})
        if isinstance(g, list):
            return g[0].get("symbol", "")
        return g.get("symbol", "")
    except Exception:
        return ""


def _extract_clinvar(hit: dict, genotype: str) -> list[dict]:
    """Extract ClinVar annotations from a MyVariant.info hit."""
    results = []
    cv = hit.get("clinvar", {})
    if not cv:
        return results

    rcvs = cv.get("rcv", [])
    if isinstance(rcvs, dict):
        rcvs = [rcvs]

    for rcv in rcvs:
        sig = rcv.get("clinical_significance", "").lower()
        if sig not in KEEP_SIG:
            continue
        conditions = ""
        cond = rcv.get("conditions", {})
        if isinstance(cond, dict):
            conditions = cond.get("name", "")
        elif isinstance(cond, list):
            conditions = "; ".join(c.get("name", "") for c in cond if isinstance(c, dict))

        results.append({
            "clinical_sig": sig,
            "conditions":   conditions,
            "source":       "clinvar",
        })
    return results


def _extract_pharmgkb(hit: dict) -> list[dict]:
    results = []
    pgkb = hit.get("pharmgkb", {})
    if not pgkb:
        return results
    anns = pgkb.get("variant_annotation", [])
    if isinstance(anns, dict):
        anns = [anns]
    for ann in anns:
        chemicals = ann.get("chemicals", {})
        if isinstance(chemicals, dict):
            chemicals = [chemicals]
        for chem in (chemicals or []):
            results.append({
                "clinical_sig": "drug response",
                "conditions":   ann.get("description", ""),
                "source":       "pharmgkb",
                "drug":         chem.get("name", ""),
            })
    return results


def annotate_via_api(snps: list[dict], progress_cb=None) -> list[dict]:
    """
    Annotate a list of SNP dicts via MyVariant.info API.
    Returns findings list.
    """
    # Index snps by rsid for fast genotype lookup
    snp_index = {s["rsid"]: s for s in snps if s["rsid"].startswith("rs")}
    rsids = list(snp_index.keys())
    chunks = [rsids[i:i+CHUNK_SIZE] for i in range(0, len(rsids), CHUNK_SIZE)]

    findings = []
    processed = 0

    def fetch_chunk(chunk):
        return _api_fetch_chunk(chunk)

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(fetch_chunk, c): c for c in chunks}
        for fut in concurrent.futures.as_completed(futures):
            chunk = futures[fut]
            try:
                data = fut.result()
            except Exception:
                data = {}

            hits = data if isinstance(data, list) else []
            for hit in hits:
                rsid = hit.get("_id", "").split(":")[0].lower()
                if not rsid.startswith("rs"):
                    continue
                snp = snp_index.get(rsid, {})
                genotype = snp.get("allele1", "") + snp.get("allele2", "")
                gene     = _extract_gene(hit)

                anns = _extract_clinvar(hit, genotype) + _extract_pharmgkb(hit)
                for ann in anns:
                    conditions = ann.get("conditions", "")
                    interp     = conditions[:300] if conditions else ann.get("drug", "")
                    findings.append({
                        "rsid":         rsid,
                        "gene":         gene,
                        "chromosome":   snp.get("chromosome", ""),
                        "position":     snp.get("position", 0),
                        "allele1":      snp.get("allele1", ""),
                        "allele2":      snp.get("allele2", ""),
                        "genotype":     genotype,
                        "clinical_sig": ann.get("clinical_sig", ""),
                        "conditions":   conditions,
                        "interpretation": interp,
                        "category":     ann.get("source", ""),
                        "silo":         classify_silo(interp),
                        "sources":      [ann.get("source", "")],
                    })

            processed += len(chunk)
            if progress_cb:
                progress_cb(processed, len(rsids))

    return findings


# ---------------------------------------------------------------------------
# Bundled reference scan (offline-first)
# ---------------------------------------------------------------------------

def annotate_bundled(snps: list[dict]) -> list[dict]:
    """
    Annotate SNPs using the bundled reference database only (no network).
    Returns findings list.
    """
    ref     = _load_bundled()
    snp_idx = {s["rsid"]: s for s in snps}
    findings = []

    for rsid, entry in ref.items():
        snp = snp_idx.get(rsid)
        if not snp:
            continue

        genotype = snp["allele1"] + snp["allele2"]
        interp   = entry.get("interpretation", entry.get("description", ""))
        gene     = entry.get("gene", "")
        category = entry.get("category", "")

        findings.append({
            "rsid":         rsid,
            "gene":         gene,
            "chromosome":   snp.get("chromosome", ""),
            "position":     snp.get("position", 0),
            "allele1":      snp["allele1"],
            "allele2":      snp["allele2"],
            "genotype":     genotype,
            "clinical_sig": entry.get("clinical_sig", "drug response"),
            "conditions":   entry.get("conditions", ""),
            "interpretation": interp,
            "category":     category,
            "silo":         classify_silo(interp),
            "sources":      ["bundled_reference"],
        })

    return findings


# ---------------------------------------------------------------------------
# Full scan pipeline
# ---------------------------------------------------------------------------

def run_scan(snps: list[dict], use_api: bool = True, progress_cb=None) -> list[dict]:
    """
    Run full annotation pipeline.
    1. Bundled reference (offline, instant)
    2. API for remaining SNPs (if use_api=True)

    Returns deduplicated findings list.
    """
    bundled_findings = annotate_bundled(snps)
    found_rsids      = {f["rsid"] for f in bundled_findings}

    all_findings = list(bundled_findings)

    if use_api:
        remaining_snps = [s for s in snps if s["rsid"] not in found_rsids]
        if remaining_snps:
            api_findings = annotate_via_api(remaining_snps, progress_cb=progress_cb)
            all_findings.extend(api_findings)

    # Deduplicate by rsid (keep highest-priority silo)
    silo_rank = {"pre_prescription": 0, "actionable": 1, "informational": 2}
    deduped   = {}
    for f in all_findings:
        rsid = f["rsid"]
        if rsid not in deduped:
            deduped[rsid] = f
        else:
            existing_rank = silo_rank.get(deduped[rsid]["silo"], 3)
            new_rank      = silo_rank.get(f["silo"], 3)
            if new_rank < existing_rank:
                deduped[rsid] = f

    return list(deduped.values())
