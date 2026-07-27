"""
routes_v2.py -- the v2 HTTP surface, as a separate blueprint.

Deliberately separate from routes.py so that every v1.2 endpoint keeps its exact
behaviour. The v1 blueprint is registered first and is untouched; this one adds
the new paths and re-serves /findings with the full filter engine.

Implements docs/API_V2.md. Scan results are cached per profile in memory and
persisted to the findings table, so filtering does not require a rescan.
"""

from __future__ import annotations

import io
import csv
import json
import threading
from datetime import datetime, timezone
from pathlib import Path

from flask import Blueprint, request, jsonify, send_file
from werkzeug.utils import secure_filename

from . import APP_VERSION
from . import database as db
from . import filters as F
from . import pipeline
from . import frequency as freq_mod
from . import scoring as scoring_mod
from .parsers import parse_dna_file, ParseError

api_v2 = Blueprint("api_v2", __name__)

BASE = Path(__file__).parent.parent
UPLOAD_DIR = BASE / "uploads"
REPORTS_DIR = BASE / "reports_output"
UPLOAD_DIR.mkdir(exist_ok=True)
REPORTS_DIR.mkdir(exist_ok=True)

ALLOWED_EXTENSIONS = {".txt", ".csv", ".tsv"}
VALID_ROLES = ("self", "mother", "father", "mate", "child", "sibling",
               "other", "ignore")

# Last full scan payload per profile. Findings are also written to the database,
# but the richer payload (genosets, traits, prs, qc, conflicts) lives here so a
# filter request never has to rebuild it.
_scan_cache: dict[int, dict] = {}
_scan_progress: dict[int, dict] = {}
_lock = threading.Lock()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _err(message: str, code: int = 400):
    return jsonify({"error": message}), code


def _sources_path(pid: int) -> Path:
    return UPLOAD_DIR / f"sources_{pid}.json"


def _load_sources(pid: int) -> list[dict]:
    path = _sources_path(pid)
    if not path.exists():
        # Fall back to the v1 single-file cache so existing profiles still work.
        legacy = UPLOAD_DIR / f"snps_{pid}.json"
        if legacy.exists():
            try:
                with open(legacy, "r", encoding="utf-8") as fh:
                    snps = json.load(fh)
                profile = db.get_profile(pid) or {}
                return [{
                    "id": 1, "label": f"{profile.get('name', 'profile')} upload",
                    "role": "self", "provider": profile.get("provider", ""),
                    "snps": snps, "uploaded_at": profile.get("created_at", ""),
                }]
            except (OSError, ValueError):
                return []
        return []
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return []


def _save_sources(pid: int, sources: list[dict]) -> None:
    with open(_sources_path(pid), "w", encoding="utf-8") as fh:
        json.dump(sources, fh)


def _cached(pid: int) -> dict | None:
    with _lock:
        return _scan_cache.get(pid)


def _unique_report_name(safe_name: str, kind: str) -> str:
    """Return a report filename that cannot collide with an existing one.

    A second-granular timestamp is not enough. Generating four reports inside
    one second, which happens routinely when a user clicks through the report
    buttons or when a test exercises them back to back, produced four database
    rows all pointing at the SAME path. Each write clobbered the last, so
    opening report 1 served report 4's content. The report ids looked fine and
    nothing errored, which is what made it worth guarding against explicitly.
    """
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    base = f"{safe_name}_{kind}_report_{stamp}"
    candidate = f"{base}.html"
    counter = 2
    while (REPORTS_DIR / candidate).exists():
        candidate = f"{base}_{counter}.html"
        counter += 1
    return candidate


# ---------------------------------------------------------------------------
# System
# ---------------------------------------------------------------------------

@api_v2.route("/api/capabilities", methods=["GET"])
def capabilities():
    """Report which subsystems have data, so the UI can hide dead controls."""
    subsystems = pipeline.available_subsystems()
    subsystems["api"] = True
    return jsonify(subsystems)


@api_v2.route("/api/populations", methods=["GET"])
def populations():
    """Every reference population, marked with whether data exists for it."""
    available = {p["code"] for p in freq_mod.available_populations()}
    rows = [{
        "code": p["code"], "label": p["label"], "brief": p["brief"],
        "superpop": p.get("superpop", ""),
        "available": p["code"] in available,
    } for p in freq_mod.POPULATIONS]
    return jsonify({
        "populations": rows,
        "default": freq_mod.DEFAULT_POPULATION,
        "aggregate_modes": list(freq_mod.AGGREGATE_MODES),
        "coverage": freq_mod.coverage_report(),
        "note": (
            "Choosing a population changes the frequency denominator only. It "
            "does not infer your ancestry and does not change which variants "
            "you carry."
        ),
    })


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

@api_v2.route("/api/profiles/<int:pid>/sources", methods=["GET"])
def list_sources(pid: int):
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    sources = _load_sources(pid)
    cached = _cached(pid) or {}
    stats = {s.get("label"): s for s in (cached.get("sources") or [])}
    out = []
    for s in sources:
        row = {
            "id": s.get("id"), "label": s.get("label", ""),
            "role": s.get("role", "other"), "provider": s.get("provider", ""),
            "snp_count": len(s.get("snps") or []),
            "uploaded_at": s.get("uploaded_at", ""),
        }
        row.update({k: v for k, v in (stats.get(s.get("label")) or {}).items()
                    if k in ("contributed", "overlapped", "conflicting")})
        out.append(row)
    return jsonify({"sources": out, "roles": list(VALID_ROLES)})


@api_v2.route("/api/profiles/<int:pid>/sources", methods=["POST"])
def add_source(pid: int):
    """Attach an additional raw file to an existing profile.

    A role of ``self`` POOLS the file into the primary genotype set, which is how
    a 23andMe file and an AncestryDNA file from the same person are combined.
    Any other role adds comparison rows only and never alters the primary calls.
    """
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    if "file" not in request.files:
        return _err("No file provided.")

    upload = request.files["file"]
    role = (request.form.get("role") or "other").strip().lower()
    if role not in VALID_ROLES:
        return _err(f"role must be one of: {', '.join(VALID_ROLES)}")

    name = secure_filename(upload.filename or "") or "dna.txt"
    if Path(name).suffix.lower() not in ALLOWED_EXTENSIONS:
        return _err("Upload an uncompressed .txt, .csv or .tsv export, not a .zip or .gz.")

    dest = UPLOAD_DIR / f"p{pid}_{len(_load_sources(pid)) + 1}_{name}"
    upload.save(str(dest))
    try:
        parsed = parse_dna_file(str(dest))
    except ParseError as exc:
        dest.unlink(missing_ok=True)
        return _err(str(exc))

    sources = _load_sources(pid)
    sources.append({
        "id": len(sources) + 1,
        "label": request.form.get("label") or name,
        "role": role,
        "provider": parsed["provider"],
        "snps": parsed["snps"],
        "uploaded_at": _now(),
    })
    _save_sources(pid, sources)
    db.record_upload(pid, name, parsed["snp_count"])

    return jsonify({
        "message": (f"Added {parsed['snp_count']:,} SNPs as role '{role}'. "
                    "Run a scan to fold it into your findings."),
        "source_count": len(sources),
        "provider": parsed["provider"],
        "snp_count": parsed["snp_count"],
    }), 201


@api_v2.route("/api/profiles/<int:pid>/sources/<int:sid>", methods=["PATCH"])
def patch_source(pid: int, sid: int):
    sources = _load_sources(pid)
    target = next((s for s in sources if s.get("id") == sid), None)
    if target is None:
        return _err("Source not found.", 404)
    body = request.get_json(silent=True) or {}
    if "role" in body:
        role = str(body["role"]).strip().lower()
        if role not in VALID_ROLES:
            return _err(f"role must be one of: {', '.join(VALID_ROLES)}")
        target["role"] = role
    if "label" in body:
        target["label"] = str(body["label"]).strip() or target["label"]
    _save_sources(pid, sources)
    return jsonify({"message": "Source updated. Re-run the scan to apply."})


@api_v2.route("/api/profiles/<int:pid>/sources/<int:sid>", methods=["DELETE"])
def delete_source(pid: int, sid: int):
    sources = _load_sources(pid)
    remaining = [s for s in sources if s.get("id") != sid]
    if len(remaining) == len(sources):
        return _err("Source not found.", 404)
    _save_sources(pid, remaining)
    return jsonify({"message": "Source removed. Re-run the scan to apply."})


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

@api_v2.route("/api/profiles/<int:pid>/scan/v2", methods=["POST"])
def scan_v2(pid: int):
    """Run the full v2 pipeline in a background thread."""
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    sources = _load_sources(pid)
    if not sources:
        return _err("No DNA data for this profile. Upload a raw file first.")

    body = request.get_json(silent=True) or {}
    options = {
        "use_api": bool(body.get("use_api", False)),
        "population": str(body.get("population") or freq_mod.DEFAULT_POPULATION),
        "include_genosets": bool(body.get("include_genosets", True)),
        "include_traits": bool(body.get("include_traits", True)),
        "include_prs": bool(body.get("include_prs", True)),
        "use_snpedia": bool(body.get("use_snpedia", False)),
    }

    with _lock:
        if _scan_progress.get(pid, {}).get("running"):
            return _err("A scan is already running for this profile.", 409)
        _scan_progress[pid] = {"running": True, "done": False, "phase": "merge",
                               "processed": 0, "total": 0, "findings": 0,
                               "error": None}

    def progress(phase: str, done: int = 0, total: int = 0) -> None:
        with _lock:
            state = _scan_progress.setdefault(pid, {})
            state["phase"] = phase
            if total:
                state["processed"], state["total"] = done, total

    def worker() -> None:
        try:
            result = pipeline.run_full_scan(sources, progress_cb=progress, **options)
            with _lock:
                _scan_cache[pid] = result
            uid = None
            conn = db.get_connection()
            row = conn.execute(
                "SELECT id FROM snp_uploads WHERE profile_id=? "
                "ORDER BY uploaded_at DESC LIMIT 1", (pid,)).fetchone()
            conn.close()
            if row:
                uid = row["id"]
            for finding in result["findings"]:
                try:
                    db.upsert_finding(pid, uid, finding)
                except Exception:
                    continue
            with _lock:
                _scan_progress[pid] = {
                    "running": False, "done": True, "phase": "complete",
                    "processed": len(result["findings"]),
                    "total": len(result["findings"]),
                    "findings": len(result["findings"]), "error": None,
                }
        except Exception as exc:
            with _lock:
                _scan_progress[pid] = {"running": False, "done": True,
                                       "phase": "error", "error": str(exc)}

    threading.Thread(target=worker, daemon=True).start()
    return jsonify({"message": "Scan started.", "profile_id": pid,
                    "options": options, "phases": list(pipeline.PHASES)})


@api_v2.route("/api/profiles/<int:pid>/scan/v2/status", methods=["GET"])
def scan_v2_status(pid: int):
    with _lock:
        return jsonify(_scan_progress.get(
            pid, {"running": False, "done": False, "phase": "", "error": None}))


# ---------------------------------------------------------------------------
# Findings, filtered
# ---------------------------------------------------------------------------

def _findings_for(pid: int, population: str | None = None) -> tuple[list[dict], dict]:
    """Return (findings, extras) from cache, falling back to the database."""
    cached = _cached(pid)
    if cached:
        findings = cached["findings"]
        if population and population != cached.get("population"):
            for f in findings:
                if f.get("entity_type") == "snp":
                    freq_mod.annotate(f, population)
            scoring_mod.score_all(findings)
            cached["population"] = population
            cached["ranges"] = pipeline.compute_ranges(findings)
        return findings, cached
    rows = db.get_findings(pid)
    for r in rows:
        r.setdefault("entity_type", "snp")
        if population:
            freq_mod.annotate(r, population)
    scoring_mod.score_all(rows)
    return rows, {"qc": {}, "ranges": pipeline.compute_ranges(rows),
                  "population": population or freq_mod.DEFAULT_POPULATION}


@api_v2.route("/api/profiles/<int:pid>/findings/v2", methods=["GET"])
def findings_v2(pid: int):
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    population = request.args.get("population") or None
    findings, extras = _findings_for(pid, population)

    params = {k: v for k, v in request.args.items()}
    result = F.filter_and_sort(findings, params)
    result["summary"] = F.summarise(findings)
    result["ranges"] = extras.get("ranges") or pipeline.compute_ranges(findings)
    result["qc"] = extras.get("qc") or {}
    result["population"] = extras.get("population") or freq_mod.DEFAULT_POPULATION
    result["sort"] = params.get("sort") or "magnitude"
    result["order"] = params.get("order") or "desc"
    return jsonify(result)


@api_v2.route("/api/profiles/<int:pid>/facets", methods=["GET"])
def facets(pid: int):
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    findings, _ = _findings_for(pid)
    return jsonify(F.build_facets(findings))


# ---------------------------------------------------------------------------
# Subsystem views
# ---------------------------------------------------------------------------

def _require_scan(pid: int):
    cached = _cached(pid)
    if cached is None:
        return None, _err("Run a v2 scan first so this view has data.", 409)
    return cached, None


@api_v2.route("/api/profiles/<int:pid>/genosets", methods=["GET"])
def genosets_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    result = dict(cached["genosets"])
    result["note"] = (
        "A genoset in the 'not testable' list requires a position your array "
        "did not genotype, so it could not be evaluated. That is different from "
        "having been checked and found absent."
    )
    return jsonify(result)


@api_v2.route("/api/profiles/<int:pid>/traits", methods=["GET"])
def traits_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    return jsonify({"traits": cached["traits"], "blood_type": cached["blood_type"]})


@api_v2.route("/api/profiles/<int:pid>/prs", methods=["GET"])
def prs_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    return jsonify({
        "results": cached["prs"],
        "disclaimer": (
            "A polygenic score is a statistical predictor, not a diagnostic "
            "test. It is computed from a small subset of known variants, it is "
            "calibrated on a European reference population and transfers poorly "
            "to others, and it ignores every environmental and lifestyle factor. "
            "A high score does not mean you will develop the condition and a low "
            "score does not mean you will not."
        ),
    })


@api_v2.route("/api/profiles/<int:pid>/pgx", methods=["GET"])
def pgx_view(pid: int):
    """Pharmacogenomic findings grouped by drug."""
    cached, error = _require_scan(pid)
    if error:
        return error
    by_drug: dict[str, list[dict]] = {}
    for f in cached["findings"]:
        for med in f.get("medicines") or []:
            by_drug.setdefault(str(med), []).append({
                "rsid": f.get("rsid"), "gene": f.get("gene"),
                "genotype": f.get("genotype"), "zygosity": f.get("zygosity"),
                "magnitude": f.get("magnitude"), "repute": f.get("repute"),
                "cpic_level": f.get("cpic_level"), "silo": f.get("silo"),
                "confidence": f.get("confidence"),
                "summary": f.get("summary") or f.get("interpretation"),
                "carrier": f.get("carrier"),
            })
    drugs = [{"drug": k, "variants": sorted(v, key=lambda x: -(x["magnitude"] or 0)),
              "highest_cpic": next((x["cpic_level"] for x in v if x["cpic_level"]), "")}
             for k, v in sorted(by_drug.items())]
    drugs.sort(key=lambda d: (d["highest_cpic"] != "A",
                              -max((x["magnitude"] or 0) for x in d["variants"])))
    return jsonify({
        "drugs": drugs,
        "count": len(drugs),
        "warning": (
            "Show this to your prescriber or pharmacist BEFORE any medication "
            "change. Do not stop, start or adjust a medication based on this "
            "page. Consumer arrays do not call star alleles completely, so a "
            "normal result here does not rule out a variant."
        ),
    })


@api_v2.route("/api/profiles/<int:pid>/conflicts", methods=["GET"])
def conflicts_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    return jsonify({
        "conflicts": cached["conflicts"],
        "count": len(cached["conflicts"]),
        "note": (
            "Both calls are kept. Nothing was reconciled and no winner was "
            "chosen, because a disagreement between two arrays is information "
            "about reliability that a silent merge would destroy."
        ),
    })


@api_v2.route("/api/profiles/<int:pid>/trio", methods=["GET"])
def trio_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    return jsonify(cached["trio"])


@api_v2.route("/api/profiles/<int:pid>/qc", methods=["GET"])
def qc_view(pid: int):
    cached, error = _require_scan(pid)
    if error:
        return error
    return jsonify({
        "qc": cached["qc"],
        "counts": cached["counts"],
        "sources": cached["sources"],
        "subsystems": cached["subsystems"],
    })


# ---------------------------------------------------------------------------
# Export honouring the active filters
# ---------------------------------------------------------------------------

_EXPORT_COLUMNS = (
    "rsid", "entity_type", "gene", "chromosome", "position", "genotype",
    "zygosity", "magnitude", "repute", "confidence", "silo", "category",
    "clinical_sig", "review_stars", "cpic_level", "evidence", "publications",
    "freq", "freq_population", "freq_band", "gmaf", "carrier",
    "variant_copies", "flipped", "ambiguous", "dubious", "count",
    "summary", "interpretation",
)


def _filtered_for_export(pid: int) -> list[dict]:
    findings, _ = _findings_for(pid, request.args.get("population") or None)
    params = {k: v for k, v in request.args.items()}
    params.setdefault("limit", 0)
    return F.filter_and_sort(findings, params)["findings"]


@api_v2.route("/api/profiles/<int:pid>/export/v2/<fmt>", methods=["GET"])
def export_v2(pid: int, fmt: str):
    profile = db.get_profile(pid)
    if not profile:
        return _err("Profile not found.", 404)
    fmt = fmt.lower()
    if fmt not in ("json", "csv", "tsv"):
        return _err("Format must be json, csv or tsv.")

    rows = _filtered_for_export(pid)
    safe = "".join(c for c in profile["name"] if c.isalnum() or c in " _-").strip().replace(" ", "_")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")

    if fmt == "json":
        payload = {
            "exported_at": _now(),
            "app_version": APP_VERSION,
            "patient": {k: profile.get(k, "") for k in ("name", "dob", "sex", "provider")},
            "filters": {k: v for k, v in request.args.items()},
            "findings_count": len(rows),
            "findings": rows,
        }
        buf = io.BytesIO(json.dumps(payload, indent=2).encode("utf-8"))
        return send_file(buf, mimetype="application/json", as_attachment=True,
                         download_name=f"{safe}_findings_{stamp}.json")

    delimiter = "\t" if fmt == "tsv" else ","
    text = io.StringIO()
    writer = csv.writer(text, delimiter=delimiter, lineterminator="\n")
    writer.writerow(_EXPORT_COLUMNS)
    for r in rows:
        writer.writerow([_flatten(r.get(c)) for c in _EXPORT_COLUMNS])
    data = io.BytesIO(text.getvalue().encode("utf-8-sig"))
    mime = "text/tab-separated-values" if fmt == "tsv" else "text/csv"
    return send_file(data, mimetype=mime, as_attachment=True,
                     download_name=f"{safe}_findings_{stamp}.{fmt}")


def _flatten(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (list, tuple)):
        return "; ".join(str(v) for v in value)
    return str(value)


# ---------------------------------------------------------------------------
# SNPedia admin, licence-gated
# ---------------------------------------------------------------------------

_harvest_state: dict = {"running": False, "done": False, "error": None}


def _snpedia():
    try:
        from . import snpedia
        return snpedia
    except Exception:
        return None


@api_v2.route("/api/profiles/<int:pid>/reports/interactive", methods=["POST"])
def interactive_report(pid: int):
    """Generate the self-contained offline interactive report.

    Honours every active filter, so the exported file contains exactly the
    selection the user was looking at. The result needs no server and makes no
    network requests, which is the whole point of it.
    """
    profile = db.get_profile(pid)
    if not profile:
        return _err("Profile not found.", 404)

    try:
        from .interactive_report import generate_interactive_report
    except Exception as exc:
        return _err(f"Interactive report support is unavailable: {exc}", 500)

    findings, extras = _findings_for(pid, request.args.get("population") or None)
    params = {k: v for k, v in request.args.items()}
    params.setdefault("limit", 0)
    selected = F.filter_and_sort(findings, params)["findings"]
    if not selected:
        return _err("No findings match the current filters, so there is nothing "
                    "to put in the report.", 400)

    cached = _cached(pid) or {}
    html = generate_interactive_report(profile, selected, {
        "population": extras.get("population") or freq_mod.DEFAULT_POPULATION,
        "counts": cached.get("counts", {}),
        "qc": cached.get("qc", {}),
        "sources": cached.get("sources", []),
        "genosets": cached.get("genosets", {}),
        "blood_type": cached.get("blood_type", {}),
    })

    safe = "".join(c for c in profile["name"] if c.isalnum() or c in " _-").strip().replace(" ", "_")
    filename = _unique_report_name(safe, "interactive")
    path = REPORTS_DIR / filename
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(html)

    rid = db.record_report(pid, "interactive", "html", str(path))
    return jsonify({
        "report_id": rid,
        "type": "interactive",
        "filename": filename,
        "view_url": f"/api/reports/{rid}/view",
        "findings_included": len(selected),
        "bytes": len(html.encode("utf-8")),
        "message": (f"Interactive report built from {len(selected)} findings. "
                    "It opens with no server and no internet."),
    }), 201


@api_v2.route("/api/admin/snpedia/status", methods=["GET"])
def snpedia_status():
    module = _snpedia()
    if module is None:
        return jsonify({"available": False, "reason": "module not installed"})
    status = dict(module.cache_status())
    status["notice"] = getattr(module, "NOTICE", "")
    status["path"] = str(module.cache_path())
    with _lock:
        status["harvest"] = dict(_harvest_state)
    return jsonify(status)


@api_v2.route("/api/admin/snpedia/harvest", methods=["POST"])
def snpedia_harvest():
    """Start a local SNPedia harvest. Requires explicit licence acceptance.

    Returns 403 with the full notice when the caller has not accepted, because
    SNPedia is CC-BY-NC-SA-3.0-US and the cache is for the user's own personal,
    non-commercial use only. Nothing harvested is ever written inside the
    repository.
    """
    module = _snpedia()
    if module is None:
        return _err("SNPedia support is not installed.", 404)
    body = request.get_json(silent=True) or {}
    if not body.get("accept_license"):
        return jsonify({
            "error": "Licence acceptance required before harvesting.",
            "notice": getattr(module, "NOTICE", ""),
            "license": "CC-BY-NC-SA-3.0-US",
            "license_url": "http://creativecommons.org/licenses/by-nc-sa/3.0/us/",
        }), 403

    with _lock:
        if _harvest_state.get("running"):
            return _err("A harvest is already running.", 409)
        _harvest_state.update({"running": True, "done": False, "error": None})

    def worker() -> None:
        try:
            module.harvest(rsids=body.get("rsids"), accept_license=True)
            with _lock:
                _harvest_state.update({"running": False, "done": True, "error": None})
        except Exception as exc:
            with _lock:
                _harvest_state.update({"running": False, "done": True,
                                       "error": str(exc)})

    threading.Thread(target=worker, daemon=True).start()
    return jsonify({"message": "Harvest started. This runs locally on your machine."}), 202


@api_v2.route("/api/admin/snpedia/harvest/status", methods=["GET"])
def snpedia_harvest_status():
    with _lock:
        return jsonify(dict(_harvest_state))


@api_v2.route("/api/admin/snpedia/cache", methods=["DELETE"])
def snpedia_purge():
    module = _snpedia()
    if module is None:
        return _err("SNPedia support is not installed.", 404)
    try:
        removed = module.purge_cache()
    except Exception as exc:
        return _err(f"Could not purge the cache: {exc}", 500)
    return jsonify({"message": "Local SNPedia cache purged.", "removed": removed})
