"""
routes.py -- Flask API endpoints for DNAInsight.

Endpoints:
  GET  /api/profiles                  -- list all profiles
  POST /api/profiles                  -- create profile + upload DNA file
  GET  /api/profiles/<id>             -- get profile detail + findings summary
  DEL  /api/profiles/<id>             -- delete profile and all data
  POST /api/profiles/<id>/scan        -- run SNP annotation scan
  GET  /api/profiles/<id>/findings    -- get all findings (optional ?silo=)
  GET  /api/profiles/<id>/reports     -- list generated reports
  POST /api/profiles/<id>/reports     -- generate report (type: genetic|doctor)
  GET  /api/reports/<id>/view         -- serve report HTML
  GET  /api/status                    -- health check
"""

import os
import io
import csv
import json
import threading
from pathlib import Path
from datetime import datetime
from flask import Blueprint, request, jsonify, send_file, current_app
from werkzeug.utils import secure_filename

from . import database as db
from . import APP_VERSION
from .parsers import parse_dna_file, ParseError
from .scanner import run_scan, lookup_bundled, zygosity_of
from .genetic_report import generate_genetic_report
from .doctor_report import generate_doctor_report
from .updater import update_bundled_reference, get_reference_metadata

# Accepted raw-DNA upload extensions and hard size cap. Consumer exports are
# typically 15-25 MB uncompressed; 64 MB leaves headroom while preventing
# memory-exhaustion via oversized uploads.
ALLOWED_EXTENSIONS = {".txt", ".csv", ".tsv"}
MAX_UPLOAD_BYTES   = 64 * 1024 * 1024

api = Blueprint("api", __name__)

UPLOAD_DIR  = Path(__file__).parent.parent / "uploads"
REPORTS_DIR = Path(__file__).parent.parent / "reports_output"
UPLOAD_DIR.mkdir(exist_ok=True)
REPORTS_DIR.mkdir(exist_ok=True)

# Progress tracker for long-running scans
_scan_progress: dict = {}
_scan_lock = threading.Lock()

# Progress tracker for database updates
_update_progress: dict = {"running": False, "done": False, "processed": 0, "total": 0}
_update_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------

@api.route("/api/status")
def status():
    return jsonify({"status": "ok", "version": APP_VERSION})


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------

@api.route("/api/profiles", methods=["GET"])
def list_profiles():
    profiles = db.list_profiles()
    for p in profiles:
        p["findings_summary"] = db.get_findings_summary(p["id"])
    return jsonify(profiles)


@api.route("/api/profiles", methods=["POST"])
def create_profile():
    if "file" not in request.files:
        return jsonify({"error": "No DNA file provided. Upload a raw DNA file (.txt)."}), 400

    file  = request.files["file"]
    name  = request.form.get("name", "").strip()
    dob   = request.form.get("dob", "")
    sex   = request.form.get("sex", "unknown")

    if not name:
        return jsonify({"error": "Name is required."}), 400
    if not file.filename:
        return jsonify({"error": "No file selected."}), 400

    # Validate extension (defence-in-depth; reject compressed/binary uploads early)
    orig_name = secure_filename(file.filename) or "dna.txt"
    ext = Path(orig_name).suffix.lower()
    if ext not in ALLOWED_EXTENSIONS:
        return jsonify({
            "error": f"Unsupported file type '{ext or 'unknown'}'. "
                     "Upload an uncompressed raw DNA export (.txt, .csv, or .tsv), not a .zip/.gz."
        }), 400

    # Enforce a hard size cap before writing to disk.
    file.seek(0, os.SEEK_END)
    size = file.tell()
    file.seek(0)
    if size == 0:
        return jsonify({"error": "Uploaded file is empty."}), 400
    if size > MAX_UPLOAD_BYTES:
        return jsonify({
            "error": f"File too large ({size // (1024*1024)} MB). Maximum is "
                     f"{MAX_UPLOAD_BYTES // (1024*1024)} MB."
        }), 413

    # Save uploaded file using a fully sanitized name (prevents path traversal).
    safe_name = "".join(c for c in name if c.isalnum() or c in " _-").strip().replace(" ", "_") or "profile"
    dest = UPLOAD_DIR / f"{safe_name}_{orig_name}"
    file.save(str(dest))

    # Parse
    try:
        parsed = parse_dna_file(str(dest))
    except ParseError as e:
        dest.unlink(missing_ok=True)
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        dest.unlink(missing_ok=True)
        return jsonify({"error": f"Parse failed: {e}"}), 500

    provider = parsed["provider"]
    snp_count = parsed["snp_count"]

    # Create profile
    pid = db.create_profile(name, dob, sex, provider)
    uid = db.record_upload(pid, orig_name, snp_count)

    # Store snps in a temp JSON for scanning
    snp_cache = UPLOAD_DIR / f"snps_{pid}.json"
    with open(snp_cache, "w") as f:
        json.dump(parsed["snps"], f)

    return jsonify({
        "profile_id": pid,
        "upload_id":  uid,
        "provider":   provider,
        "snp_count":  snp_count,
        "message":    f"Uploaded {snp_count:,} SNPs from {provider} array.",
    }), 201


@api.route("/api/profiles/<int:pid>", methods=["GET"])
def get_profile(pid: int):
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404
    profile["findings_summary"] = db.get_findings_summary(pid)
    profile["reports"] = db.get_reports(pid)
    return jsonify(profile)


@api.route("/api/profiles/<int:pid>", methods=["DELETE"])
def delete_profile(pid: int):
    db.delete_profile(pid)
    snp_cache = UPLOAD_DIR / f"snps_{pid}.json"
    snp_cache.unlink(missing_ok=True)
    return jsonify({"message": "Profile deleted."})


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

@api.route("/api/profiles/<int:pid>/scan", methods=["POST"])
def start_scan(pid: int):
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404

    snp_cache = UPLOAD_DIR / f"snps_{pid}.json"
    if not snp_cache.exists():
        return jsonify({"error": "DNA data not found. Re-upload the file."}), 400

    body      = request.get_json(silent=True) or {}
    use_api   = body.get("use_api", True)

    with _scan_lock:
        if _scan_progress.get(pid, {}).get("running"):
            return jsonify({"error": "Scan already in progress for this profile."}), 409
        _scan_progress[pid] = {"running": True, "processed": 0, "total": 0, "findings": 0}

    def _run():
        try:
            with open(snp_cache, "r") as f:
                snps = json.load(f)

            def progress_cb(processed, total):
                with _scan_lock:
                    _scan_progress[pid]["processed"] = processed
                    _scan_progress[pid]["total"]      = total

            findings = run_scan(snps, use_api=use_api, progress_cb=progress_cb)

            # Get last upload id
            conn = db.get_connection()
            row  = conn.execute(
                "SELECT id FROM snp_uploads WHERE profile_id=? ORDER BY uploaded_at DESC LIMIT 1",
                (pid,)
            ).fetchone()
            conn.close()
            uid = row["id"] if row else None

            for finding in findings:
                db.upsert_finding(pid, uid, finding)

            with _scan_lock:
                _scan_progress[pid] = {
                    "running":  False,
                    "processed": len(snps),
                    "total":     len(snps),
                    "findings":  len(findings),
                    "done":      True,
                }
        except Exception as e:
            with _scan_lock:
                _scan_progress[pid] = {"running": False, "error": str(e), "done": True}

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()

    return jsonify({"message": "Scan started.", "profile_id": pid})


@api.route("/api/profiles/<int:pid>/scan/status", methods=["GET"])
def scan_status(pid: int):
    with _scan_lock:
        status_data = _scan_progress.get(pid, {"running": False, "done": False})
    return jsonify(status_data)


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------

@api.route("/api/profiles/<int:pid>/findings", methods=["GET"])
def get_findings(pid: int):
    silo = request.args.get("silo")
    findings = db.get_findings(pid, silo=silo)
    summary  = db.get_findings_summary(pid)
    return jsonify({"findings": findings, "summary": summary, "total": len(findings)})


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

@api.route("/api/profiles/<int:pid>/reports", methods=["GET"])
def list_reports(pid: int):
    return jsonify(db.get_reports(pid))


@api.route("/api/profiles/<int:pid>/reports", methods=["POST"])
def generate_report(pid: int):
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404

    body        = request.get_json(silent=True) or {}
    report_type = body.get("type", "genetic")

    if report_type not in ("genetic", "doctor"):
        return jsonify({"error": "type must be 'genetic' or 'doctor'"}), 400

    findings = db.get_findings(pid)
    if not findings:
        return jsonify({"error": "No findings yet. Run a scan first."}), 400

    safe_name = "".join(c for c in profile["name"] if c.isalnum() or c in " _-").strip().replace(" ", "_")
    from datetime import datetime
    ts = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    filename = f"{safe_name}_{report_type}_report_{ts}.html"
    filepath = REPORTS_DIR / filename

    if report_type == "genetic":
        html = generate_genetic_report(profile, findings)
    else:
        html = generate_doctor_report(profile, findings)

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(html)

    rid = db.record_report(pid, report_type, "html", str(filepath))

    return jsonify({
        "report_id":   rid,
        "type":        report_type,
        "filename":    filename,
        "view_url":    f"/api/reports/{rid}/view",
        "message":     f"{report_type.title()} report generated.",
    }), 201


@api.route("/api/reports/<int:rid>/view", methods=["GET"])
def view_report(rid: int):
    conn = db.get_connection()
    row  = conn.execute("SELECT * FROM reports WHERE id=?", (rid,)).fetchone()
    conn.close()
    if not row:
        return jsonify({"error": "Report not found."}), 404
    filepath = Path(dict(row)["filepath"])
    if not filepath.exists():
        return jsonify({"error": "Report file missing from disk."}), 404
    return send_file(str(filepath), mimetype="text/html")


# ---------------------------------------------------------------------------
# Admin: Database management
# ---------------------------------------------------------------------------

@api.route("/api/admin/db-status", methods=["GET"])
def db_status():
    """Return current bundled reference metadata (version, last updated, SNP count)."""
    meta = get_reference_metadata()
    with _update_lock:
        meta["update_running"] = _update_progress.get("running", False)
    return jsonify(meta)


@api.route("/api/admin/update-databases", methods=["POST"])
def start_db_update():
    """
    Trigger a background re-fetch of all bundled SNP annotations from MyVariant.info.
    Returns 409 if an update is already in progress.
    """
    with _update_lock:
        if _update_progress.get("running"):
            return jsonify({"error": "Database update already in progress."}), 409
        _update_progress.clear()
        _update_progress.update({
            "running":   True,
            "done":      False,
            "processed": 0,
            "total":     0,
            "error":     None,
        })

    def _run_update():
        def progress_cb(processed, total):
            with _update_lock:
                _update_progress["processed"] = processed
                _update_progress["total"]      = total

        try:
            result = update_bundled_reference(progress_cb=progress_cb)
            with _update_lock:
                _update_progress.update({
                    "running":   False,
                    "done":      True,
                    "error":     None,
                    "updated":   result.get("updated", 0),
                    "skipped":   result.get("skipped", 0),
                    "errors":    result.get("errors", 0),
                    "timestamp": result.get("timestamp"),
                    "snp_count": result.get("snp_count", 0),
                })
        except Exception as e:
            with _update_lock:
                _update_progress.update({
                    "running": False,
                    "done":    True,
                    "error":   str(e),
                })

    thread = threading.Thread(target=_run_update, daemon=True)
    thread.start()

    return jsonify({"message": "Database update started."}), 202


@api.route("/api/admin/update-databases/status", methods=["GET"])
def db_update_status():
    """Poll the status of a running database update."""
    with _update_lock:
        return jsonify(dict(_update_progress))


# ---------------------------------------------------------------------------
# Export endpoints
# ---------------------------------------------------------------------------

@api.route("/api/profiles/<int:pid>/export/json", methods=["GET"])
def export_findings_json(pid: int):
    """Export all findings for a profile as a downloadable JSON file."""
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404

    findings = db.get_findings(pid)
    payload = {
        "exported_at": datetime.utcnow().isoformat() + "Z",
        "app_version": APP_VERSION,
        "patient": {
            "name": profile["name"],
            "dob":  profile.get("dob", ""),
            "sex":  profile.get("sex", ""),
            "provider": profile.get("provider", ""),
        },
        "findings_count": len(findings),
        "findings": findings,
    }

    buf = io.BytesIO(json.dumps(payload, indent=2).encode("utf-8"))
    safe = profile["name"].replace(" ", "_")
    ts   = datetime.utcnow().strftime("%Y%m%d")
    return send_file(
        buf,
        mimetype="application/json",
        as_attachment=True,
        download_name=f"{safe}_findings_{ts}.json",
    )


@api.route("/api/profiles/<int:pid>/export/csv", methods=["GET"])
def export_findings_csv(pid: int):
    """Export all findings for a profile as a downloadable CSV file."""
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404

    findings = db.get_findings(pid)

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([
        "rsid", "gene", "genotype", "zygosity", "category", "silo",
        "clinical_sig", "interpretation", "discovered_at",
    ])
    for f in findings:
        gt = f.get("genotype") or f"{f.get('allele1','')}{f.get('allele2','')}"
        writer.writerow([
            f.get("rsid", ""),
            f.get("gene", ""),
            gt,
            f.get("zygosity") or zygosity_of(f.get("allele1", ""), f.get("allele2", "")),
            f.get("category", ""),
            f.get("silo", ""),
            f.get("clinical_sig", ""),
            f.get("interpretation", ""),
            f.get("discovered_at", ""),
        ])

    bytes_buf = io.BytesIO(buf.getvalue().encode("utf-8-sig"))  # utf-8-sig for Excel BOM
    safe = profile["name"].replace(" ", "_")
    ts   = datetime.utcnow().strftime("%Y%m%d")
    return send_file(
        bytes_buf,
        mimetype="text/csv",
        as_attachment=True,
        download_name=f"{safe}_findings_{ts}.csv",
    )


# ---------------------------------------------------------------------------
# Version check endpoint
# ---------------------------------------------------------------------------

@api.route("/api/version", methods=["GET"])
def version_info():
    """Return current app version. Frontend uses this for update-available banners."""
    return jsonify({
        "version": APP_VERSION,
        "snp_count": get_reference_metadata().get("snp_count", 0),
    })


# ---------------------------------------------------------------------------
# Single-SNP lookup
# ---------------------------------------------------------------------------

@api.route("/api/profiles/<int:pid>/lookup/<rsid>", methods=["GET"])
def lookup_snp(pid: int, rsid: str):
    """
    Look up a single rsID for a profile.

    Returns the person's genotype (from their uploaded raw data), the zygosity,
    and — if present — the bundled reference annotation for that variant. Enables
    ad-hoc questions like "what's my rs1801133 genotype?" without a full scan.
    """
    profile = db.get_profile(pid)
    if not profile:
        return jsonify({"error": "Profile not found."}), 404

    rsid = rsid.strip().lower()
    if not rsid.startswith("rs"):
        return jsonify({"error": "rsID must start with 'rs' (e.g. rs1801133)."}), 400

    snp_cache = UPLOAD_DIR / f"snps_{pid}.json"
    if not snp_cache.exists():
        return jsonify({"error": "DNA data not found. Re-upload the file."}), 400

    with open(snp_cache, "r") as fh:
        snps = json.load(fh)

    match = next((s for s in snps if s.get("rsid", "").lower() == rsid), None)
    if not match:
        return jsonify({
            "rsid": rsid,
            "in_your_data": False,
            "message": "This rsID was not genotyped on your array.",
        })

    a1, a2 = match.get("allele1", ""), match.get("allele2", "")
    ref = lookup_bundled(rsid) or {}
    return jsonify({
        "rsid":         rsid,
        "in_your_data": True,
        "genotype":     f"{a1}{a2}",
        "zygosity":     zygosity_of(a1, a2),
        "chromosome":   match.get("chromosome", ""),
        "position":     match.get("position", 0),
        "in_reference": bool(ref),
        "gene":         ref.get("gene", ""),
        "category":     ref.get("category", ""),
        "clinical_sig": ref.get("clinical_sig", ""),
        "interpretation": ref.get("interpretation", ""),
    })
