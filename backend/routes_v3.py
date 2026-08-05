"""
routes_v3.py -- the v3 HTTP surface, as a third blueprint.

Deliberately separate from routes.py and routes_v2.py for the same reason those
two are separate from each other: every v1 and v2 endpoint keeps its exact
behaviour, and a v3 subsystem that fails to import must not take the rest of the
application down with it. app.py registers the blueprints in order and wraps
each optional one in a try, so a checkout with no v3 data still boots.

Implements docs/API_V3.md.

THE DEGRADATION CONTRACT, WHICH IS THE POINT OF THIS FILE
---------------------------------------------------------
Six of the ten v3 capabilities need a third-party tool the user installs
themselves. None of them may raise when that tool is absent. Every one of them
returns the standard payload from ``backend.external.unavailable``, which
carries ``available: False``, ``not_attempted: True`` and a plain-English reason.

That last flag is the load-bearing one. "We looked and found nothing" and "we
could not look at all" are different claims, and collapsing them is the exact
failure this project already refuses in three other places: a genoset over
positions the array never read is reported as not testable rather than absent, a
strand that cannot be verified is badged rather than guessed, and a no-call
scores zero rather than counting as a negative finding. Ancestry, haplogroups,
imputation and IBD now obey the same rule.

LICENCE GATE
------------
``POST /api/v3/tools/<id>/licence`` returns HTTP 403 with the full licence
notice unless the body carries ``accept_license: true``. That is the same shape
as the SNPedia harvest gate in routes_v2.py, on purpose: one rule to learn.
Tools whose licence forbids redistribution or commercial use return 409 and can
never be accepted at all.
"""

from __future__ import annotations

import json
import threading
from datetime import datetime, timezone
from pathlib import Path

from flask import Blueprint, request, jsonify
from werkzeug.utils import secure_filename

from . import APP_VERSION
from . import database as db
from . import external
from . import merge as merge_mod
from . import pipeline
from .routes import _bounded_upload_path
from .routes_v2 import (
    UPLOAD_DIR, REPORTS_DIR, VALID_ROLES,
    _load_sources, _save_sources, _cached, _err, _now,
)

api_v3 = Blueprint("api_v3", __name__)

_lock = threading.Lock()

# Long-running v3 jobs keep their state here. Imputation over a whole genome is
# minutes, not milliseconds, so it cannot block a request thread.
_jobs: dict[str, dict] = {}


# ---------------------------------------------------------------------------
# Optional subsystem imports.
#
# Same defensive pattern as pipeline._try. A module that fails to import is a
# missing capability, not a broken server.
# ---------------------------------------------------------------------------

def _mod(name: str):
    try:
        return __import__(f"backend.{name}", fromlist=[name])
    except Exception:
        return None


_ledger = _mod("ledger")
_provenance = _mod("provenance")
_sequencing = _mod("sequencing")
_haplogroups = _mod("haplogroups")
_relatedness = _mod("relatedness")
_imputation = _mod("imputation")
_ancestry = _mod("ancestry")
_diplotype = _mod("diplotype")
_carrier = _mod("carrier")
_assistant = _mod("assistant")


def _need(module, name: str):
    """Return a 501 payload when a subsystem module is absent.

    501 rather than 404, because the path exists and the capability is real; it
    is this installation that cannot serve it. A 404 would tell the frontend the
    endpoint was never part of the API, which is a different and wrong story.
    """
    if module is None:
        return jsonify({
            "available": False,
            "not_attempted": True,
            "error": f"The {name} subsystem is not installed in this build.",
            "capability": name,
        }), 501
    return None


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _merged_for(pid: int) -> dict | None:
    """Pool this profile's sources without requiring a completed scan.

    Ancestry, haplogroups and IBD all work off raw genotypes rather than
    annotated findings, so gating them behind a full scan would be an
    artificial dependency.
    """
    sources = _load_sources(pid)
    if not sources:
        return None
    try:
        return merge_mod.merge_sources(sources)
    except Exception:
        return None


def _genotypes_for(pid: int) -> dict | None:
    merged = _merged_for(pid)
    if merged is None:
        return None
    return pipeline.build_genotype_map(merged.get("genotypes") or {})


def _findings_for(pid: int) -> list[dict]:
    """Findings from the in-memory scan cache, falling back to the database."""
    cached = _cached(pid)
    if cached and cached.get("findings"):
        return list(cached["findings"])
    try:
        return db.get_findings(pid)
    except Exception:
        return []


def _profile_or_404(pid: int):
    if not db.get_profile(pid):
        return _err("Profile not found.", 404)
    return None


_NO_DNA = "No DNA data for this profile. Upload a raw file first."


def _guard_loaded(module, name: str, pid: int, loader):
    """Subsystem, profile and DNA checks in the one order that tells the truth.

    Returns ``(error_response, loaded)`` with exactly one of the two set to
    None. The ORDER is why this is a helper rather than nine lines copied into
    every endpoint that reads a genome: a missing subsystem is 501 and a missing
    profile is 404, so testing for genotypes first would answer "no DNA data,
    upload a raw file" to somebody whose profile does not exist or whose build
    cannot serve the capability at all. That is the same class of wrong answer
    the degradation contract at the top of this file exists to prevent, and it
    is not a mistake that should be re-avoidable once per endpoint.
    """
    guard = _need(module, name)
    if guard:
        return guard, None
    missing = _profile_or_404(pid)
    if missing:
        return missing, None
    loaded = loader(pid)
    if loaded is None:
        return _err(_NO_DNA), None
    return None, loaded


def _guard_genotypes(module, name: str, pid: int):
    """Guarded flat rsID to alleles map, for endpoints that score positions."""
    return _guard_loaded(module, name, pid, _genotypes_for)


def _guard_merged(module, name: str, pid: int):
    """Guarded merge result, for endpoints that need more than a flat map.

    Haplogroups needs the mtDNA position view and IBD needs the per-source
    split, and both are gone once the sources are flattened into a genotype
    map, so these endpoints cannot use ``_guard_genotypes``.
    """
    return _guard_loaded(module, name, pid, _merged_for)


def _job(kind: str, pid: int) -> str:
    return f"{kind}:{pid}"


# ---------------------------------------------------------------------------
# 1. External tools, panels and the licence gate
# ---------------------------------------------------------------------------

@api_v3.route("/api/v3/tools", methods=["GET"])
def tools_index():
    """Every external tool, every permanently excluded tool, every panel."""
    return jsonify(external.status_all())


@api_v3.route("/api/v3/tools/<tool_id>", methods=["GET"])
def tool_detail(tool_id: str):
    state = external.status(tool_id)
    state["licence_notice"] = external.licence_notice(tool_id)
    state["install"] = external.install_hint(tool_id)
    return jsonify(state)


@api_v3.route("/api/v3/tools/<tool_id>/licence", methods=["POST"])
def tool_accept_licence(tool_id: str):
    """Record acceptance of an external tool's licence.

    Refuses with 403 and the full notice until ``accept_license`` is true, so a
    user cannot end up having run something they never agreed to. Permanently
    excluded tools return 409: that state is not a missing consent and no body
    clears it.
    """
    body = request.get_json(silent=True) or {}
    if external.is_blocked(tool_id):
        return jsonify({
            "error": "This tool is permanently excluded on licence grounds.",
            "notice": external.licence_notice(tool_id),
            "status": external.status(tool_id),
        }), 409
    if not body.get("accept_license"):
        return jsonify({
            "error": "Licence acceptance required.",
            "notice": external.licence_notice(tool_id),
            "hint": "Repeat this request with {\"accept_license\": true}.",
        }), 403
    try:
        record = external.accept_licence(tool_id, accept=True)
    except external.ToolBlocked as exc:
        return jsonify({"error": str(exc)}), 409
    except external.ExternalError as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"accepted": record, "status": external.status(tool_id)})


@api_v3.route("/api/v3/tools/<tool_id>/licence", methods=["DELETE"])
def tool_revoke_licence(tool_id: str):
    """Withdraw acceptance. The capability reports unavailable again at once."""
    removed = external.revoke_licence(tool_id)
    return jsonify({"revoked": removed, "status": external.status(tool_id)})


@api_v3.route("/api/v3/panels/<panel_id>", methods=["GET"])
def panel_detail(panel_id: str):
    payload = external.panel_status(panel_id)
    if _ancestry is not None:
        try:
            payload["manifest"] = _ancestry.panel_manifest(panel_id)
        except Exception:
            payload["manifest"] = None
    return jsonify(payload)


@api_v3.route("/api/v3/licence-audit", methods=["GET"])
def licence_audit():
    """Prove the bundling rule in data/DATA_SOURCES.md still holds.

    This is the runtime half of a document that would otherwise drift. A
    non-empty ``violations`` list means something non-redistributable reached a
    bundled artefact, which is the failure DATA_SOURCES.md section 9 exists to
    prevent.
    """
    guard = _need(_provenance, "provenance")
    if guard:
        return guard
    return jsonify(_provenance.licence_audit())


# ---------------------------------------------------------------------------
# 2. Sequencing ingest
# ---------------------------------------------------------------------------

_SEQ_EXTENSIONS = {".vcf", ".gz", ".bgz", ".bam", ".cram"}


@api_v3.route("/api/profiles/<int:pid>/sources/sequencing", methods=["POST"])
def add_sequencing_source(pid: int):
    """Add a VCF, gVCF, BAM or CRAM source to a profile.

    Genome build is detected from contig lengths, which is the only header field
    that cannot quietly disagree with the coordinates in the body of the file. A
    build that does not match the bundled GRCh37 reference is REFUSED rather
    than silently accepted. Mixing builds is the most common way this class of
    tool produces confidently wrong answers, and it deserves the same loudness
    the project already gives strand ambiguity.
    """
    guard = _need(_sequencing, "sequencing")
    if guard:
        return guard
    missing = _profile_or_404(pid)
    if missing:
        return missing

    upload = request.files.get("file")
    if upload is None or not upload.filename:
        return _err("No file uploaded.")
    name = secure_filename(upload.filename)
    suffix = Path(name).suffix.lower()
    if suffix not in _SEQ_EXTENSIONS:
        return _err(
            "Expected a .vcf, .vcf.gz, .bam or .cram file. For array exports "
            "use the standard source endpoint instead."
        )

    role = str(request.form.get("role") or "self").strip().lower()
    if role not in VALID_ROLES:
        return _err(f"Role must be one of: {', '.join(VALID_ROLES)}.")

    dest = _bounded_upload_path(UPLOAD_DIR, f"{pid}_{name}")
    upload.save(dest)

    try:
        parsed = _sequencing.parse_sequencing_file(
            str(dest), sample=request.form.get("sample") or None
        )
    except _sequencing.BuildMismatch as exc:
        return jsonify({
            "error": str(exc),
            "detected_build": getattr(exc, "detected", None),
            "expected_build": getattr(exc, "expected", None),
            "hint": (
                "DNAInsight's bundled reference is GRCh37 plus strand. Lift the "
                "file over before uploading, or rebuild the reference against "
                "your build. Coordinates are not translated silently."
            ),
        }), 422
    except Exception as exc:
        return _err(f"Could not read this file: {exc}")

    if not parsed.get("snps"):
        return jsonify({
            "error": "No usable genotypes were read from this file.",
            "detail": parsed,
        }), 422

    sources = _load_sources(pid)
    next_id = max([s.get("id", 0) for s in sources], default=0) + 1
    sources.append({
        "id": next_id,
        "label": request.form.get("label") or name,
        "role": role,
        "provider": parsed.get("provider") or "vcf",
        "snps": parsed["snps"],
        "uploaded_at": _now(),
        "build": parsed.get("build"),
        "build_confidence": parsed.get("build_confidence"),
        "skipped": parsed.get("skipped"),
    })
    _save_sources(pid, sources)
    try:
        db.record_upload(pid, name, len(parsed["snps"]))
    except Exception:
        pass

    return jsonify({
        "message": "Sequencing source added.",
        "source_id": next_id,
        "format": parsed.get("format"),
        "build": parsed.get("build"),
        "build_confidence": parsed.get("build_confidence"),
        "build_evidence": parsed.get("build_evidence"),
        "sample": parsed.get("sample"),
        "sample_count": parsed.get("sample_count"),
        "snp_count": parsed.get("snp_count"),
        # Skipped records are reported, never swallowed. A user who uploads a
        # 4-million-record gVCF and gets 600,000 calls is owed the arithmetic.
        "skipped": parsed.get("skipped"),
        "warnings": parsed.get("warnings") or [],
    })


@api_v3.route("/api/v3/sequencing/inspect", methods=["POST"])
def sequencing_inspect():
    """Report format and build for an already-uploaded file without ingesting it."""
    guard = _need(_sequencing, "sequencing")
    if guard:
        return guard
    body = request.get_json(silent=True) or {}
    path = _bounded_upload_path(UPLOAD_DIR, str(body.get("filename") or ""))
    if not Path(path).exists():
        return _err("File not found in the upload directory.", 404)
    try:
        header = _sequencing.read_header(str(path))
        return jsonify({
            "format": _sequencing.detect_format(str(path)),
            "build": _sequencing.detect_build(header),
            "sample_count": header.sample_count,
            "liftover_available": _sequencing.liftover_available("GRCh38", "GRCh37"),
        })
    except Exception as exc:
        return _err(f"Could not inspect this file: {exc}")


# ---------------------------------------------------------------------------
# 3. Ancestry
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/ancestry", methods=["GET"])
def ancestry_view(pid: int):
    """Global ancestry with per-population marker coverage and intervals.

    A population the array cannot resolve is returned in ``not_resolvable`` with
    a null proportion, never as zero percent. Zero percent is a measurement.
    Unmeasurable is not, and reporting one as the other is how every incumbent
    turns a model artefact into an apparent fact about a person.
    """
    error, genotypes = _guard_genotypes(_ancestry, "ancestry", pid)
    if error:
        return error

    panel = request.args.get("panel") or "onekg_sgdp"
    mode = request.args.get("mode") or "projection"
    try:
        result = _ancestry.global_ancestry(genotypes, panel=panel, mode=mode)
    except Exception as exc:
        return _err(f"Ancestry estimation failed: {exc}", 500)
    result.setdefault("caveats", _ancestry.ancestry_caveats(panel=panel))
    return jsonify(result)


@api_v3.route("/api/profiles/<int:pid>/ancestry/painting", methods=["GET"])
def ancestry_painting(pid: int):
    """Chromosome painting segments, ready for the offline report renderer."""
    guard = _need(_ancestry, "ancestry")
    if guard:
        return guard
    missing = _profile_or_404(pid)
    if missing:
        return missing
    body_path = request.args.get("phased_vcf") or ""
    if not body_path:
        return jsonify(external.unavailable(
            "flare", "ancestry_local",
            detail=("Local ancestry needs a phased VCF. Run imputation first, "
                    "which phases as a side effect, then request painting "
                    "against its output."),
        ))
    try:
        local = _ancestry.local_ancestry(
            str(_bounded_upload_path(UPLOAD_DIR, body_path)),
            panel=request.args.get("panel") or "onekg_sgdp",
        )
    except Exception as exc:
        return _err(f"Local ancestry failed: {exc}", 500)
    if not local.get("available", True):
        return jsonify(local)
    return jsonify(_ancestry.chromosome_painting(local))


# ---------------------------------------------------------------------------
# 4. Haplogroups
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/haplogroups", methods=["GET"])
def haplogroups_view(pid: int):
    """Y and mtDNA haplogroups with an explicit resolution ceiling.

    The ceiling is the differentiator and it is computed, not decorative. A
    consumer array carries roughly 2,000 Y-SNPs against a tree with far more,
    so telling the user where their data stops is more honest than any
    incumbent manages.
    """
    error, merged = _guard_merged(_haplogroups, "haplogroups", pid)
    if error:
        return error

    genotypes = pipeline.build_genotype_map(merged.get("genotypes") or {})
    profile = db.get_profile(pid) or {}
    sex_hint = str(profile.get("sex") or "").strip().lower() or None
    try:
        mt_genotypes = _haplogroups.mt_positions_from_merged(merged)
    except Exception:
        mt_genotypes = None
    try:
        result = _haplogroups.analyse(
            genotypes, sex_hint=sex_hint, mt_genotypes=mt_genotypes,
            verified_only=str(request.args.get("verified_only", "")).lower()
            in ("1", "true", "yes"),
        )
    except Exception as exc:
        return _err(f"Haplogroup analysis failed: {exc}", 500)

    # The backbone tree bundled with v3.0 is provisional and says so. Surfacing
    # the unverified marker count on every response is the difference between a
    # caveat somebody wrote once and a caveat the user actually sees.
    try:
        result["unverified_markers"] = _haplogroups.unverified_markers()
    except Exception:
        pass
    result["tree"] = _haplogroups.tree_stamp()
    return jsonify(result)


# ---------------------------------------------------------------------------
# 5. Household genomics
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/household", methods=["GET"])
def household_view(pid: int):
    """Private IBD across the kits loaded on this machine, and nothing else.

    This is deliberately not a matching service. GEDmatch's profile count is a
    network effect, not a software moat, and a local tool cannot copy it. What a
    local tool CAN do is the family case, which needs no database and is immune
    to the opt-out circumvention documented at GEDmatch in 2023 and at
    MyHeritage in 2025. The response says so in ``scope`` so the limitation is
    never mistaken for a bug.
    """
    error, merged = _guard_merged(_relatedness, "relatedness", pid)
    if error:
        return error
    try:
        result = _relatedness.analyse_household(merged)
    except Exception as exc:
        return _err(f"Household analysis failed: {exc}", 500)
    result["scope"] = (
        "Only the DNA files loaded into this profile are compared. DNAInsight "
        "has no matching database and cannot find relatives you have not "
        "already uploaded. That is a design decision, not a limitation to be "
        "fixed."
    )
    return jsonify(result)


@api_v3.route("/api/profiles/<int:pid>/household/browser", methods=["GET"])
def household_browser(pid: int):
    """Segment data for one pair, shaped for the SVG chromosome browser."""
    guard = _need(_relatedness, "relatedness")
    if guard:
        return guard
    merged = _merged_for(pid)
    if merged is None:
        return _err("No DNA data for this profile.")
    label_a = request.args.get("a") or ""
    label_b = request.args.get("b") or ""
    try:
        household = _relatedness.analyse_household(merged)
    except Exception as exc:
        return _err(f"Household analysis failed: {exc}", 500)
    for pair in household.get("pairs") or []:
        if {str(pair.get("a")), str(pair.get("b"))} == {label_a, label_b}:
            return jsonify(_relatedness.chromosome_browser_data(pair))
    return _err("No comparison found for that pair of labels.", 404)


@api_v3.route("/api/profiles/<int:pid>/phasing", methods=["GET"])
def phasing_view(pid: int):
    """Which allele came from which parent, where the trio makes it resolvable."""
    guard = _need(_relatedness, "relatedness")
    if guard:
        return guard
    merged = _merged_for(pid)
    if merged is None:
        return _err("No DNA data for this profile.")
    try:
        return jsonify(_relatedness.phase_by_parents(merged))
    except Exception as exc:
        return _err(f"Phasing failed: {exc}", 500)


# ---------------------------------------------------------------------------
# 6. Imputation
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/imputation", methods=["POST"])
def imputation_start(pid: int):
    """Impute untyped genotypes, carrying DR2 as a first-class field.

    Runs in a background thread because whole-genome imputation is minutes.
    Every imputed call is hard-capped in magnitude below any typed call and the
    cap is written into the audit trail as a named step, because an unexplained
    score change is precisely what this project refuses to ship.
    """
    error, genotypes = _guard_genotypes(_imputation, "imputation", pid)
    if error:
        return error

    body = request.get_json(silent=True) or {}
    panel = str(body.get("panel") or "onekg_sgdp")
    threshold = float(body.get("dr2_threshold") or 0.8)

    key = _job("imputation", pid)
    with _lock:
        if _jobs.get(key, {}).get("running"):
            return _err("Imputation is already running for this profile.", 409)
        _jobs[key] = {"running": True, "done": False, "error": None, "result": None}

    def worker() -> None:
        try:
            result = _imputation.impute(genotypes, panel=panel,
                                        dr2_threshold=threshold)
            with _lock:
                _jobs[key] = {"running": False, "done": True,
                              "error": None, "result": result}
        except Exception as exc:
            with _lock:
                _jobs[key] = {"running": False, "done": True,
                              "error": str(exc), "result": None}

    threading.Thread(target=worker, daemon=True).start()
    return jsonify({"message": "Imputation started.", "profile_id": pid,
                    "panel": panel, "dr2_threshold": threshold,
                    "caveats": _imputation.build_caveats(panel=panel)})


@api_v3.route("/api/profiles/<int:pid>/imputation/status", methods=["GET"])
def imputation_status(pid: int):
    guard = _need(_imputation, "imputation")
    if guard:
        return guard
    with _lock:
        state = dict(_jobs.get(_job("imputation", pid))
                     or {"running": False, "done": False, "error": None})
    return jsonify(state)


@api_v3.route("/api/profiles/<int:pid>/imputation/coverage", methods=["GET"])
def imputation_coverage(pid: int):
    guard = _need(_imputation, "imputation")
    if guard:
        return guard
    with _lock:
        state = _jobs.get(_job("imputation", pid)) or {}
    result = state.get("result") or {}
    variants = result.get("variants") or []
    try:
        report = _imputation.coverage_report(variants)
    except Exception as exc:
        return _err(f"Coverage report failed: {exc}", 500)
    report["caveats"] = _imputation.build_caveats(report.get("coverage"))
    return jsonify(report)


@api_v3.route("/api/profiles/<int:pid>/imputation/safety", methods=["GET"])
def imputation_safety(pid: int):
    """Assert no imputed pathogenic call escaped without a quality figure.

    This is invariant 1 of the project, "do not alarm a non-carrier", restated
    for imputation. An imputed pathogenic call the user does not actually carry
    is the worst failure available here, so it is a checked guarantee with an
    endpoint rather than a comment somebody wrote once.
    """
    guard = _need(_imputation, "imputation")
    if guard:
        return guard
    findings = _findings_for(pid)
    violations = _imputation.assert_no_imputed_pathogenic_without_quality(
        findings, strict=False
    )
    return jsonify({"ok": not violations, "violations": violations,
                    "checked": len(findings)})


# ---------------------------------------------------------------------------
# 7. Pharmacogenomics: diplotypes and the prescription guard
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/pgx/diplotypes", methods=["GET"])
def diplotypes_view(pid: int):
    """Star-allele diplotypes with CPIC phenotype translation.

    Indeterminate is the default when the evidence is insufficient, never
    Normal. Defaulting to Normal is how this class of tool tells somebody they
    metabolise a drug fine when nobody actually checked.
    """
    error, genotypes = _guard_genotypes(_diplotype, "diplotype", pid)
    if error:
        return error

    results = []
    for gene in _diplotype.GENES:
        try:
            call = _diplotype.call_diplotype(gene, genotypes)
            call["phenotype"] = _diplotype.translate_phenotype(gene, call)
        except Exception as exc:
            call = {"gene": gene, "error": str(exc)}
        results.append(call)
    return jsonify({
        "diplotypes": results,
        "provisional_genes": sorted(_diplotype.PROVISIONAL_GENES),
        "unverified": _diplotype.unverified_entries(),
        "disclaimer": _diplotype.DISCLAIMER,
    })


@api_v3.route("/api/profiles/<int:pid>/pgx/prescription-guard", methods=["POST"])
def prescription_guard(pid: int):
    """Flag only the gene-drug pairs that apply to the medications supplied.

    Output is framed for a conversation with a prescriber. It never instructs
    anyone to start, stop or change a medication, and the module carries a test
    asserting no output string contains imperative dosing language.
    """
    guard = _need(_diplotype, "diplotype")
    if guard:
        return guard
    genotypes = _genotypes_for(pid)
    if genotypes is None:
        return _err("No DNA data for this profile.")
    body = request.get_json(silent=True) or {}
    medications = body.get("medications") or []
    if not isinstance(medications, list):
        return _err("medications must be a list of drug names.")

    diplotypes = {}
    for gene in _diplotype.GENES:
        try:
            diplotypes[gene] = _diplotype.call_diplotype(gene, genotypes)
        except Exception:
            continue
    try:
        result = _diplotype.prescription_guard(medications, diplotypes)
    except Exception as exc:
        return _err(f"Prescription guard failed: {exc}", 500)
    result["disclaimer"] = _diplotype.DISCLAIMER
    return jsonify(result)


# ---------------------------------------------------------------------------
# 8. Carrier screening
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/carrier", methods=["GET"])
def carrier_view(pid: int):
    """Carrier panel with residual risk.

    The wording is enforced in the module, not merely intended: no output path
    is allowed to emit the bare phrase "not a carrier". It is always "not a
    carrier for the N variants tested", because CFTR alone has over two
    thousand known pathogenic variants and a consumer array reads a few dozen of
    them. The residual risk arithmetic is the whole product here.
    """
    error, genotypes = _guard_genotypes(_carrier, "carrier", pid)
    if error:
        return error
    population = request.args.get("population") or ""
    try:
        report = _carrier.carrier_report(genotypes, population=population)
    except Exception as exc:
        return _err(f"Carrier report failed: {exc}", 500)
    report["unverified"] = _carrier.unverified_figures()
    report["disclaimer"] = _carrier.DISCLAIMER
    return jsonify(report)


@api_v3.route("/api/profiles/<int:pid>/carrier/joint", methods=["POST"])
def carrier_joint(pid: int):
    """Joint reproductive risk using a loaded partner file.

    The ``mate`` role already exists in the v2 source API, so this needs no new
    upload path. Results are ranges rather than single numbers, because the
    inputs are population estimates and false precision here would be worse than
    an honest interval.
    """
    guard = _need(_carrier, "carrier")
    if guard:
        return guard
    merged = _merged_for(pid)
    if merged is None:
        return _err("No DNA data for this profile.")
    body = request.get_json(silent=True) or {}
    gene = str(body.get("gene") or "").strip().upper()
    if not gene:
        return _err("gene is required.")

    comparison = (merged.get("comparison") or {})
    mate = comparison.get("mate")
    if not mate:
        return _err(
            "No partner file loaded. Add a source with role 'mate' first.", 409
        )

    genotypes_self = pipeline.build_genotype_map(merged.get("genotypes") or {})
    genotypes_mate = {
        rsid: (row.get("allele1", "N"), row.get("allele2", "N"))
        for rsid, row in mate.items()
    }
    population = str(body.get("population") or "")
    try:
        a = _carrier.carrier_status(gene, genotypes_self)
        b = _carrier.carrier_status(gene, genotypes_mate)
        risk_a = _carrier.residual_risk(gene, a.get("tested_negative", []), population)
        risk_b = _carrier.residual_risk(gene, b.get("tested_negative", []), population)
        joint = _carrier.joint_reproductive_risk(gene, risk_a, risk_b)
    except Exception as exc:
        return _err(f"Joint risk failed: {exc}", 500)
    return jsonify({"gene": gene, "self": a, "mate": b,
                    "residual_self": risk_a, "residual_mate": risk_b,
                    "joint": joint, "disclaimer": _carrier.DISCLAIMER})


@api_v3.route("/api/profiles/<int:pid>/acmg", methods=["GET"])
def acmg_view(pid: int):
    """ACMG secondary-findings coverage, which for an array is close to zero.

    Reporting that plainly is the point. A panel that covers a handful of
    positions in BRCA1 is not a BRCA1 test, and a report that implies otherwise
    is worse than no report.
    """
    guard = _need(_carrier, "carrier")
    if guard:
        return guard
    genotypes = _genotypes_for(pid)
    if genotypes is None:
        return _err("No DNA data for this profile.")
    try:
        return jsonify(_carrier.acmg_coverage_report(genotypes))
    except Exception as exc:
        return _err(f"ACMG coverage failed: {exc}", 500)


# ---------------------------------------------------------------------------
# 9. Reclassification ledger
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/snapshots", methods=["GET"])
def snapshots_view(pid: int):
    guard = _need(_ledger, "ledger")
    if guard:
        return guard
    _ledger.init_ledger()
    return jsonify({"snapshots": _ledger.list_snapshots(pid),
                    "latest": _ledger.latest_snapshot(pid)})


@api_v3.route("/api/profiles/<int:pid>/changes", methods=["GET"])
def changes_view(pid: int):
    """What changed for this person since the last scan.

    Not what changed in the databases. A VUS that became Pathogenic in ClinVar
    is a news item; a VUS that became Pathogenic AND that this person carries is
    the single highest-value event in personal genomics, and no consumer product
    surfaces it as an event.
    """
    guard = _need(_ledger, "ledger")
    if guard:
        return guard
    _ledger.init_ledger()
    since = request.args.get("since")
    limit = request.args.get("limit")
    try:
        return jsonify(_ledger.changes_for(
            pid, since=since, limit=int(limit) if limit else None
        ))
    except Exception as exc:
        return _err(f"Change query failed: {exc}", 500)


@api_v3.route("/api/profiles/<int:pid>/addendum", methods=["POST"])
def addendum_view(pid: int):
    """Generate a dated report addendum.

    Additive by construction. It never rewrites or supersedes the original
    report, because a report the user already printed and took to a clinician
    must stay exactly what it was.
    """
    guard = _need(_ledger, "ledger")
    if guard:
        return guard
    _ledger.init_ledger()
    body = request.get_json(silent=True) or {}
    try:
        payload = _ledger.addendum(pid, old_id=body.get("old_id"),
                                   new_id=body.get("new_id"))
    except Exception as exc:
        return _err(f"Addendum failed: {exc}", 500)

    if body.get("format") == "html":
        html = _ledger.render_addendum_html(payload)
        profile = db.get_profile(pid) or {}
        safe = secure_filename(str(profile.get("name") or f"profile_{pid}"))
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        path = REPORTS_DIR / f"{safe}_addendum_{stamp}.html"
        path.write_text(html, encoding="utf-8")
        try:
            rid = db.record_report(pid, "addendum", "html", str(path))
        except Exception:
            rid = None
        return jsonify({"report_id": rid, "path": str(path),
                        "addendum": payload})
    return jsonify(payload)


# ---------------------------------------------------------------------------
# 10. Provenance and the signed manifest
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/manifest", methods=["POST"])
def manifest_build(pid: int):
    """Emit a signed, reproducible manifest for this profile's current state.

    A clinician holding a DNAInsight report can ask which ClinVar release said
    this, and when. Without a manifest that question has no answer. With one,
    the report can be reproduced exactly or shown to have drifted.

    The signature proves the manifest was not altered after generation ON THIS
    MACHINE. It is an HMAC over a local key, not a public-key attestation, so it
    does not prove authorship to a third party. Overclaiming that would be worse
    than not signing at all.
    """
    guard = _need(_provenance, "provenance")
    if guard:
        return guard
    missing = _profile_or_404(pid)
    if missing:
        return missing
    _provenance.init_provenance()
    body = request.get_json(silent=True) or {}

    input_files = []
    try:
        conn = db.get_connection()
        rows = conn.execute(
            "SELECT filename FROM snp_uploads WHERE profile_id=?", (pid,)
        ).fetchall()
        conn.close()
        for row in rows:
            candidate = UPLOAD_DIR / str(row["filename"])
            if candidate.exists():
                input_files.append(str(candidate))
    except Exception:
        pass

    findings = _findings_for(pid)
    manifest = _provenance.build_manifest(
        profile_id=pid, findings=findings, input_files=input_files,
        report_type=str(body.get("report_type") or ""),
        extra={"scan_parameters": body.get("scan_parameters") or {},
               "app_version": APP_VERSION},
    )
    signed = _provenance.sign_manifest(manifest)
    return jsonify({"manifest": signed,
                    "text": _provenance.render_manifest_text(signed)})


@api_v3.route("/api/v3/verify-manifest", methods=["POST"])
def manifest_verify():
    """Verify a manifest. Returns a structured verdict, never a bare boolean."""
    guard = _need(_provenance, "provenance")
    if guard:
        return guard
    body = request.get_json(silent=True) or {}
    payload = body.get("manifest") if "manifest" in body else body
    try:
        return jsonify(_provenance.verify_manifest(payload))
    except Exception as exc:
        return jsonify({"ok": False, "error": str(exc),
                        "field": "manifest"}), 400


@api_v3.route("/api/profiles/<int:pid>/conflicts/sources", methods=["GET"])
def source_conflicts(pid: int):
    """Where ClinVar, the GWAS Catalog and CPIC disagree about this person.

    Disagreements are displayed, not resolved. That mirrors how pooled DNA files
    are handled when two arrays read the same position differently: both calls
    stay, because a disagreement is information about reliability that a silent
    reconciliation would destroy.
    """
    guard = _need(_provenance, "provenance")
    if guard:
        return guard
    out = []
    for finding in _findings_for(pid):
        try:
            conflicts = _provenance.detect_conflicts(finding)
        except Exception:
            continue
        if conflicts:
            out.append({"rsid": finding.get("rsid"),
                        "gene": finding.get("gene"),
                        "conflicts": conflicts})
    return jsonify({"count": len(out), "findings": out,
                    "policy": "Conflicts are surfaced, never resolved."})


# ---------------------------------------------------------------------------
# 11. Grounded local assistant
# ---------------------------------------------------------------------------

@api_v3.route("/api/profiles/<int:pid>/assistant", methods=["POST"])
def assistant_ask(pid: int):
    """Answer a question strictly from this profile's own findings.

    Refusal is the default, not the exception. Retrieval is limited to the local
    evidence store, every claim must cite a finding id that was actually in the
    context, and a response citing anything else is rejected and replaced with
    the refusal. Genotypes are stripped before anything leaves this process.

    The alternatives in this market are cloud services that require handing over
    the genome. This one talks to a model on loopback or it does not answer.
    """
    guard = _need(_assistant, "assistant")
    if guard:
        return guard
    missing = _profile_or_404(pid)
    if missing:
        return missing
    body = request.get_json(silent=True) or {}
    question = str(body.get("question") or "").strip()
    if not question:
        return _err("question is required.")
    findings = _findings_for(pid)
    try:
        result = _assistant.ask(question, findings,
                                model=body.get("model") or None)
    except Exception as exc:
        return _err(f"Assistant failed: {exc}", 500)
    return jsonify(result)


@api_v3.route("/api/v3/assistant/contract", methods=["GET"])
def assistant_contract():
    """Publish the grounding contract the assistant runs under.

    A safety rule nobody can read is a safety rule nobody can check.
    """
    guard = _need(_assistant, "assistant")
    if guard:
        return guard
    return jsonify({
        "contract": _assistant.GROUNDING_CONTRACT,
        "refusal": _assistant.REFUSAL,
        "redacted_fields": sorted(_assistant.REDACTED_FIELDS),
        "banned_phrases": sorted(_assistant.BANNED_ADVICE_PHRASES),
    })


# ---------------------------------------------------------------------------
# 12. Capability map
# ---------------------------------------------------------------------------

@api_v3.route("/api/v3/capabilities", methods=["GET"])
def capabilities_v3():
    """Everything the UI needs to decide which controls to render.

    Three separate axes, kept separate on purpose. ``subsystems`` is code that
    ships with DNAInsight. ``external`` is third-party tools the user installed.
    ``panels`` is reference data the user built. A control needs all three of
    its dependencies before it is worth showing, and collapsing them into one
    boolean would make an unbuilt panel look like a broken feature.
    """
    payload = {
        "version": APP_VERSION,
        "subsystems": pipeline.available_subsystems(),
        "external": external.capability_report(),
        "panels": {pid: external.panel_status(pid)["available"]
                   for pid in external.PANELS},
        "offline": {
            "core": True,
            "note": (
                "Every core capability runs with no network. Builders and tool "
                "installers may download, but only when you run them yourself."
            ),
        },
    }
    return jsonify(payload)
