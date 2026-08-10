"""
pipeline.py -- the v2 scan orchestrator.

Turns a set of uploaded raw DNA files into the complete finding set described by
docs/API_V2.md. This is the module that makes the v2 subsystems reachable:
without it, orientation, frequency, genosets, traits, prs, merge, scoring and
the SNPedia cache are seven independent libraries that never meet.

ORDER MATTERS, and this is why:

  1. merge      Pool every ``self`` file first, so everything downstream sees
                one genotype set. Conflicting calls are RETAINED, never voted on.
  2. bundled    Annotate against the curated offline reference.
  3. api        Optionally refine against MyVariant (rsIDs only, no genotypes).
  4. orientation Resolve strand BEFORE frequency and BEFORE scoring, because a
                flipped genotype looks up the wrong frequency and then scores
                off it. Fixing strand after scoring would be too late.
  5. frequency  Population genotype frequency and the rarity band.
  6. snpedia    Optional local cache overlay, only if the user opted in.
  7. scoring    Magnitude, repute and confidence LAST, because it consumes
                carrier status, frequency band, publication count and the
                ambiguity flag that the earlier stages produce.
  8. genosets, traits, prs  Independent entity types appended to the same list.

Nothing here talks to the network except the optional MyVariant step, which is
inherited unchanged from scanner.py and sends rsIDs only.
"""

from __future__ import annotations

from typing import Any, Callable, Iterable

from . import scanner
from . import merge as merge_mod
from . import orientation as orient_mod
from . import frequency as freq_mod
from . import scoring as scoring_mod

__all__ = [
    "PHASES", "DEFAULT_OPTIONS", "build_genotype_map", "snps_from_merged",
    "attach_provenance", "apply_orientation", "enrich_findings",
    "collect_genosets", "collect_traits", "collect_prs",
    "compute_qc", "compute_ranges", "compute_summary", "run_full_scan",
]

PHASES: tuple[str, ...] = (
    "merge", "bundled", "orientation", "frequency", "snpedia",
    "scoring", "genosets", "traits", "prs", "api", "writing", "complete",
)

DEFAULT_OPTIONS: dict[str, Any] = {
    "use_api": False,
    "population": "CEU",
    "include_genosets": True,
    "include_traits": True,
    "include_prs": True,
    "use_snpedia": False,
}

_NOCALL = {"", "N", "-", "--", "0", "00", "NN", "?"}


# ---------------------------------------------------------------------------
# Optional subsystems
#
# Every one of these is imported defensively. A fresh checkout that has not run
# the data builders yet must still be able to scan, just with fewer entity
# types, rather than failing at import time.
# ---------------------------------------------------------------------------

def _try(name: str):
    try:
        module = __import__(f"backend.{name}", fromlist=[name])
        return module
    except Exception:
        return None


_genosets = _try("genosets")
_traits = _try("traits")
_prs = _try("prs")
_snpedia = _try("snpedia")


def available_subsystems() -> dict[str, bool]:
    """Report which optional subsystems have both code and data present."""
    out = {
        "frequency": bool(freq_mod.load_frequencies()),
        "genosets": False,
        "traits": False,
        "prs": False,
        "snpedia": False,
    }
    if _genosets is not None:
        try:
            out["genosets"] = bool(_genosets.load_genosets())
        except Exception:
            out["genosets"] = False
    if _traits is not None:
        out["traits"] = bool(getattr(_traits, "TRAITS", None))
    if _prs is not None:
        try:
            out["prs"] = bool(_prs.load_models())
        except Exception:
            out["prs"] = False
    if _snpedia is not None:
        try:
            out["snpedia"] = bool(_snpedia.cache_status().get("available"))
        except Exception:
            out["snpedia"] = False

    # v3.0 subsystems.
    #
    # Two different kinds of availability are reported here and conflating them
    # would mislead the UI. A module flag means the CODE is present and works
    # offline with no external help. A capability flag from backend.external
    # means a THIRD-PARTY TOOL the user installed themselves is present and its
    # licence has been accepted. Ancestry has both: the code always loads, but
    # it can only produce an answer when fastmixture and a panel exist.
    for name in ("ledger", "provenance", "sequencing", "haplogroups",
                 "relatedness", "imputation", "ancestry", "diplotype",
                 "carrier", "assistant", "concordance"):
        out[name] = _try(name) is not None

    try:
        from . import external as _external
        # PREFIXED, and this is not cosmetic. Beagle's capability is literally
        # named "imputation" and Ollama's is "assistant", which are also the
        # names of the DNAInsight modules that drive them. An unprefixed update
        # silently overwrote the module flag with the tool flag, so a user
        # without Beagle installed was told DNAInsight's imputation module was
        # missing. That is the exact conflation this map exists to prevent, and
        # it was caught by the test that asserts the two namespaces stay apart.
        for name, ready in _external.capability_report().items():
            out[f"tool_{name}"] = ready
    except Exception:
        # A missing or broken external registry must leave the offline
        # subsystems reporting truthfully rather than blanking the whole map.
        pass
    return out


# ---------------------------------------------------------------------------
# Genotype plumbing
# ---------------------------------------------------------------------------

def build_genotype_map(merged_genotypes: dict) -> dict:
    """Return {rsid: (allele1, allele2)} for the genoset, trait and PRS engines.

    Those three modules take a plain mapping rather than the richer merged
    entry, and they expect a no-call to be a genuine no-call so that a missing
    or failed position evaluates false instead of being imputed.
    """
    out: dict[str, tuple[str, str]] = {}
    for rsid, entry in (merged_genotypes or {}).items():
        a1 = str(entry.get("allele1") or "").strip().upper()
        a2 = str(entry.get("allele2") or "").strip().upper()
        out[str(rsid).strip().lower()] = (a1, a2)
    return out


def snps_from_merged(merged_genotypes: dict) -> list[dict]:
    """Flatten merged genotypes back into the SNP-dict list the scanner takes."""
    snps: list[dict] = []
    for rsid, entry in (merged_genotypes or {}).items():
        snps.append({
            "rsid": rsid,
            "chromosome": entry.get("chromosome", ""),
            "position": entry.get("position", 0),
            "allele1": entry.get("allele1", ""),
            "allele2": entry.get("allele2", ""),
        })
    return snps


# ---------------------------------------------------------------------------
# Enrichment stages
# ---------------------------------------------------------------------------

def attach_provenance(finding: dict, merged: dict) -> dict:
    """Copy multi-file provenance from the merged genotype set onto a finding.

    Sets count, labels, conflict, calls and comparison. When two pooled files
    disagree, BOTH calls are carried through in ``calls`` and ``conflict`` is
    true. Nothing is reconciled, because silently choosing a winner would hide
    exactly the information the user needs in order to distrust that position.
    """
    rsid = str(finding.get("rsid") or "").strip().lower()
    entry = (merged.get("genotypes") or {}).get(rsid) or {}

    finding["count"] = int(entry.get("count") or 0)
    finding["labels"] = list(entry.get("labels") or [])
    finding["conflict"] = bool(entry.get("conflict"))
    finding["calls"] = list(entry.get("calls") or [])

    rows: list[dict] = []
    own = str(finding.get("genotype") or "").strip().upper()
    for role, table in (merged.get("comparison") or {}).items():
        other = table.get(rsid)
        if not other:
            continue
        their = str(other.get("genotype") or "").strip().upper()
        rows.append({
            "label": other.get("label", ""),
            "role": role,
            "genotype": their,
            "shared": bool(own) and bool(their) and sorted(own) == sorted(their),
        })
    finding["comparison"] = rows
    return finding


def apply_orientation(finding: dict) -> dict:
    """Resolve strand for one finding and record what happened.

    Sets orientation, stabilized_orientation, flipped, ambiguous and token.

    The bundled reference stores alleles on the GRCh37 plus strand, which is
    what 23andMe and AncestryDNA report, so no flip is needed for the offline
    path. The flags still have to be set, because a palindromic A/T or C/G
    heterozygote cannot be strand-verified at all and must be marked so it is
    not treated as settled fact downstream.
    """
    a1 = str(finding.get("allele1") or "").strip().upper()
    a2 = str(finding.get("allele2") or "").strip().upper()

    stabilized = str(finding.get("stabilized_orientation") or "").strip().lower()
    finding.setdefault("orientation", "")
    finding["stabilized_orientation"] = stabilized

    try:
        ambiguous = bool(orient_mod.is_ambiguous_pair(a1, a2))
    except Exception:
        ambiguous = False

    flipped = False
    if stabilized == "minus":
        try:
            result = orient_mod.orient_to_snpedia(
                a1, a2, stabilized_orientation="minus")
            flipped = bool(result.get("flipped"))
            finding["snpedia_token"] = result.get("token", "")
        except Exception:
            flipped = False

    finding["flipped"] = flipped
    finding["ambiguous"] = ambiguous
    if not finding.get("token"):
        pair = f"{a1};{a2}" if a1 and a2 else ""
        finding["token"] = f"({pair})" if pair else ""
    return finding


def _ensure_contract_keys(finding: dict) -> dict:
    """Guarantee every key docs/API_V2.md promises exists on a finding."""
    defaults = {
        "entity_type": "snp", "gene": "", "chromosome": "", "position": 0,
        "allele1": "", "allele2": "", "genotype": "", "token": "",
        "zygosity": "", "magnitude": None, "magnitude_source": "",
        "repute": "", "summary": "", "interpretation": "", "confidence": "none",
        "clinical_sig": "", "clinvar_sig_code": None, "review_status": "",
        "review_stars": 0, "cpic_level": "", "pgx_level": "", "evidence": "",
        "publications": 0, "conditions": "", "conditions_list": [],
        "sources": [], "orientation": "", "stabilized_orientation": "",
        "flipped": False, "ambiguous": False, "dubious": False,
        "variant_allele": "", "variant_copies": None, "carrier": None,
        "count": 0, "labels": [], "calls": [], "comparison": [],
        "probability": None, "mendelian_ok": None,
        "topics": [], "medicines": [],
        "criteria": "", "matched_rsids": [], "coverage": None,
        "percentile": None, "band": "", "reliable": None, "caveats": [],

        # Identity. A finding with no rsid at all is malformed, but the contract
        # still promises the key exists so template code cannot raise.
        "rsid": "",

        # Section 2.4, the frequency block. These used to be set only by
        # frequency.annotate, which enrich_findings calls for entity_type "snp"
        # alone. That left every genoset, trait and polygenic score in the
        # payload missing all eleven keys, so any consumer that read them
        # uniformly (a report template, an export, the table view) hit a
        # KeyError on exactly the rows that are not SNPs.
        "freq": None, "freq_population": "", "freq_band": "unknown",
        "freq_color": "", "freq_derived": False, "freq_method": "unavailable",
        "freq_flipped": False, "freq_ambiguous": False, "freq_queried": "",
        "gmaf": None, "minor_allele": "", "population_series": [],

        # Section 2.6. Only attach_provenance set this, and the appended
        # entity types never pass through it.
        "conflict": False,
    }

    # Two distinct jobs here, and they must not be conflated.
    #
    #   1. A MISSING key gets its default.
    #   2. A key that is present but explicitly None gets coerced, but only for
    #      the three fields where None would break string handling downstream.
    #
    # setdefault cannot do job 2, because the key already exists. The previous
    # version called setdefault for both cases, so the coercion silently never
    # happened and entity_type could reach the UI as None.
    COERCE_NONE = ("entity_type", "gene", "genotype", "rsid", "token",
                   "summary", "interpretation", "repute", "clinical_sig",
                   "conditions", "silo", "category")
    for key, value in defaults.items():
        if key not in finding:
            finding[key] = value
        elif finding[key] is None and key in COERCE_NONE:
            finding[key] = value

    if not finding.get("conditions_list") and finding.get("conditions"):
        parts = [p.strip() for p in str(finding["conditions"]).split(";")]
        finding["conditions_list"] = [p for p in parts if p]

    copies = finding.get("variant_copies")
    if isinstance(copies, int):
        finding["carrier"] = copies > 0
    return finding


def enrich_findings(findings: list[dict], merged: dict, *,
                    population: str = "CEU",
                    use_snpedia: bool = False,
                    progress_cb: Callable[[str, int, int], None] | None = None
                    ) -> list[dict]:
    """Run the full per-finding enrichment chain in the correct order."""
    total = len(findings)

    for index, finding in enumerate(findings):
        _ensure_contract_keys(finding)
        attach_provenance(finding, merged)
        apply_orientation(finding)

        # Frequency for positional entities only. A genoset, a trait or a
        # polygenic score has no single position, so it has no frequency.
        if finding.get("entity_type") == "snp":
            try:
                freq_mod.annotate(finding, population)
            except Exception:
                finding.setdefault("freq", None)
                finding.setdefault("freq_band", "unknown")

            if use_snpedia and _snpedia is not None:
                try:
                    _snpedia.annotate(finding)
                except Exception:
                    pass

        if progress_cb and total and index % 200 == 0:
            progress_cb("scoring", index, total)

    scoring_mod.score_all(findings, prefer_snpedia=use_snpedia)
    return findings


# ---------------------------------------------------------------------------
# Non-SNP entity types
# ---------------------------------------------------------------------------

def collect_genosets(genotypes: dict) -> dict:
    """Evaluate the genoset corpus, returning matched, unmatched, incomplete.

    A genoset whose required rsIDs were not all genotyped is INCOMPLETE, which
    means not testable on this array. That is a different fact from absent, and
    conflating the two is how a report tells someone they do not have something
    it never actually checked. A genoset that matched on the positions that were
    available is reported as matched with ``partial_coverage`` set, and is
    removed from the incomplete list so it is never double counted.
    """
    empty = {"matched": [], "unmatched": [], "incomplete": []}
    if _genosets is None:
        return empty
    try:
        corpus = _genosets.load_genosets()
        if not corpus:
            return empty
        verbose = _genosets.evaluate_all_verbose(genotypes, corpus)
    except Exception:
        return empty

    matched = list(verbose.get("matched") or [])
    matched_ids = {g.get("rsid") for g in matched}

    for g in matched:
        coverage = g.get("coverage")
        g["partial_coverage"] = bool(
            isinstance(coverage, (int, float)) and coverage < 1.0)

    incomplete = [g for g in (verbose.get("incomplete") or [])
                  if g.get("rsid") not in matched_ids]
    unmatched = [g for g in (verbose.get("unmatched") or [])
                 if g.get("rsid") not in matched_ids]

    return {"matched": matched, "unmatched": unmatched, "incomplete": incomplete}


def collect_traits(genotypes: dict) -> dict:
    """Predict traits and blood type, returning both plus renderable findings."""
    empty = {"traits": [], "blood_type": {}, "findings": []}
    if _traits is None:
        return empty
    try:
        called = _traits.predict_traits(genotypes)
        blood = _traits.predict_blood_type(genotypes)
        findings = _traits.to_findings(called, blood)
    except Exception:
        return empty
    # A trait is never good or bad. Enforce it here as well as in scoring, so a
    # future change to either module cannot start colouring traits.
    for f in findings:
        f["repute"] = ""
        f.setdefault("entity_type", "trait")
    return {"traits": called, "blood_type": blood, "findings": findings}


def collect_prs(genotypes: dict) -> dict:
    """Compute polygenic scores, returning results plus renderable findings."""
    empty = {"results": [], "findings": []}
    if _prs is None:
        return empty
    try:
        results = _prs.compute_all(genotypes)
        findings = _prs.to_findings(results)
    except Exception:
        return empty
    for f in findings:
        f["repute"] = ""
        f.setdefault("entity_type", "prs")
    return {"results": results, "findings": findings}


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

def compute_qc(findings: Iterable[dict], merged: dict) -> dict:
    """Summarise strand and call quality across the finding set."""
    total = flipped = ambiguous = no_call = unknown_orientation = dubious = 0
    non_carrier = 0
    for f in findings:
        if f.get("entity_type") != "snp":
            continue
        total += 1
        if f.get("flipped"):
            flipped += 1
        if f.get("ambiguous") or f.get("freq_ambiguous"):
            ambiguous += 1
        if f.get("zygosity") == "no_call":
            no_call += 1
        if not f.get("stabilized_orientation"):
            unknown_orientation += 1
        if f.get("dubious"):
            dubious += 1
        if f.get("carrier") is False:
            non_carrier += 1

    counts = merged.get("counts") or {}
    return {
        "total": total,
        "flipped": flipped,
        "ambiguous": ambiguous,
        "no_call": no_call,
        "unknown_orientation": unknown_orientation,
        "dubious": dubious,
        "non_carrier": non_carrier,
        "conflicts": int(counts.get("conflicts") or 0),
        "pooled_sources": int(counts.get("pooled_sources") or 0),
        "positions": int(counts.get("total_positions") or 0),
        "note": (
            "Palindromic sites (an A/T or C/G heterozygote) cannot be strand "
            "verified, because complementing the reported alleles yields the "
            "other observed allele. Those calls are capped and flagged rather "
            "than trusted."
        ),
    }


def compute_ranges(findings: Iterable[dict]) -> dict:
    """Data-derived slider bounds, so the UI never hard-codes them."""
    mags: list[float] = []
    pubs: list[int] = []
    freqs: list[float] = []
    for f in findings:
        m = f.get("magnitude")
        if isinstance(m, (int, float)):
            mags.append(float(m))
        p = f.get("publications")
        if isinstance(p, int):
            pubs.append(p)
        q = f.get("freq")
        if isinstance(q, (int, float)):
            freqs.append(float(q))
    return {
        "magnitude": [0.0, round(max(mags), 2) if mags else 10.0],
        "publications": [0, max(pubs) if pubs else 0],
        "frequency": [0.0, round(max(freqs), 2) if freqs else 100.0],
    }


def compute_summary(findings: Iterable[dict]) -> dict:
    """Count findings per silo, plus per entity type and per repute."""
    silos: dict[str, int] = {}
    entities: dict[str, int] = {}
    reputes = {"Good": 0, "Bad": 0, "unset": 0}
    for f in findings:
        silo = f.get("silo") or "informational"
        silos[silo] = silos.get(silo, 0) + 1
        ent = f.get("entity_type") or "snp"
        entities[ent] = entities.get(ent, 0) + 1
        rep = f.get("repute") or "unset"
        reputes[rep if rep in ("Good", "Bad") else "unset"] += 1
    return {"silos": silos, "entity_types": entities, "reputes": reputes}


# ---------------------------------------------------------------------------
# The scan
# ---------------------------------------------------------------------------

def run_full_scan(sources: list[dict], *,
                  use_api: bool = False,
                  population: str = "CEU",
                  include_genosets: bool = True,
                  include_traits: bool = True,
                  include_prs: bool = True,
                  use_snpedia: bool = False,
                  progress_cb: Callable[[str, int, int], None] | None = None
                  ) -> dict:
    """Run the complete v2 scan over one or more raw DNA sources.

    ``sources`` is a list of ``{label, role, provider, snps}`` dicts. Roles are
    normalised by ``backend.merge``; every ``self`` source is pooled, every other
    role becomes comparison rows, and ``ignore`` is dropped.

    Returns the full payload described in docs/API_V2.md section 3.3 and 3.4.
    """
    def report(phase: str, done: int = 0, total: int = 0) -> None:
        if progress_cb:
            progress_cb(phase, done, total)

    # 1. Pool the files.
    report("merge")
    merged = merge_mod.merge_sources(sources)
    genotype_entries = merged.get("genotypes") or {}
    genotypes = build_genotype_map(genotype_entries)
    snps = snps_from_merged(genotype_entries)

    # 2. Offline reference, then optionally the API.
    report("bundled", 0, len(snps))
    if use_api:
        def api_progress(done: int, total: int) -> None:
            report("api", done, total)
        findings = scanner.run_scan(snps, use_api=True, progress_cb=api_progress)
    else:
        findings = scanner.annotate_bundled(snps)

    # 3. Strand, frequency, optional SNPedia overlay, then scoring.
    report("orientation", 0, len(findings))
    enrich_findings(findings, merged, population=population,
                    use_snpedia=use_snpedia, progress_cb=progress_cb)

    # 4. The other entity types.
    genoset_result = {"matched": [], "unmatched": [], "incomplete": []}
    if include_genosets:
        report("genosets")
        genoset_result = collect_genosets(genotypes)
        for g in genoset_result["matched"]:
            _ensure_contract_keys(g)
            g["entity_type"] = "genoset"
            findings.append(g)

    trait_result = {"traits": [], "blood_type": {}, "findings": []}
    if include_traits:
        report("traits")
        trait_result = collect_traits(genotypes)
        for t in trait_result["findings"]:
            _ensure_contract_keys(t)
            findings.append(t)

    prs_result = {"results": [], "findings": []}
    if include_prs:
        report("prs")
        prs_result = collect_prs(genotypes)
        for p in prs_result["findings"]:
            _ensure_contract_keys(p)
            findings.append(p)

    # 5. Score the appended entity types too, then finalise.
    report("scoring", len(findings), len(findings))
    scoring_mod.score_all(findings, prefer_snpedia=use_snpedia)

    trio = {}
    try:
        trio = merge_mod.trio_annotate(merged)
    except Exception:
        trio = {"trio_available": False}

    report("complete", len(findings), len(findings))
    return {
        "findings": findings,
        "summary": compute_summary(findings),
        "ranges": compute_ranges(findings),
        "qc": compute_qc(findings, merged),
        "genosets": genoset_result,
        "traits": trait_result["traits"],
        "blood_type": trait_result["blood_type"],
        "prs": prs_result["results"],
        "conflicts": merged.get("conflicts") or [],
        "sources": merged.get("sources") or [],
        "counts": merged.get("counts") or {},
        "trio": trio,
        "population": population,
        "subsystems": available_subsystems(),
    }
