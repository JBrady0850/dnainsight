"""vcomplete.py -- every file the v2 contract and docs promise actually exists."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()

REQUIRED = {
    "app.py": "entry point",
    "requirements.txt": "runtime dependencies",
    "requirements-dev.txt": "test dependencies",
    "README.md": "user documentation",
    "CHANGELOG.md": "release history",
    "LICENSE": "MIT licence",
    ".gitignore": "hazard exclusions",
    "backend/__init__.py": "version single source",
    "backend/parsers.py": "raw file parsing",
    "backend/database.py": "SQLite layer",
    "backend/scanner.py": "annotation engine",
    "backend/orientation.py": "strand engine",
    "backend/genosets.py": "genoset criteria engine",
    "backend/frequency.py": "population frequency",
    "backend/prs.py": "polygenic scores",
    "backend/merge.py": "multi-file pooling",
    "backend/traits.py": "traits and blood type",
    "backend/snpedia.py": "opt-in local cache",
    "backend/scoring.py": "magnitude and repute",
    "backend/pipeline.py": "scan orchestrator",
    "backend/filters.py": "filter engine",
    "backend/routes.py": "v1 endpoints",
    "backend/routes_v2.py": "v2 endpoints",
    "backend/updater.py": "reference refresh",
    "backend/genetic_report.py": "report for the user",
    "backend/doctor_report.py": "report for a clinician",
    "backend/interactive_report.py": "self-contained offline report",
    "frontend/index.html": "single page app",
    "data/build_reference.py": "curated table builder",
    "data/evidence_overlay.py": "evidence layer",
    "data/build_genosets.py": "genoset corpus builder",
    "data/build_frequencies.py": "frequency builder",
    "data/build_prs.py": "polygenic model builder",
    "data/build_full_reference.py": "optional Tier 2 builder",
    "data/DATA_SOURCES.md": "per-source licence record",
    "data/snp_reference.json": "bundled reference",
    "data/genosets.json": "genoset corpus",
    "data/prs_models.json": "polygenic models",
    "data/frequencies.json": "population frequencies",
    "docs/API_V2.md": "the authoritative contract",
    "CONTRIBUTING.md": "contributor guide and architecture notes",
    "tools/golive.py": "release gate",
    "tests/test_parsers.py": "parser tests",
    "tests/test_database.py": "database path resolution tests",
    "tests/test_uploads.py": "upload path bounding tests",
    "tests/test_snpedia_ratelimit.py": "rate limiter tests",
    "tests/test_scanner.py": "scanner tests",
    "tests/test_reference.py": "reference tests",
    "tests/test_orientation.py": "orientation tests",
    "tests/test_genosets.py": "genoset tests",
    "tests/test_frequency.py": "frequency tests",
    "tests/test_prs.py": "polygenic tests",
    "tests/test_merge.py": "merge tests",
    "tests/test_traits.py": "trait tests",
    "tests/test_snpedia.py": "snpedia tests",
    "tests/test_scoring.py": "scoring tests",
    "tests/test_filters.py": "filter tests",
    "tests/test_pipeline.py": "pipeline tests",
}

print("=" * 74)
print("COMPLETENESS")
print("=" * 74)

missing = []
total_lines = 0
for rel, what in sorted(REQUIRED.items()):
    p = ROOT / rel
    if p.exists():
        try:
            n = len(p.read_text(encoding="utf-8", errors="replace").splitlines())
        except OSError:
            n = 0
        total_lines += n
        print(f"  [ ok ] {rel:<36} {n:>6} lines   {what}")
    else:
        missing.append(rel)
        print(f"  [FAIL] {rel:<36} {'':>6}          MISSING ({what})")

# The frontend must actually be the v2 build, not the v1 file.
fe = ROOT / "frontend" / "index.html"
print()
if fe.exists():
    body = fe.read_text(encoding="utf-8", errors="replace")
    markers = ["findings/v2", "SLIDER_DEFS", "renderGenosets", "renderQcBanner",
               "CLINVAR_CODES", "onPopulationChange"]
    hit = [m for m in markers if m in body]
    if len(hit) == len(markers):
        print(f"  [ ok ] frontend is the v2 build "
              f"({len(body.splitlines())} lines, {len(hit)}/{len(markers)} markers)")
    else:
        missing.append("frontend/index.html is not the v2 build")
        print(f"  [FAIL] frontend has only {len(hit)}/{len(markers)} v2 markers, "
              "it may still be the v1 file")
    if "localStorage" in body.replace("NO localStorage", ""):
        missing.append("frontend uses localStorage")
        print("  [FAIL] frontend appears to use localStorage, which is prohibited")

print(f"\n  {len(REQUIRED) - len(missing)} of {len(REQUIRED)} required files present, "
      f"{total_lines:,} lines total")
print("=" * 74)
print("COMPLETE" if not missing else f"MISSING: {len(missing)}")
for m in missing:
    print("  -", m)
sys.exit(0 if not missing else 1)
