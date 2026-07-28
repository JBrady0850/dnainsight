"""vpipe.py -- end-to-end verification of the v2 pipeline on real upload files."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

from backend import pipeline, parsers
from sample import ensure_fixture

UPLOADS = ROOT / "uploads"

print("=" * 76)
print("DNAInsight v2 PIPELINE VERIFICATION")
print("=" * 76)

print("\nsubsystems available:", json.dumps(pipeline.available_subsystems()))

ensure_fixture(ROOT)
files = sorted(p for p in UPLOADS.glob("*.txt"))
print(f"\nupload files found: {[p.name for p in files]}")
if not files:
    print("NO UPLOAD FILES, cannot run end to end")
    sys.exit(1)

sources = []
for i, p in enumerate(files):
    parsed = parsers.parse_dna_file(str(p))
    sources.append({
        "label": p.name,
        "role": "self" if i == 0 else "other",
        "provider": parsed["provider"],
        "snps": parsed["snps"],
    })
    print(f"  {p.name:<34} provider={parsed['provider']:<12} snps={parsed['snp_count']:,}")

phases = []
def progress(phase, done=0, total=0):
    if not phases or phases[-1][0] != phase:
        phases.append((phase, done, total))

result = pipeline.run_full_scan(
    sources,
    use_api=False,
    population="CEU",
    include_genosets=True,
    include_traits=True,
    include_prs=True,
    use_snpedia=False,
    progress_cb=progress,
)

print("\nphases reported:", [p[0] for p in phases])

f = result["findings"]
print(f"\nTOTAL FINDINGS: {len(f)}")
print("summary:", json.dumps(result["summary"], indent=2))
print("ranges :", json.dumps(result["ranges"]))
print("counts :", json.dumps(result["counts"], indent=2))
print("qc     :", json.dumps({k: v for k, v in result["qc"].items() if k != "note"}, indent=2))

print("\ngenosets: matched=%d unmatched=%d incomplete=%d" % (
    len(result["genosets"]["matched"]),
    len(result["genosets"]["unmatched"]),
    len(result["genosets"]["incomplete"])))
overlap = ({g["rsid"] for g in result["genosets"]["matched"]}
           & {g["rsid"] for g in result["genosets"]["incomplete"]})
print("matched/incomplete overlap (must be empty):", overlap or "none")

print("traits called:", len(result["traits"]))
bt = result["blood_type"] or {}
print("blood type   :", bt.get("blood_type"), "confidence", bt.get("confidence"))
print("prs models   :", len(result["prs"]))
print("conflicts    :", len(result["conflicts"]))
print("trio         :", json.dumps({k: v for k, v in (result["trio"] or {}).items()
                                    if k != "note"})[:200])

# Contract compliance: every promised key present on every finding.
REQUIRED = [
    "rsid", "entity_type", "gene", "chromosome", "position", "allele1",
    "allele2", "genotype", "token", "zygosity", "magnitude",
    "magnitude_source", "repute", "summary", "interpretation", "confidence",
    "clinical_sig", "clinvar_sig_code", "review_status", "review_stars",
    "cpic_level", "pgx_level", "evidence", "publications", "conditions",
    "conditions_list", "sources", "orientation", "stabilized_orientation",
    "flipped", "ambiguous", "dubious", "variant_allele", "variant_copies",
    "carrier", "count", "labels", "calls", "comparison", "topics", "medicines",
    "criteria", "matched_rsids", "coverage",
]
missing = {}
for item in f:
    for key in REQUIRED:
        if key not in item:
            missing.setdefault(key, 0)
            missing[key] += 1
print("\nCONTRACT KEY CHECK:", "all present" if not missing else missing)

# Repute discipline.
bad_repute = [x["rsid"] for x in f
              if x.get("entity_type") in ("trait", "prs") and x.get("repute")]
print("traits/prs with a repute (must be none):", bad_repute or "none")

# Magnitude sanity.
mags = [x["magnitude"] for x in f if isinstance(x.get("magnitude"), (int, float))]
print(f"magnitude: n={len(mags)} min={min(mags):.2f} max={max(mags):.2f}"
      if mags else "magnitude: none")
out_of_range = [x["rsid"] for x in f
                if isinstance(x.get("magnitude"), (int, float))
                and not (0.0 <= x["magnitude"] <= 10.0)]
print("out of 0-10 range (must be none):", out_of_range or "none")

nocall_scored = [x["rsid"] for x in f
                 if x.get("zygosity") == "no_call" and (x.get("magnitude") or 0) > 0]
print("no-calls with magnitude > 0 (must be none):", nocall_scored or "none")

print("\n--- TOP 12 BY MAGNITUDE ---")
ranked = sorted(f, key=lambda x: (x.get("magnitude") or 1.0), reverse=True)[:12]
for x in ranked:
    print(f"  {x['magnitude']:>5} {x.get('repute') or '-':<5} "
          f"{x['entity_type']:<8} {x['rsid']:<14} {(x.get('gene') or ''):<10} "
          f"{(x.get('summary') or x.get('interpretation') or '')[:58]}")

print("\n--- A SCORED SNP WITH ITS AUDIT TRAIL ---")
sample = next((x for x in ranked if x["entity_type"] == "snp"), None)
if sample:
    for k in ("rsid", "gene", "genotype", "zygosity", "magnitude",
              "magnitude_base", "magnitude_source", "repute", "confidence",
              "evidence", "review_stars", "clinvar_sig_code", "freq",
              "freq_band", "freq_flipped", "freq_ambiguous", "carrier",
              "variant_copies", "count", "labels", "dubious"):
        print(f"  {k:<20} {sample.get(k)!r}")
    print("  magnitude_factors:")
    for line in sample.get("magnitude_factors", []):
        print("     -", line)

print("\n" + "=" * 76)
ok = (not missing) and (not bad_repute) and (not out_of_range) \
     and (not nocall_scored) and (not overlap) and len(f) > 0
print("PIPELINE OK" if ok else "PIPELINE PROBLEMS FOUND")
sys.exit(0 if ok else 1)
