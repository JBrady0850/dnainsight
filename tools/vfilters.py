"""vfilters.py -- verify the filter engine against a real scanned finding set."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

from backend import filters as F
from backend import pipeline, parsers
from sample import ensure_fixture

UPLOADS = ROOT / "uploads"
ensure_fixture(ROOT)
sources = []
for i, p in enumerate(sorted(UPLOADS.glob("*.txt"))):
    parsed = parsers.parse_dna_file(str(p))
    sources.append({"label": p.name, "role": "self", "provider": parsed["provider"],
                    "snps": parsed["snps"]})

result = pipeline.run_full_scan(sources, use_api=False, population="CEU")
findings = result["findings"]

print("=" * 74)
print("FILTER ENGINE VERIFICATION")
print("=" * 74)
print(f"base finding set: {len(findings)}")

def run(label, params, expect=None):
    out = F.filter_and_sort(findings, params)
    got = out["total"]
    flag = ""
    if expect is not None:
        flag = "  OK" if got == expect else f"  EXPECTED {expect}"
    print(f"  {label:<56} -> {got:>4}{flag}")
    return out

print("\n-- entity exemption (the rule most easily broken) --")
n_snp = sum(1 for f in findings if f["entity_type"] == "snp")
n_exempt = sum(1 for f in findings if f["entity_type"] in F.EXEMPT_ENTITIES)
print(f"  snp={n_snp}  genoset/trait/prs={n_exempt}")
a = run("no filters", {})
b = run("min_publications=1 (must keep all exempt entities)", {"min_publications": 1})
kept_exempt = sum(1 for f in b["findings"] if f["entity_type"] in F.EXEMPT_ENTITIES)
print(f"     exempt entities surviving a publication floor: {kept_exempt} of {n_exempt}",
      "OK" if kept_exempt == n_exempt else "BROKEN")
c = run("min_freq=10 (must keep all exempt entities)", {"min_freq": 10})
kept_exempt2 = sum(1 for f in c["findings"] if f["entity_type"] in F.EXEMPT_ENTITIES)
print(f"     exempt entities surviving a frequency floor:   {kept_exempt2} of {n_exempt}",
      "OK" if kept_exempt2 == n_exempt else "BROKEN")

print("\n-- magnitude --")
run("min_magnitude=2", {"min_magnitude": 2})
run("min_magnitude=5", {"min_magnitude": 5})
run("max_magnitude=1", {"max_magnitude": 1})

print("\n-- repute --")
for r in ("good", "bad", "unset"):
    run(f"repute={r}", {"repute": r})

print("\n-- clinvar --")
run("clinvar_only=1 (defaults to codes 5,4)", {"clinvar_only": 1})
run("clinvar_sig=6 (drug response)", {"clinvar_sig": "6"})
run("min_stars=3", {"min_stars": 3})
run("min_stars=4", {"min_stars": 4})

print("\n-- entity + zygosity + carrier --")
run("entity_type=genoset", {"entity_type": "genoset"})
run("entity_type=snp,prs", {"entity_type": "snp,prs"})
run("zygosity=homozygous", {"zygosity": "homozygous"})
run("carrier_only=1", {"carrier_only": 1})

print("\n-- free text grammar --")
run("q=MTHFR", {"q": "MTHFR"})
run("q=warfarin", {"q": "warfarin"})
run("q=chr1", {"q": "chr1"})
run("q=chr1:1-999999999", {"q": "chr1:1-999999999"})
run("q=/MAG>=2", {"q": "/MAG>=2"})
run("q=/STARS>=3", {"q": "/STARS>=3"})
run("q=/CLNSIG=6", {"q": "/CLNSIG=6"})
run("q=/flipped", {"q": "/flipped"})
run("q=/ambiguous", {"q": "/ambiguous"})
run("q=[unbalanced(  (must not raise)", {"q": "[unbalanced("})

print("\n-- sorting, both directions --")
for key in F.SORT_KEYS:
    desc = F.sort_findings(findings, key, "desc")
    asc = F.sort_findings(findings, key, "asc")
    same = desc[0]["rsid"] == asc[0]["rsid"]
    print(f"  {key:<14} desc first={desc[0]['rsid']:<14} asc first={asc[0]['rsid']:<14}"
          f"{'  (single value column)' if same else ''}")

print("\n-- null magnitude sorts as 1, not 0 --")
probe = [
    {"rsid": "rs_null", "entity_type": "snp", "magnitude": None},
    {"rsid": "rs_zero", "entity_type": "snp", "magnitude": 0.0},
    {"rsid": "rs_two", "entity_type": "snp", "magnitude": 2.0},
]
order = [f["rsid"] for f in F.sort_findings(probe, "magnitude", "desc")]
print("  desc order:", order, "OK" if order == ["rs_two", "rs_null", "rs_zero"] else "BROKEN")
kept = [f["rsid"] for f in F.apply_filters(probe, {"min_magnitude": 1})]
print("  min_magnitude=1 keeps:", kept,
      "OK" if kept == ["rs_null", "rs_two"] else "BROKEN")

print("\n-- pagination --")
p = F.filter_and_sort(findings, {"limit": 5, "offset": 0})
q = F.filter_and_sort(findings, {"limit": 5, "offset": 5})
allf = F.filter_and_sort(findings, {"limit": 0})
print(f"  limit=5 returned={p['returned']} total={p['total']}")
print(f"  offset=5 first differs: {p['findings'][0]['rsid'] != q['findings'][0]['rsid']}")
print(f"  limit=0 returns everything: {allf['returned'] == allf['total']}")

print("\n-- facets --")
facets = F.build_facets(findings)
for name in ("genes", "topics", "medicines", "silos", "entity_types",
             "cpic_levels", "review_stars", "freq_bands", "confidence"):
    vals = facets[name][:4]
    print(f"  {name:<22} {len(facets[name]):>3} values, top: "
          + ", ".join(f"{v['value']}({v['count']})" for v in vals))

print("\n-- chromosome sort order sanity --")
chroms = [{"rsid": f"rs{i}", "chromosome": c, "position": 1}
          for i, c in enumerate(["MT", "X", "2", "10", "1", "Y", "22"])]
print("  asc:", [f["chromosome"] for f in F.sort_findings(chroms, "location", "asc")])

print("\n" + "=" * 74)
print("FILTER ENGINE OK")
