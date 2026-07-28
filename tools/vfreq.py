"""vfreq.py -- verify strand-tolerant frequency lookups after the patch."""
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

from backend import frequency as f

print("=" * 74)
print("STRAND-TOLERANT FREQUENCY VERIFICATION")
print("=" * 74)

# rs1801133: array reports C/T, Ensembl stores G/A. This is the regression case.
print("\n-- rs1801133 (MTHFR C677T): array C/T vs table G/A --")
print("observed in table  :", sorted(f.observed_alleles("rs1801133")))
r = f.resolve_strand("rs1801133", "C", "T")
print("resolve_strand     :", r)
d = f.genotype_frequency_detail("rs1801133", "C", "T", "CEU")
print("CEU C;T detail     :", d)
print("CEU allele T       :", f.allele_frequency("rs1801133", "T", "CEU"))
print("CEU allele T literal:", f.allele_frequency("rs1801133", "T", "CEU", strand_tolerant=False))
print("gmaf               :", f.gmaf("rs1801133"))
print("aggregate MAX      :", f.aggregate_frequency("rs1801133", "C", "T", "MAX"))
print("aggregate MIN      :", f.aggregate_frequency("rs1801133", "C", "T", "MIN"))
print("band / colour      :", f.rarity_band(d["frequency"]), f.rarity_color(d["frequency"]))

print("\n-- population series (first 8) --")
for s in f.population_series("rs1801133", "C", "T")[:8]:
    print(f"   {s['code']:<5} {str(s.get('frequency')):>7}   yours={s.get('yours')}   {s.get('label','')[:40]}")

print("\n-- annotate() on a bare finding --")
fin = f.annotate({"rsid": "rs1801133", "allele1": "C", "allele2": "T"})
for k in sorted(k for k in fin if k.startswith("freq") or k in ("gmaf", "minor_allele")):
    v = fin[k]
    if isinstance(v, list):
        v = f"<{len(v)} populations>"
    print(f"   {k:<20} {v}")

print("\n-- plus-strand control: rs671 (ALDH2), array G/A --")
print("observed           :", sorted(f.observed_alleles("rs671")))
print("resolve            :", f.resolve_strand("rs671", "G", "A"))
print("CEU G;A            :", f.genotype_frequency_detail("rs671", "G", "A", "CEU"))
print("JPT G;A            :", f.genotype_frequency_detail("rs671", "G", "A", "JPT"))

print("\n-- palindromic case: is it flagged rather than guessed? --")
for rs in ("rs9939609", "rs4988235", "rs6025", "rs4680"):
    obs = sorted(f.observed_alleles(rs))
    if len(obs) == 2 and f.is_palindromic(obs[0], obs[1]):
        rr = f.resolve_strand(rs, obs[0], obs[1])
        print(f"   {rs} observed={obs} PALINDROMIC ambiguous={rr['ambiguous']} flipped={rr['flipped']}")
    else:
        print(f"   {rs} observed={obs} not palindromic")

print("\n-- coverage across the whole bundled set --")
sys.path.insert(0, str(ROOT / "data"))
from data.build_reference import REFERENCE

resolved = flipped = ambiguous = missing = 0
for row in REFERENCE:
    rs = row[0]
    obs = f.observed_alleles(rs)
    if not obs:
        missing += 1
        continue
    obs_l = sorted(obs)
    a1 = obs_l[0]
    a2 = obs_l[-1] if len(obs_l) > 1 else obs_l[0]
    rr = f.resolve_strand(rs, a1, a2)
    if rr["resolved"]:
        resolved += 1
    if rr["flipped"]:
        flipped += 1
    if rr["ambiguous"]:
        ambiguous += 1

print(f"   bundled rsIDs      : {len(REFERENCE)}")
print(f"   with frequency data: {len(REFERENCE) - missing}")
print(f"   no data            : {missing}")
print(f"   palindromic sites  : {ambiguous}  (irreducibly ambiguous, must be flagged in UI)")
print(f"   coverage_report    : {f.coverage_report()}")

print("\n" + "=" * 74)
ok = d["frequency"] is not None and r["flipped"] is True
print("REGRESSION FIXED" if ok else "STILL BROKEN")
sys.exit(0 if ok else 1)
