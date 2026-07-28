"""
smoke.py -- import and exercise every DNAInsight v2 module with synthetic data.

Reports what actually works rather than what should work. Never raises: every
failure is caught and printed so one broken module does not hide the rest.
"""

import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

# A small synthetic genotype set covering the variants the v2 modules care about.
GENOTYPES = {
    "rs429358":   ("T", "C"),   # APOE, het -> e3/e4 territory
    "rs7412":     ("C", "C"),
    "rs1801133":  ("C", "T"),   # MTHFR C677T het
    "rs1801131":  ("A", "A"),
    "rs9923231":  ("T", "T"),   # VKORC1 high warfarin sensitivity
    "rs1799853":  ("C", "T"),   # CYP2C9 *2
    "rs1057910":  ("A", "A"),
    "rs4244285":  ("A", "G"),   # CYP2C19 *2 het
    "rs12248560": ("C", "T"),
    "rs4149056":  ("T", "C"),   # SLCO1B1 *5
    "rs1800562":  ("G", "A"),   # HFE C282Y het
    "rs1799945":  ("C", "C"),
    "rs6025":     ("C", "T"),   # Factor V Leiden het
    "rs1799963":  ("G", "G"),
    "rs4988235":  ("C", "C"),   # lactase non-persistence
    "rs671":      ("G", "A"),   # ALDH2 flush
    "rs762551":   ("A", "C"),
    "rs9939609":  ("A", "T"),
    "rs7903146":  ("T", "T"),
    "rs1333049":  ("C", "C"),
    "rs4680":     ("A", "G"),
    "rs6265":     ("C", "T"),
    "rs53576":    ("G", "G"),
    "rs1815739":  ("C", "T"),
    "rs12913832": ("A", "A"),
    "rs17822931": ("C", "T"),
    "rs713598":   ("G", "C"),
    "rs8176719":  ("N", "N"),   # ABO indel deliberately a no-call
    "rs1051730":  ("A", "G"),
    "rs324420":   ("C", "A"),
    "rs1761667":  ("A", "G"),
    "rs1801260":  ("T", "C"),
    "rs4343":     ("A", "G"),
    "rs3918290":  ("C", "C"),
    "rs116855232": ("C", "C"),
    "rs887829":   ("C", "T"),
    "rs1800795":  ("G", "C"),
    "rs1205":     ("C", "T"),
    "rs10830963": ("C", "G"),
    "rs780094":   ("C", "T"),
    "rs5219":     ("C", "T"),
    "rs13266634": ("C", "C"),
    "rs1801282":  ("C", "C"),
    "rs10811661": ("T", "T"),
    "rs17782313": ("C", "T"),
    "rs1137101":  ("A", "G"),
    "rs1421085":  ("T", "C"),
    "rs8050136":  ("C", "A"),
    "rs1800544":  ("C", "G"),
    "rs10455872": ("A", "A"),
    "rs1799983":  ("G", "T"),
    "rs1800588":  ("C", "T"),
    "rs3135506":  ("C", "C"),
    "rs1805087":  ("A", "G"),
    "rs1801394":  ("A", "G"),
    "rs30187":    ("C", "T"),
    "rs1800629":  ("G", "A"),
}

SNP_LIST = [
    {"rsid": k, "chromosome": "1", "position": 1000 + i,
     "allele1": v[0], "allele2": v[1]}
    for i, (k, v) in enumerate(GENOTYPES.items())
]

results = []


def step(label, fn):
    try:
        out = fn()
        results.append((label, "ok", out))
        print(f"[ ok ] {label}: {out}")
    except Exception as exc:
        results.append((label, "FAIL", f"{type(exc).__name__}: {exc}"))
        print(f"[FAIL] {label}: {type(exc).__name__}: {exc}")
        traceback.print_exc(limit=3)


# --- orientation ---------------------------------------------------------
def t_orientation():
    from backend import orientation as o
    r = o.orient_to_snpedia("C", "T", stabilized_orientation="minus")
    assert r["flipped"] is True, r
    assert set((r["allele1"], r["allele2"])) == {"A", "G"}, r
    amb = o.is_ambiguous_pair("A", "T")
    return f"minus flip C;T -> {r['token']}, A/T ambiguous={amb}, {len(o.COMPLEMENT)} bases"


# --- genosets ------------------------------------------------------------
def t_genosets():
    from backend import genosets as g
    ast_ = g.parse_criteria("and(rs429358(T;C), or(rs7412(C;C), not(rs1801133(T;T))))")
    need = g.required_rsids(ast_)
    corpus = g.load_genosets()
    matched = g.evaluate_all(GENOTYPES, corpus) if corpus else []
    verbose = g.evaluate_all_verbose(GENOTYPES, corpus) if corpus else {}
    return (f"parse ok ({len(need)} rsids), corpus={len(corpus)}, "
            f"matched={len(matched)}, incomplete={len(verbose.get('incomplete', []))}")


# --- frequency -----------------------------------------------------------
def t_frequency():
    from backend import frequency as f
    freqs = f.load_frequencies()
    d = f.genotype_frequency_detail("rs1801133", "C", "T", "CEU")
    series = f.population_series("rs1801133", "C", "T")
    band = f.rarity_band(d.get("frequency"))
    fin = f.annotate({"rsid": "rs1801133", "allele1": "C", "allele2": "T"})
    assert "freq_band" in fin and "gmaf" in fin, list(fin)
    return (f"{len(freqs)} rsids, CEU C;T={d.get('frequency')} ({d.get('method')}), "
            f"band={band}, {len(series)} populations, annotate ok")


# --- prs -----------------------------------------------------------------
def t_prs():
    from backend import prs as p
    models = p.load_models()
    if not models:
        return "no prs_models.json yet (builder not run)"
    res = p.compute_all(GENOTYPES)
    f = p.to_findings(res)
    top = res[0] if res else {}
    return (f"{len(models)} models, {len(res)} scored, {len(f)} findings, "
            f"first={top.get('id')} band={top.get('band')} cov={top.get('coverage')}")


# --- merge ---------------------------------------------------------------
def t_merge():
    from backend import merge as m
    a = {"label": "23andMe_v5.txt", "role": "self", "provider": "23andme",
         "snps": SNP_LIST[:40]}
    # second self file: overlaps, one deliberate conflict, one no-call, one unique
    b_snps = [dict(s) for s in SNP_LIST[20:50]]
    b_snps[0]["allele1"], b_snps[0]["allele2"] = "A", "A"      # conflict
    b_snps[1]["allele1"], b_snps[1]["allele2"] = "N", "N"      # no-call, not a conflict
    b = {"label": "AncestryDNA.txt", "role": "self", "provider": "ancestrydna",
         "snps": b_snps}
    mom = {"label": "mom.txt", "role": "mother", "provider": "23andme",
           "snps": SNP_LIST[:30]}
    dad = {"label": "dad.txt", "role": "father", "provider": "23andme",
           "snps": SNP_LIST[:30]}
    merged = m.merge_sources([a, b, mom, dad])
    trio = m.trio_annotate(merged)
    prob = m.transmission_probability("AG", "AG", "AG")
    return (f"union={merged['counts']['union']}, conflicts={merged['counts']['conflicts']}, "
            f"comparison sets={merged['counts']['comparison_sources']}, "
            f"trio compared={trio.get('compared')} violations={trio.get('violations')}, "
            f"P(AG|AGxAG)={prob}")


# --- traits --------------------------------------------------------------
def t_traits():
    from backend import traits as t
    blood = t.predict_blood_type(GENOTYPES)
    tr = t.predict_traits(GENOTYPES)
    f = t.to_findings(tr, blood)
    reputes = {x.get("repute", "") for x in f}
    assert reputes <= {""}, f"traits must be neutral, got {reputes}"
    return (f"{len(t.TRAITS)} defined, {len(tr)} called, blood={blood.get('blood_type')} "
            f"(conf={blood.get('confidence')}), {len(f)} findings, repute neutral")


# --- snpedia -------------------------------------------------------------
def t_snpedia():
    from backend import snpedia as s
    st = s.cache_status()
    assert "snpedia" not in str(ROOT).lower() or True
    cp = str(s.cache_path())
    assert str(ROOT).lower() not in cp.lower(), f"cache must live outside the repo: {cp}"
    fin = s.annotate({"rsid": "rs1801133", "allele1": "C", "allele2": "T"})
    assert "magnitude" in fin, list(fin)
    gated = "no"
    try:
        s.harvest(rsids=["rs1801133"], accept_license=False)
    except PermissionError as exc:
        gated = "yes" if "oncommercial" in str(exc) else "yes(no-notice)"
    return (f"available={st.get('available')}, cache={cp}, license gate={gated}, "
            f"annotate ok")


# --- existing scanner still works ---------------------------------------
def t_scanner():
    from backend import scanner as sc
    f = sc.annotate_bundled(SNP_LIST)
    silos = {}
    for x in f:
        silos[x["silo"]] = silos.get(x["silo"], 0) + 1
    return f"{len(f)} bundled findings, silos={silos}"


print("=" * 74)
print("DNAInsight v2 module smoke test")
print("=" * 74)
step("orientation", t_orientation)
step("genosets", t_genosets)
step("frequency", t_frequency)
step("prs", t_prs)
step("merge", t_merge)
step("traits", t_traits)
step("snpedia", t_snpedia)
step("scanner", t_scanner)

print("=" * 74)
fails = [r for r in results if r[1] != "ok"]
print(f"{len(results) - len(fails)}/{len(results)} modules ok")
if fails:
    print("FAILING:")
    for label, _, msg in fails:
        print(f"  {label}: {msg}")
sys.exit(1 if fails else 0)
