"""
build_prs.py -- builds data/prs_models.json for the DNAInsight PRS layer.

Two modes.

Default mode writes the bundled models: seven original additive scores whose
weights are natural-log odds ratios taken from published per-allele effect
sizes, restricted to rsIDs that data/build_reference.py already bundles so
every model is scorable from a consumer array file. No PGS Catalog weight
file is copied in this mode. The trailing comment on each weight row gives
the odds ratio the beta was taken from, so a reviewer can check ln(OR)
against the stored weight by hand. The builder checks it too and refuses to
write a model whose beta and odds ratio disagree.

``--from-pgs PGS000123`` imports one score from the PGS Catalog REST API. The
license string is printed verbatim and then gated: anything non-commercial,
no-derivatives, academic-only or research-only is refused with a non-zero
exit, because bundling it would silently relicense this project.

``--validate`` reloads the built file and recomputes every reference mean and
standard deviation, failing when a stored value has drifted.

Reference distribution
----------------------
The reference mean and sd are computed analytically from the effect allele
frequencies under Hardy-Weinberg rather than simulated, so the build is
deterministic and reviewable. For an additive dosage D in {0, 1, 2} with
effect allele frequency f, Hardy-Weinberg gives E[D] = 2f and
Var[D] = 2f(1 - f). Summing independent variants:

    mean = sum over variants of 2 * f * w
    var  = sum over variants of 2 * f * (1 - f) * w * w
    sd   = sqrt(var)

Independence is the assumption doing the work. Where two variants of a model
sit in one LD block (rs7903146 with rs12255372, and the three FTO markers)
the real variance is larger than that formula gives, so the sd is a lower
bound and the percentiles it feeds are over-dispersed at the tails. The
affected models say so in their own description field rather than leaving it
to this comment.

Effect allele frequencies
-------------------------
Frequencies come from data/frequencies.json, averaged over the European
panels CEU, GBR, IBS, TSI and FIN. That file is in dbSNP orientation while
the effect alleles below are quoted on the array plus strand, so every
lookup goes through backend.frequency.resolve_strand and
backend.frequency.allele_frequency rather than indexing the JSON directly.
A rsID absent from the bundled data falls back to a literature estimate.
Every variant records which of the two it used in ``af_source``.

Usage:
    python data/build_prs.py
    python data/build_prs.py --validate
    python data/build_prs.py --from-pgs PGS000123
"""

import argparse
import gzip
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_BASE = Path(__file__).resolve().parent.parent
if str(_BASE) not in sys.path:
    sys.path.insert(0, str(_BASE))

from backend import frequency  # noqa: E402

try:  # importable both as data.build_prs and as a script inside data/
    from data.build_reference import REFERENCE  # noqa: E402
except ImportError:  # pragma: no cover
    from build_reference import REFERENCE  # noqa: E402


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
OUT_FILE = _BASE / "data" / "prs_models.json"

VERSION = "2.0.0"
GENOME_BUILD = "GRCh37"
LICENSE = "MIT (DNAInsight-authored weights from published effect sizes)"
SOURCE = "DNAInsight additive model from published per-allele effect sizes"
REFERENCE_POPULATION = "EUR"

# The EUR-ish panels averaged to get an effect allele frequency. Deliberately
# five panels rather than one, so a single small panel cannot swing a weight.
EUR_PANELS: tuple[str, ...] = ("CEU", "GBR", "IBS", "TSI", "FIN")

# ln(OR) is stored to six decimals, so the beta / odds ratio cross-check
# cannot be tighter than half a unit in the last place.
WEIGHT_TOLERANCE = 1e-5

# --validate fails above this absolute drift in a stored mean or sd.
DRIFT_TOLERANCE = 1e-6

# A palindromic site (A/T or C/G heterozygote) cannot have its strand checked
# against the table, so the resolved frequency is compared with the literature
# estimate instead and a warning is printed when they disagree this much.
PALINDROME_WARN = 0.15

META_NOTE = (
    "Additive scores only. Weights are natural-log odds ratios from published "
    "per-allele effect sizes, restricted to rsIDs bundled in "
    "data/build_reference.py. Reference mean and sd are analytic under "
    "Hardy-Weinberg from European effect allele frequencies, so they assume "
    "the variants of a model are independent. Percentiles are meaningful only "
    "for European-ancestry users."
)


# ---------------------------------------------------------------------------
# Bundled models
#
# Row format: (rsid, effect allele, other allele, weight, odds ratio)
#
# The effect allele is on the GRCh37 PLUS strand, which is what consumer
# arrays report and what backend.prs.compute_model expects to compare
# against. It is NOT necessarily the orientation data/frequencies.json
# stores; MTHFR rs1801133 is the classic case, array C/T against dbSNP G/A.
# ---------------------------------------------------------------------------
MODEL_SPECS: tuple[dict[str, Any], ...] = (
    {
        "id": "t2d",
        "trait": "Type 2 diabetes",
        "efo": "EFO_0001360",
        "description": (
            "Ten common variants with replicated per-allele associations with "
            "type 2 diabetes, covering insulin secretion (TCF7L2, KCNJ11, "
            "SLC30A8, MTNR1B, CDKN2A/B), insulin sensitivity (PPARG, IRS1), "
            "appetite (MC4R) and hepatic glucose handling (GCKR). rs7903146 "
            "and rs12255372 are both TCF7L2 markers in strong linkage "
            "disequilibrium, so their contributions overlap and the analytic "
            "sd understates the true spread."
        ),
        "citation": "Zeggini E, Nature Genetics, 2008",
        "pmid": "18372903",
        "rows": (
            ("rs7903146",  "T", "C", 0.314811, 1.37),  # OR 1.37  TCF7L2 intron 3
            ("rs1801282",  "C", "G", 0.131028, 1.14),  # OR 1.14  PPARG Pro12 versus Ala12
            ("rs10811661", "T", "C", 0.182322, 1.20),  # OR 1.20  CDKN2A/B beta cell reserve
            ("rs13266634", "C", "T", 0.113329, 1.12),  # OR 1.12  SLC30A8 Arg325
            ("rs5219",     "T", "C", 0.131028, 1.14),  # OR 1.14  KCNJ11 Lys23
            ("rs780094",   "C", "T", 0.058269, 1.06),  # OR 1.06  GCKR, T lowers fasting glucose
            ("rs10830963", "G", "C", 0.139762, 1.15),  # OR 1.15  MTNR1B fasting glucose
            ("rs2943641",  "C", "T", 0.173953, 1.19),  # OR 1.19  IRS1 muscle insulin resistance
            ("rs17782313", "C", "T", 0.076961, 1.08),  # OR 1.08  MC4R
            ("rs12255372", "T", "G", 0.270027, 1.31),  # OR 1.31  TCF7L2, LD with rs7903146
        ),
    },
    {
        "id": "cad",
        "trait": "Coronary artery disease",
        "efo": "EFO_0001645",
        "description": (
            "Seven variants spanning the dominant 9p21 susceptibility locus, "
            "lipoprotein(a) concentration, endothelial nitric oxide "
            "synthesis, the renin-angiotensin axis and two thrombophilias "
            "with replicated but modest coronary effects. rs4149056 is "
            "carried at a deliberately small weight: it is a statin "
            "transporter variant that acts on tolerability and adherence "
            "rather than a primary susceptibility locus."
        ),
        "citation": "Schunkert H, Nature Genetics, 2011",
        "pmid": "21378990",
        "rows": (
            ("rs1333049",  "C", "G", 0.254642, 1.29),  # OR 1.29  9p21 CDKN2B-AS1
            ("rs10455872", "G", "A", 0.530628, 1.70),  # OR 1.70  LPA, Clarke NEJM 2009
            ("rs1799983",  "T", "G", 0.104360, 1.11),  # OR 1.11  NOS3 Asp298
            ("rs6025",     "T", "C", 0.157004, 1.17),  # OR 1.17  factor V Leiden, Ye Lancet 2006
            ("rs1799963",  "A", "G", 0.270027, 1.31),  # OR 1.31  prothrombin G20210A
            ("rs4149056",  "C", "T", 0.048790, 1.05),  # OR 1.05  SLCO1B1 *5 statin transport
            ("rs5186",     "C", "A", 0.086178, 1.09),  # OR 1.09  AGTR1 A1166C
        ),
    },
    {
        "id": "bmi",
        "trait": "Body mass index",
        "efo": "EFO_0004340",
        "description": (
            "Seven common variants associated with body mass index and "
            "obesity risk. Weights are obesity odds ratios rather than "
            "kg/m2 betas, which keeps every model in this file on one scale. "
            "rs9939609, rs1421085 and rs8050136 are three markers of the same "
            "FTO intron 1 haplotype and are not independent, so the analytic "
            "sd is a lower bound. rs1801282 enters with the Ala12 allele as "
            "the effect allele, the opposite direction to the type 2 diabetes "
            "model, because Ala12 improves insulin sensitivity while still "
            "associating with higher measured BMI."
        ),
        "citation": "Speliotes EK, Nature Genetics, 2010",
        "pmid": "20935630",
        "rows": (
            ("rs9939609",  "A", "T", 0.277632, 1.32),  # OR 1.32  FTO, Frayling Science 2007
            ("rs1421085",  "C", "T", 0.262364, 1.30),  # OR 1.30  FTO, LD with rs9939609
            ("rs8050136",  "A", "C", 0.239017, 1.27),  # OR 1.27  FTO, LD with rs9939609
            ("rs17782313", "C", "T", 0.113329, 1.12),  # OR 1.12  MC4R, Loos 2008
            ("rs1137101",  "G", "A", 0.076961, 1.08),  # OR 1.08  LEPR Arg223
            ("rs1800544",  "G", "C", 0.058269, 1.06),  # OR 1.06  ADRA2A -1291
            ("rs1801282",  "G", "C", 0.058269, 1.06),  # OR 1.06  PPARG Ala12
        ),
    },
    {
        "id": "vte",
        "trait": "Venous thromboembolism",
        "efo": "",
        "description": (
            "The two established European thrombophilias plus one endothelial "
            "modifier. Weights are heterozygote odds ratios applied "
            "per-allele, which models a heterozygote correctly and "
            "understates the homozygote: factor V Leiden homozygotes carry "
            "far more than the squared heterozygote risk this additive form "
            "implies. A high score here reflects one or two large-effect "
            "alleles rather than a polygenic burden, so read the underlying "
            "genotypes, not the percentile alone."
        ),
        "citation": "Emmerich J, Thrombosis and Haemostasis, 2001",
        "pmid": "",
        "rows": (
            ("rs6025",     "T", "C", 1.589235, 4.90),  # OR 4.90  factor V Leiden heterozygote
            ("rs1799963",  "A", "G", 1.029619, 2.80),  # OR 2.80  prothrombin G20210A heterozygote
            ("rs1799983",  "T", "G", 0.198851, 1.22),  # OR 1.22  NOS3 Asp298
        ),
    },
    {
        "id": "ldl",
        "trait": "Elevated LDL cholesterol",
        "efo": "EFO_0004611",
        "description": (
            "Five variants acting on circulating LDL and remnant "
            "lipoproteins. The two APOE sites pull in opposite directions "
            "and rs7412 therefore carries a negative weight: the e2 allele "
            "lowers LDL. A carrier of both raising and lowering alleles can "
            "land mid-distribution while still having an informative "
            "genotype, which is the case for reading the variant list rather "
            "than the band alone."
        ),
        "citation": "Willer CJ, Nature Genetics, 2013",
        "pmid": "24097068",
        "rows": (
            ("rs429358",  "C", "T",  0.262364, 1.30),  # OR 1.30  APOE e4 defining site
            ("rs7412",    "T", "C", -0.328504, 0.72),  # OR 0.72  APOE e2 defining site, protective
            ("rs1800588", "T", "C",  0.067659, 1.07),  # OR 1.07  LIPC -514 hepatic lipase
            ("rs3135506", "C", "G",  0.350657, 1.42),  # OR 1.42  APOA5 Trp19
            ("rs1799983", "T", "G",  0.039221, 1.04),  # OR 1.04  NOS3, small endothelial modifier
        ),
    },
    {
        "id": "homocysteine",
        "trait": "Elevated homocysteine",
        "efo": "",
        "description": (
            "Five variants in the folate and cobalamin remethylation "
            "pathway. Four of the five are stored in data/frequencies.json "
            "on the opposite strand to the array orientation used here, so "
            "this model is the sharpest test of the strand layer in the "
            "build: MTHFR rs1801133 is quoted C/T on an array and A/G by "
            "dbSNP. Homocysteine is a modifiable biochemical intermediate, "
            "so a high score is a prompt to measure the analyte rather than "
            "a finding in itself."
        ),
        "citation": "Klerk M, JAMA, 2002",
        "pmid": "12387654",
        "rows": (
            ("rs1801133", "T", "C", 0.336472, 1.40),  # OR 1.40  MTHFR 677T, dbSNP stores A/G
            ("rs1801131", "C", "A", 0.139762, 1.15),  # OR 1.15  MTHFR 1298C, dbSNP stores G/T
            ("rs1805087", "G", "A", 0.113329, 1.12),  # OR 1.12  MTR 2756G
            ("rs1801394", "G", "A", 0.086178, 1.09),  # OR 1.09  MTRR 66G
            ("rs601338",  "A", "G", 0.095310, 1.10),  # OR 1.10  FUT2 non-secretor
        ),
    },
    {
        "id": "inflammation",
        "trait": "Systemic inflammation (hsCRP)",
        "efo": "EFO_0004458",
        "description": (
            "Five regulatory variants influencing the baseline inflammatory "
            "set point read out by high-sensitivity CRP. Two sit in or near "
            "CRP itself and three in upstream cytokine genes. Weights are "
            "odds ratios for the upper tertile of hsCRP rather than betas on "
            "mg/L, so the score orders people and does not predict a "
            "concentration. hsCRP moves with acute infection, adiposity and "
            "smoking, all of which swamp this genetic component."
        ),
        "citation": "Dehghan A, Circulation, 2011",
        "pmid": "21300955",
        "rows": (
            ("rs1205",    "C", "T", 0.190620, 1.21),  # OR 1.21  CRP 3 prime UTR
            ("rs30187",   "T", "C", 0.165514, 1.18),  # OR 1.18  ERAP1 Lys528
            ("rs1800795", "G", "C", 0.139762, 1.15),  # OR 1.15  IL6 -174
            ("rs1800629", "A", "G", 0.182322, 1.20),  # OR 1.20  TNF -308
            ("rs1143627", "C", "T", 0.113329, 1.12),  # OR 1.12  IL1B -31, dbSNP stores A/G
        ),
    },
)


# ---------------------------------------------------------------------------
# Literature fallback frequencies
#
# Used only when data/frequencies.json has nothing for an rsID, and as the
# cross-check for palindromic sites whose strand cannot be verified against
# the table. Keys are "rsid:effect_allele" on the array plus strand. Values
# are European effect allele frequencies as published.
# ---------------------------------------------------------------------------
LITERATURE_AF: dict[str, float] = {
    "rs7903146:T":  0.30,
    "rs1801282:C":  0.88,
    "rs1801282:G":  0.12,
    "rs10811661:T": 0.83,
    "rs13266634:C": 0.72,
    "rs5219:T":     0.36,
    "rs780094:C":   0.60,
    "rs10830963:G": 0.29,
    "rs2943641:C":  0.63,
    "rs17782313:C": 0.24,
    "rs12255372:T": 0.29,
    "rs1333049:C":  0.47,
    "rs10455872:G": 0.07,
    "rs1799983:T":  0.34,
    "rs6025:T":     0.02,
    "rs1799963:A":  0.02,
    "rs4149056:C":  0.16,
    "rs5186:C":     0.28,
    "rs9939609:A":  0.41,
    "rs1421085:C":  0.42,
    "rs8050136:A":  0.41,
    "rs1137101:G":  0.47,
    "rs1800544:G":  0.26,
    "rs429358:C":   0.15,
    "rs7412:T":     0.07,
    "rs1800588:T":  0.21,
    "rs3135506:C":  0.06,
    "rs1801133:T":  0.34,
    "rs1801131:C":  0.31,
    "rs1805087:G":  0.18,
    "rs1801394:G":  0.53,
    "rs601338:A":   0.45,
    "rs1205:C":     0.68,
    "rs30187:T":    0.35,
    "rs1800795:G":  0.58,
    "rs1800629:A":  0.13,
    "rs1143627:C":  0.35,
}


# ---------------------------------------------------------------------------
# Bundled rsID gate
# ---------------------------------------------------------------------------

def bundled_rsids() -> set[str]:
    """Return every rsID present in data/build_reference.py.

    This is the authoritative list of positions the app can annotate, so a
    model variant outside it would be unusable and the build refuses it.
    """
    return {str(row[0]).strip().lower() for row in REFERENCE}


# ---------------------------------------------------------------------------
# Effect allele frequency resolution
# ---------------------------------------------------------------------------

def literature_frequency(rsid: str, effect_allele: str) -> float | None:
    """Return the published European effect allele frequency, or None."""
    return LITERATURE_AF.get(f"{rsid}:{effect_allele.upper()}")


def panel_frequency(rsid: str, table_effect: str,
                    table_other: str, population: str) -> float | None:
    """Return one panel's frequency for the effect allele, or None.

    Queries ``backend.frequency.allele_frequency``, which is strand tolerant.
    When the effect allele is genuinely absent from a panel's table but the
    other allele is present, the biallelic complement ``1 - f(other)`` is
    used: for a two-allele site an absent key means the allele was not seen,
    which is a frequency of zero rather than unknown. Skipping those panels
    instead would bias a rare effect allele upward, because only the panels
    that happened to observe it would count.
    """
    direct = frequency.allele_frequency(rsid, table_effect, population)
    if direct is not None:
        return direct
    other = frequency.allele_frequency(rsid, table_other, population)
    if other is None:
        return None
    return 1.0 - other


def effect_allele_frequency(rsid: str, effect_allele: str,
                            other_allele: str) -> tuple[float, str, dict]:
    """Resolve one variant's European effect allele frequency.

    Returns ``(frequency, af_source, strand)``. The effect allele given is on
    the array plus strand; ``backend.frequency.resolve_strand`` maps it into
    the orientation data/frequencies.json actually stores before any lookup,
    so a minus-strand variant such as MTHFR rs1801133 resolves instead of
    silently reporting no data. The value is the mean over the EUR_PANELS
    that carry the site. A site with no bundled data falls back to
    :func:`literature_frequency`, and the returned ``af_source`` says which
    path was taken and flags an unverifiable palindromic strand.

    Raises SystemExit when neither the bundled data nor the literature table
    can supply a frequency, because a variant with no frequency cannot be
    mean-imputed by backend.prs and would silently reduce coverage.
    """
    strand = frequency.resolve_strand(rsid, effect_allele, other_allele)
    table_effect = strand["allele1"] or effect_allele.upper()
    table_other = strand["allele2"] or other_allele.upper()

    values = [
        f for f in (
            panel_frequency(rsid, table_effect, table_other, pop)
            for pop in EUR_PANELS
        ) if f is not None
    ]

    literature = literature_frequency(rsid, effect_allele)

    if not values:
        if literature is None:
            raise SystemExit(
                f"FATAL: no effect allele frequency for {rsid} {effect_allele}. "
                "Add a literature value to LITERATURE_AF before building."
            )
        return literature, "literature EUR estimate (rsID absent from bundled data)", strand

    value = round(sum(values) / len(values), 6)

    if strand["ambiguous"]:
        source = ("frequencies.json EUR panel average, palindromic site so the "
                  "strand could not be verified against the table")
        if literature is not None and abs(value - literature) > PALINDROME_WARN:
            print(f"  WARNING: {rsid} {effect_allele} palindromic, bundled "
                  f"{value:.4f} versus literature {literature:.4f}. "
                  "Check the orientation by hand.")
    elif strand["flipped"]:
        source = ("frequencies.json EUR panel average, complemented from dbSNP "
                  "orientation to the array plus strand")
    else:
        source = "frequencies.json EUR panel average"
    return value, source, strand


# ---------------------------------------------------------------------------
# Analytic reference distribution
# ---------------------------------------------------------------------------

def reference_moments(variants: list[dict]) -> tuple[float, float]:
    """Return the analytic ``(mean, sd)`` of a model under Hardy-Weinberg.

        mean = sum over variants of 2 * f * w
        var  = sum over variants of 2 * f * (1 - f) * w * w
        sd   = sqrt(var)

    with f the effect allele frequency and w the log odds ratio. Assumes the
    variants are independent, which overstates precision for a model carrying
    two markers of one LD block.
    """
    mean = 0.0
    variance = 0.0
    for variant in variants:
        f = float(variant["effect_allele_frequency"])
        w = float(variant["weight"])
        mean += 2.0 * f * w
        variance += 2.0 * f * (1.0 - f) * w * w
    return mean, math.sqrt(variance)


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def build_variants(spec: dict) -> list[dict]:
    """Turn one spec's weight rows into scorable variant dicts.

    Verifies that each stored beta really is ln of the odds ratio in its
    trailing comment, and that every rsID is one data/build_reference.py
    bundles, before resolving the effect allele frequency.
    """
    bundled = bundled_rsids()
    variants: list[dict] = []
    for rsid, effect, other, weight, odds in spec["rows"]:
        if rsid.lower() not in bundled:
            raise SystemExit(
                f"FATAL: {spec['id']} uses {rsid}, which data/build_reference.py "
                "does not bundle. Drop the variant or add it to REFERENCE."
            )
        expected = math.log(float(odds))
        if abs(float(weight) - expected) > WEIGHT_TOLERANCE:
            raise SystemExit(
                f"FATAL: {spec['id']} {rsid} weight {weight} is not ln({odds}) "
                f"= {expected:.6f}. Fix the row before building."
            )
        eaf, af_source, strand = effect_allele_frequency(rsid, effect, other)
        variants.append({
            "rsid":                    rsid,
            "effect_allele":           effect.upper(),
            "other_allele":            other.upper(),
            "weight":                  float(weight),
            "effect_allele_frequency": eaf,
            "odds_ratio":              float(odds),
            "af_source":               af_source,
            "strand_flipped":          bool(strand["flipped"]),
            "strand_ambiguous":        bool(strand["ambiguous"]),
        })
    return variants


def build_models() -> dict[str, dict]:
    """Build every bundled model, keyed by model id."""
    models: dict[str, dict] = {}
    for spec in MODEL_SPECS:
        variants = build_variants(spec)
        mean, sd = reference_moments(variants)
        models[spec["id"]] = {
            "id":            spec["id"],
            "trait":         spec["trait"],
            "efo":           spec["efo"],
            "description":   spec["description"],
            "build":         GENOME_BUILD,
            "source":        SOURCE,
            "citation":      spec["citation"],
            "pmid":          spec["pmid"],
            "license":       LICENSE,
            "variant_count": len(variants),
            "variants":      variants,
            "reference": {
                "population": REFERENCE_POPULATION,
                "mean":       round(mean, 9),
                "sd":         round(sd, 9),
            },
        }
        flipped = sum(1 for v in variants if v["strand_flipped"])
        ambiguous = sum(1 for v in variants if v["strand_ambiguous"])
        print(f"  {spec['id']:<14} {len(variants):>2} variants  "
              f"mean {mean:+.4f}  sd {sd:.4f}  "
              f"strand flipped {flipped}, palindromic {ambiguous}")
    return models


def write_models(models: dict[str, dict], out_path: Path) -> None:
    """Write the ``{_meta, models}`` document to disk as UTF-8 with LF endings."""
    document = {
        "_meta": {
            "version":     VERSION,
            "built_at":    datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "model_count": len(models),
            "license":     LICENSE,
            "note":        META_NOTE,
        },
        "models": models,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(document, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def run_build(out_path: Path) -> int:
    """Build the bundled models and write them. Returns a process exit code."""
    print(f"Building {len(MODEL_SPECS)} PRS models")
    models = build_models()
    write_models(models, out_path)
    total = sum(m["variant_count"] for m in models.values())
    print(f"WROTE {out_path}  ({len(models)} models, {total} variant rows)")
    return 0


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

def run_validate(out_path: Path) -> int:
    """Recompute every stored mean and sd and report any drift.

    Returns 0 when every model is within DRIFT_TOLERANCE of its recomputed
    moments, 1 otherwise. Also re-checks variant_count and the license gate,
    since both are cheap and both are things a hand edit could break.
    """
    if not out_path.exists():
        print(f"FAIL: {out_path} does not exist. Run the builder first.")
        return 1
    with open(out_path, "r", encoding="utf-8") as handle:
        document = json.load(handle)

    models = document.get("models") or {}
    if not models:
        print(f"FAIL: {out_path} carries no models.")
        return 1

    failures = 0
    for model_id, model in sorted(models.items()):
        variants = model.get("variants") or []
        mean, sd = reference_moments(variants)
        stored = model.get("reference") or {}
        d_mean = abs(float(stored.get("mean", 0.0)) - mean)
        d_sd = abs(float(stored.get("sd", 0.0)) - sd)
        ok = d_mean <= DRIFT_TOLERANCE and d_sd <= DRIFT_TOLERANCE

        if int(model.get("variant_count") or -1) != len(variants):
            print(f"FAIL {model_id}: variant_count {model.get('variant_count')} "
                  f"does not match {len(variants)} stored variants")
            ok = False
        allowed, term = license_is_allowed(str(model.get("license") or ""))
        if not allowed:
            print(f"FAIL {model_id}: license contains the forbidden term {term!r}")
            ok = False

        verdict = "OK  " if ok else "FAIL"
        print(f"{verdict} {model_id:<14} mean drift {d_mean:.3e}  "
              f"sd drift {d_sd:.3e}  ({len(variants)} variants)")
        failures += 0 if ok else 1

    if failures:
        print(f"VALIDATE FAILED: {failures} of {len(models)} models drifted")
        return 1
    print(f"VALIDATE CLEAN: {len(models)} models, all moments within "
          f"{DRIFT_TOLERANCE:g}")
    return 0


# ---------------------------------------------------------------------------
# PGS Catalog import
#
# The license gate is the point of this code path. A PGS Catalog score may
# carry author terms that forbid commercial or derivative use, and bundling
# one of those into an MIT-licensed application would relicense the
# application by accident. The gate is a refusal by default: only CC0, CC BY
# 4.0 and the default EMBL-EBI terms pass, and the raw license string is
# printed verbatim first so a human can overrule the machine.
# ---------------------------------------------------------------------------
PGS_API = "https://www.pgscatalog.org/rest/score/{pgs_id}"
HTTP_TIMEOUT = 30
USER_AGENT = (
    "DNAInsight/2.0 (local polygenic score builder; contact: "
    "https://github.com/dnainsight/dnainsight)"
)

FORBIDDEN_LICENSE_TERMS: tuple[str, ...] = (
    "NonCommercial",
    "NoDerivatives",
    "ND 4.0",
    "academic",
    "not-for-profit",
    "research purposes only",
    "commercial purposes should contact",
)

RSID_COLUMNS: tuple[str, ...] = ("rsid", "rsids", "snp", "snpid", "variant_id")
EFFECT_COLUMNS: tuple[str, ...] = ("effect_allele", "a1")
OTHER_COLUMNS: tuple[str, ...] = (
    "other_allele", "reference_allele", "noneffect_allele", "a2",
    "hm_inferotherallele",
)
WEIGHT_COLUMNS: tuple[str, ...] = (
    "effect_weight", "beta", "weight", "or", "odds_ratio", "hazard_ratio",
)
NON_ADDITIVE_COLUMNS: tuple[str, ...] = (
    "is_dominant", "is_recessive", "is_interaction",
)
TRUE_TOKENS = {"true", "1", "yes", "t"}


def license_is_allowed(license_text: str) -> tuple[bool, str]:
    """Return ``(allowed, matched_term)`` for a license string.

    Matching is case-insensitive substring matching against
    FORBIDDEN_LICENSE_TERMS. An empty license is treated as allowed, because
    the PGS Catalog leaves the field blank when the default EMBL-EBI terms
    apply, but the caller still prints it for a human to read.
    """
    lowered = str(license_text or "").lower()
    for term in FORBIDDEN_LICENSE_TERMS:
        if term.lower() in lowered:
            return False, term
    return True, ""


def fetch_score_metadata(pgs_id: str) -> dict:
    """Fetch one score's metadata from the PGS Catalog REST API."""
    import requests

    url = PGS_API.format(pgs_id=pgs_id)
    response = requests.get(
        url, timeout=HTTP_TIMEOUT,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict) or not payload.get("id"):
        raise SystemExit(f"FATAL: {url} returned no score record for {pgs_id}")
    return payload


def pick_scoring_url(meta: dict) -> str:
    """Return the best scoring-file URL from a score record.

    Prefers the GRCh37 harmonized file, because the bundled reference and
    every consumer array file in this project are GRCh37. Falls back to
    GRCh38 and then to the author-submitted file.
    """
    harmonized = meta.get("ftp_harmonized_scoring_files") or {}
    for build in ("GRCh37", "GRCh38"):
        block = harmonized.get(build) or {}
        for key in ("positions", "additional"):
            url = block.get(key)
            if url:
                return str(url)
    url = meta.get("ftp_scoring_file")
    if not url:
        raise SystemExit("FATAL: score record carries no scoring file URL")
    return str(url)


def download_scoring_text(url: str) -> str:
    """Download a scoring file and return its text, transparently gunzipped."""
    import requests

    response = requests.get(
        url, timeout=HTTP_TIMEOUT, headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    body = response.content
    if url.endswith(".gz") or body[:2] == b"\x1f\x8b":
        body = gzip.decompress(body)
    return body.decode("utf-8", errors="replace")


def parse_header_block(text: str) -> dict[str, str]:
    """Parse a scoring file's ``#key=value`` header block, keyed by name.

    Header lines are read by KEY rather than by position, because the number
    of header lines and their order both vary between Catalog releases.
    ``##`` banner lines carry no ``=`` and are skipped.
    """
    header: dict[str, str] = {}
    for line in text.splitlines():
        if not line.startswith("#"):
            break
        stripped = line.lstrip("#").strip()
        if "=" not in stripped:
            continue
        key, _, value = stripped.partition("=")
        header[key.strip().lower()] = value.strip()
    return header


def sniff_columns(header_line: str) -> dict[str, int]:
    """Map logical column names onto indexes in a scoring file's header row.

    Returns a dict that always carries ``rsid``, ``effect_allele``, ``weight``
    and optionally ``other_allele`` plus any non-additive flag columns.
    Necessary because the column set genuinely varies: some files ship
    ``reference_allele``, others ``other_allele``, others neither, and
    harmonized files add ``hm_`` prefixed duplicates.
    """
    fields = [f.strip().lower() for f in header_line.split("\t")]
    index = {name: position for position, name in enumerate(fields)}
    out: dict[str, int] = {}

    def first_of(candidates: tuple[str, ...]) -> int | None:
        for candidate in candidates:
            if candidate in index:
                return index[candidate]
        return None

    for logical, candidates in (
        ("rsid", RSID_COLUMNS),
        ("effect_allele", EFFECT_COLUMNS),
        ("other_allele", OTHER_COLUMNS),
        ("weight", WEIGHT_COLUMNS),
    ):
        position = first_of(candidates)
        if position is not None:
            out[logical] = position

    for flag in NON_ADDITIVE_COLUMNS:
        if flag in index:
            out[flag] = index[flag]

    missing = [name for name in ("rsid", "effect_allele", "weight")
               if name not in out]
    if missing:
        raise SystemExit(
            f"FATAL: scoring file has no usable column for {', '.join(missing)}. "
            f"Columns seen: {', '.join(fields)}"
        )
    return out


def to_natural_log_beta(raw_weight: float, weight_type: str) -> float:
    """Convert a scoring file weight into a natural-log beta.

    backend.prs sums dosage times a natural-log odds ratio, so a file quoting
    plain odds ratios, log base 2 or log base 10 has to be rescaled first.
    Real Catalog files use all of these: PGS000123 ships ``Log2(OR)``, and
    treating that as a beta would understate every weight by a factor of
    ln 2. An unrecognised or absent type is taken to be a beta already, which
    is what ``beta``, ``log(OR)`` and the common ``NR`` placeholder mean.
    """
    kind = str(weight_type or "").strip().lower().replace(" ", "")
    if kind in ("or", "odds_ratio", "oddsratio", "hr", "hazard_ratio",
                "hazardratio", "rr", "risk_ratio"):
        return math.log(raw_weight)
    if kind in ("log2(or)", "log2or", "log2(hr)"):
        return raw_weight * math.log(2.0)
    if kind in ("log10(or)", "log10or", "log10(hr)"):
        return raw_weight * math.log(10.0)
    return raw_weight


def parse_scoring_rows(text: str, columns: dict[str, int],
                       weight_type: str) -> tuple[list[dict], int, int]:
    """Parse the data rows of a scoring file into additive variant dicts.

    Returns ``(variants, non_additive_refused, unusable_skipped)``. Rows
    flagged ``is_dominant``, ``is_recessive`` or ``is_interaction`` are
    refused rather than approximated, because backend.prs implements a simple
    additive model and scoring a dominant term additively would be wrong
    without saying so. Rows with no rsID or no parsable weight are skipped.
    """
    variants: list[dict] = []
    refused = 0
    skipped = 0
    in_body = False

    for line in text.splitlines():
        if line.startswith("#"):
            continue
        if not in_body:
            in_body = True  # this is the column header row
            continue
        if not line.strip():
            continue
        fields = line.split("\t")

        if any(
            len(fields) > columns[flag]
            and fields[columns[flag]].strip().lower() in TRUE_TOKENS
            for flag in NON_ADDITIVE_COLUMNS if flag in columns
        ):
            refused += 1
            continue

        try:
            rsid = fields[columns["rsid"]].strip().lower()
            effect = fields[columns["effect_allele"]].strip().upper()
            raw_weight = float(fields[columns["weight"]].strip())
        except (IndexError, ValueError):
            skipped += 1
            continue
        if not rsid.startswith("rs") or not effect:
            skipped += 1
            continue

        other_index = columns.get("other_allele")
        other = ""
        if other_index is not None and len(fields) > other_index:
            other = fields[other_index].strip().upper()

        weight = to_natural_log_beta(raw_weight, weight_type)
        variants.append({
            "rsid":                    rsid,
            "effect_allele":           effect,
            "other_allele":            other,
            "weight":                  weight,
            "effect_allele_frequency": None,
            "af_source":               "not supplied by the scoring file",
        })
    return variants, refused, skipped


def run_from_pgs(pgs_id: str) -> int:
    """Inspect one PGS Catalog score and report whether it is usable here.

    Prints the license verbatim, applies the license gate, parses the
    harmonized scoring file and reports how many of its variants intersect
    the bundled rsID list. Nothing is written to data/prs_models.json: a
    network fetch must never overwrite a reviewed, hand-checked bundle, so
    importing is an inspection step and adding a model stays a deliberate
    edit to MODEL_SPECS.
    """
    print(f"Fetching metadata for {pgs_id}")
    meta = fetch_score_metadata(pgs_id)

    license_text = str(meta.get("license") or "")
    print("---- license string, verbatim ----")
    print(license_text if license_text else "(empty: default EMBL-EBI terms)")
    print("---- end license string ----")

    allowed, term = license_is_allowed(license_text)
    if not allowed:
        print(f"REFUSED: license contains {term!r}. Only CC0, CC BY 4.0 and "
              "the default EMBL-EBI terms may be bundled with this project.")
        return 2
    print("License gate: accepted")

    print(f"  trait:    {meta.get('trait_reported', '')}")
    print(f"  variants: {meta.get('variants_number', 'unknown')}")
    print(f"  build:    {meta.get('variants_genomebuild', 'unknown')}")

    url = pick_scoring_url(meta)
    print(f"Downloading {url}")
    text = download_scoring_text(url)

    header = parse_header_block(text)
    weight_type = header.get("weight_type", "")
    print(f"  header keys:  {len(header)} parsed by key")
    print(f"  pgs_id:       {header.get('pgs_id', '(absent)')}")
    print(f"  genome_build: {header.get('genome_build', '(absent)')}")
    print(f"  weight_type:  {weight_type or '(absent, assuming beta)'}")

    header_line = ""
    for line in text.splitlines():
        if not line.startswith("#"):
            header_line = line
            break
    columns = sniff_columns(header_line)
    print(f"  columns:      {columns}")

    variants, refused, skipped = parse_scoring_rows(text, columns, weight_type)

    bundled = bundled_rsids()
    usable = [v for v in variants if v["rsid"] in bundled]
    print(f"  additive rows:        {len(variants)}")
    print(f"  non-additive refused: {refused}")
    print(f"  unparsable skipped:   {skipped}")
    print(f"  INTERSECTION with the bundled rsID list: {len(usable)} of "
          f"{len(variants)}")
    if not usable:
        print("  This score shares no variants with the bundled reference, so "
              "it is not scorable here.")
    else:
        share = 100.0 * len(usable) / max(len(variants), 1)
        print(f"  {share:.2f} percent of the score is covered. Below roughly "
              "90 percent, backend.prs will flag results reliable=False.")
        print(f"  rsIDs: {', '.join(v['rsid'] for v in usable[:20])}")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch to build, validate or PGS import."""
    parser = argparse.ArgumentParser(
        prog="build_prs.py",
        description="Build data/prs_models.json, or inspect a PGS Catalog score.",
    )
    parser.add_argument("--from-pgs", dest="from_pgs", metavar="PGS_ID",
                        help="inspect one PGS Catalog score, for example PGS000123")
    parser.add_argument("--validate", action="store_true",
                        help="recompute every stored reference mean and sd")
    parser.add_argument("--out", default=str(OUT_FILE), metavar="PATH",
                        help="output path (default data/prs_models.json)")
    args = parser.parse_args(argv)

    out_path = Path(args.out)
    if args.from_pgs:
        return run_from_pgs(args.from_pgs)
    if args.validate:
        return run_validate(out_path)
    return run_build(out_path)


if __name__ == "__main__":
    sys.exit(main())
