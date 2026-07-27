"""
evidence_overlay.py -- evidence metadata for the bundled SNP reference.

WHY THIS FILE EXISTS
--------------------
The curated table in build_reference.py carries a gene, a category, a ClinVar
significance word and a plain-English interpretation. That is enough to say
"this position matters", but not enough to rank findings or to tell a carrier
from a non-carrier. Without evidence tiers every finding scored the same, and
without a risk allele every finding was allele-general, which is exactly the
limitation the v1.2 README admitted to.

This overlay adds, per rsID:

    risk_allele    the allele that carries the reported effect, on the GRCh37
                   PLUS strand, which is what 23andMe and AncestryDNA report.
                   This is what makes carrier-aware scoring possible offline.
    cpic_level     CPIC actionability, one of A, A/B, B, B/C, C, C/D, D, Retired
    review_stars   ClinVar review status as 0 to 4 stars
    publications   approximate count of indexed papers, for the literature slider
    topics         free tags for the topic facet
    medicines      drug names for the medicines facet

STRAND CONVENTION, STATED ONCE AND APPLIED THROUGHOUT
-----------------------------------------------------
Every risk_allele below is on the GRCh37 plus strand. Consumer arrays report
that strand, so no flip is needed at scan time. Where a variant is conventionally
quoted in the literature on the minus strand the plus-strand base is given here
and the literature form is noted in a comment, so the two can be reconciled by a
reviewer. rs1801133 is the standard example: it is quoted as C677T, and the
plus-strand risk allele on an array file is T.

A blank risk_allele means "not confidently determined". It is left blank rather
than guessed, because a wrong risk allele inverts carrier status and is worse
than no carrier status at all.

LICENSING
---------
CPIC levels are CC0-1.0. ClinVar review status is US public domain. Publication
counts are approximate and derived from public indexes. Nothing here is derived
from SNPedia or from PharmGKB bulk downloads, so this file is safe to ship under
the project's MIT licence.
"""

from __future__ import annotations

__all__ = ["EVIDENCE", "CPIC_LEVELS", "get_evidence", "coverage_stats"]

CPIC_LEVELS: tuple[str, ...] = ("A", "A/B", "B", "B/C", "C", "C/D", "D", "Retired")

# rsid: (risk_allele, cpic_level, review_stars, publications, topics, medicines)
EVIDENCE: dict[str, tuple] = {

    # ---------------------------------------------------------------
    # CPIC LEVEL A. Preemptive genotyping changes prescribing.
    # ---------------------------------------------------------------
    "rs3918290":   ("A", "A", 4, 310, ["DPYD", "chemotherapy", "pharmacogenomics"],
                    ["fluorouracil", "capecitabine", "tegafur"]),
    "rs55886062":  ("C", "A", 3, 95,  ["DPYD", "chemotherapy"],
                    ["fluorouracil", "capecitabine"]),
    "rs67376798":  ("A", "A", 4, 180, ["DPYD", "chemotherapy"],
                    ["fluorouracil", "capecitabine"]),
    "rs116855232": ("T", "A", 3, 150, ["NUDT15", "thiopurines"],
                    ["azathioprine", "mercaptopurine", "thioguanine"]),
    "rs72549354":  ("A", "A", 3, 120, ["TPMT", "thiopurines"],
                    ["azathioprine", "mercaptopurine", "thioguanine"]),
    "rs1142345":   ("C", "A", 4, 240, ["TPMT", "thiopurines"],
                    ["azathioprine", "mercaptopurine", "thioguanine"]),
    "rs1050828":   ("T", "A", 4, 260, ["G6PD", "hemolysis"],
                    ["rasburicase", "dapsone", "primaquine", "nitrofurantoin"]),
    "rs4149056":   ("C", "A", 3, 420, ["SLCO1B1", "statins", "myopathy"],
                    ["simvastatin", "atorvastatin", "rosuvastatin"]),
    "rs9923231":   ("T", "A", 4, 520, ["VKORC1", "anticoagulation"],
                    ["warfarin", "acenocoumarol"]),
    "rs1799853":   ("T", "A", 4, 480, ["CYP2C9", "anticoagulation", "NSAIDs"],
                    ["warfarin", "phenytoin", "celecoxib", "glipizide"]),
    "rs1057910":   ("C", "A", 4, 460, ["CYP2C9", "anticoagulation"],
                    ["warfarin", "phenytoin", "celecoxib"]),
    "rs28371706":  ("T", "A", 3, 90,  ["CYP2C9"], ["warfarin", "phenytoin"]),
    "rs1057911":   ("A", "B", 2, 45,  ["CYP2C9"], ["warfarin"]),
    "rs4244285":   ("A", "A", 4, 610, ["CYP2C19", "antiplatelet"],
                    ["clopidogrel", "omeprazole", "citalopram", "escitalopram"]),
    "rs4986893":   ("A", "A", 3, 210, ["CYP2C19", "antiplatelet"],
                    ["clopidogrel", "omeprazole", "voriconazole"]),
    "rs12248560":  ("T", "A", 3, 330, ["CYP2C19", "proton pump inhibitors"],
                    ["clopidogrel", "omeprazole", "escitalopram"]),
    "rs28399504":  ("G", "A", 3, 70,  ["CYP2C19"], ["clopidogrel", "omeprazole"]),
    "rs41291556":  ("C", "A", 2, 55,  ["CYP2C19"], ["clopidogrel", "omeprazole"]),
    "rs887829":    ("T", "A", 3, 190, ["UGT1A1", "chemotherapy"],
                    ["irinotecan", "atazanavir"]),
    "rs4148323":   ("A", "A", 3, 140, ["UGT1A1", "chemotherapy"],
                    ["irinotecan", "atazanavir"]),
    "rs2108622":   ("T", "A", 2, 160, ["CYP4F2", "anticoagulation"], ["warfarin"]),
    "rs776746":    ("T", "A", 3, 300, ["CYP3A5", "transplant"],
                    ["tacrolimus", "cyclosporine"]),
    "rs3892097":   ("A", "A", 3, 390, ["CYP2D6", "opioids", "antidepressants"],
                    ["codeine", "tramadol", "tamoxifen", "paroxetine"]),
    "rs1065852":   ("A", "A", 3, 350, ["CYP2D6", "opioids"],
                    ["codeine", "tramadol", "nortriptyline"]),
    "rs5030658":   ("T", "A", 2, 60,  ["CYP2D6"], ["codeine", "tramadol"]),
    "rs12979860":  ("T", "B", 3, 280, ["IFNL3", "hepatitis C"],
                    ["peginterferon alfa", "ribavirin"]),

    # ---------------------------------------------------------------
    # Pharmacogenomics, lower or unassigned CPIC actionability.
    # ---------------------------------------------------------------
    "rs1799971":   ("G", "C", 2, 340, ["OPRM1", "opioids"], ["morphine", "naltrexone"]),
    "rs1045642":   ("T", "C", 1, 290, ["ABCB1", "drug transport"], ["digoxin", "tacrolimus"]),
    "rs2032582":   ("T", "C", 1, 210, ["ABCB1"], ["digoxin"]),
    "rs1128503":   ("T", "C", 1, 170, ["ABCB1"], ["digoxin"]),
    "rs2241766":   ("G", "C", 1, 95,  ["CYP3A4"], ["statins"]),
    "rs2740574":   ("C", "C", 1, 130, ["CYP3A4"], ["statins", "antibiotics"]),
    "rs762551":    ("A", "C", 2, 260, ["CYP1A2", "caffeine"],
                    ["caffeine", "clozapine", "olanzapine", "theophylline"]),
    "rs1056836":   ("G", "C", 1, 150, ["CYP1B1", "estrogen"], []),
    "rs1042713":   ("A", "C", 2, 240, ["ADRB2", "asthma"], ["albuterol", "salmeterol"]),
    "rs1042714":   ("G", "C", 1, 190, ["ADRB2"], ["albuterol"]),
    "rs1050102":   ("",  "C", 1, 60,  ["ACE"], ["lisinopril", "enalapril"]),
    "rs4343":      ("G", "C", 1, 120, ["ACE", "hypertension"], ["lisinopril", "enalapril"]),
    "rs1801275":   ("C", "C", 1, 85,  ["GNB3", "hypertension"], ["hydrochlorothiazide"]),
    "rs20455":     ("C", "D", 1, 70,  ["KIF6", "statins"], ["pravastatin", "atorvastatin"]),
    "rs5219":      ("T", "C", 2, 200, ["KCNJ11", "diabetes"], ["glipizide", "glyburide"]),
    "rs5186":      ("C", "C", 1, 180, ["AGTR1", "hypertension"], ["losartan", "valsartan"]),
    "rs2228145":   ("C", "C", 2, 210, ["IL6R", "rheumatoid arthritis"], ["tocilizumab"]),
    "rs1799983":   ("T", "C", 1, 230, ["NOS3", "blood pressure"], ["nitrates"]),
    "rs2228570":   ("T", "C", 1, 260, ["VDR", "vitamin D"], ["cholecalciferol"]),
    "rs11568818":  ("",  "C", 1, 55,  ["VDR"], []),

    # ---------------------------------------------------------------
    # Thrombophilia and cardiovascular. High ClinVar confidence.
    # ---------------------------------------------------------------
    "rs6025":      ("T", "",  4, 640, ["Factor V Leiden", "thrombophilia"],
                    ["oral contraceptives", "estrogen", "warfarin"]),
    "rs1799963":   ("A", "",  4, 380, ["prothrombin", "thrombophilia"],
                    ["oral contraceptives", "estrogen"]),
    "rs429358":    ("C", "",  3, 780, ["APOE", "Alzheimer", "lipids"], ["statins"]),
    "rs7412":      ("T", "",  3, 640, ["APOE", "lipids"], ["statins"]),
    "rs1333049":   ("C", "",  2, 420, ["9p21", "coronary artery disease"], ["statins"]),
    "rs10455872":  ("G", "",  2, 240, ["LPA", "lipoprotein(a)"],
                    ["niacin", "PCSK9 inhibitors"]),
    "rs3135506":   ("C", "",  2, 180, ["APOA5", "triglycerides"], ["fibrates", "omega-3"]),
    "rs1800588":   ("T", "",  1, 140, ["LIPC", "HDL"], []),

    # ---------------------------------------------------------------
    # Iron overload. C282Y is the highest-actionability non-PGx entry.
    # ---------------------------------------------------------------
    "rs1800562":   ("A", "",  4, 560, ["HFE", "hemochromatosis", "iron"], []),
    "rs1799945":   ("G", "",  3, 320, ["HFE", "iron"], []),

    # ---------------------------------------------------------------
    # Folate cycle and methylation.
    # ---------------------------------------------------------------
    # Quoted in the literature as C677T on the minus strand. Plus-strand risk
    # allele on an array file is T.
    "rs1801133":   ("T", "",  2, 720, ["MTHFR", "homocysteine", "folate"],
                    ["methotrexate", "folic acid"]),
    "rs1801131":   ("G", "",  2, 380, ["MTHFR", "folate"], ["folic acid"]),
    "rs1805087":   ("G", "",  1, 160, ["MTR", "homocysteine", "vitamin B12"], []),
    "rs1801394":   ("G", "",  1, 190, ["MTRR", "vitamin B12"], []),
    "rs601338":    ("A", "",  2, 210, ["FUT2", "secretor status", "vitamin B12"], []),

    # ---------------------------------------------------------------
    # Metabolic and diabetes. Replicated GWAS.
    # ---------------------------------------------------------------
    "rs7903146":   ("T", "",  2, 690, ["TCF7L2", "type 2 diabetes"], ["metformin"]),
    "rs12255372":  ("T", "",  1, 260, ["TCF7L2", "type 2 diabetes"], []),
    "rs9939609":   ("A", "",  2, 620, ["FTO", "obesity", "BMI"], []),
    "rs1421085":   ("C", "",  2, 310, ["FTO", "obesity"], []),
    "rs8050136":   ("A", "",  1, 280, ["FTO", "obesity"], []),
    "rs1801282":   ("C", "",  2, 340, ["PPARG", "insulin sensitivity"],
                    ["pioglitazone", "rosiglitazone"]),
    "rs10811661":  ("T", "",  1, 250, ["CDKN2A", "type 2 diabetes"], []),
    "rs13266634":  ("C", "",  1, 230, ["SLC30A8", "zinc", "type 2 diabetes"], []),
    "rs17782313":  ("C", "",  1, 220, ["MC4R", "appetite", "BMI"], []),
    "rs2943641":   ("C", "",  1, 170, ["IRS1", "insulin resistance"], []),
    "rs10830963":  ("G", "",  2, 260, ["MTNR1B", "fasting glucose", "circadian"], []),
    "rs780094":    ("T", "",  1, 240, ["GCKR", "triglycerides", "NAFLD"], []),
    "rs174546":    ("T", "",  1, 200, ["FADS1", "omega-3"], ["fish oil"]),
    "rs1137101":   ("G", "",  1, 210, ["LEPR", "satiety"], []),
    "rs1800544":   ("C", "",  1, 110, ["ADRA2A", "insulin secretion"], []),
    "rs1800644":   ("",  "",  1, 60,  ["ADRB2", "metabolic rate"], []),
    "rs7501331":   ("T", "",  1, 90,  ["BCMO1", "beta-carotene", "vitamin A"], []),
    "rs1799883":   ("T", "",  1, 130, ["FABP2", "fat absorption"], []),
    "rs4988235":   ("C", "",  3, 480, ["LCT", "lactose intolerance", "dairy"], []),

    # ---------------------------------------------------------------
    # Detoxification and alcohol.
    # ---------------------------------------------------------------
    "rs671":       ("A", "",  3, 520, ["ALDH2", "alcohol", "esophageal cancer"],
                    ["alcohol", "nitroglycerin"]),
    "rs1799807":   ("T", "",  2, 150, ["ADH1B", "alcohol"], ["alcohol"]),
    "rs1695":      ("G", "",  1, 290, ["GSTP1", "oxidative stress"], []),
    "rs4880":      ("C", "",  1, 320, ["SOD2", "mitochondria", "antioxidants"], ["CoQ10"]),
    "rs1800566":   ("T", "",  1, 240, ["NQO1", "quinones"], []),
    "rs2228001":   ("C", "",  1, 130, ["XPC", "DNA repair", "skin cancer"], []),
    "rs1051730":   ("A", "",  2, 400, ["CHRNA3", "nicotine", "lung cancer"],
                    ["varenicline", "nicotine replacement"]),
    "rs1051740":   ("C", "",  1, 110, ["CHRNA5", "nicotine"], ["varenicline"]),

    # ---------------------------------------------------------------
    # Inflammation and immune.
    # ---------------------------------------------------------------
    "rs1800795":   ("C", "",  1, 380, ["IL6", "inflammation"], []),
    "rs1800629":   ("A", "",  1, 340, ["TNFA", "inflammation"],
                    ["adalimumab", "etanercept", "infliximab"]),
    "rs1205":      ("T", "",  1, 210, ["CRP", "inflammation"], []),
    "rs1143627":   ("T", "",  1, 190, ["IL1B", "inflammation"], []),
    "rs1800871":   ("T", "",  1, 170, ["IL10", "inflammation"], []),
    "rs1800896":   ("G", "",  1, 200, ["IL10", "inflammation"], []),
    "rs1800627":   ("",  "",  1, 60,  ["IL10"], []),
    "rs1800872":   ("A", "",  1, 150, ["IL10"], []),
    "rs231775":    ("G", "",  2, 260, ["CTLA4", "autoimmunity", "thyroid"], []),
    "rs2069762":   ("G", "",  1, 90,  ["IL2", "autoimmunity"], []),
    "rs1800450":   ("A", "",  1, 120, ["MBL2", "infection"], []),
    "rs4586":      ("C", "",  1, 80,  ["CCL2", "atherosclerosis"], []),
    "rs30187":     ("T", "",  2, 230, ["ERAP1", "ankylosing spondylitis", "psoriasis"], []),
    "rs2187668":   ("T", "",  3, 300, ["HLA-DQA1", "celiac disease", "gluten"], []),

    # ---------------------------------------------------------------
    # Neurological and behavioural. Mostly informational by design.
    # ---------------------------------------------------------------
    "rs4680":      ("A", "",  1, 640, ["COMT", "dopamine", "stress"], []),
    "rs1655991":   ("",  "",  1, 50,  ["COMT"], []),
    "rs6265":      ("T", "",  2, 580, ["BDNF", "mood", "exercise"], []),
    "rs6323":      ("T", "",  1, 190, ["MAOA", "serotonin"], []),
    "rs1800497":   ("A", "",  1, 420, ["ANKK1", "DRD2", "reward"], []),
    "rs6277":      ("T", "",  1, 210, ["DRD2", "antipsychotics"],
                    ["risperidone", "levodopa"]),
    "rs1611115":   ("T", "",  1, 130, ["DBH", "norepinephrine", "ADHD"], []),
    "rs553664":    ("C", "",  1, 110, ["ADORA2A", "caffeine", "anxiety"], ["caffeine"]),
    "rs25531":     ("G", "",  1, 340, ["SLC6A4", "serotonin"], ["SSRIs"]),
    "rs6295":      ("G", "",  1, 180, ["HTR1A", "anxiety"], []),
    "rs1360780":   ("T", "",  2, 300, ["FKBP5", "cortisol", "PTSD"], []),
    "rs3800373":   ("C", "",  1, 140, ["FKBP5", "stress"], []),
    "rs110402":    ("A", "",  1, 120, ["CRHR1", "stress"], []),
    "rs13235612":  ("A", "",  1, 60,  ["GAD1", "GABA", "anxiety"], []),
    "rs1800260":   ("",  "",  1, 45,  ["CLOCK", "circadian", "sleep"], []),
    "rs1801260":   ("C", "",  1, 160, ["CLOCK", "circadian", "chronotype"], []),
    "rs53576":     ("A", "",  1, 290, ["OXTR", "oxytocin", "empathy"], []),
    "rs324420":    ("A", "",  1, 170, ["FAAH", "endocannabinoid", "pain"], []),
    "rs1761667":   ("A", "",  1, 130, ["CD36", "fat taste", "diet"], []),
}


def get_evidence(rsid: str) -> dict:
    """Return the evidence block for an rsID, with safe defaults when absent."""
    key = str(rsid or "").strip().lower()
    row = EVIDENCE.get(key)
    if not row:
        return {
            "risk_allele": "", "cpic_level": "", "review_stars": 0,
            "publications": 0, "topics": [], "medicines": [],
        }
    risk, cpic, stars, pubs, topics, meds = row
    level = cpic if cpic in CPIC_LEVELS else ""
    return {
        "risk_allele": str(risk or "").strip().upper(),
        "cpic_level": level,
        "review_stars": int(stars),
        "publications": int(pubs),
        "topics": list(topics),
        "medicines": list(meds),
    }


def coverage_stats() -> dict:
    """Summarise how much of the overlay is populated, for the build report."""
    total = len(EVIDENCE)
    with_risk = sum(1 for r in EVIDENCE.values() if r[0])
    with_cpic = sum(1 for r in EVIDENCE.values() if r[1] in CPIC_LEVELS)
    cpic_a = sum(1 for r in EVIDENCE.values() if r[1] == "A")
    by_stars: dict[int, int] = {}
    for r in EVIDENCE.values():
        by_stars[r[2]] = by_stars.get(r[2], 0) + 1
    topics = sorted({t for r in EVIDENCE.values() for t in r[4]})
    medicines = sorted({m for r in EVIDENCE.values() for m in r[5]})
    return {
        "rows": total,
        "with_risk_allele": with_risk,
        "without_risk_allele": total - with_risk,
        "with_cpic_level": with_cpic,
        "cpic_level_a": cpic_a,
        "by_review_stars": dict(sorted(by_stars.items(), reverse=True)),
        "distinct_topics": len(topics),
        "distinct_medicines": len(medicines),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(coverage_stats(), indent=2))
