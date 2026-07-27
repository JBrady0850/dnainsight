"""
build_genosets.py -- Generates the bundled genoset corpus (data/genosets.json).

A genoset is a named boolean rule over several SNP genotypes: it is either
present or absent, never scored and never quantitative. See
backend/genosets.py for the criteria language and the evaluator.

Usage:
    python data/build_genosets.py

This ships as a pre-built file (genosets.json), so end users do NOT need to
run this script. Maintainers run it to update the bundled corpus.

LICENSING
---------
Every summary and interpretation below is written from scratch for DNAInsight
out of well established pharmacogenomic and clinical-genetics knowledge. No
text is copied or paraphrased from SNPedia (CC-BY-NC-SA), so this repository
stays MIT-clean. Genosets use a DNAInsight-native "dgsNNN" namespace rather
than SNPedia "gsNNN" identifiers. The optional "aka" field carries the
standard clinical name of a combination where one exists.

ALLELE ORIENTATION
------------------
All genotypes are written on the dbSNP plus (forward) strand, which is the
orientation used by 23andMe and AncestryDNA raw data files.

  ORIENTATION NOTES -- variants whose plus-strand variant allele is worth
  re-verifying against dbSNP when the reference table is next refreshed, and
  which are therefore only ever used as one branch of an or()/atleast() rule
  so a mis-set allele cannot silently create a false positive on its own:
    rs72549354 (TPMT reduced function, paired with rs1142345)
    rs1799807  (paired with rs671)
    rs1051740  (paired with rs1051730)
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.genosets import CriteriaError, parse_criteria, required_rsids, referenced_genosets, topological_order
from data.build_reference import REFERENCE

GENOSET_VERSION = "2.0.0"

VALID_CATEGORIES = {"PHARM", "METAB", "INFLAM", "NEURO", "DETOX", "CARDIO"}
VALID_SILOS = {"pre_prescription", "actionable", "informational"}
VALID_REPUTES = {"Good", "Bad", ""}


def _gs(name, criteria, magnitude, repute, summary, interpretation,
        category, silo, topics, medicines, evidence, aka=""):
    """Build one corpus entry. Keeps the table below readable."""
    return {
        "name":           name,
        "aka":            aka,
        "criteria":       criteria,
        "magnitude":      float(magnitude),
        "repute":         repute,
        "summary":        summary,
        "interpretation": interpretation,
        "category":       category,
        "silo":           silo,
        "topics":         list(topics),
        "medicines":      list(medicines),
        "evidence":       evidence,
    }


# ---------------------------------------------------------------------------
# The corpus
#
# Magnitude scale: 0 common and uninteresting, 1 unset or semi-plausible,
# 2 worth reading, 3 probably worth your time, 4+ definitely worth attention,
# 10 really significant.
# ---------------------------------------------------------------------------

GENOSETS = [

    # -- APOE diplotypes -------------------------------------------------
    # rs429358: T carries e2 or e3, C carries e4.
    # rs7412:   T carries e2, C carries e3 or e4.
    _gs("dgs001", "and(rs429358(T;T), rs7412(T;T))", 4.0, "",
        "You carry two copies of APOE e2, the rarest of the three common APOE forms.",
        "e2/e2 usually means low LDL cholesterol and the lowest genetic Alzheimer risk of any APOE pair. "
        "The trade-off is a small chance of type III hyperlipoproteinaemia, which shows up as high "
        "triglycerides and remnant particles rather than high LDL. Ask for a full lipid panel that "
        "includes triglycerides rather than LDL alone.",
        "CARDIO", "actionable", ["APOE", "cholesterol", "Alzheimer"], [], "OMIM / replicated GWAS",
        aka="APOE e2/e2"),

    _gs("dgs002", "and(rs429358(T;T), rs7412(C;T))", 2.0, "Good",
        "You carry one APOE e2 and one APOE e3, a combination associated with favourable cholesterol.",
        "e2/e3 is the most protective common APOE pairing: lower average LDL and slightly lower "
        "Alzheimer risk than the e3/e3 majority. No specific action is needed beyond ordinary "
        "cardiovascular care, though triglycerides are still worth watching because e2 raises them.",
        "CARDIO", "informational", ["APOE", "cholesterol"], [], "replicated GWAS",
        aka="APOE e2/e3"),

    _gs("dgs003", "and(rs429358(C;T), rs7412(C;T))", 3.0, "",
        "You carry one APOE e2 and one APOE e4, so the two alleles pull your risk in opposite directions.",
        "e2/e4 combines the triglyceride-raising e2 allele with the Alzheimer-risk and LDL-raising e4 "
        "allele, and the net effect is roughly average overall risk with an unusual lipid pattern. "
        "Because e4 is present, treat lipid and blood-pressure control as worthwhile, and interpret any "
        "single-number cardiovascular risk score with caution.",
        "CARDIO", "actionable", ["APOE", "cholesterol", "Alzheimer"], [], "replicated GWAS",
        aka="APOE e2/e4"),

    _gs("dgs004", "and(rs429358(T;T), rs7412(C;C))", 1.0, "",
        "You carry the most common APOE pairing, e3/e3, which is the population reference.",
        "e3/e3 is the baseline against which other APOE results are compared, so this result neither "
        "raises nor lowers your Alzheimer or lipid risk. Ordinary preventive care applies.",
        "CARDIO", "informational", ["APOE"], [], "replicated GWAS",
        aka="APOE e3/e3"),

    _gs("dgs005", "and(rs429358(C;T), rs7412(C;C))", 4.5, "Bad",
        "You carry one copy of APOE e4, which raises Alzheimer disease and cardiovascular risk.",
        "A single e4 allele roughly two- to three-fold raises lifetime Alzheimer risk relative to e3/e3 "
        "and modestly raises LDL cholesterol. This is a risk factor, not a diagnosis, and most single-e4 "
        "carriers never develop dementia. The evidence-supported levers are blood-pressure control, "
        "aerobic exercise, sleep quality, hearing correction, and not smoking.",
        "NEURO", "actionable", ["APOE", "Alzheimer", "cholesterol"], [], "OMIM / replicated GWAS",
        aka="APOE e3/e4"),

    _gs("dgs006", "and(rs429358(C;C), rs7412(C;C))", 8.0, "Bad",
        "You carry two copies of APOE e4, the strongest common genetic risk factor for Alzheimer disease.",
        "e4/e4 carries a substantially higher lifetime Alzheimer risk and an earlier average age of onset "
        "than e3/e3, along with higher LDL cholesterol. It is still a probability rather than a certainty. "
        "This result is worth discussing with a clinician: it affects anti-amyloid therapy eligibility and "
        "ARIA monitoring, and it makes aggressive vascular risk-factor control especially worthwhile. "
        "Genetic counselling is reasonable before sharing the result with family.",
        "NEURO", "actionable", ["APOE", "Alzheimer", "cholesterol"], [], "OMIM / replicated GWAS",
        aka="APOE e4/e4"),
]

# -- Folate / methylation cycle ------------------------------------------
# rs1801133 C677T: plus-strand variant allele is A. rs1801131 A1298C: G.
GENOSETS += [

    _gs("dgs007", "rs1801133(A;A)", 3.0, "Bad",
        "You carry two copies of the MTHFR C677T variant, which slows folate activation.",
        "The enzyme runs at roughly a third of normal speed, so homocysteine can drift upward when folate "
        "or B12 intake is low. Ask for a plasma homocysteine measurement; if it is raised, supplement with "
        "methylfolate (5-MTHF) rather than plain folic acid, alongside B12 and B6. This variant is common "
        "and is not on its own a reason for alarm.",
        "NEURO", "actionable", ["folate", "homocysteine", "methylation"], [], "OMIM / replicated GWAS",
        aka="MTHFR 677TT"),

    _gs("dgs008", "rs1801131(G;G)", 2.5, "Bad",
        "You carry two copies of the MTHFR A1298C variant, a milder brake on folate metabolism.",
        "A1298C reduces enzyme activity less than C677T does and rarely raises homocysteine by itself. "
        "It matters most when combined with C677T. A single homocysteine measurement is enough to decide "
        "whether methylated B-vitamin support is worth taking.",
        "NEURO", "informational", ["folate", "methylation"], [], "replicated GWAS",
        aka="MTHFR 1298CC"),

    _gs("dgs009", "and(rs1801133(A;G), rs1801131(G;T))", 3.5, "Bad",
        "You carry one copy each of MTHFR C677T and A1298C, the classic compound heterozygote pattern.",
        "Carrying one of each variant reduces enzyme activity to roughly the same degree as two copies of "
        "C677T, so the two results should be read together rather than separately. Measure homocysteine, "
        "and if it is elevated use methylfolate with B12 rather than folic acid. Relevant in pregnancy "
        "planning and before any high-dose folate regimen.",
        "NEURO", "actionable", ["folate", "homocysteine", "methylation", "pregnancy"], [],
        "replicated GWAS", aka="MTHFR compound heterozygote"),
]

# -- Warfarin dosing -----------------------------------------------------
# rs9923231 T = VKORC1 low-expression (sensitive). rs1799853 T = CYP2C9*2,
# rs1057910 C = *3, rs28371706 A = *8. rs2108622 T = CYP4F2*3.
GENOSETS += [

    _gs("dgs010", "or(and(rs9923231(T;T), or(rs1799853(T), rs1057910(C), rs28371706(A))), rs1057910(C;C), rs1799853(T;T))",
        7.0, "Bad",
        "Your VKORC1 and CYP2C9 genotypes together predict high warfarin sensitivity.",
        "You have both a low-expression VKORC1 promoter and at least one reduced-function CYP2C9 allele, "
        "or two reduced-function CYP2C9 alleles. That combination means a therapeutic INR is reached at a "
        "much lower daily dose than average, and a standard 5 mg start carries real bleeding risk. Give "
        "this result to the prescriber before the first dose so a genotype-guided algorithm can be used, "
        "and expect more frequent INR checks during induction. Direct oral anticoagulants are unaffected.",
        "PHARM", "pre_prescription", ["warfarin", "anticoagulation", "bleeding"],
        ["warfarin", "acenocoumarol", "phenprocoumon"], "CPIC Level A"),

    _gs("dgs011", "and(not(dgs010), or(rs9923231(T), rs1799853(T), rs1057910(C), rs28371706(A)))",
        4.5, "Bad",
        "You carry one warfarin-sensitivity factor, predicting an intermediate dose requirement.",
        "A single low-expression VKORC1 allele or a single reduced-function CYP2C9 allele typically shifts "
        "the maintenance dose modestly downward rather than dramatically. Genotype-guided dosing still "
        "shortens the time spent outside the therapeutic range, so pass the result to the prescriber. "
        "Routine INR monitoring during induction remains the safety net.",
        "PHARM", "pre_prescription", ["warfarin", "anticoagulation"],
        ["warfarin", "acenocoumarol"], "CPIC Level A"),

    _gs("dgs012", "and(rs9923231(C;C), rs1799853(C;C), rs1057910(A;A), rs28371706(G;G), not(rs2108622(T;T)))",
        1.5, "Good",
        "None of the tested VKORC1, CYP2C9 or CYP4F2 warfarin-dose variants are present.",
        "You carry the reference genotype at every warfarin variant on this panel, so a standard "
        "clinical starting dose is appropriate and no genotype-driven reduction is indicated. This does "
        "not remove the need for INR monitoring, and rarer CYP2C9 alleles are not covered by this array.",
        "PHARM", "pre_prescription", ["warfarin", "anticoagulation"], ["warfarin"], "CPIC Level A"),

    _gs("dgs013", "rs2108622(T;T)", 2.5, "",
        "You carry two copies of CYP4F2 *3, which shifts warfarin dose requirement slightly upward.",
        "CYP4F2 clears vitamin K1, and the *3 allele slows that clearance, leaving more vitamin K "
        "available. Homozygotes therefore need a modestly higher warfarin dose, on the order of an extra "
        "milligram a day, to reach the same INR. This works in the opposite direction to VKORC1 and "
        "CYP2C9 sensitivity, which is why dosing algorithms combine all three rather than reading any one "
        "in isolation.",
        "PHARM", "pre_prescription", ["warfarin", "vitamin K"], ["warfarin"], "CPIC Level A",
        aka="CYP4F2 *3/*3"),
]

# -- CYP2C19 metaboliser phenotypes --------------------------------------
# Loss of function: rs4244285 A (*2), rs4986893 A (*3), rs28399504 G (*4),
# rs41291556 C (*8). Increased function: rs12248560 T (*17).
GENOSETS += [

    _gs("dgs014",
        "or(rs4244285(A;A), rs4986893(A;A), rs28399504(G;G), rs41291556(C;C),"
        " and(rs4244285(A;G), rs4986893(A;G)), and(rs4244285(A;G), rs28399504(A;G)),"
        " and(rs4244285(A;G), rs41291556(C;T)), and(rs4986893(A;G), rs28399504(A;G)),"
        " and(rs4986893(A;G), rs41291556(C;T)), and(rs28399504(A;G), rs41291556(C;T)))",
        6.5, "Bad",
        "You carry two non-functional CYP2C19 alleles, making you a predicted poor metaboliser.",
        "With no working copies of the enzyme, drugs that CYP2C19 activates fail and drugs it clears "
        "accumulate. Clopidogrel is the headline problem because it is a prodrug that will not be "
        "converted, so prasugrel or ticagrelor is preferred where the indication allows. Proton pump "
        "inhibitors, escitalopram, sertraline, voriconazole and some tricyclics need lower doses. Tell "
        "every prescriber; this is one of the most consequential pharmacogenomic results on the panel.",
        "PHARM", "pre_prescription", ["CYP2C19", "antiplatelet", "PPI", "antidepressant"],
        ["clopidogrel", "omeprazole", "escitalopram", "voriconazole", "amitriptyline"], "CPIC Level A",
        aka="CYP2C19 poor metaboliser"),

    _gs("dgs015",
        "and(not(dgs014), or(rs4244285(A), rs4986893(A), rs28399504(G), rs41291556(C)))",
        4.5, "Bad",
        "You carry one non-functional CYP2C19 allele, making you a predicted intermediate metaboliser.",
        "One working copy gives roughly half-normal enzyme activity. Clopidogrel activation is measurably "
        "reduced and guidelines favour an alternative antiplatelet after coronary stenting or stroke. "
        "Proton pump inhibitor and SSRI exposure runs higher than average, which is often clinically "
        "helpful rather than harmful. Worth flagging before antiplatelet or antifungal therapy.",
        "PHARM", "pre_prescription", ["CYP2C19", "antiplatelet", "PPI"],
        ["clopidogrel", "omeprazole", "escitalopram", "voriconazole"], "CPIC Level A",
        aka="CYP2C19 intermediate metaboliser"),

    _gs("dgs016",
        "and(rs4244285(G;G), rs4986893(G;G), rs28399504(A;A), rs41291556(T;T), rs12248560(C;C))",
        1.0, "Good",
        "You carry two fully functional CYP2C19 alleles across every variant tested.",
        "This is the normal metaboliser phenotype, so standard doses of clopidogrel, proton pump "
        "inhibitors, escitalopram and voriconazole are expected to behave as intended. Rarer CYP2C19 "
        "alleles are not covered by consumer arrays, so a clinically important discrepancy should still "
        "prompt confirmatory testing.",
        "PHARM", "pre_prescription", ["CYP2C19"], ["clopidogrel", "omeprazole"], "CPIC Level A",
        aka="CYP2C19 normal metaboliser"),

    _gs("dgs017",
        "and(rs12248560(C;T), not(rs4244285(A), rs4986893(A), rs28399504(G), rs41291556(C)))",
        2.5, "",
        "You carry one CYP2C19 *17 increased-function allele and no loss-of-function alleles.",
        "The *17 promoter variant raises transcription, so this predicts a rapid metaboliser. Proton pump "
        "inhibitors clear faster and may under-treat reflux at standard doses; escitalopram levels run "
        "lower than expected. Clopidogrel activation is slightly increased, which marginally raises "
        "bleeding risk rather than reducing efficacy.",
        "PHARM", "pre_prescription", ["CYP2C19", "PPI"],
        ["omeprazole", "escitalopram", "clopidogrel"], "CPIC Level A",
        aka="CYP2C19 rapid metaboliser"),

    _gs("dgs018",
        "and(rs12248560(T;T), not(rs4244285(A), rs4986893(A), rs28399504(G), rs41291556(C)))",
        3.5, "",
        "You carry two CYP2C19 *17 alleles and no loss-of-function alleles: an ultrarapid metaboliser.",
        "Two increased-function copies mean CYP2C19 substrates are cleared unusually fast. Proton pump "
        "inhibitor failure at standard dose is the most common practical consequence, and rabeprazole or "
        "a dose increase is the usual answer. Escitalopram and sertraline may under-perform. Mention this "
        "if a reflux or antidepressant regimen appears to be doing nothing.",
        "PHARM", "pre_prescription", ["CYP2C19", "PPI", "antidepressant"],
        ["omeprazole", "rabeprazole", "escitalopram", "sertraline"], "CPIC Level A",
        aka="CYP2C19 ultrarapid metaboliser"),

    _gs("dgs019", "or(dgs014, dgs015)", 5.0, "Bad",
        "Your CYP2C19 genotype predicts reduced clopidogrel activation.",
        "Clopidogrel is inactive until CYP2C19 converts it, so carrying any loss-of-function allele blunts "
        "platelet inhibition and raises the rate of stent thrombosis and recurrent ischaemic events. "
        "Prasugrel and ticagrelor do not depend on CYP2C19 and are the guideline-preferred substitutes "
        "after percutaneous coronary intervention. Do not stop clopidogrel on the strength of this result "
        "alone; take it to the cardiologist who prescribed it.",
        "PHARM", "pre_prescription", ["clopidogrel", "antiplatelet", "stent"],
        ["clopidogrel", "prasugrel", "ticagrelor"], "CPIC Level A"),
]

# -- Statin myopathy (SLCO1B1, with and without APOE e4) -----------------
GENOSETS += [

    _gs("dgs020", "rs4149056(C;C)", 5.0, "Bad",
        "You carry two copies of SLCO1B1 *5, which strongly raises statin muscle-toxicity risk.",
        "SLCO1B1 pumps statins out of the blood and into the liver. Two reduced-function copies leave "
        "several times more drug circulating in muscle, and the risk of myopathy on high-dose simvastatin "
        "is substantial. Simvastatin above 20 mg should be avoided; pravastatin, rosuvastatin at low dose, "
        "or fluvastatin are better tolerated. Report new muscle aching or dark urine promptly.",
        "PHARM", "pre_prescription", ["statin", "myopathy", "cholesterol"],
        ["simvastatin", "atorvastatin", "pravastatin", "rosuvastatin"], "CPIC Level A",
        aka="SLCO1B1 *5/*5"),

    _gs("dgs021", "and(rs4149056(C), or(dgs005, dgs006))", 5.5, "Bad",
        "You carry SLCO1B1 *5 together with APOE e4, combining statin intolerance risk with a strong reason to treat.",
        "This is the awkward combination: APOE e4 raises LDL cholesterol and cardiovascular risk, so lipid "
        "lowering is clearly worthwhile, while SLCO1B1 *5 makes the most commonly prescribed statin more "
        "likely to cause muscle symptoms. The practical route is a hydrophilic statin such as pravastatin "
        "or low-dose rosuvastatin, escalated slowly, with ezetimibe added rather than pushing the statin "
        "dose. Do not simply avoid treatment.",
        "PHARM", "pre_prescription", ["statin", "myopathy", "APOE", "cholesterol"],
        ["pravastatin", "rosuvastatin", "ezetimibe", "simvastatin"], "CPIC Level A"),

    _gs("dgs022", "and(rs4149056(C), not(dgs005, dgs006))", 4.0, "Bad",
        "You carry SLCO1B1 *5 without an APOE e4 allele.",
        "Statin exposure in muscle runs higher than average, so muscle aching is more likely, but you lack "
        "the extra APOE-driven lipid and vascular risk. If a statin is indicated, start low and prefer "
        "pravastatin or rosuvastatin over simvastatin. Creatine kinase measurement is worthwhile if "
        "symptoms appear rather than routinely.",
        "PHARM", "pre_prescription", ["statin", "myopathy"],
        ["simvastatin", "pravastatin", "rosuvastatin"], "CPIC Level A"),
]

# -- Thiopurines (TPMT, NUDT15) ------------------------------------------
GENOSETS += [

    _gs("dgs023", "or(rs1142345(C;C), and(rs72549354(T), rs1142345(C)), and(rs72549354(T;T), rs1142345(T)))",
        8.0, "Bad",
        "Your TPMT genotype predicts markedly reduced thiopurine methyltransferase activity.",
        "Thiopurines are broken down by TPMT, and with little or no enzyme the drug accumulates as toxic "
        "nucleotides in bone marrow. Conventional azathioprine or mercaptopurine dosing can cause "
        "life-threatening neutropenia within weeks. Guidelines call for a drastic dose reduction, often to "
        "around a tenth of standard, or an alternative immunosuppressant entirely, with weekly blood "
        "counts. This must reach the prescriber before the first tablet.",
        "PHARM", "pre_prescription", ["thiopurine", "myelosuppression", "IBD", "leukaemia"],
        ["azathioprine", "mercaptopurine", "thioguanine"], "CPIC Level A"),

    _gs("dgs024", "or(dgs023, rs116855232(T))", 7.0, "Bad",
        "You carry a TPMT or NUDT15 variant that raises thiopurine myelosuppression risk.",
        "TPMT and NUDT15 sit on independent arms of thiopurine metabolism, and a defect in either produces "
        "the same clinical picture of severe, sometimes abrupt bone-marrow suppression. NUDT15 deficiency "
        "is the dominant cause in people of East Asian, South Asian and Hispanic ancestry, where TPMT "
        "testing alone misses it. Dose reduction plus close blood-count monitoring is required; "
        "homozygotes generally need a different drug class.",
        "PHARM", "pre_prescription", ["thiopurine", "NUDT15", "TPMT", "myelosuppression"],
        ["azathioprine", "mercaptopurine", "thioguanine"], "CPIC Level A"),
]

# -- Fluoropyrimidines (DPYD) --------------------------------------------
GENOSETS += [

    _gs("dgs025", "or(rs3918290(T;T), rs55886062(C;C), rs67376798(A;A))", 9.0, "Bad",
        "You carry two defective DPYD alleles, predicting near-complete loss of fluoropyrimidine clearance.",
        "DPD is the enzyme that disposes of about 80 percent of a 5-fluorouracil or capecitabine dose. "
        "With essentially no activity, standard chemotherapy dosing is life-threatening and produces "
        "severe mucositis, neutropenia and neurotoxicity. Fluoropyrimidines are contraindicated, and an "
        "alternative regimen should be chosen. This result belongs in the oncology record before any "
        "colorectal, breast or gastric chemotherapy decision.",
        "PHARM", "pre_prescription", ["chemotherapy", "DPYD", "toxicity"],
        ["fluorouracil", "capecitabine", "tegafur"], "CPIC Level A"),

    _gs("dgs026", "and(not(dgs025), or(rs3918290(T), rs55886062(C), rs67376798(A)))", 7.0, "Bad",
        "You carry one defective DPYD allele, predicting partial loss of fluoropyrimidine clearance.",
        "Heterozygotes have roughly half-normal DPD activity and a several-fold higher rate of grade 3 or 4 "
        "toxicity on standard 5-FU or capecitabine dosing. Guidelines recommend starting at 25 to 50 "
        "percent of the usual dose with escalation guided by tolerance, which preserves efficacy while "
        "avoiding the worst reactions. Give this to the oncology team before treatment planning.",
        "PHARM", "pre_prescription", ["chemotherapy", "DPYD", "toxicity"],
        ["fluorouracil", "capecitabine"], "CPIC Level A"),
]

# -- Irinotecan (UGT1A1) -------------------------------------------------
GENOSETS += [

    _gs("dgs027", "or(rs887829(T;T), rs4148323(A;A), and(rs887829(C;T), rs4148323(A;G)))", 5.0, "Bad",
        "Your UGT1A1 genotype predicts substantially reduced glucuronidation capacity.",
        "UGT1A1 attaches sugar groups to bilirubin and to the active irinotecan metabolite so both can be "
        "excreted. Two reduced-function copies produce mild lifelong unconjugated hyperbilirubinaemia, the "
        "benign pattern usually labelled Gilbert syndrome, and they slow clearance of several drugs. The "
        "practical significance is mostly oncological, but it also explains harmless jaundice during "
        "fasting or illness.",
        "PHARM", "pre_prescription", ["UGT1A1", "bilirubin", "Gilbert"],
        ["irinotecan", "atazanavir", "nilotinib"], "CPIC Level A",
        aka="UGT1A1 poor metaboliser"),

    _gs("dgs028", "or(dgs027, rs887829(T), rs4148323(A))", 4.5, "Bad",
        "You carry a UGT1A1 reduced-function allele that raises irinotecan neutropenia risk.",
        "Irinotecan's active metabolite SN-38 is cleared by UGT1A1, so reduced-function carriers reach "
        "higher SN-38 exposure and have more grade 3 or 4 neutropenia and diarrhoea on FOLFIRI. Guidelines "
        "support a lower irinotecan starting dose in homozygotes, with escalation if tolerated. Bring this "
        "to the oncology team; it changes starting dose rather than drug choice.",
        "PHARM", "pre_prescription", ["irinotecan", "neutropenia", "chemotherapy"],
        ["irinotecan"], "CPIC Level A"),
]

# -- Hereditary haemochromatosis (HFE) -----------------------------------
GENOSETS += [

    _gs("dgs029", "rs1800562(A;A)", 8.0, "Bad",
        "You carry two copies of HFE C282Y, the classic hereditary haemochromatosis genotype.",
        "C282Y homozygosity removes the brake on intestinal iron absorption, and iron gradually loads into "
        "liver, pancreas, heart and joints. Penetrance is incomplete, and many homozygotes never develop "
        "organ damage, but the condition is easy to detect and easy to treat if found early. Ask for serum "
        "ferritin and transferrin saturation now; if they are raised, therapeutic phlebotomy is curative "
        "for the iron loading. Avoid iron and high-dose vitamin C supplements and inform first-degree "
        "relatives, who each have a one-in-four chance of the same genotype.",
        "METAB", "actionable", ["iron", "haemochromatosis", "liver", "ferritin"], [], "OMIM",
        aka="HFE C282Y/C282Y"),

    _gs("dgs030", "and(rs1800562(A;G), rs1799945(C;G))", 4.0, "Bad",
        "You are a compound heterozygote for HFE C282Y and H63D.",
        "One C282Y plus one H63D allele gives a modest increase in iron absorption. Most people with this "
        "combination keep normal iron studies for life, but a minority load iron, particularly alongside "
        "alcohol use, fatty liver disease or another liver insult. A one-off ferritin and transferrin "
        "saturation measurement settles the question, repeated every few years if either is borderline.",
        "METAB", "actionable", ["iron", "haemochromatosis", "ferritin"], [], "OMIM",
        aka="HFE C282Y/H63D compound heterozygote"),

    _gs("dgs031", "and(rs1799945(G;G), rs1800562(G;G))", 2.0, "",
        "You carry two copies of HFE H63D and no C282Y allele.",
        "H63D is a mild variant, and H63D homozygosity is common in many populations without causing "
        "clinical iron overload. On its own it is not a reason for surveillance. A single ferritin "
        "measurement is reasonable for reassurance, especially if there is unexplained fatigue, joint pain "
        "or abnormal liver enzymes.",
        "METAB", "informational", ["iron", "haemochromatosis"], [], "OMIM",
        aka="HFE H63D/H63D"),
]

# -- Thrombophilia (F5, F2) ----------------------------------------------
GENOSETS += [

    _gs("dgs032", "rs6025(C;T)", 5.0, "Bad",
        "You carry one copy of Factor V Leiden, the most common inherited clotting-risk variant.",
        "The variant makes factor V resistant to inactivation, so clotting continues longer than it should. "
        "Heterozygotes have roughly four to eight times the background risk of deep vein thrombosis, which "
        "is still a low absolute risk in ordinary life but climbs during pregnancy, after surgery, on long "
        "flights and on combined oestrogen contraception. Discuss contraception choices and surgical clot "
        "prophylaxis with a clinician; no treatment is needed while you are well.",
        "CARDIO", "actionable", ["thrombosis", "clotting", "pregnancy", "contraception"], [], "OMIM",
        aka="Factor V Leiden heterozygote"),

    _gs("dgs033", "rs6025(T;T)", 8.0, "Bad",
        "You carry two copies of Factor V Leiden.",
        "Homozygosity raises venous thrombosis risk far more than carrying one copy, into the range where "
        "unprovoked clots become genuinely likely over a lifetime. Combined oestrogen contraception is "
        "generally avoided, pregnancy warrants haematology input, and any surgery or immobility needs "
        "active clot prophylaxis. This result deserves a referral rather than self-management, and "
        "first-degree relatives should be told they may be carriers.",
        "CARDIO", "actionable", ["thrombosis", "clotting", "pregnancy", "surgery"], [], "OMIM",
        aka="Factor V Leiden homozygote"),

    _gs("dgs034", "rs1799963(A)", 4.0, "Bad",
        "You carry the prothrombin G20210A variant.",
        "This variant raises circulating prothrombin, giving carriers roughly two to three times the "
        "background venous thrombosis risk. As with Factor V Leiden the absolute everyday risk stays low, "
        "and the variant matters mainly around pregnancy, oestrogen use, surgery and prolonged immobility. "
        "Mention it before any of those, and know the warning signs of calf swelling and unexplained "
        "breathlessness.",
        "CARDIO", "actionable", ["thrombosis", "clotting", "pregnancy"], [], "OMIM",
        aka="Prothrombin G20210A"),

    _gs("dgs035", "and(rs6025(T), rs1799963(A))", 7.0, "Bad",
        "You carry both Factor V Leiden and prothrombin G20210A.",
        "The two variants act on different steps of the clotting cascade, and their risks multiply rather "
        "than simply add, putting this combination well above either variant alone. Combined oestrogen "
        "contraception is best avoided, pregnancy should be managed with haematology involvement, and "
        "thromboprophylaxis is warranted for surgery, immobility and long-haul travel. A haematology "
        "review is appropriate even without a personal clot history.",
        "CARDIO", "actionable", ["thrombosis", "clotting", "pregnancy", "contraception"], [], "OMIM",
        aka="Factor V Leiden compound"),
]

# -- Lactase persistence (LCT/MCM6) --------------------------------------
GENOSETS += [

    _gs("dgs036", "rs4988235(A)", 1.5, "Good",
        "You carry the lactase persistence allele, so you most likely digest lactose as an adult.",
        "A single copy of the persistence allele is usually enough to keep lactase switched on after "
        "weaning, which is why most people of northern European ancestry drink milk without trouble. If "
        "dairy still causes symptoms the cause is more likely a milk protein sensitivity, irritable bowel "
        "syndrome or small-intestinal bacterial overgrowth than lactose itself.",
        "METAB", "informational", ["lactose", "dairy", "digestion"], [], "OMIM / replicated GWAS",
        aka="Lactase persistence"),

    _gs("dgs037", "rs4988235(G;G)", 2.0, "",
        "You carry two ancestral LCT alleles, the genotype for adult lactase non-persistence.",
        "Lactase production falls after childhood in this genotype, which is the global majority pattern "
        "rather than a disease. Undigested lactose ferments in the colon and causes bloating, wind and "
        "loose stools in proportion to the dose. Fermented dairy such as hard cheese and yoghurt is "
        "usually well tolerated, and lactase tablets cover occasional exposure. Keep calcium and vitamin D "
        "intake up if you cut dairy substantially.",
        "METAB", "informational", ["lactose", "dairy", "digestion", "calcium"], [],
        "OMIM / replicated GWAS", aka="Lactase non-persistence"),
]

# -- Alcohol and acetaldehyde (ALDH2, ADH1B) -----------------------------
GENOSETS += [

    _gs("dgs038", "rs671(A)", 4.5, "Bad",
        "You carry the ALDH2 variant that causes the alcohol flushing reaction.",
        "ALDH2 clears acetaldehyde, the first breakdown product of alcohol. The variant enzyme is close to "
        "inactive and behaves dominantly, so even one copy causes acetaldehyde to accumulate: facial "
        "flushing, a racing heart, nausea and headache after modest drinking. Acetaldehyde is a recognised "
        "carcinogen, and carriers who drink regularly have a markedly higher risk of oesophageal and "
        "head-and-neck cancer. Reducing or avoiding alcohol is the single most effective response.",
        "DETOX", "actionable", ["alcohol", "acetaldehyde", "cancer risk"], [], "OMIM / replicated GWAS",
        aka="ALDH2 alcohol flush"),

    _gs("dgs039", "and(rs671(A), rs1799807(C))", 5.5, "Bad",
        "You carry the ALDH2 flushing variant together with the alcohol dehydrogenase variant on this panel.",
        "Acetaldehyde is produced quickly and then cleared slowly, so exposure per drink is higher than "
        "with either variant alone. Epidemiological studies of this combination show the steepest "
        "alcohol-related risk of upper gastrointestinal cancer. Practical advice is unambiguous: treat "
        "alcohol as something to avoid rather than moderate, and take persistent swallowing difficulty or "
        "hoarseness seriously.",
        "DETOX", "actionable", ["alcohol", "acetaldehyde", "oesophageal cancer"], [],
        "replicated GWAS"),
]

# -- Type 2 diabetes polygenic burden ------------------------------------
_T2D_LOCI = ("rs7903146(T), rs1801282(C;C), rs10811661(T), rs13266634(C), "
             "rs5219(T), rs780094(T), rs10830963(G)")

GENOSETS += [

    _gs("dgs040", f"atleast(3, {_T2D_LOCI})", 2.5, "Bad",
        "You carry at least three of the seven type 2 diabetes risk variants on this panel.",
        "Each of these variants nudges insulin secretion, insulin sensitivity or fasting glucose in the "
        "wrong direction, and their effects add up. Three or more puts you above average genetic "
        "susceptibility, which shifts risk by a modest amount rather than making diabetes inevitable. The "
        "response is ordinary but effective: keep an eye on HbA1c every few years, keep weight stable, and "
        "do regular resistance and aerobic exercise, which blunts most of this genetic effect.",
        "METAB", "actionable", ["type 2 diabetes", "insulin", "HbA1c", "polygenic"], [],
        "replicated GWAS"),

    _gs("dgs041", f"atleast(5, {_T2D_LOCI})", 4.0, "Bad",
        "You carry at least five of the seven type 2 diabetes risk variants on this panel.",
        "This is a high common-variant burden, concentrated on beta-cell function rather than obesity, so "
        "risk can be raised even at a normal body weight. Annual or biennial HbA1c or fasting glucose is "
        "worthwhile from mid-life, and earlier if there is a family history or gestational diabetes. "
        "Structured lifestyle intervention has been shown to work at least as well in people with high "
        "genetic risk as in those without.",
        "METAB", "actionable", ["type 2 diabetes", "insulin", "HbA1c", "polygenic"],
        ["metformin"], "replicated GWAS"),
]

# -- Obesity polygenic burden (FTO, MC4R, LEPR) --------------------------
GENOSETS += [

    _gs("dgs042",
        "atleast(3, rs9939609(A), rs1421085(C), rs8050136(A), rs17782313(C), rs1137101(G))",
        2.5, "Bad",
        "You carry at least three of the FTO, MC4R and leptin-receptor variants associated with higher body weight.",
        "These variants act mainly through appetite and satiety signalling rather than metabolic rate, so "
        "the lived experience is feeling less full after a meal rather than burning fewer calories. Average "
        "effect sizes are small, a few kilograms across a population, and they are strongly modifiable. "
        "Higher-protein meals, adequate sleep and regular activity reliably reduce the measured effect of "
        "the FTO variants in particular.",
        "METAB", "informational", ["obesity", "appetite", "FTO", "MC4R", "polygenic"], [],
        "replicated GWAS"),
]

# -- Coronary artery disease burden --------------------------------------
GENOSETS += [

    _gs("dgs043", "atleast(2, rs1333049(C), rs10455872(G), rs6025(T), rs1799983(T))", 3.5, "Bad",
        "You carry at least two of the coronary artery disease risk variants on this panel.",
        "The 9p21 locus, lipoprotein(a), factor V Leiden and the eNOS variant raise coronary risk through "
        "independent mechanisms, so carrying more than one compounds the effect. None of them is "
        "deterministic, and conventional risk factors still dominate. Worth acting on: know your blood "
        "pressure and lipid numbers, ask whether a lipoprotein(a) measurement is indicated, and treat "
        "smoking cessation as the highest-value change available.",
        "CARDIO", "actionable", ["coronary artery disease", "cholesterol", "polygenic"], [],
        "replicated GWAS"),

    _gs("dgs044", "and(rs10455872(G), dgs043)", 4.5, "Bad",
        "You carry the lipoprotein(a)-raising variant on top of a broader coronary risk burden.",
        "Lipoprotein(a) is largely genetically fixed, is not lowered by statins or by diet, and is an "
        "independent driver of atherosclerosis and aortic valve calcification. Combined with other "
        "coronary risk variants it justifies measuring Lp(a) once, since a single result guides lifetime "
        "management. If it is high, the response is tighter LDL targets rather than any Lp(a)-specific "
        "drug, though targeted agents are in late-stage trials.",
        "CARDIO", "actionable", ["lipoprotein(a)", "coronary artery disease", "aortic valve"],
        ["statins", "PCSK9 inhibitors"], "replicated GWAS"),
]

# -- Caffeine ------------------------------------------------------------
# NOTE: the ADORA2A anxiety marker rs5751876 is not in build_reference.py,
# so the ADORA2A variant carried by this project (rs553664) is used instead.
GENOSETS += [

    _gs("dgs045", "and(rs762551(C;C), rs553664(C;C))", 3.0, "Bad",
        "You are a slow caffeine metaboliser who also carries the caffeine-sensitive adenosine receptor genotype.",
        "CYP1A2 clears caffeine, and the slow genotype leaves it circulating for many hours; the ADORA2A "
        "genotype makes the receptor it acts on more likely to produce jitteriness and anxiety. Together "
        "they explain poor sleep from afternoon coffee and unpleasant rather than pleasant stimulation. "
        "Keeping intake below roughly 200 mg a day and stopping by early afternoon usually resolves both "
        "problems. Slow metabolisers also clear clozapine, olanzapine and theophylline more slowly.",
        "DETOX", "actionable", ["caffeine", "sleep", "anxiety"],
        ["clozapine", "olanzapine", "theophylline"], "replicated GWAS"),
]

# -- COMT / BDNF stress-response combinations ----------------------------
GENOSETS += [

    _gs("dgs046", "and(rs4680(A;A), rs6265(T))", 3.0, "",
        "You carry slow COMT together with at least one BDNF Met allele, a stress-sensitive combination.",
        "Slow COMT leaves more dopamine in the prefrontal cortex, which helps cognition when calm but "
        "tips into overload under pressure. The BDNF Met allele reduces activity-dependent release of the "
        "growth factor that supports adapting to that pressure. The combination is associated with sharper "
        "performance in low-stress conditions and worse performance and mood under sustained stress. "
        "Aerobic exercise raises BDNF and is the best-evidenced counterweight.",
        "NEURO", "informational", ["stress", "dopamine", "BDNF", "cognition"], [],
        "replicated candidate-gene studies"),

    _gs("dgs047", "and(rs4680(G;G), rs6265(C;C))", 2.0, "Good",
        "You carry fast COMT together with two BDNF Val alleles, a stress-resilient combination.",
        "Fast dopamine clearance plus normal activity-dependent BDNF release is associated with steadier "
        "performance under acute stress, at the cost of a slightly lower working-memory ceiling when "
        "conditions are calm. This is a trait pattern rather than a health risk, and effect sizes in "
        "individuals are small.",
        "NEURO", "informational", ["stress", "dopamine", "BDNF", "resilience"], [],
        "replicated candidate-gene studies"),

    _gs("dgs048", "and(rs4680(A;A), rs1360780(T))", 2.5, "",
        "You carry slow COMT together with the FKBP5 stress-reactivity variant.",
        "FKBP5 regulates how quickly the cortisol response switches off, and the risk allele slows that "
        "shutdown; slow COMT independently amplifies the catecholamine side of the stress response. "
        "Together they are associated with longer physiological recovery after stressors and, in people "
        "with significant early-life adversity, higher rates of anxiety and post-traumatic symptoms. This "
        "is a susceptibility pattern, not a diagnosis, and it responds to the usual levers of sleep, "
        "exercise and trauma-informed therapy where relevant.",
        "NEURO", "informational", ["stress", "cortisol", "HPA axis", "anxiety"], [],
        "replicated candidate-gene studies"),
]

# -- Celiac disease ------------------------------------------------------
GENOSETS += [

    _gs("dgs049", "rs2187668(T)", 4.0, "Bad",
        "You carry the HLA-DQ2.5 tag allele associated with celiac disease susceptibility.",
        "Almost everyone with celiac disease carries HLA-DQ2 or DQ8, so this marker is the main genetic "
        "gate on the condition. Carrying it is common and most carriers never develop celiac disease, but "
        "the negative predictive value is what makes the test useful: without it, celiac disease is very "
        "unlikely. If you have unexplained anaemia, diarrhoea, weight loss or persistent bloating, ask for "
        "tissue transglutaminase IgA testing while still eating gluten, since a pre-emptive gluten-free "
        "diet invalidates the result.",
        "INFLAM", "actionable", ["celiac", "HLA", "gluten", "autoimmunity"], [], "OMIM / HLA typing"),

    _gs("dgs050", "and(dgs049, rs231775(G))", 4.0, "Bad",
        "You carry the celiac HLA risk marker together with the CTLA4 autoimmunity variant.",
        "CTLA4 sets the threshold for T-cell activation, and its risk allele is associated with a cluster "
        "of autoimmune conditions including autoimmune thyroid disease, type 1 diabetes and celiac disease. "
        "Combined with HLA-DQ2.5 this suggests a general autoimmune predisposition rather than celiac "
        "disease specifically. Reasonable action is a low threshold for checking thyroid function and "
        "celiac serology if symptoms appear, not routine screening while well.",
        "INFLAM", "actionable", ["celiac", "autoimmunity", "thyroid", "HLA"], [],
        "replicated GWAS"),
]

# -- Nicotine dependence -------------------------------------------------
# NOTE: rs16969968 (CHRNA5) is not in build_reference.py; the CHRNA3 and
# CHRNA5 markers this project does carry are used instead.
GENOSETS += [

    _gs("dgs051", "atleast(2, rs1051730(A), rs1051740(C), rs1800497(A))", 3.0, "Bad",
        "You carry at least two variants linked to nicotine dependence and reward sensitivity.",
        "The nicotinic receptor variants in this cluster are among the most reproducible genetic "
        "predictors of how heavily someone smokes once they start, and the dopamine receptor variant adds "
        "general reward-seeking. This does not make cessation impossible, but it does predict a harder "
        "withdrawal and a higher relapse rate on willpower alone. Carriers do measurably better with "
        "pharmacotherapy, and varenicline outperforms nicotine replacement in this group.",
        "NEURO", "actionable", ["nicotine", "smoking", "addiction", "reward"],
        ["varenicline", "bupropion", "nicotine replacement"], "replicated GWAS"),

    _gs("dgs052", "and(rs1051730(A;A), rs1051740(C))", 3.5, "Bad",
        "You carry the homozygous CHRNA3 risk genotype together with the CHRNA5 variant on this panel.",
        "This is the heavier end of the nicotinic receptor risk spectrum: on average more cigarettes per "
        "day, earlier morning craving, and a higher smoking-attributable lung cancer risk than other "
        "smokers with the same pack-year history. The variant has no effect at all in people who never "
        "start. If you do smoke, combination pharmacotherapy plus behavioural support is the "
        "evidence-based route, and lung cancer screening eligibility is worth checking.",
        "NEURO", "actionable", ["nicotine", "smoking", "lung cancer", "addiction"],
        ["varenicline", "bupropion"], "replicated GWAS"),
]

# -- DNAInsight-authored combinations ------------------------------------
GENOSETS += [

    _gs("dgs053", "and(or(rs1801133(A), rs1801131(G)), or(rs1805087(G), rs1801394(G)))", 3.0, "Bad",
        "You carry reduced-function variants on both the folate and the B12 arms of the methylation cycle.",
        "MTHFR supplies the methyl group and MTR with MTRR hand it on using vitamin B12. A variant on only "
        "one arm is usually well compensated; variants on both make the pathway more dependent on adequate "
        "intake of both nutrients. Check homocysteine and B12 rather than guessing, and if either is off, "
        "methylfolate plus methylcobalamin is the sensible combination. High-dose folic acid alone can mask "
        "a B12 deficiency and is best avoided here.",
        "NEURO", "actionable", ["folate", "B12", "homocysteine", "methylation"], [],
        "OMIM / replicated GWAS"),

    _gs("dgs054", "and(rs1801133(A;A), rs601338(A))", 3.5, "Bad",
        "You combine homozygous MTHFR C677T with the FUT2 non-secretor allele.",
        "Slow folate activation raises the demand for B12-dependent remethylation, while FUT2 non-secretor "
        "status is associated with lower B12 absorption and a different gut microbiome. Together they make "
        "raised homocysteine more likely than either variant alone. Measure homocysteine and B12 once; if "
        "they are abnormal, methylcobalamin taken separately from meals plus methylfolate is the "
        "straightforward correction.",
        "METAB", "actionable", ["B12", "folate", "homocysteine", "microbiome"], [],
        "replicated GWAS"),

    _gs("dgs055", "atleast(2, rs1695(G;G), rs4880(T;T), rs1800566(T))", 2.5, "Bad",
        "You carry at least two reduced-capacity variants in the antioxidant and detoxification enzymes tested.",
        "GSTP1, SOD2 and NQO1 handle different steps of neutralising reactive oxygen species and quinones. "
        "Reduced capacity across more than one of them is associated with greater sensitivity to smoking, "
        "air pollution and some occupational exposures rather than with any specific disease. The useful "
        "response is exposure reduction and a diet genuinely rich in vegetables, not high-dose antioxidant "
        "supplements, which have repeatedly failed in trials.",
        "DETOX", "informational", ["oxidative stress", "detoxification", "pollution"], [],
        "replicated candidate-gene studies"),

    _gs("dgs056", "rs1050828(T)", 6.0, "Bad",
        "You carry the G6PD A- deficiency allele.",
        "G6PD protects red cells from oxidative stress, and deficiency means certain drugs and foods can "
        "trigger acute haemolysis: rasburicase, dapsone, primaquine and tafenoquine, nitrofurantoin, "
        "methylene blue, high-dose aspirin and fava beans are the classic triggers. Because the gene is on "
        "the X chromosome, males with one copy are fully affected while females vary. Add G6PD deficiency "
        "to your medical record and mention it before antimalarials, urate-lowering therapy for tumour "
        "lysis, or treatment for a urinary infection.",
        "PHARM", "pre_prescription", ["G6PD", "haemolysis", "antimalarial"],
        ["rasburicase", "dapsone", "primaquine", "nitrofurantoin", "methylene blue"], "CPIC Level A"),

    _gs("dgs057", "and(rs3135506(G), rs780094(T))", 3.0, "Bad",
        "You carry both the APOA5 and GCKR variants associated with higher fasting triglycerides.",
        "APOA5 slows triglyceride clearance and the GCKR variant pushes the liver toward converting sugar "
        "into fat, so the two act in the same direction. Carriers tend toward higher fasting triglycerides "
        "and a modestly higher chance of fatty liver, with severe hypertriglyceridaemia and pancreatitis "
        "possible at the extreme. Reducing alcohol, refined carbohydrate and fructose has a much larger "
        "effect here than reducing dietary fat; omega-3 at pharmacological dose is the usual add-on.",
        "METAB", "actionable", ["triglycerides", "fatty liver", "pancreatitis"],
        ["omega-3 ethyl esters", "fibrates"], "replicated GWAS"),

    _gs("dgs058", "and(rs174546(T;T), rs1799883(A))", 2.5, "",
        "You carry the low-conversion FADS1 genotype together with the FABP2 fat-absorption variant.",
        "FADS1 controls how efficiently plant-derived alpha-linolenic acid becomes EPA and DHA, and the "
        "low-activity genotype converts poorly; the FABP2 variant increases absorption of dietary fat. "
        "Practically, flaxseed and other plant omega-3 sources will not raise your EPA and DHA much, so "
        "oily fish or an algal or fish oil supplement is the more reliable route. Nothing here is harmful "
        "in itself.",
        "METAB", "informational", ["omega-3", "EPA", "DHA", "fat absorption"],
        ["fish oil", "algal DHA"], "replicated GWAS"),

    _gs("dgs059", "and(rs2228570(T), rs11568818(G))", 2.0, "",
        "You carry two vitamin D receptor variants associated with reduced receptor signalling.",
        "These variants change how strongly tissues respond to a given blood level of vitamin D rather "
        "than how much you make or absorb. Carriers may need a level toward the upper part of the normal "
        "range to get the same effect on bone and immune function. Measure 25-hydroxyvitamin D rather than "
        "supplementing blind, and re-check after three months on a fixed dose.",
        "METAB", "informational", ["vitamin D", "bone health", "immunity"], ["cholecalciferol"],
        "replicated candidate-gene studies"),

    _gs("dgs060",
        "atleast(3, rs1800795(G;G), rs1800629(A), rs1143627(C), rs1800896(A;A), rs2228145(C))",
        3.0, "Bad",
        "You carry at least three variants that shift your inflammatory set-point upward.",
        "These variants raise pro-inflammatory cytokine output or reduce the anti-inflammatory brake, and "
        "carrying several together is associated with a higher baseline hsCRP. Chronic low-grade "
        "inflammation is a shared driver of cardiovascular, metabolic and depressive risk, and it is "
        "measurable and modifiable. Ask for an hsCRP; if it is raised without an obvious cause, the "
        "highest-yield changes are treating sleep apnoea and periodontal disease, regular exercise, and a "
        "Mediterranean-pattern diet.",
        "INFLAM", "actionable", ["inflammation", "hsCRP", "cytokines"], [], "replicated GWAS"),

    _gs("dgs061", "and(rs4343(G;G), rs5186(C))", 2.0, "",
        "Your ACE and angiotensin receptor genotypes together describe a renin-angiotensin response profile.",
        "The ACE genotype is associated with higher circulating enzyme activity and the AGTR1 variant with "
        "a more reactive receptor, so this combination is associated with salt-sensitive blood pressure and "
        "a good response to drugs that block the pathway. Evidence is not strong enough to choose a drug on "
        "genotype alone, but it is a reasonable tiebreaker if an ACE inhibitor or an angiotensin receptor "
        "blocker is already under consideration, and it strengthens the case for reducing dietary sodium.",
        "CARDIO", "informational", ["blood pressure", "renin-angiotensin", "sodium"],
        ["lisinopril", "ramipril", "losartan", "valsartan"], "replicated candidate-gene studies"),

    _gs("dgs062", "and(rs1801260(C), rs10830963(G))", 2.5, "Bad",
        "You carry an evening-chronotype CLOCK variant together with the MTNR1B glucose variant.",
        "The MTNR1B variant prolongs melatonin's suppression of insulin release, so eating while melatonin "
        "is high produces an exaggerated glucose rise; the CLOCK variant pushes your natural body clock "
        "later, which makes late eating more likely. The combination is a genuine gene-behaviour "
        "interaction and it is fixable by timing rather than by diet composition: finish the last "
        "substantial meal at least three hours before sleep and keep a consistent wake time.",
        "METAB", "actionable", ["circadian rhythm", "glucose", "melatonin", "meal timing"], [],
        "replicated GWAS"),

    _gs("dgs063", "and(rs1799971(G), or(rs3892097(A), rs1065852(A)))", 3.5, "Bad",
        "You carry the opioid receptor variant together with a reduced-function CYP2D6 allele.",
        "Codeine and tramadol are prodrugs that CYP2D6 must convert into the active opioid, so a "
        "reduced-function allele can leave them ineffective for pain; the opioid receptor variant is "
        "separately associated with needing higher doses for the same analgesia. Together they predict poor "
        "response to codeine-type analgesics specifically. Morphine, oxycodone and non-opioid options are "
        "unaffected by the conversion step. Raise this before elective surgery or dental work.",
        "PHARM", "pre_prescription", ["opioid", "analgesia", "CYP2D6", "surgery"],
        ["codeine", "tramadol", "morphine", "oxycodone"], "CPIC Level A"),

    _gs("dgs064", "and(or(dgs032, dgs033), dgs043)", 5.5, "Bad",
        "You carry Factor V Leiden on top of a broader coronary and vascular risk burden.",
        "Venous and arterial thrombosis have different mechanisms, but carrying a strong venous "
        "thrombophilia alongside multiple arterial risk variants raises the overall vascular burden and "
        "makes situational risks such as surgery, immobility, oestrogen therapy and pregnancy more "
        "consequential. This combination is worth a single consolidated conversation with a clinician "
        "covering contraception, perioperative prophylaxis and cardiovascular prevention rather than "
        "treating each result separately.",
        "CARDIO", "actionable", ["thrombosis", "coronary artery disease", "surgery"], [],
        "replicated GWAS"),

    _gs("dgs065", "not(dgs010, dgs014, dgs020, dgs023, dgs025, dgs027, dgs056)", 1.5, "Good",
        "None of the high-impact pharmacogenomic risk genosets on this panel are present.",
        "You do not carry the tested variants for warfarin hypersensitivity, CYP2C19 poor metabolism, "
        "SLCO1B1 statin myopathy, TPMT deficiency, DPYD deficiency, UGT1A1 poor glucuronidation or G6PD "
        "deficiency. Standard dosing of the drugs affected by those genes is expected to behave normally. "
        "This is reassurance about a specific list, not a clean bill of pharmacogenomic health: consumer "
        "arrays miss rare alleles and do not cover CYP2D6 copy number, so clinical judgement still leads.",
        "PHARM", "pre_prescription", ["pharmacogenomics", "drug safety"], [], "CPIC Level A"),
]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

REQUIRED_KEYS = ("name", "aka", "criteria", "magnitude", "repute", "summary",
                 "interpretation", "category", "silo", "topics", "medicines",
                 "evidence")

REFERENCE_RSIDS = {row[0].lower() for row in REFERENCE}


def _number_of(name: str) -> int:
    digits = "".join(ch for ch in name if ch.isdigit())
    return int(digits) if digits else -1


def find_duplicates() -> list[str]:
    """Return genoset names that appear more than once in GENOSETS."""
    seen, dupes = set(), []
    for entry in GENOSETS:
        name = entry.get("name")
        if name in seen:
            dupes.append(name)
        seen.add(name)
    return dupes


def build_corpus() -> dict:
    """Return {name: entry} with the redundant 'name' key removed from values."""
    corpus = {}
    for entry in GENOSETS:
        body = {k: v for k, v in entry.items() if k != "name"}
        corpus.setdefault(entry["name"], body)
    return corpus


def validate() -> list[str]:
    """Check the corpus end to end. Prints and returns a list of error strings."""
    errors: list[str] = []

    for name in find_duplicates():
        errors.append(f"duplicate genoset name: {name}")

    names = {e["name"] for e in GENOSETS}

    for entry in GENOSETS:
        name = entry.get("name", "<unnamed>")
        for key in REQUIRED_KEYS:
            if key not in entry:
                errors.append(f"{name}: missing key {key!r}")
        if entry.get("category") not in VALID_CATEGORIES:
            errors.append(f"{name}: invalid category {entry.get('category')!r}")
        if entry.get("silo") not in VALID_SILOS:
            errors.append(f"{name}: invalid silo {entry.get('silo')!r}")
        if entry.get("repute") not in VALID_REPUTES:
            errors.append(f"{name}: invalid repute {entry.get('repute')!r}")
        mag = entry.get("magnitude")
        if not isinstance(mag, (int, float)) or not 0.0 <= float(mag) <= 10.0:
            errors.append(f"{name}: magnitude {mag!r} outside 0-10")
        if not (entry.get("summary") or "").strip():
            errors.append(f"{name}: empty summary")
        if not (entry.get("interpretation") or "").strip():
            errors.append(f"{name}: empty interpretation")
        if not (entry.get("evidence") or "").strip():
            errors.append(f"{name}: empty evidence")

        # Parse the criteria and check every rsID and every reference.
        try:
            node = parse_criteria(entry.get("criteria", ""))
        except CriteriaError as exc:
            errors.append(f"{name}: unparseable criteria: {exc}")
            continue

        for rsid in sorted(required_rsids(node)):
            if rsid not in REFERENCE_RSIDS:
                errors.append(f"{name}: rsID {rsid} is not in build_reference.REFERENCE")

        for ref in sorted(referenced_genosets(node)):
            if ref not in names:
                errors.append(f"{name}: references unknown genoset {ref}")
            elif _number_of(ref) >= _number_of(name):
                errors.append(f"{name}: may only reference lower-numbered genosets, got {ref}")

    # Cycle detection over the whole corpus.
    if not errors:
        try:
            order = topological_order(build_corpus())
        except CriteriaError as exc:
            errors.append(f"dependency ordering failed: {exc}")
        else:
            if len(order) != len(names):
                errors.append(f"topological_order returned {len(order)} of {len(names)} genosets")

    if errors:
        print(f"VALIDATION FAILED: {len(errors)} problem(s)")
        for err in errors:
            print(f"  - {err}")
    return errors


if __name__ == "__main__":
    if validate():
        sys.exit(1)

    corpus = build_corpus()
    payload = {
        "_meta": {
            "version":       GENOSET_VERSION,
            "genoset_count": len(corpus),
            "built_at":      datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "license":       "MIT",
            "note":          "Authored for DNAInsight. Not derived from SNPedia.",
        },
        "genosets": corpus,
    }
    out = Path(__file__).parent / "genosets.json"
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    print(f"Built genoset corpus v{GENOSET_VERSION}: {len(corpus)} genosets -> {out}")
