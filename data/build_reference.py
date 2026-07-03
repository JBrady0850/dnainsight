"""
build_reference.py -- Generates the bundled SNP reference database.

Run this script once to create data/snp_reference.json from the embedded
SNP table below. The output is a JSON dict keyed by rsID that the scanner
uses for offline-first annotation.

This ships as a pre-built file (snp_reference.json), so end users do NOT
need to run this script. Maintainers run it to update the bundled reference.

Usage:
    python data/build_reference.py
"""

import json
import datetime
from datetime import timezone
from pathlib import Path

REFERENCE_VERSION = "1.2.0"

# ---------------------------------------------------------------------------
# Master SNP reference table
# Format: rsid | gene | category | clinical_sig | interpretation
#
# Categories: PHARM, METAB, INFLAM, NEURO, DETOX, CARDIO
# Silos: pre_prescription (drug-safety), actionable (health action needed),
#         informational (background)
# ---------------------------------------------------------------------------

REFERENCE = [
    # PHARMACOGENOMICS -- Prescription-critical
    ("rs429358",   "APOE",    "PHARM",  "risk factor",    "APOE e4 allele: increased Alzheimer disease risk and altered lipid metabolism. Statins may have reduced efficacy. Discuss with cardiologist before lipid-lowering therapy."),
    ("rs7412",     "APOE",    "PHARM",  "risk factor",    "APOE e2 allele marker: associated with lower LDL but elevated triglycerides. Modifies cardiovascular and neurological risk profile."),
    ("rs1799971",  "CYP2D6",  "PHARM",  "drug response",  "CYP2D6 reduced function: impacts metabolism of codeine, tramadol, oxycodone, tamoxifen, antidepressants (fluoxetine, paroxetine), antipsychotics, and beta-blockers. Alert prescriber before use."),
    ("rs1065852",  "CYP2D6",  "PHARM",  "drug response",  "CYP2D6 *10 reduced function allele: decreased enzyme activity for pain medications, antidepressants, and opioids. Dose adjustment may be required."),
    ("rs3892097",  "CYP2D6",  "PHARM",  "drug response",  "CYP2D6 *4 non-functional allele (most common European PM allele): zero enzyme contribution from this allele. Significant impact on opioid efficacy and antidepressant dosing."),
    ("rs5030658",  "CYP2D6",  "PHARM",  "drug response",  "CYP2D6 *6 null allele: no functional enzyme from this allele. Combined with other reduced-function alleles may result in poor metabolizer status."),
    ("rs1057910",  "CYP2C9",  "PHARM",  "drug response",  "CYP2C9 *3 allele: significantly reduced metabolism of warfarin, NSAIDs, phenytoin, and sulfonylureas. Requires dose reduction for warfarin; risk of drug toxicity."),
    ("rs1799853",  "CYP2C9",  "PHARM",  "drug response",  "CYP2C9 *2 allele: intermediate metabolism of warfarin and NSAIDs. Warfarin dose adjustment recommended. Monitor INR closely."),
    ("rs1057911",  "CYP2C9",  "PHARM",  "drug response",  "CYP2C9 variant: slower clearance of anti-inflammatory drugs and oral hypoglycemics. Dose reduction may be needed."),
    ("rs4244285",  "CYP2C19", "PHARM",  "drug response",  "CYP2C19 *2 loss-of-function allele: poor metabolism of clopidogrel (antiplatelet), PPIs, SSRIs, and TCAs. Clopidogrel may be ineffective -- consider alternative antiplatelet agent."),
    ("rs12248560", "CYP2C19", "PHARM",  "drug response",  "CYP2C19 *17 ultrarapid metabolizer allele: accelerated clearance of PPIs, SSRIs, and clopidogrel. PPIs may fail at standard dose; dose increase or switch to rabeprazole warranted."),
    ("rs1058204",  "CYP2C19", "PHARM",  "drug response",  "CYP2C19 variant: influences proton pump inhibitor effectiveness and psychiatric medication response."),
    ("rs28399504", "CYP2C19", "PHARM",  "drug response",  "CYP2C19 *4 loss-of-function allele: poor metabolizer for clopidogrel, PPIs, SSRIs, and TCAs. CPIC Level A: consider alternative antiplatelet therapy; adjust GI and psychiatric medication dosing accordingly."),
    ("rs41291556", "CYP2C19", "PHARM",  "drug response",  "CYP2C19 *8 reduced-function allele: partial loss of enzyme activity. Contributes to intermediate or poor metabolizer phenotype when combined with other CYP2C19 loss-of-function alleles."),
    ("rs72549354", "TPMT",    "PHARM",  "drug response",  "TPMT *3A variant: severely reduced thiopurine methyltransferase activity. Thiopurine drugs (azathioprine, 6-mercaptopurine, thioguanine) carry life-threatening myelosuppression risk. CPIC Level A: 90% dose reduction or switch to alternative therapy is mandatory."),
    ("rs1142345",  "TPMT",    "PHARM",  "drug response",  "TPMT *3C allele: reduced-function variant for thiopurine metabolism; most common reduced-function allele in East Asian and African populations. CPIC requires preemptive genotyping before thiopurine initiation; dose reduction is warranted."),
    ("rs2108622",  "CYP4F2",  "PHARM",  "drug response",  "CYP4F2 *3 allele: reduced vitamin K metabolism, leading to higher warfarin sensitivity. Requires lower initial warfarin dose."),
    ("rs9923231",  "VKORC1",  "PHARM",  "drug response",  "VKORC1 variant (rs9923231): lower enzyme activity, significantly increased warfarin and anticoagulant sensitivity. Reduce starting dose; use CPIC warfarin dosing algorithm."),
    ("rs4149056",  "SLCO1B1", "PHARM",  "drug response",  "SLCO1B1 *5 allele: reduced statin transport into hepatocytes, leading to higher plasma statin concentrations and myopathy/rhabdomyolysis risk. Use lowest effective statin dose; consider pravastatin."),
    ("rs2241766",  "CYP3A4",  "PHARM",  "drug response",  "CYP3A4 *1B variant: altered metabolism of immunosuppressants, statins, antibiotics, and hormonal medications."),
    ("rs2740574",  "CYP3A4",  "PHARM",  "drug response",  "CYP3A4 variant: impacts metabolism of many common antibiotics and statins. Review drug-drug interactions carefully."),
    ("rs776746",   "CYP3A5",  "PHARM",  "drug response",  "CYP3A5 *3 non-expresser allele: standard tacrolimus and cyclosporine starting dose for non-expressers. Transplant patients: use clinical dosing protocol."),
    ("rs1045642",  "ABCB1",   "PHARM",  "drug response",  "ABCB1 (P-gp) variant: altered multidrug transporter efficiency affecting absorption of digoxin, antiretrovirals, and chemotherapy agents."),
    ("rs2032582",  "ABCB1",   "PHARM",  "drug response",  "ABCB1 variant: affects blood-brain barrier drug transport, impacting CNS drug concentrations."),
    ("rs1128503",  "ABCB1",   "PHARM",  "drug response",  "ABCB1 variant: impacts absorption of multiple clinical medications including antiepileptics and opioids."),
    ("rs762551",   "CYP1A2",  "DETOX",  "drug response",  "CYP1A2 fast/slow metabolizer: affects caffeine, clozapine, olanzapine, theophylline metabolism. Slow metabolizers: higher caffeine and drug exposure."),
    ("rs10455872", "LPA",     "CARDIO", "risk factor",    "LPA variant: elevated lipoprotein(a) levels, a significant independent cardiovascular risk factor. Discuss with cardiologist; niacin or PCSK9 inhibitor may be indicated."),
    ("rs1799983",  "NOS3",    "CARDIO", "drug response",  "NOS3 (eNOS) variant: reduced nitric oxide production, impacting blood pressure regulation and cardiovascular drug response (ACE inhibitors, nitrates)."),
    ("rs1050102",  "ACE",     "PHARM",  "drug response",  "ACE gene variant: influences response to ACE inhibitor medications (lisinopril, enalapril) and cardiovascular risk."),
    ("rs4343",     "ACE",     "PHARM",  "drug response",  "ACE D/I polymorphism: linked to ACE inhibitor drug response and hypertension management. DD genotype associated with higher ACE activity."),
    ("rs1801275",  "GNB3",    "PHARM",  "drug response",  "GNB3 C825T variant: influences blood pressure response to thiazide diuretics and other antihypertensives."),
    ("rs20455",    "KIF6",    "PHARM",  "drug response",  "KIF6 variant: associated with differential statin benefit for coronary heart disease risk reduction."),

    # PHARMACOGENOMICS -- CPIC Level A: Chemotherapy & Specialty Drugs
    ("rs3918290",  "DPYD",    "PHARM",  "drug response",  "DPYD *2A (IVS14+1G>A): complete loss of dihydropyrimidine dehydrogenase function. Fluoropyrimidine chemotherapy (5-fluorouracil, capecitabine) is contraindicated or requires 50%+ dose reduction. CPIC Level A: preemptive genotyping mandatory before treatment."),
    ("rs55886062", "DPYD",    "PHARM",  "drug response",  "DPYD *13 (c.1679T>G): severely reduced DPD activity. Fluoropyrimidine-induced life-threatening toxicity (mucositis, myelosuppression, neurotoxicity) risk. CPIC Level A: alternative therapy or 50% dose reduction required."),
    ("rs67376798", "DPYD",    "PHARM",  "drug response",  "DPYD c.2846A>T: reduced DPD enzyme activity. Increases risk of severe 5-FU and capecitabine toxicity. CPIC Level A: 25% starting dose reduction recommended for heterozygous carriers."),
    ("rs887829",   "UGT1A1",  "PHARM",  "drug response",  "UGT1A1 *80 (*28 tag SNP): reduced bilirubin and irinotecan glucuronidation. Homozygotes at elevated risk for irinotecan-induced severe neutropenia and diarrhea. CPIC Level A: dose reduction warranted in *28/*28 patients on FOLFIRI regimens."),
    ("rs4148323",  "UGT1A1",  "PHARM",  "drug response",  "UGT1A1 *6 (G71R): loss-of-function variant predominant in East Asian populations (~13%). Causes unconjugated hyperbilirubinemia and increases irinotecan toxicity risk. CPIC Level A: dose reduction recommended."),
    ("rs116855232","NUDT15",  "PHARM",  "drug response",  "NUDT15 *3 (Arg139Cys): dramatically reduced nucleoside diphosphate enzyme activity. Thiopurine drugs (azathioprine, 6-mercaptopurine) cause severe myelosuppression in heterozygotes and are contraindicated in homozygotes. CPIC Level A; especially prevalent in East Asian and Hispanic populations."),
    ("rs1050828",  "G6PD",    "PHARM",  "drug response",  "G6PD A- variant: glucose-6-phosphate dehydrogenase deficiency. Oxidant drugs (rasburicase, dapsone, primaquine, high-dose aspirin) and fava beans can trigger acute hemolytic anemia. CPIC Level A: these agents are contraindicated; inform all prescribers of G6PD status."),
    ("rs4986893",  "CYP2C19", "PHARM",  "drug response",  "CYP2C19 *3 null allele: premature stop codon, zero enzyme activity. Common in East Asian populations (~5%). Contributes to poor metabolizer phenotype for clopidogrel, PPIs, and SSRIs when combined with other loss-of-function alleles."),
    ("rs28371706", "CYP2C9",  "PHARM",  "drug response",  "CYP2C9 *5 allele: reduced-function variant more prevalent in individuals of African ancestry. Reduces clearance of warfarin, NSAIDs, and sulfonylureas; requires dose reduction and close INR monitoring when initiating anticoagulation therapy."),
    ("rs12979860", "IFNL3",   "PHARM",  "drug response",  "IFNL3 (IL28B) CC genotype: favorable interferon-lambda response for hepatitis C treatment. CC carriers achieve significantly higher sustained virologic response rates with pegIFN/ribavirin; also influences COVID-19 severity outcomes per replicated GWAS data."),

    ("rs6025",     "F5",      "CARDIO", "risk factor",    "Factor V Leiden (rs6025, c.1691G>A): most common inherited thrombophilia. Heterozygotes: 4-8x elevated venous thrombosis risk; homozygotes: 50-80x risk. Clinically significant with estrogen use, pregnancy, or surgery. Discuss anticoagulation strategy with physician."),
    ("rs1799963",  "F2",      "CARDIO", "risk factor",    "Prothrombin G20210A: elevated prothrombin levels and venous clotting risk. Heterozygotes: 2-3x increased thrombosis risk. Risk compounds when co-present with Factor V Leiden (rs6025). Review contraceptive options and surgical anticoagulation prophylaxis with a physician."),

    # PHARMACOGENOMICS -- Informational
    ("rs1056836",  "CYP1B1",  "PHARM",  "drug response",  "CYP1B1 variant: phase-1 enzyme for estrogen and environmental toxin activation. Relevant for estrogen-based therapies and cancer risk screening."),
    ("rs2228570",  "VDR",     "METAB",  "drug response",  "VDR (vitamin D receptor) variant: altered vitamin D sensitivity. May require higher supplementation to achieve target 25-OH-D levels."),
    ("rs11568818", "VDR",     "METAB",  "drug response",  "VDR secondary variant: secondary marker for vitamin D receptor signaling. Consider vitamin D monitoring."),

    # METABOLIC HEALTH
    ("rs9939609",  "FTO",     "METAB",  "risk factor",    "FTO A allele (obesity risk): master regulator of appetite and fat storage. Associated with 20-30% higher obesity risk. High-protein, lower-carb diet and regular aerobic exercise shown to attenuate genetic effect."),
    ("rs1421085",  "FTO",     "METAB",  "risk factor",    "FTO rs1421085 T>C: 'master switch' for thermogenesis and fat browning. C/C genotype associated with reduced beige fat activation and higher obesity susceptibility."),
    ("rs8050136",  "FTO",     "METAB",  "risk factor",    "FTO secondary obesity marker: preference for high-calorie foods and reduced satiety signaling. Caloric density awareness and mindful eating strategies recommended."),
    ("rs7903146",  "TCF7L2",  "METAB",  "risk factor",    "TCF7L2 T allele: highest-impact common genetic marker for type 2 diabetes risk. Associated with impaired insulin secretion. Prioritize low-glycemic diet, regular exercise, and HbA1c monitoring."),
    ("rs12255372", "TCF7L2",  "METAB",  "risk factor",    "TCF7L2 secondary diabetes marker: correlated with reduced insulin production capacity."),
    ("rs1801282",  "PPARG",   "METAB",  "risk factor",    "PPARG Pro12Ala: the Ala allele is associated with improved insulin sensitivity and reduced T2DM risk. Pro/Pro genotype linked to higher adipose insulin resistance."),
    ("rs10811661", "CDKN2A",  "METAB",  "risk factor",    "CDKN2A/B variant: reduced pancreatic beta-cell regeneration capacity and insulin production. Associated with T2DM risk."),
    ("rs13266634", "SLC30A8", "METAB",  "risk factor",    "SLC30A8 variant: altered zinc transport in insulin-secreting beta cells. Low zinc intake may amplify risk; zinc supplementation worth discussing with provider."),
    ("rs17782313", "MC4R",    "METAB",  "risk factor",    "MC4R variant: reduced melanocortin receptor signaling, impairing energy homeostasis and satiety. Associated with higher BMI and binge-eating tendency."),
    ("rs2943641",  "IRS1",    "METAB",  "risk factor",    "IRS1 variant: reduced insulin receptor substrate signaling, contributing to muscle insulin resistance."),
    ("rs5219",     "KCNJ11",  "METAB",  "drug response",  "KCNJ11 E23K variant: linked to differential response to sulfonylurea diabetes medications and T2DM susceptibility."),
    ("rs1800644",  "ADRB2",   "METAB",  "drug response",  "ADRB2 Arg16Gly: beta-2 receptor variant influencing fat mobilization, asthma medication response, and metabolic rate."),
    ("rs1042713",  "ADRB2",   "PHARM",  "drug response",  "ADRB2 variant: influences response to asthma medications (albuterol, salmeterol) and cardiovascular beta-agonists."),
    ("rs1042714",  "ADRB2",   "PHARM",  "drug response",  "ADRB2 Glu27Gln: affects metabolic rate and exercise-induced fat burning. Common in obesity resistance."),
    ("rs1800544",  "ADRB3",   "METAB",  "risk factor",    "ADRB3 Trp64Arg variant: associated with abdominal fat storage and reduced basal metabolic rate in some ethnic groups."),
    ("rs7501331",  "BCMO1",   "METAB",  "risk factor",    "BCMO1 variant: reduced beta-carotene to vitamin A conversion. Increased dietary preformed vitamin A (retinol) may be needed."),
    ("rs601338",   "FUT2",    "METAB",  "risk factor",    "FUT2 secretor status: non-secretors have reduced gut B12 absorption and altered microbiome composition. Monitor B12; consider methylcobalamin supplementation."),
    ("rs1799883",  "FABP2",   "METAB",  "risk factor",    "FABP2 Ala54Thr: Thr allele associated with higher postprandial fat absorption and insulin resistance."),
    ("rs1800588",  "LIPC",    "METAB",  "risk factor",    "LIPC variant: hepatic lipase activity affecting HDL cholesterol metabolism and cardiovascular risk."),
    ("rs1137101",  "LEPR",    "METAB",  "risk factor",    "Leptin receptor variant: reduced satiety signaling efficiency. Protein-rich diet and resistance training may partially compensate."),
    ("rs1800562",  "HFE",     "METAB",  "risk factor",    "HFE C282Y: primary hereditary hemochromatosis mutation. Homozygotes face high risk of iron overload, liver cirrhosis, and cardiomyopathy if untreated. Obtain serum ferritin and transferrin saturation testing; therapeutic phlebotomy is indicated if iron is elevated."),
    ("rs1799945",  "HFE",     "METAB",  "risk factor",    "HFE H63D: milder hemochromatosis variant. Compound heterozygotes (H63D + C282Y) carry elevated iron overload risk. Monitor ferritin and iron studies annually; avoid vitamin C supplements taken with iron-rich meals."),
    ("rs10830963", "MTNR1B",  "METAB",  "risk factor",    "MTNR1B (melatonin receptor 1B) risk allele: impairs glucose-stimulated insulin secretion and elevates fasting glucose. Carriers show higher T2D risk, especially with late eating patterns. Aligning meal timing with circadian rhythms (earlier dinner, consistent sleep schedule) reduces the phenotypic impact."),
    ("rs780094",   "GCKR",    "METAB",  "risk factor",    "GCKR (glucokinase regulatory protein) variant: increases hepatic glucose uptake but raises fasting triglycerides. Carriers show elevated plasma triglycerides and modestly increased NAFLD risk. Reduce refined carbohydrate and fructose intake; monitor fasting lipid panel."),
    ("rs174546",   "FADS1",   "METAB",  "risk factor",    "FADS1 (fatty acid desaturase 1) variant: reduced delta-5-desaturase efficiency, impairing conversion of ALA to EPA/DHA. Carriers may benefit from preformed omega-3 supplementation (fish oil) rather than flaxseed ALA sources alone."),
    ("rs3135506",  "APOA5",   "METAB",  "risk factor",    "APOA5 variant: associated with elevated fasting triglycerides and reduced lipoprotein lipase activity. Carriers are at higher risk for hypertriglyceridemia and pancreatitis with high-fat or high-carbohydrate diets. Omega-3 supplementation and alcohol limitation are evidence-based interventions."),
    ("rs1333049",  "CDKN2B",  "CARDIO", "risk factor",    "9p21 locus (CDKN2B-AS1 region): the strongest and most replicated common genetic risk locus for coronary artery disease (CAD) in genome-wide association studies. Risk allele approximately doubles lifetime CAD risk independent of traditional risk factors. Prioritize aggressive primary prevention: statin therapy discussion, blood pressure control, smoking cessation."),
    ("rs5186",     "AGTR1",   "CARDIO", "drug response",  "AGTR1 A1166C (angiotensin II receptor type 1): C allele associated with higher vascular reactivity, hypertension susceptibility, and enhanced response to ARB medications (losartan, valsartan). CC homozygotes may achieve greater blood pressure reduction with ARBs versus ACE inhibitors."),
    ("rs4988235",  "LCT",     "METAB",  "informational",  "LCT -13910C>T (lactase persistence): T allele confers continued lactase expression into adulthood. C/C homozygotes are genetically lactose intolerant and may benefit from lactase enzyme supplements or lactose-reduced dietary substitutions to avoid GI symptoms."),

    # INFLAMMATION
    ("rs1800795",  "IL6",     "INFLAM", "risk factor",    "IL-6 promoter variant: elevated systemic inflammation and cytokine production. Monitor hsCRP; anti-inflammatory diet (Mediterranean), omega-3 supplementation, and exercise recommended."),
    ("rs1800629",  "TNFA",    "INFLAM", "risk factor",    "TNF-alpha -308G>A: increased TNF-alpha production, promoting chronic low-grade inflammation. Higher risk of autoimmune and metabolic conditions."),
    ("rs1205",     "CRP",     "INFLAM", "risk factor",    "CRP genetic variant: influences baseline C-reactive protein levels. Elevated baseline CRP predicts cardiovascular and metabolic risk."),
    ("rs1143627",  "IL1B",    "INFLAM", "risk factor",    "IL-1beta variant: pro-inflammatory cytokine production. Associated with joint inflammation, metabolic syndrome, and mood dysregulation."),
    ("rs1800871",  "IL10",    "INFLAM", "risk factor",    "IL-10 -592A>C: variant affecting anti-inflammatory cytokine balance. Lower IL-10 associated with higher inflammatory set-point."),
    ("rs1800896",  "IL10",    "INFLAM", "risk factor",    "IL-10 -1082A>G: regulatory variant for immune system dampening. Impacts autoimmune susceptibility."),
    ("rs1800627",  "IL10",    "INFLAM", "risk factor",    "IL-10 variant: promotes balanced inflammatory response; altered alleles linked to chronic inflammatory conditions."),
    ("rs1800872",  "IL10",    "INFLAM", "risk factor",    "IL-10 secondary variant: genetic predisposition to higher inflammatory set-point."),
    ("rs231775",   "CTLA4",   "INFLAM", "risk factor",    "CTLA4 Thr17Ala: T cell activation regulator. Associated with autoimmune susceptibility (thyroid disease, type 1 diabetes, rheumatoid arthritis, Celiac disease)."),
    ("rs2069762",  "IL2",     "INFLAM", "risk factor",    "IL-2 variant: T cell growth and immune signaling. Relevant for autoimmune condition risk and immune modulation."),
    ("rs1800450",  "MBL2",    "INFLAM", "risk factor",    "MBL2 (mannose-binding lectin): innate immune defense variant. Low MBL2 associated with increased infection susceptibility."),
    ("rs4586",     "CCL2",    "INFLAM", "risk factor",    "CCL2 (MCP-1) variant: monocyte chemoattractant protein influencing vessel inflammation and atherosclerosis risk."),
    ("rs2228145",  "IL6R",    "INFLAM", "drug response",  "IL-6 receptor Asp358Ala: reduced IL-6 receptor shedding. Predictive marker for tocilizumab (IL-6 blocker) response in rheumatoid arthritis."),
    ("rs30187",    "CRP",     "INFLAM", "risk factor",    "CRP variant: baseline systemic inflammation predictor. Elevated values correlate with cardiovascular and metabolic risk."),
    ("rs2187668",  "HLA-DQA1","INFLAM", "risk factor",    "HLA-DQA1 rs2187668 (DQ2.5 haplotype marker): associated with celiac disease susceptibility. HLA-DQ2 haplotype carriers have ~3-5% lifetime celiac risk; risk rises with first-degree relatives affected. Negative result substantially reduces celiac likelihood. If symptomatic: serum tTG-IgA testing and gastroenterology referral recommended."),

    # NEUROLOGICAL
    ("rs1801133",  "MTHFR",   "NEURO",  "risk factor",    "MTHFR C677T: reduced methylenetetrahydrofolate reductase activity. Impairs folate conversion and homocysteine clearance. Elevated homocysteine raises cardiovascular and neurological risk. Supplementation: methylfolate (5-MTHF), NOT folic acid."),
    ("rs1801131",  "MTHFR",   "NEURO",  "risk factor",    "MTHFR A1298C: secondary folate metabolism variant. Combined with C677T creates compound heterozygosity. Check homocysteine level; methylated B-vitamin supplementation recommended."),
    ("rs4680",     "COMT",    "NEURO",  "informational",  "COMT Val158Met: dictates dopamine and estrogen clearance. Val/Val (fast COMT): lower dopamine in prefrontal cortex, better stress resilience but lower working memory ceiling. Met/Met (slow COMT): higher dopamine, better cognition under low stress, worse under high stress."),
    ("rs1805087",  "MTR",     "NEURO",  "risk factor",    "MTR (methionine synthase): converts homocysteine to methionine. Variant reduces B12-dependent remethylation; monitor homocysteine and B12."),
    ("rs1801394",  "MTRR",    "NEURO",  "risk factor",    "MTRR (methionine synthase reductase): recycles vitamin B12 cofactor. Reduced efficiency compounds MTHFR-related methylation concerns."),
    ("rs6265",     "BDNF",    "NEURO",  "informational",  "BDNF Val66Met: brain-derived neurotrophic factor variant. Met carriers show reduced activity-dependent BDNF secretion, linked to lower resilience to stress, depression susceptibility, and impaired fear extinction. Regular aerobic exercise is the most evidence-supported intervention."),
    ("rs6323",     "MAOA",    "NEURO",  "informational",  "MAOA (monoamine oxidase A): regulates serotonin and norepinephrine clearance. High-activity allele: faster neurotransmitter breakdown; lower-activity: slower clearance with implications for mood regulation."),
    ("rs1800497",  "ANKK1",   "NEURO",  "informational",  "DRD2/ANKK1 Taq1A: reduces dopamine receptor density in reward pathways. Associated with reward sensitivity, addiction vulnerability, and preference for high-stimulation behaviors."),
    ("rs1611115",  "DBH",     "NEURO",  "informational",  "DBH (dopamine beta-hydroxylase): converts dopamine to norepinephrine. Variant associated with altered dopamine-to-norepinephrine ratio and ADHD-related traits."),
    ("rs553664",   "ADORA2A", "NEURO",  "informational",  "ADORA2A (adenosine A2A receptor): caffeine sensitivity and anxiety. C/C genotype: increased caffeine-induced anxiety; consider limiting caffeine to <200mg/day."),
    ("rs25531",    "SLC6A4",  "NEURO",  "informational",  "5-HTTLPR (serotonin transporter): short allele associated with reduced serotonin transporter expression, increased stress reactivity, and depression risk under adverse conditions. Note: SSRI response is NOT reliably predicted by this SNP."),
    ("rs6295",     "HTR1A",   "NEURO",  "informational",  "HTR1A (serotonin 1A receptor): G allele associated with altered receptor sensitivity. Impacts emotional regulation and social anxiety."),
    ("rs1360780",  "FKBP5",   "NEURO",  "risk factor",    "FKBP5 variant: high-impact cortisol regulation marker. Associated with PTSD susceptibility, glucocorticoid resistance, and stress response dysregulation."),
    ("rs1800260",  "CLOCK",   "NEURO",  "informational",  "CLOCK gene variant: circadian rhythm and sleep cycle timing. Associated with eveningness (night-owl chronotype) and increased risk of sleep disorders and mood dysregulation."),
    ("rs1801260",  "CLOCK",   "NEURO",  "informational",  "CLOCK circadian rhythm variant: determines peak performance timing, sleep architecture, and metabolic timing of meals."),
    ("rs110402",   "CRHR1",   "NEURO",  "risk factor",    "CRHR1 (corticotropin-releasing hormone receptor): influences HPA axis stress response and resilience to trauma."),
    ("rs3800373",  "FKBP5",   "NEURO",  "risk factor",    "FKBP5 secondary variant: linked to mood stability under chronic stress and glucocorticoid sensitivity."),
    ("rs13235612", "GAD1",    "NEURO",  "risk factor",    "GAD1 (glutamate decarboxylase): GABA synthesis enzyme. Variant associated with anxiety vulnerability and panic disorder risk."),
    ("rs1051740",  "CHRNA5",  "NEURO",  "risk factor",    "CHRNA5 (nicotinic acetylcholine receptor): strongly associated with nicotine dependence and lung cancer risk in smokers."),
    ("rs1655991",  "COMT",    "NEURO",  "informational",  "COMT secondary variant: influences cognitive processing speed and emotional regulation; modulates primary COMT effect."),
    ("rs53576",    "OXTR",    "NEURO",  "informational",  "Oxytocin receptor (OXTR) rs53576: GG homozygotes demonstrate higher empathy, social sensitivity, and prosocial behavior. AA carriers show lower social bonding and modestly elevated risk for autism spectrum traits. Context-dependent; not diagnostic."),
    ("rs324420",   "FAAH",    "NEURO",  "informational",  "FAAH (fatty acid amide hydrolase) C385A: the A allele reduces FAAH enzyme activity, increasing endocannabinoid tone (anandamide levels). Associated with lower anxiety, pain sensitivity, and stress reactivity. May influence response to cannabinoid-based therapies and analgesics."),
    ("rs6277",     "DRD2",    "NEURO",  "informational",  "DRD2 C957T (rs6277): synonymous variant affecting dopamine D2 receptor mRNA stability and translation. T allele associated with lower striatal D2 receptor density, reduced reward sensitivity, and altered response to antipsychotic medications and dopaminergic therapies for Parkinson disease."),
    ("rs1761667",  "CD36",    "NEURO",  "informational",  "CD36 rs1761667: fatty acid translocase variant affecting oral fat taste perception. AA genotype reduces sensitivity to dietary fat; associated with higher fat intake and preference for fatty foods. Awareness of this preference may support personalized dietary counseling."),

    # DETOX / OXIDATIVE STRESS
    ("rs1695",     "GSTP1",   "DETOX",  "risk factor",    "GSTP1 Ile105Val: reduced glutathione-S-transferase activity. Decreased detoxification of carcinogens and chemotherapy agents. Increased sensitivity to oxidative stress and environmental pollutants."),
    ("rs4880",     "SOD2",    "DETOX",  "risk factor",    "SOD2 Ala16Val: mitochondrial superoxide dismutase. Val/Val (homozygous T on forward strand): reduced mitochondrial import and antioxidant capacity. Consider CoQ10 and mitochondrial antioxidant support."),
    ("rs1800566",  "NQO1",    "DETOX",  "risk factor",    "NQO1 Pro187Ser: near-zero enzyme activity for homozygous Ser/Ser. Impaired quinone reduction and CoQ10 recycling. Avoid benzene exposure; may increase sensitivity to certain chemotherapy drugs."),
    ("rs1799807",  "ADH1B",   "DETOX",  "informational",  "ADH1B (alcohol dehydrogenase): determines speed of alcohol-to-acetaldehyde conversion. Fast convertors accumulate acetaldehyde more rapidly."),
    ("rs671",      "ALDH2",   "DETOX",  "risk factor",    "ALDH2 *2 allele: severely impaired acetaldehyde clearance. Alcohol causes rapid acetaldehyde buildup (flushing, nausea). Significantly elevated esophageal cancer risk with alcohol use. Alcohol avoidance is the evidence-based recommendation."),
    ("rs2228001",  "XPC",     "DETOX",  "risk factor",    "XPC (DNA repair): nucleotide excision repair capacity. Variant associated with reduced UV-induced DNA repair and elevated skin cancer risk."),
    ("rs1051730",  "CHRNA3",  "DETOX",  "risk factor",    "CHRNA3 rs1051730: nicotinic acetylcholine receptor subunit alpha-3. Risk allele strongly associated with nicotine dependence, heavier daily smoking, and lung cancer risk. Carriers typically require higher-intensity cessation support; varenicline is the preferred pharmacotherapy per clinical evidence."),
]

# ---------------------------------------------------------------------------
# Build the reference dict
# ---------------------------------------------------------------------------

def find_duplicates() -> list[str]:
    """Return rsIDs that appear more than once in the REFERENCE table."""
    seen, dupes = set(), []
    for row in REFERENCE:
        rsid = row[0]
        if rsid in seen:
            dupes.append(rsid)
        seen.add(rsid)
    return dupes


def build_reference() -> dict:
    ref = {}
    for rsid, gene, category, clinical_sig, interpretation in REFERENCE:
        if rsid not in ref:
            ref[rsid] = {
                "gene":          gene,
                "category":      category,
                "clinical_sig":  clinical_sig,
                "interpretation": interpretation,
                "conditions":    interpretation,
            }
    return ref


if __name__ == "__main__":
    out = Path(__file__).parent / "snp_reference.json"

    dupes = find_duplicates()
    if dupes:
        print(f"WARNING: {len(dupes)} duplicate rsID(s) in REFERENCE table "
              f"(first occurrence kept): {', '.join(sorted(set(dupes)))}")

    ref = build_reference()
    payload = {
        "_meta": {
            "version":    REFERENCE_VERSION,
            "snp_count":  len(ref),
            "built_at":   datetime.datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sources": [
                "CPIC guidelines — https://cpicpgx.org/guidelines/",
                "PharmGKB — https://www.pharmgkb.org/",
                "ClinVar — https://www.ncbi.nlm.nih.gov/clinvar/",
                "GWAS Catalog — https://www.ebi.ac.uk/gwas/",
            ],
        },
        "snps": ref,
    }
    with open(out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"Built bundled reference v{REFERENCE_VERSION}: {len(ref)} SNPs -> {out}")
