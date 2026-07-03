"""
doctor_report.py -- Generate the Doctor Discussion Report.

Produces a clinical-style HTML document designed to be brought to a
physician, pharmacist, or genetic counselor appointment.

Content:
  - Patient summary block
  - Prescription-critical findings table (flagged for prescriber review)
  - Actionable findings (monitoring, supplements, lifestyle)
  - Drug interaction warnings organized by gene/drug class
  - Recommended lab follow-ups
  - AI analysis prompt block (for use with Grok or other LLMs)
"""

from datetime import datetime


DRUG_CLASS_MAP = {
    "CYP2D6":  "Antidepressants, antipsychotics, opioids, beta-blockers, tamoxifen",
    "CYP2C19": "PPIs, clopidogrel, SSRIs, tricyclic antidepressants, voriconazole",
    "CYP2C9":  "Warfarin, NSAIDs, phenytoin, sulfonylureas",
    "CYP3A4":  "Immunosuppressants, statins, HIV medications, benzodiazepines",
    "CYP3A5":  "Tacrolimus, cyclosporine, sirolimus",
    "CYP2B6":  "Bupropion, efavirenz, methadone, ketamine",
    "CYP1A2":  "Caffeine, clozapine, olanzapine, theophylline, fluvoxamine",
    "VKORC1":  "Warfarin, acenocoumarol, phenprocoumon (vitamin K antagonists)",
    "SLCO1B1": "Statins — simvastatin, atorvastatin, rosuvastatin (myopathy risk)",
    "TPMT":    "Azathioprine, 6-mercaptopurine, thioguanine (hematologic toxicity)",
    "DPYD":    "Fluorouracil, capecitabine (severe 5-FU toxicity)",
    "NUDT15":  "Thiopurines — azathioprine, 6-MP (marrow suppression)",
    "MTHFR":   "Methotrexate, nitrous oxide sensitivity, folate metabolism",
    "COMT":    "Catecholamine drugs, methylphenidate, levodopa",
    "NAT2":    "Isoniazid, hydralazine, procainamide, caffeine (slow/fast acetylation)",
    "APOE":    "Statins (response variability), Alzheimer risk stratification",
    "HLA-B":   "Abacavir (HIV), carbamazepine, allopurinol — severe hypersensitivity",
    "G6PD":    "Rasburicase, primaquine, dapsone — hemolytic anemia risk",
}

LAB_RECOMMENDATIONS = {
    "CYP2C19": ["Helicobacter pylori breath test if on PPI long-term", "CBC if on clopidogrel"],
    "VKORC1":  ["INR monitoring (warfarin)", "Vitamin K status"],
    "SLCO1B1": ["CK (creatine kinase) if on statin", "Liver function tests"],
    "TPMT":    ["CBC before starting thiopurine therapy", "6-TGN/6-MMP metabolite levels"],
    "MTHFR":   ["Homocysteine level", "Folate/B12 panel", "RBC folate"],
    "APOE":    ["Fasting lipid panel", "ApoB", "Lp(a)"],
    "DPYD":    ["Consider pre-treatment DPYD genotyping if fluoropyrimidine therapy planned"],
}


def _build_drug_warning_table(findings: list[dict]) -> str:
    # Gather genes with pre_prescription or actionable findings
    genes = {}
    for f in findings:
        gene = f.get("gene", "")
        if not gene:
            continue
        if f["silo"] in ("pre_prescription", "actionable") and gene in DRUG_CLASS_MAP:
            genes[gene] = {
                "genotype":  f.get("genotype", ""),
                "rsid":      f.get("rsid", ""),
                "interp":    f.get("interpretation", "")[:200],
                "drug_classes": DRUG_CLASS_MAP[gene],
            }

    if not genes:
        return "<p style='color:#666;font-style:italic;'>No prescription-critical gene variants identified in this scan.</p>"

    rows = ""
    for gene, info in sorted(genes.items()):
        rows += f"""
      <tr>
        <td><strong>{gene}</strong><br><code style='font-size:0.8em;'>{info['rsid']}</code></td>
        <td>{info['genotype']}</td>
        <td style='font-size:0.85em;'>{info['drug_classes']}</td>
        <td style='font-size:0.85em;color:#c0392b;'>{info['interp']}</td>
      </tr>"""

    return f"""
    <table style="width:100%;border-collapse:collapse;font-size:0.9em;">
      <thead>
        <tr style="background:#2c3e50;color:#fff;">
          <th style="padding:8px;text-align:left;">Gene</th>
          <th style="padding:8px;text-align:left;">Genotype</th>
          <th style="padding:8px;text-align:left;">Affected Drug Classes</th>
          <th style="padding:8px;text-align:left;">Clinical Note</th>
        </tr>
      </thead>
      <tbody>
        {rows}
      </tbody>
    </table>"""


def _build_lab_section(findings: list[dict]) -> str:
    genes_found = {f.get("gene", "") for f in findings if f["silo"] in ("pre_prescription", "actionable")}
    labs = set()
    for gene in genes_found:
        for lab in LAB_RECOMMENDATIONS.get(gene, []):
            labs.add(lab)
    if not labs:
        return "<p style='color:#666;font-style:italic;'>No specific lab follow-ups indicated by current findings.</p>"
    items = "".join(f"<li>{lab}</li>" for lab in sorted(labs))
    return f"<ul style='margin:0;padding-left:20px;'>{items}</ul>"


def _pre_rx_table(findings: list[dict]) -> str:
    items = [f for f in findings if f["silo"] == "pre_prescription"]
    if not items:
        return "<p style='color:#27ae60;'>No prescription-critical variants identified in this scan. Always verify with clinical testing.</p>"

    rows = ""
    for f in items:
        gene   = f.get("gene", "N/A")
        rsid   = f.get("rsid", "")
        gt     = f.get("genotype", "")
        zyg    = (f.get("zygosity") or "").replace("_", " ")
        clsig  = f.get("clinical_sig", "").title()
        interp = f.get("interpretation", f.get("conditions", ""))[:250]
        gt_cell = gt + (f"<br><span style='font-size:0.78em;color:#888;'>{zyg}</span>" if zyg else "")
        rows += f"""
      <tr style="border-bottom:1px solid #eee;">
        <td style="padding:8px;"><strong>{gene}</strong></td>
        <td style="padding:8px;"><code>{rsid}</code></td>
        <td style="padding:8px;font-weight:bold;">{gt_cell}</td>
        <td style="padding:8px;color:#c0392b;">{clsig}</td>
        <td style="padding:8px;font-size:0.85em;">{interp}</td>
      </tr>"""

    return f"""
    <table style="width:100%;border-collapse:collapse;font-size:0.9em;">
      <thead>
        <tr style="background:#c0392b;color:#fff;">
          <th style="padding:8px;text-align:left;">Gene</th>
          <th style="padding:8px;text-align:left;">rsID</th>
          <th style="padding:8px;text-align:left;">Genotype</th>
          <th style="padding:8px;text-align:left;">Classification</th>
          <th style="padding:8px;text-align:left;">Clinical Note</th>
        </tr>
      </thead>
      <tbody>{rows}</tbody>
    </table>"""


def generate_doctor_report(profile: dict, findings: list[dict]) -> str:
    name     = profile.get("name", "Unknown")
    dob      = profile.get("dob", "N/A")
    sex      = profile.get("sex", "N/A").title()
    provider = profile.get("provider", "Unknown").title()
    gen_date = datetime.utcnow().strftime("%Y-%m-%d")

    pre_rx    = [f for f in findings if f["silo"] == "pre_prescription"]
    actionable = [f for f in findings if f["silo"] == "actionable"]

    genes_list = sorted({f.get("gene","") for f in pre_rx + actionable if f.get("gene")})
    genes_str  = ", ".join(genes_list) if genes_list else "None"

    # Build JSON summary for AI prompt block
    import json
    ai_findings = []
    for f in (pre_rx + actionable)[:30]:
        ai_findings.append({
            "rsid": f.get("rsid"), "gene": f.get("gene"),
            "genotype": f.get("genotype"), "zygosity": f.get("zygosity"),
            "silo": f.get("silo"),
            "clinical_sig": f.get("clinical_sig"),
            "interpretation": f.get("interpretation","")[:150],
        })
    ai_json = json.dumps({"patient": name, "sex": sex, "findings": ai_findings}, indent=2)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Doctor Discussion Report — {name}</title>
<style>
  body {{font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:0;background:#f4f6f9;color:#2c3e50;}}
  .container {{max-width:900px;margin:0 auto;padding:24px;}}
  .header {{background:#1a3a6b;color:#fff;padding:28px;border-radius:8px;margin-bottom:20px;}}
  .header h1 {{margin:0 0 8px;font-size:1.6em;}}
  .section {{background:#fff;border-radius:8px;padding:20px;margin-bottom:20px;box-shadow:0 1px 4px rgba(0,0,0,0.1);}}
  .section h2 {{margin:0 0 14px;font-size:1.1em;color:#1a3a6b;border-bottom:2px solid #1a3a6b;padding-bottom:6px;}}
  .alert-box {{background:#fdf3f3;border-left:4px solid #c0392b;padding:12px 16px;border-radius:0 6px 6px 0;margin-bottom:16px;}}
  .info-grid {{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-bottom:0;}}
  .info-cell {{background:#f4f6f9;border-radius:6px;padding:12px;}}
  .info-cell label {{font-size:0.75em;color:#666;text-transform:uppercase;}}
  .info-cell value {{display:block;font-weight:bold;font-size:1em;}}
  .code-block {{background:#1e1e1e;color:#d4d4d4;border-radius:6px;padding:16px;font-family:monospace;font-size:0.8em;white-space:pre-wrap;overflow-x:auto;max-height:300px;overflow-y:auto;}}
  .disclaimer {{background:#fff8e1;border:1px solid #ffc107;border-radius:6px;padding:12px 16px;font-size:0.82em;line-height:1.5;}}
  @media print {{body{{background:#fff;}} .container{{padding:0;max-width:100%;}}}}
</style>
</head>
<body>
<div class="container">

  <div class="header">
    <h1>Doctor Discussion Report</h1>
    <p style="margin:0;opacity:0.85;font-size:0.9em;">
      Genetic summary for clinical review &mdash; {gen_date}
    </p>
  </div>

  <div class="disclaimer">
    <strong>FOR HEALTHCARE PROVIDER USE:</strong> This document summarizes consumer DNA array findings for discussion with a licensed provider. Consumer arrays are NOT clinical-grade tests. Results require clinical validation before any prescribing or diagnostic decision. This document does not constitute a medical record.
  </div>
  <br>

  <div class="section">
    <h2>Patient Information</h2>
    <div class="info-grid">
      <div class="info-cell"><label>Full Name</label><value>{name}</value></div>
      <div class="info-cell"><label>Date of Birth</label><value>{dob}</value></div>
      <div class="info-cell"><label>Biological Sex</label><value>{sex}</value></div>
      <div class="info-cell"><label>DNA Source</label><value>{provider}</value></div>
      <div class="info-cell"><label>Prescription-Critical Variants</label><value style="color:#c0392b;">{len(pre_rx)}</value></div>
      <div class="info-cell"><label>Actionable Findings</label><value style="color:#e67e22;">{len(actionable)}</value></div>
    </div>
    <p style="margin-top:14px;font-size:0.9em;">
      <strong>Genes with variants of clinical interest:</strong> {genes_str}
    </p>
  </div>

  <div class="section">
    <h2>Prescription-Critical Variants</h2>
    <div class="alert-box">
      These variants may require dose adjustment, alternative drug selection, or enhanced monitoring before prescribing affected drug classes.
    </div>
    {_pre_rx_table(findings)}
  </div>

  <div class="section">
    <h2>Drug Class Interaction Summary</h2>
    <p style="font-size:0.85em;color:#666;margin-bottom:12px;">
      Drug classes affected by this patient's identified gene variants:
    </p>
    {_build_drug_warning_table(findings)}
  </div>

  <div class="section">
    <h2>Recommended Laboratory Follow-Ups</h2>
    {_build_lab_section(findings)}
  </div>

  <div class="section">
    <h2>AI-Assisted Analysis Prompt (Grok / Claude / ChatGPT)</h2>
    <p style="font-size:0.85em;color:#555;margin-bottom:10px;">
      Copy the text below and paste it into Grok, Claude, or ChatGPT for a detailed AI-assisted pharmacogenomics interpretation. Do not share this with services you do not trust.
    </p>
    <div class="code-block">You are a clinical pharmacogenomics specialist. Analyze the following patient DNA findings and provide:
1. A plain-language summary of each finding and its clinical significance.
2. Specific drug classes the patient should discuss with their prescriber.
3. Any lifestyle or supplement considerations based on the variants.
4. Recommended monitoring or follow-up labs.
5. Questions the patient should bring to their next appointment.

IMPORTANT: Separate confirmed clinical findings from consumer-array limitations. Flag any finding where chip-based detection is unreliable.

Patient findings (JSON):
{ai_json}
</div>
  </div>

  <div class="section">
    <h2>Notes for Prescriber</h2>
    <p style="font-size:0.9em;line-height:1.6;">
      The patient presents consumer genetic data from a <strong>{provider}</strong> DNA array.
      Consumer arrays genotype approximately 600,000&ndash;700,000 SNPs using a chip-based method.
      They <strong>cannot</strong> reliably detect: copy number variants (e.g., CYP2D6 duplications),
      star allele haplotypes requiring phasing, structural variants, or rare variants not on the array.
    </p>
    <p style="font-size:0.9em;line-height:1.6;">
      For clinical prescribing decisions, order a <strong>validated clinical PGx panel</strong>
      (e.g., GeneSight, Genomind, Invitae PGx, or equivalent) before making medication changes
      solely based on this report.
    </p>
  </div>

</div>
</body>
</html>"""

    return html
