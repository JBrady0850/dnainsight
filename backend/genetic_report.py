"""
genetic_report.py -- Generate the Genetic Health Report.

Produces an HTML report covering:
  - Pharmacogenomics (drug metabolism)
  - Metabolic health markers
  - Inflammation / immune markers
  - Neurological / cognitive markers
  - Detox / toxin clearance
  - Cardiovascular risk

Output: HTML string (saved to disk by caller)
"""

from datetime import datetime
from pathlib import Path

from . import APP_VERSION


CATEGORY_LABELS = {
    "PHARM":  "Pharmacogenomics & Drug Metabolism",
    "METAB":  "Metabolic Health & Diabetes Risk",
    "INFLAM": "Inflammation & Immune Response",
    "NEURO":  "Neurological & Cognitive Markers",
    "DETOX":  "Detox & Toxin Clearance",
    "CARDIO": "Cardiovascular Risk",
    "OTHER":  "Other Notable Variants",
}

SILO_LABELS = {
    "pre_prescription": ("Prescription-Critical", "#c0392b", "Requires prescriber review before any medication in this pathway"),
    "actionable":       ("Actionable Health Finding", "#e67e22", "Lifestyle, supplement, or monitoring action recommended"),
    "informational":    ("Informational", "#2980b9", "Background information; no immediate action required"),
}


def _badge(silo: str) -> str:
    label, color, _ = SILO_LABELS.get(silo, ("Unknown", "#7f8c8d", ""))
    return (
        f'<span style="background:{color};color:#fff;padding:2px 8px;'
        f'border-radius:3px;font-size:0.75em;font-weight:bold;">{label}</span>'
    )


def _silo_section(findings: list[dict], silo: str) -> str:
    items = [f for f in findings if f["silo"] == silo]
    if not items:
        return ""

    label, color, desc = SILO_LABELS.get(silo, ("Unknown", "#7f8c8d", ""))
    html = f"""
    <div style="margin-bottom:32px;">
      <div style="background:{color};color:#fff;padding:10px 16px;border-radius:6px 6px 0 0;">
        <h2 style="margin:0;font-size:1.1em;">{label}</h2>
        <div style="font-size:0.85em;opacity:0.9;">{desc}</div>
      </div>
      <div style="border:1px solid {color};border-top:none;border-radius:0 0 6px 6px;overflow:hidden;">
    """

    # Group by gene
    by_gene: dict[str, list[dict]] = {}
    for f in items:
        gene = f.get("gene") or "Unknown"
        by_gene.setdefault(gene, []).append(f)

    for gene, gfindings in sorted(by_gene.items()):
        html += f"""
        <div style="background:#f9f9f9;padding:8px 16px;border-bottom:1px solid #e0e0e0;">
          <strong style="font-size:0.95em;color:#2c3e50;">Gene: {gene}</strong>
        </div>
        """
        for f in gfindings:
            rsid    = f.get("rsid", "")
            gt      = f.get("genotype", f.get("allele1","") + f.get("allele2",""))
            zyg     = (f.get("zygosity") or "").replace("_", " ")
            cond    = f.get("conditions", "")
            interp  = f.get("interpretation", cond)[:400]
            clsig   = f.get("clinical_sig", "").title()
            sources = ", ".join(f.get("sources") or []).upper()
            html += f"""
        <div style="padding:12px 16px;border-bottom:1px solid #f0f0f0;">
          <div style="display:flex;justify-content:space-between;align-items:flex-start;flex-wrap:wrap;gap:4px;">
            <div>
              <code style="background:#eee;padding:1px 5px;border-radius:3px;font-size:0.85em;">{rsid}</code>
              &nbsp; Genotype: <strong>{gt}</strong>
              {f'&nbsp; <span style="color:#888;font-size:0.78em;">{zyg}</span>' if zyg else ''}
              &nbsp; <span style="color:#666;font-size:0.85em;">({clsig})</span>
            </div>
            <div style="font-size:0.8em;color:#888;">{sources}</div>
          </div>
          <p style="margin:6px 0 0 0;font-size:0.9em;color:#333;">{interp}</p>
        </div>
            """

    html += "</div></div>"
    return html


def generate_genetic_report(profile: dict, findings: list[dict]) -> str:
    """
    Generate full genetic HTML report.

    Args:
        profile:  dict from database.get_profile()
        findings: list of finding dicts from database.get_findings()

    Returns: HTML string
    """
    name     = profile.get("name", "Unknown")
    dob      = profile.get("dob", "N/A")
    sex      = profile.get("sex", "N/A").title()
    provider = profile.get("provider", "Unknown").title()
    gen_date = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    total     = len(findings)
    pre_rx    = len([f for f in findings if f["silo"] == "pre_prescription"])
    actionable = len([f for f in findings if f["silo"] == "actionable"])
    info      = len([f for f in findings if f["silo"] == "informational"])

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genetic Health Report — {name}</title>
<style>
  body {{font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:0;background:#f4f6f9;color:#2c3e50;}}
  .container {{max-width:900px;margin:0 auto;padding:24px;}}
  .header {{background:linear-gradient(135deg,#1a3a6b,#2980b9);color:#fff;padding:32px;border-radius:10px;margin-bottom:24px;}}
  .header h1 {{margin:0 0 8px 0;font-size:1.8em;}}
  .header .meta {{font-size:0.9em;opacity:0.85;line-height:1.8;}}
  .summary-grid {{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:24px;}}
  .summary-card {{background:#fff;border-radius:8px;padding:16px;text-align:center;box-shadow:0 1px 4px rgba(0,0,0,0.1);}}
  .summary-card .num {{font-size:2em;font-weight:bold;}}
  .disclaimer {{background:#fff3cd;border:1px solid #ffc107;border-radius:6px;padding:14px 16px;margin-bottom:24px;font-size:0.85em;line-height:1.5;}}
  .footer {{text-align:center;font-size:0.8em;color:#999;margin-top:32px;padding-top:16px;border-top:1px solid #e0e0e0;}}
  @media print {{body{{background:#fff;}} .container{{padding:0;}}}}
</style>
</head>
<body>
<div class="container">

  <div class="header">
    <h1>Genetic Health Report</h1>
    <div class="meta">
      <strong>Patient:</strong> {name} &nbsp;|&nbsp;
      <strong>DOB:</strong> {dob} &nbsp;|&nbsp;
      <strong>Sex:</strong> {sex}<br>
      <strong>DNA Source:</strong> {provider} Array &nbsp;|&nbsp;
      <strong>Generated:</strong> {gen_date}
    </div>
  </div>

  <div class="disclaimer">
    <strong>IMPORTANT DISCLAIMER:</strong> This report is for personal informational use only. It is not a substitute for professional medical advice, diagnosis, or treatment. Consumer DNA arrays do not perform clinical-grade sequencing. Always consult a licensed healthcare provider or clinical geneticist before making any medical decisions based on genetic data.
  </div>

  <div class="summary-grid">
    <div class="summary-card">
      <div class="num" style="color:#2c3e50;">{total}</div>
      <div>Total Findings</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#c0392b;">{pre_rx}</div>
      <div>Prescription-Critical</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#e67e22;">{actionable}</div>
      <div>Actionable</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#2980b9;">{info}</div>
      <div>Informational</div>
    </div>
  </div>

  {_silo_section(findings, "pre_prescription")}
  {_silo_section(findings, "actionable")}
  {_silo_section(findings, "informational")}

  <div class="footer">
    DNAInsight v{APP_VERSION} &nbsp;|&nbsp; Open-source personal DNA analysis tool &nbsp;|&nbsp;
    Data sourced from ClinVar, PharmGKB, and MyVariant.info<br>
    This report does not constitute medical advice. Not for clinical use.
  </div>

</div>
</body>
</html>"""

    return html
