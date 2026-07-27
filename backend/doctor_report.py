"""
doctor_report.py -- the Doctor Discussion Report, v2 field set.

WHAT THIS IS
------------
The document a clinician reads. It is meant to be printed and carried into an
appointment with a physician, pharmacist or genetic counsellor. One
self-contained HTML file, inline CSS, no external request of any kind.

WHY IT LOOKS LIKE THIS
----------------------
A clinician reads top down and stops when they have what they need, so the
prescription-critical table leads, sorted by CPIC assignment level and then by
magnitude. Everything that would change a prescribing decision is above
everything that would not.

interactive_report.py is the reference for how a v2 finding is presented, and
this file matches its standards: the same repute colours, the same word
"unscored" for a null magnitude, caveats inline rather than folded away, and the
same escaping rule. Every value that reaches the document passes through _esc(),
because a raw DNA file is attacker controllable in principle and an
interpretation string is never trusted as markup.

The v1.2 content is preserved: the same three silos, the same disclaimers, the
drug class table, the lab follow-ups and the AI analysis prompt block.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

from . import APP_VERSION

__all__ = ["generate_doctor_report"]


DRUG_CLASS_MAP: dict[str, str] = {
    "CYP2D6":  "Antidepressants, antipsychotics, opioids, beta-blockers, tamoxifen",
    "CYP2C19": "PPIs, clopidogrel, SSRIs, tricyclic antidepressants, voriconazole",
    "CYP2C9":  "Warfarin, NSAIDs, phenytoin, sulfonylureas",
    "CYP3A4":  "Immunosuppressants, statins, HIV medications, benzodiazepines",
    "CYP3A5":  "Tacrolimus, cyclosporine, sirolimus",
    "CYP2B6":  "Bupropion, efavirenz, methadone, ketamine",
    "CYP1A2":  "Caffeine, clozapine, olanzapine, theophylline, fluvoxamine",
    "VKORC1":  "Warfarin, acenocoumarol, phenprocoumon (vitamin K antagonists)",
    "SLCO1B1": "Statins: simvastatin, atorvastatin, rosuvastatin (myopathy risk)",
    "TPMT":    "Azathioprine, 6-mercaptopurine, thioguanine (hematologic toxicity)",
    "DPYD":    "Fluorouracil, capecitabine (severe 5-FU toxicity)",
    "NUDT15":  "Thiopurines: azathioprine, 6-MP (marrow suppression)",
    "MTHFR":   "Methotrexate, nitrous oxide sensitivity, folate metabolism",
    "COMT":    "Catecholamine drugs, methylphenidate, levodopa",
    "NAT2":    "Isoniazid, hydralazine, procainamide, caffeine (slow/fast acetylation)",
    "APOE":    "Statins (response variability), Alzheimer risk stratification",
    "HLA-B":   "Abacavir (HIV), carbamazepine, allopurinol: severe hypersensitivity",
    "G6PD":    "Rasburicase, primaquine, dapsone: hemolytic anemia risk",
}

LAB_RECOMMENDATIONS: dict[str, list[str]] = {
    "CYP2C19": ["Helicobacter pylori breath test if on PPI long-term",
                "CBC if on clopidogrel"],
    "VKORC1":  ["INR monitoring (warfarin)", "Vitamin K status"],
    "SLCO1B1": ["CK (creatine kinase) if on statin", "Liver function tests"],
    "TPMT":    ["CBC before starting thiopurine therapy",
                "6-TGN/6-MMP metabolite levels"],
    "MTHFR":   ["Homocysteine level", "Folate/B12 panel", "RBC folate"],
    "APOE":    ["Fasting lipid panel", "ApoB", "Lp(a)"],
    "DPYD":    ["Consider pre-treatment DPYD genotyping if fluoropyrimidine "
                "therapy planned"],
}

# CPIC assignment levels in clinical priority order. Eight values, not four:
# there are three split levels plus a non-letter "Retired", per API_V2 section
# 2.3. Anything unrecognised or empty sorts last.
CPIC_ORDER: dict[str, int] = {
    "A": 0, "A/B": 1, "B": 2, "B/C": 3, "C": 4, "C/D": 5, "D": 6, "Retired": 7,
}

# Repute colours are fixed by docs/API_V2.md section 5 and shared with the
# interactive report, so a clinician reading both sees one colour language.
REPUTE_GOOD = "#60B060"
REPUTE_BAD = "#FF9090"
REPUTE_UNSET = "#C0C0C0"

# Never coloured: a trait is not good or bad, and a polygenic score is a
# statistic. Colouring either would be an editorial claim about the patient.
NEUTRAL_ENTITIES = ("trait", "prs")

# Magnitude at or above which a finding is flagged for confirmatory testing.
CONFIRM_THRESHOLD = 6.0

PRS_DISCLAIMER = (
    "A polygenic score is a statistical predictor, not a diagnostic test. It is "
    "computed from a small subset of known variants, it is calibrated on a "
    "European reference population and transfers poorly to others, and it "
    "ignores every environmental and lifestyle factor. A high score is not a "
    "diagnosis and a low score is not a clearance."
)

# The three statements the report is required to make in its own voice. Kept as
# data so they cannot drift apart from the section that renders them.
LIMITATIONS = (
    ("Coverage",
     "Consumer arrays do not call star alleles completely. They genotype "
     "individual positions, not haplotypes, and they miss copy number variants, "
     "phased star allele definitions, structural variants and rare variants that "
     "are not on the chip. A normal result here does not rule out a variant, and "
     "must not be read as a negative test."),
    ("Provenance of the scores",
     "Magnitude and repute in this document are computed by DNAInsight from CC0 "
     "and public-domain evidence: CPIC assignment level, ClinVar review status, "
     "FDA label tier, population frequency and publication counts. They are not "
     "SNPedia values of the same name, they are not clinical severity scores, and "
     "they carry no regulatory standing."),
    ("Strand ambiguity",
     "A finding marked strand-ambiguous is an A/T or C/G genotype whose strand "
     "cannot be verified from the data, so the reading may be the complement. "
     "Confirm any strand-ambiguous finding on an orthogonal assay before acting "
     "on it."),
    ("Pooled disagreements",
     "Where two of the patient's own files disagree at a position, no winner is "
     "chosen. Both calls are kept and neither reading was selected, because a "
     "disagreement between two arrays is itself information about reliability."),
    ("Non-carrier findings",
     "A classification describes an allele, not a position. Where this document "
     "says the patient does not carry the reported variant, the classification "
     "attached to that variant does not apply to them."),
)


_CSS = """
:root{--good:#60B060;--bad:#FF9090;--unset:#C0C0C0;--blue:#1a3a6b;
--mid:#2980b9;--rx:#c0392b;--act:#e67e22;--bg:#f4f6f9;--text:#2c3e50;
--muted:#7f8c8d;--line:#e3e9ef}
*{box-sizing:border-box}
body{font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:0;
background:var(--bg);color:var(--text);line-height:1.5}
.container{max-width:980px;margin:0 auto;padding:24px}
.header{background:#1a3a6b;color:#fff;padding:28px;border-radius:8px;
margin-bottom:18px}
.header h1{margin:0 0 8px;font-size:1.6em}
.section{background:#fff;border-radius:8px;padding:20px;margin-bottom:18px;
box-shadow:0 1px 4px rgba(0,0,0,.1)}
.section h2{margin:0 0 14px;font-size:1.1em;color:var(--blue);
border-bottom:2px solid var(--blue);padding-bottom:6px}
.section h3{margin:14px 0 6px;font-size:.95em;color:var(--blue)}
.section p{margin:0 0 9px}
.alert-box{background:#fdf3f3;border-left:4px solid var(--rx);padding:12px 16px;
border-radius:0 6px 6px 0;margin-bottom:14px;font-size:.88em}
.info-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}
.info-cell{background:var(--bg);border-radius:6px;padding:12px}
.info-cell label{font-size:.74em;color:#666;text-transform:uppercase;
display:block}
.info-cell .value{display:block;font-weight:bold;font-size:1em}
.code-block{background:#1e1e1e;color:#d4d4d4;border-radius:6px;padding:16px;
font-family:Consolas,monospace;font-size:.8em;white-space:pre-wrap;
overflow-x:auto;max-height:340px;overflow-y:auto}
.disclaimer{background:#fff8e1;border:1px solid #ffc107;border-radius:6px;
padding:12px 16px;font-size:.84em;margin-bottom:18px}
table.grid{width:100%;border-collapse:collapse;font-size:.85em}
table.grid th{background:var(--rx);color:#fff;padding:7px 8px;text-align:left;
font-weight:600;vertical-align:bottom}
table.grid td{padding:7px 8px;border-bottom:1px solid var(--line);
vertical-align:top}
table.grid tr:nth-child(even) td{background:#fafbfc}
table.med{width:100%;border-collapse:collapse;font-size:.83em;margin-top:5px}
table.med th{background:#eef2f7;color:var(--blue);padding:5px 7px;
text-align:left}
table.med td{padding:5px 7px;border-bottom:1px solid var(--line)}
.mono{font-family:Consolas,monospace}
.stars{color:#7a5c00;letter-spacing:1px}
.silo{margin-bottom:22px}
.silo-head{color:#fff;padding:9px 15px;border-radius:6px 6px 0 0}
.silo-head h3{margin:0;color:#fff;font-size:1.02em}
.silo-head .sub{font-size:.83em;opacity:.9}
.silo-body{border:1px solid var(--line);border-top:none;
border-radius:0 0 6px 6px;background:#fff;padding:11px}
.finding{border:1px solid var(--line);border-left:5px solid var(--unset);
background:#fff;border-radius:0 6px 6px 0;padding:11px 13px;margin-bottom:9px}
.finding.good{border-left-color:var(--good)}
.finding.bad{border-left-color:var(--bad)}
.finding.dub{border-left-style:dashed}
.fh{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:4px}
.rsid{font-family:Consolas,monospace;font-weight:700}
.tok{font-family:Consolas,monospace;color:var(--mid);font-size:.9em}
.gene{color:var(--muted)}
.mag{background:var(--blue);color:#fff;border-radius:11px;padding:1px 9px;
font-size:.75em;font-weight:700}
.mag.zero{background:#b9c1c9}.mag.low{background:#8fa6bd}
.mag.mid{background:var(--mid)}.mag.high{background:#b5341f}
.tag{font-size:.71em;padding:2px 7px;border-radius:9px;background:#ecf0f1;
color:#556}
.tag.cpic{background:#e8f0fe;color:#17408b}
.tag.star{background:#fef6d8;color:#7a5c00}
.tag.flag{background:#fff8e1;color:#7d5900}
.sum{font-weight:600;font-size:.93em}
.body{font-size:.87em;color:#445}
.crit{font-family:Consolas,monospace;font-size:.79em;background:#f7f9fb;
border:1px solid var(--line);border-radius:4px;padding:6px 8px;margin-top:6px;
white-space:pre-wrap;word-break:break-word}
.nocarry{background:#eef7ee;border:2px solid var(--good);border-radius:6px;
padding:9px 12px;margin:8px 0;font-size:.88em}
.warn{background:#fff8e1;border:1px solid #ffc107;border-radius:5px;
padding:8px 11px;margin-top:8px;font-size:.84em}
.quiet{color:var(--muted);font-size:.81em;margin-top:6px}
.cav{margin-top:8px;background:#f7f9fb;border-left:3px solid var(--muted);
padding:7px 10px;font-size:.84em}
.cav ul{margin:4px 0 0 16px;padding:0}
.cfl{margin-top:8px;background:#fff8e1;border:1px solid #ffc107;
border-radius:5px;padding:8px 11px;font-size:.84em}
.calls{display:flex;flex-wrap:wrap;gap:10px;margin-top:6px}
.callcell{border:1px solid var(--line);background:#fff;border-radius:5px;
padding:5px 9px;min-width:118px}
.calllabel{font-size:.77em;color:var(--muted)}
.callgt{font-family:Consolas,monospace;font-weight:700}
table.kv{width:100%;border-collapse:collapse;font-size:.81em;margin-top:8px}
table.kv td{padding:3px 6px;vertical-align:top;border-top:1px solid var(--line)}
table.kv td:first-child{color:var(--muted);width:32%}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:.85em;margin-bottom:8px}
.dot{display:inline-block;width:11px;height:11px;border-radius:3px;
margin-right:5px}
.none{color:var(--muted);font-style:italic;font-size:.89em}
.lab-driver{color:var(--muted);font-size:.85em}
.footer{text-align:center;font-size:.79em;color:var(--muted);margin-top:26px;
padding-top:14px;border-top:1px solid var(--line);line-height:1.6}
@media(max-width:820px){.info-grid{grid-template-columns:repeat(2,1fr)}}
@media print{body{background:#fff}.container{padding:0;max-width:100%}
.section{box-shadow:none;border:1px solid var(--line)}
.finding,table.grid tr{break-inside:avoid}}
"""


# ---------------------------------------------------------------------------
# Escaping. Identical discipline to interactive_report.py.
# ---------------------------------------------------------------------------

_HTML_ESCAPES = {
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
}


def _esc(value: Any) -> str:
    """HTML-escape any value. Never trust a string that came from a data file."""
    if value is None:
        return ""
    return "".join(_HTML_ESCAPES.get(c, c) for c in str(value))


# ---------------------------------------------------------------------------
# Value helpers
# ---------------------------------------------------------------------------

def _num(value: Any) -> float | None:
    """Return ``value`` as a float, or None when it is not a usable number."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _trim(value: float) -> str:
    """Format a float without trailing zeros, so 4.0 prints as 4."""
    return f"{value:.2f}".rstrip("0").rstrip(".") or "0"


def _mag_sort(finding: dict) -> float:
    """Magnitude sort key. Null counts as 1, per API_V2 section 2.2."""
    value = _num(finding.get("magnitude"))
    return 1.0 if value is None else value


def _mag_text(finding: dict) -> str:
    """Magnitude with the scale stated, or the word unscored.

    A null is never printed as 0. Zero is an earned score, a no-call gets one;
    unscored means no evidence was available to score it at all.
    """
    value = _num(finding.get("magnitude"))
    if value is None:
        return "unscored"
    return f"{_trim(value)} out of 10"


def _mag_class(finding: dict) -> str:
    """Colour band for the magnitude pill."""
    value = _num(finding.get("magnitude"))
    if value is None:
        return "low"
    if value == 0:
        return "zero"
    if value < 2:
        return "low"
    if value < CONFIRM_THRESHOLD:
        return "mid"
    return "high"


def _repute_class(finding: dict) -> str:
    """CSS class for the coloured left border."""
    if str(finding.get("entity_type") or "snp") in NEUTRAL_ENTITIES:
        return "unset"
    repute = str(finding.get("repute") or "")
    if repute == "Good":
        return "good"
    if repute == "Bad":
        return "bad"
    return "unset"


def _repute_text(finding: dict) -> str:
    """Repute in words, including why it is unset when it is unset."""
    if str(finding.get("entity_type") or "snp") in NEUTRAL_ENTITIES:
        return "unset by design, a trait or score has no direction of effect"
    repute = str(finding.get("repute") or "")
    if repute in ("Good", "Bad"):
        return repute
    return "unset, direction of effect neutral, conflicting or unknown"


def _star_count(finding: dict) -> int:
    """ClinVar review stars as an integer, clamped to the documented 0 to 4."""
    raw = finding.get("review_stars")
    return raw if isinstance(raw, int) and 0 <= raw <= 4 else 0


def _stars_html(finding: dict) -> str:
    """Review stars as filled and hollow star characters, 0 to 4."""
    count = _star_count(finding)
    return ('<span class="stars">' + "&#9733;" * count
            + "&#9734;" * (4 - count) + "</span>")


def _confidence_text(finding: dict) -> str:
    """Confidence in words, with 'none' stated rather than left blank."""
    value = str(finding.get("confidence") or "").strip().lower()
    if value in ("high", "moderate", "low"):
        return value
    return "none, no usable evidence tier"


def _cpic_text(finding: dict) -> str:
    """CPIC assignment level, or an explicit statement that there is none."""
    level = str(finding.get("cpic_level") or "").strip()
    return level if level else "none assigned"


def _cpic_rank(finding: dict) -> int:
    """Clinical priority rank for sorting. Unassigned sorts last."""
    level = str(finding.get("cpic_level") or "").strip()
    return CPIC_ORDER.get(level, len(CPIC_ORDER))


def _frequency_html(finding: dict) -> str:
    """Frequency: value, band, derivation method and source population.

    A null and a 0.0 are different facts and are worded differently. Null means
    the panel has no data for this genotype. Zero means the panel was checked and
    the genotype was not seen in it. Collapsing the two would overstate what is
    known, which in a prescribing context is the expensive kind of error.
    """
    value = _num(finding.get("freq"))
    if value is None:
        return '<span class="none">no data</span>'

    population = _esc(finding.get("freq_population") or "")
    where = f" in {population}" if population else ""
    method = str(finding.get("freq_method") or "")
    if finding.get("freq_derived") or method == "hardy_weinberg":
        how = "derived under Hardy-Weinberg, not directly counted"
    elif method == "observed":
        how = "observed in the panel"
    else:
        how = "derivation not recorded"

    if value == 0.0:
        return (f"not observed in this panel{where} "
                f'<span class="none">({how})</span>')

    band = str(finding.get("freq_band") or "unknown").replace("_", " ")
    band_html = "" if band == "unknown" else f", {_esc(band)}"
    return (f"{_trim(value)}%{where}{band_html} "
            f'<span class="none">({how})</span>')


def _drugs_text(finding: dict) -> str:
    """Affected drugs for a finding: its own medicines, else the gene's classes."""
    medicines = [str(m).strip() for m in (finding.get("medicines") or [])
                 if str(m).strip()]
    if medicines:
        return ", ".join(sorted(set(medicines))[:12])
    return DRUG_CLASS_MAP.get(str(finding.get("gene") or "").upper(), "")


# ---------------------------------------------------------------------------
# Honesty blocks, in clinician voice
# ---------------------------------------------------------------------------

def _carrier_html(finding: dict) -> str:
    """Non-carrier banner, or nothing.

    A banner rather than a badge, and above the interpretation rather than below
    it. A ClinVar classification describes an allele; asserting it against a
    patient who holds two reference copies is the single most consequential error
    this document could make, so it is called out where it cannot be skimmed past.
    """
    if finding.get("carrier") is not False:
        return ""
    return (
        '<div class="nocarry"><strong>The patient is NOT a carrier of this '
        "variant.</strong> The genotype called here does not carry the reported "
        "variant, so the classification below does not apply to this patient. It "
        "is retained only to show that the position was examined.</div>"
    )


def _strand_html(finding: dict) -> str:
    """Strand notes: visible for ambiguity, quiet for a routine complement."""
    parts: list[str] = []
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        parts.append(
            '<div class="warn"><strong>Strand ambiguous.</strong> This genotype '
            "is A/T or C/G, so the strand cannot be verified from the data and "
            "the reading shown may be the complement. Confirm on an orthogonal "
            "assay before acting on it.</div>")
    if finding.get("flipped") or finding.get("freq_flipped"):
        # Deliberately quiet. A complement applied during annotation is routine
        # bookkeeping; presenting it as a warning devalues the real warnings.
        parts.append(
            '<div class="quiet">Alleles were complemented to match the reference '
            "strand during annotation. This is routine.</div>")
    return "".join(parts)


def _confirm_html(finding: dict) -> str:
    """Confirmatory-testing warning for a high magnitude call."""
    value = _num(finding.get("magnitude"))
    if value is None or value < CONFIRM_THRESHOLD:
        return ""
    return (
        '<div class="warn"><strong>Confirm before acting.</strong> A DNAInsight '
        f"magnitude at or above {_trim(CONFIRM_THRESHOLD)} out of 10 is a strong "
        "claim from a chip-based assay, where a rare high impact call is "
        "sometimes a false positive. Confirm with a clinically validated test "
        "before this finding informs prescribing.</div>")


def _caveats_html(finding: dict) -> str:
    """Every caveat string, inline. Nothing is collapsed behind a toggle."""
    items = [c for c in (finding.get("caveats") or []) if str(c).strip()]
    single = finding.get("caveat")
    if single and str(single).strip():
        items.append(single)
    if not items:
        return ""
    body = "".join(f"<li>{_esc(c)}</li>" for c in items)
    return f'<div class="cav"><strong>Caveats</strong><ul>{body}</ul></div>'


def _provenance_html(finding: dict) -> str:
    """Multi-file provenance and any unresolved disagreement."""
    parts: list[str] = []
    count = finding.get("count")
    if isinstance(count, int) and count > 1:
        labels = [str(x) for x in (finding.get("labels") or []) if str(x).strip()]
        named = f' ({_esc(", ".join(labels))})' if labels else ""
        parts.append(f'<div class="quiet">Called by {count} pooled source '
                     f"files{named}.</div>")
    if finding.get("conflict"):
        cells = []
        for call in finding.get("calls") or []:
            genotype = call.get("genotype") or (
                str(call.get("allele1") or "") + str(call.get("allele2") or ""))
            cells.append('<div class="callcell"><div class="calllabel">'
                         + _esc(call.get("label"))
                         + '</div><div class="callgt">' + _esc(genotype)
                         + "</div></div>")
        parts.append('<div class="cfl"><strong>Pooled sources disagree at this '
                     "position.</strong> Both calls are kept and neither was "
                     "chosen. Treat the genotype here as unresolved until it is "
                     'repeated.<div class="calls">' + "".join(cells)
                     + "</div></div>")
    return "".join(parts)


def _reliability_html(finding: dict) -> str:
    """Flag a polygenic score whose coverage is below the reliability threshold."""
    if str(finding.get("entity_type") or "") != "prs":
        return ""
    if finding.get("reliable") is not False:
        return ""
    return ('<div class="warn"><strong>Score not reliable.</strong> Array '
            "coverage of this model is below the 0.90 threshold, so the "
            "percentile below is not dependable and should not be quoted to the "
            "patient as a risk figure.</div>")


# ---------------------------------------------------------------------------
# One finding
# ---------------------------------------------------------------------------

def _badges_html(finding: dict) -> str:
    """Tags after the identifier, ordered so the caveats are read first."""
    tags: list[tuple[str, str]] = []
    if finding.get("carrier") is False:
        tags.append(("flag", "not a carrier"))
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        tags.append(("flag", "strand ambiguous"))
    if finding.get("conflict"):
        tags.append(("flag", "sources disagree"))
    if str(finding.get("zygosity") or "") == "no_call":
        tags.append(("flag", "no call"))
    if finding.get("flipped"):
        tags.append(("", "strand flipped"))
    if _star_count(finding):
        tags.append(("star", f"{_stars_html(finding)} {_star_count(finding)} of 4"))
    level = str(finding.get("cpic_level") or "").strip()
    if level:
        tags.append(("cpic", f"CPIC {_esc(level)}"))
    entity = str(finding.get("entity_type") or "snp")
    if entity != "snp":
        tags.append(("", _esc(entity)))
    zygosity = str(finding.get("zygosity") or "")
    if zygosity and zygosity != "no_call":
        tags.append(("", _esc(zygosity.replace("_", " "))))
    return "".join(f'<span class="tag {cls}">{text}</span>' for cls, text in tags)


def _detail_rows(finding: dict) -> str:
    """Per-finding detail table. Absent values are omitted, never faked."""
    entity = str(finding.get("entity_type") or "snp")
    rows: list[tuple[str, str]] = [
        ("DNAInsight magnitude", _esc(_mag_text(finding))),
        ("Repute", _esc(_repute_text(finding))),
        ("Confidence", _esc(_confidence_text(finding))),
        ("ClinVar review stars",
         f"{_stars_html(finding)} {_star_count(finding)} of 4"),
        ("CPIC level", _esc(_cpic_text(finding))),
    ]

    def add(label: str, value: str) -> None:
        if value:
            rows.append((label, value))

    add("Evidence label", _esc(finding.get("evidence")))
    add("ClinVar classification", _esc(finding.get("clinical_sig")))
    if entity == "snp":
        rows.append(("Population frequency", _frequency_html(finding)))
        gmaf = _num(finding.get("gmaf"))
        add("GMAF", "" if gmaf is None else _trim(gmaf))
        if finding.get("chromosome"):
            add("Locus", _esc(f"chr{finding.get('chromosome')}:"
                              f"{finding.get('position')}"))
        pubs = finding.get("publications")
        add("Publications", str(pubs) if isinstance(pubs, int) and pubs else "")

    carrier = finding.get("carrier")
    if carrier is True:
        add("Carrier", "yes, the reported variant is present")
    elif carrier is False:
        add("Carrier", "no, the patient does not carry the reported variant")
    copies = finding.get("variant_copies")
    add("Variant copies", f"{copies} of 2" if isinstance(copies, int) else "")
    count = finding.get("count")
    add("Pooled calls", str(count) if isinstance(count, int) and count > 1 else "")
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        add("Strand", "cannot be verified, palindromic site")
    elif finding.get("flipped"):
        add("Strand", "complemented during annotation, routine")

    coverage = _num(finding.get("coverage"))
    add("Model coverage", "" if coverage is None else f"{_trim(coverage * 100)}%")
    percentile = _num(finding.get("percentile"))
    add("Percentile", "" if percentile is None else _trim(percentile))
    add("Band", _esc(str(finding.get("band") or "").replace("_", " ")))
    if entity == "prs":
        add("Reliable", "yes" if finding.get("reliable")
            else "no, coverage below 0.90")
    matched = finding.get("matched_rsids") or []
    add("Positions matched", str(len(matched)) if matched else "")
    add("Conditions", _esc(finding.get("conditions")))
    add("Affected drugs", _esc(_drugs_text(finding)))
    topics = [str(t) for t in (finding.get("topics") or []) if str(t).strip()]
    add("Topics", _esc(", ".join(topics[:14])))

    body = "".join(f"<tr><td>{label}</td><td>{value}</td></tr>"
                   for label, value in rows)
    return f'<table class="kv"><tbody>{body}</tbody></table>'


def _finding_html(finding: dict) -> str:
    """Render one finding block. Caveats above interpretation, detail below."""
    genotype = finding.get("genotype") or (
        str(finding.get("allele1") or "") + str(finding.get("allele2") or ""))
    header = ['<div class="fh"><span class="rsid">',
              _esc(finding.get("rsid")), "</span>"]
    if finding.get("token"):
        header.append(f'<span class="tok">{_esc(finding.get("token"))}</span>')
    elif genotype:
        header.append(f'<span class="tok">{_esc(genotype)}</span>')
    if finding.get("gene"):
        header.append(f'<span class="gene">{_esc(finding.get("gene"))}</span>')
    header.append(f'<span class="mag {_mag_class(finding)}">'
                  f"{_esc(_mag_text(finding))}</span>")
    header.append(_badges_html(finding))
    header.append("</div>")

    summary = _esc(finding.get("summary"))
    interpretation = _esc(finding.get("interpretation") or finding.get("conditions"))
    text = ""
    if summary:
        text += f'<div class="sum">{summary}</div>'
    if interpretation and interpretation != summary:
        text += f'<div class="body">{interpretation}</div>'
    criteria = _esc(finding.get("criteria"))
    if criteria:
        text += f'<div class="crit">{criteria}</div>'

    classes = f"finding {_repute_class(finding)}"
    if finding.get("dubious"):
        classes += " dub"
    return (f'<div class="{classes}">' + "".join(header)
            + _carrier_html(finding) + _reliability_html(finding) + text
            + _strand_html(finding) + _confirm_html(finding)
            + _provenance_html(finding) + _caveats_html(finding)
            + _detail_rows(finding) + "</div>")


# ---------------------------------------------------------------------------
# Selection helpers
# ---------------------------------------------------------------------------

def _positional(findings: list[dict]) -> list[dict]:
    """Findings that describe a single position."""
    return [f for f in findings
            if str(f.get("entity_type") or "snp") not in ("genoset", "trait", "prs")]


def _by_magnitude(items: list[dict]) -> list[dict]:
    """Highest magnitude first, a null counting as 1."""
    return sorted(items, key=_mag_sort, reverse=True)


def _prescription_relevant(findings: list[dict]) -> list[dict]:
    """Findings that could change a prescribing decision.

    Judgement call: the silo is not the only gate. A finding carrying a CPIC
    assignment level is prescription relevant whatever silo the scanner put it
    in, because CPIC levels exist precisely to say "this pair changes
    prescribing". Taking the union avoids a genuinely actionable pair being
    demoted into the informational section and read past.
    """
    chosen: dict[str, dict] = {}
    for finding in _positional(findings):
        silo = str(finding.get("silo") or "")
        level = str(finding.get("cpic_level") or "").strip()
        if silo != "pre_prescription" and not level:
            continue
        key = str(finding.get("rsid") or "") + "|" + str(finding.get("genotype") or "")
        chosen.setdefault(key, finding)
    return sorted(chosen.values(), key=lambda f: (_cpic_rank(f), -_mag_sort(f)))


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def _prescription_table(findings: list[dict]) -> str:
    """The lead table, sorted by CPIC level then magnitude descending."""
    items = _prescription_relevant(findings)
    if not items:
        return ('<p class="none">No prescription-critical variant was identified '
                "in this scan. That is not a negative test: consumer arrays do "
                "not call star alleles completely, so always verify with "
                "clinical testing before relying on it.</p>")

    rows = []
    for finding in items:
        genotype = _esc(finding.get("genotype") or "")
        if finding.get("carrier") is False:
            genotype += ('<br><strong style="color:#1c7c2c;">NOT A CARRIER'
                         "</strong>")
        implication = _esc(str(finding.get("summary")
                               or finding.get("interpretation")
                               or finding.get("conditions") or "")[:260])
        if finding.get("carrier") is False:
            implication = ("<em>The patient does not carry the reported variant, "
                           "so this classification does not apply to them.</em> "
                           + implication)
        if finding.get("ambiguous") or finding.get("freq_ambiguous"):
            implication += ('<br><strong>Strand ambiguous, the reading cannot be '
                            "verified.</strong>")
        rows.append(
            "<tr>"
            f"<td><strong>{_esc(finding.get('gene') or 'not mapped')}</strong></td>"
            f'<td class="mono">{_esc(finding.get("rsid"))}'
            f'<br><span class="none">{_esc(_mag_text(finding))}</span></td>'
            f'<td class="mono">{genotype}</td>'
            f"<td>{_esc(str(finding.get('zygosity') or '').replace('_', ' '))}</td>"
            f"<td>{_esc(_cpic_text(finding))}</td>"
            f"<td>{_stars_html(finding)} {_star_count(finding)}</td>"
            f"<td>{implication}</td>"
            f"<td>{_esc(_drugs_text(finding))}</td>"
            "</tr>")

    return ('<table class="grid"><thead><tr>'
            "<th>Gene</th><th>Variant</th><th>Genotype</th><th>Zygosity</th>"
            "<th>CPIC level</th><th>Review stars</th>"
            "<th>Phenotype implication</th><th>Affected drugs</th>"
            "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _medicine_summary(findings: list[dict]) -> str:
    """Drug interaction summary grouped by medicine, from the medicines field."""
    by_medicine: dict[str, list[dict]] = {}
    for finding in _positional(findings):
        for medicine in finding.get("medicines") or []:
            name = str(medicine).strip()
            if name:
                by_medicine.setdefault(name, []).append(finding)
    if not by_medicine:
        return ('<p class="none">No finding in this scan carries a drug '
                "annotation, so there is nothing to group by medicine. The gene "
                "level table below still applies.</p>")

    blocks = []
    for name in sorted(by_medicine, key=str.lower):
        group = sorted(by_medicine[name], key=lambda f: (_cpic_rank(f), -_mag_sort(f)))
        lines = []
        for finding in group:
            note = ""
            if finding.get("carrier") is False:
                note = ' <strong>the patient does not carry this variant</strong>'
            lines.append(
                f"<tr><td><strong>{_esc(finding.get('gene') or '')}</strong></td>"
                f'<td class="mono">{_esc(finding.get("rsid"))}</td>'
                f'<td class="mono">{_esc(finding.get("genotype"))}</td>'
                f"<td>{_esc(_cpic_text(finding))}</td>"
                f"<td>{_stars_html(finding)}</td>"
                f"<td>{_esc(_mag_text(finding))}</td>"
                f"<td>{_esc(_repute_text(finding))}{note}</td></tr>")
        top = group[0]
        blocks.append(
            f"<h3>{_esc(name)}</h3>"
            f'<table class="med"><thead><tr><th>Gene</th><th>Variant</th>'
            "<th>Genotype</th><th>CPIC</th><th>Stars</th><th>Magnitude</th>"
            "<th>Direction</th></tr></thead><tbody>"
            + "".join(lines) + "</tbody></table>"
            f'<p class="quiet">Highest priority for this medicine: '
            f"{_esc(finding_label(top))}.</p>")
    return "".join(blocks)


def finding_label(finding: dict) -> str:
    """Short human label for a finding, used in cross references."""
    gene = str(finding.get("gene") or "").strip()
    rsid = str(finding.get("rsid") or "").strip()
    level = str(finding.get("cpic_level") or "").strip()
    parts = [p for p in (gene, rsid) if p]
    label = " ".join(parts) if parts else "unlabelled finding"
    return f"{label} (CPIC {level})" if level else label


def _drug_class_table(findings: list[dict]) -> str:
    """The v1.2 gene to drug class table, kept because it is still useful."""
    genes: dict[str, dict] = {}
    for finding in _positional(findings):
        gene = str(finding.get("gene") or "").upper()
        if not gene or gene not in DRUG_CLASS_MAP:
            continue
        if str(finding.get("silo") or "") not in ("pre_prescription", "actionable"):
            continue
        current = genes.get(gene)
        if current is None or _mag_sort(finding) > _mag_sort(current):
            genes[gene] = finding

    if not genes:
        return ('<p class="none">No gene in this scan maps to a known affected '
                "drug class.</p>")
    rows = []
    for gene in sorted(genes):
        finding = genes[gene]
        note = _esc(str(finding.get("interpretation")
                        or finding.get("summary") or "")[:200])
        if finding.get("carrier") is False:
            note = ("<em>Not a carrier, the classification does not apply.</em> "
                    + note)
        rows.append(
            f"<tr><td><strong>{_esc(gene)}</strong>"
            f'<br><span class="mono none">{_esc(finding.get("rsid"))}</span></td>'
            f'<td class="mono">{_esc(finding.get("genotype"))}</td>'
            f"<td>{_esc(DRUG_CLASS_MAP[gene])}</td>"
            f"<td>{note}</td></tr>")
    return ('<table class="grid"><thead><tr><th>Gene</th><th>Genotype</th>'
            "<th>Affected drug classes</th><th>Clinical note</th></tr></thead>"
            "<tbody>" + "".join(rows) + "</tbody></table>")


def _lab_section(findings: list[dict]) -> str:
    """Laboratory follow-ups derived from the findings actually present.

    Every line here is produced by a specific finding and names it, rather than
    being a static list keyed off a gene panel. A recommendation a clinician
    cannot trace back to a result is a recommendation they cannot weigh.
    """
    recommendations: dict[str, set[str]] = {}

    def add(text: str, driver: str) -> None:
        recommendations.setdefault(text, set()).add(driver)

    for finding in _positional(findings):
        silo = str(finding.get("silo") or "")
        gene = str(finding.get("gene") or "").upper()
        rsid = str(finding.get("rsid") or "")
        label = finding_label(finding)
        magnitude = _num(finding.get("magnitude"))
        level = str(finding.get("cpic_level") or "").strip()
        carrier = finding.get("carrier")

        # A non-carrier drives nothing. Ordering a lab because of a variant the
        # patient does not have is the failure mode this whole report guards.
        if carrier is False:
            continue

        if silo in ("pre_prescription", "actionable"):
            for lab in LAB_RECOMMENDATIONS.get(gene, []):
                add(lab, label)
        if level in ("A", "A/B", "B"):
            add(f"Validated pharmacogenomic panel covering {gene or rsid} before "
                "prescribing in this pathway", label)
        if magnitude is not None and magnitude >= CONFIRM_THRESHOLD:
            add(f"Orthogonal confirmation of {rsid} on a clinically validated "
                "assay before this result informs care", label)
        if finding.get("ambiguous") or finding.get("freq_ambiguous"):
            add(f"Strand-resolving genotyping at {rsid}, because the array call "
                "is palindromic and cannot be verified", label)
        if str(finding.get("zygosity") or "") == "no_call":
            add(f"Repeat genotyping at {rsid}, the array returned no call at "
                "this position", label)
        if finding.get("conflict"):
            add(f"Repeat genotyping at {rsid}, the patient's own source files "
                "disagree there", label)

    if not recommendations:
        return ('<p class="none">No laboratory follow-up is indicated by the '
                "findings in this scan.</p>")
    items = []
    for text in sorted(recommendations):
        drivers = ", ".join(sorted(recommendations[text]))
        items.append(f"<li>{_esc(text)}<br>"
                     f'<span class="lab-driver">Driven by: {_esc(drivers)}'
                     "</span></li>")
    return ('<ul style="margin:0;padding-left:20px;font-size:.9em;">'
            + "".join(items) + "</ul>")


def _limitations_html() -> str:
    """Confidence and limitations. The statements this document must make."""
    blocks = "".join(f"<h3>{_esc(title)}</h3><p>{_esc(text)}</p>"
                     for title, text in LIMITATIONS)
    return ('<div class="legend">'
            f'<span><span class="dot" style="background:{REPUTE_GOOD};"></span>'
            "Green border: favourable direction</span>"
            f'<span><span class="dot" style="background:{REPUTE_BAD};"></span>'
            "Red border: unfavourable direction</span>"
            f'<span><span class="dot" style="background:{REPUTE_UNSET};"></span>'
            "Grey border: unset, neutral or not applicable</span></div>"
            "<p>Traits and polygenic scores are always grey, whatever their "
            "content, because neither has a direction of effect to report. A "
            "magnitude of 0 is an earned score, typically a no-call; a finding "
            "with no scoreable evidence is shown as unscored instead.</p>"
            + blocks)


def _silo_section(findings: list[dict], silo: str, label: str, colour: str,
                  description: str) -> str:
    """One silo of positional findings, magnitude descending."""
    items = [f for f in _positional(findings) if str(f.get("silo") or "") == silo]
    if not items:
        return ""
    cards = "".join(_finding_html(f) for f in _by_magnitude(items))
    return (f'<div class="silo"><div class="silo-head" style="background:{colour};">'
            f"<h3>{_esc(label)} ({len(items)})</h3>"
            f'<div class="sub">{_esc(description)}</div></div>'
            f'<div class="silo-body">{cards}</div></div>')


def _all_silos(findings: list[dict]) -> str:
    """The three silos in fixed clinical order."""
    order = (
        ("pre_prescription", "Prescription-Critical", "#c0392b",
         "Requires prescriber review before any medication in this pathway"),
        ("actionable", "Actionable Health Finding", "#e67e22",
         "Lifestyle, supplement, or monitoring action recommended"),
        ("informational", "Informational", "#2980b9",
         "Background information; no immediate action required"),
    )
    return "".join(_silo_section(findings, silo, label, colour, description)
                   for silo, label, colour, description in order)


def _non_snp_section(findings: list[dict]) -> str:
    """Genosets, traits and polygenic scores, grouped and labelled as such.

    Held apart from the silos because none of them is a single position: a
    genoset is a rule, a trait is a phenotype call and a polygenic score is a
    statistic. Listing them beside single variants would imply they can be read
    the same way, which they cannot.
    """
    groups = (
        ("genoset", "Genosets (rule-based findings)",
         "Each entry is a rule over several positions. The rule text is shown "
         "verbatim so the logic can be checked.", ""),
        ("trait", "Traits",
         "Phenotype calls. No direction of effect is assigned, so these are "
         "never coloured.", ""),
        ("prs", "Polygenic scores",
         "Read the limitation first, then the numbers.", PRS_DISCLAIMER),
    )
    blocks: list[str] = []
    for entity, title, intro, disclaimer in groups:
        items = [f for f in findings if str(f.get("entity_type") or "") == entity]
        if not items:
            continue
        # Disclaimer above the results, deliberately. A caveat placed under a
        # percentile is read after the number has already been believed.
        warning = f'<div class="warn">{_esc(disclaimer)}</div>' if disclaimer else ""
        cards = "".join(_finding_html(f) for f in _by_magnitude(items))
        blocks.append(f"<h3>{_esc(title)} ({len(items)})</h3><p>{_esc(intro)}</p>"
                      f"{warning}{cards}")
    if not blocks:
        return ""
    return ('<div class="section"><h2>Genosets, Traits and Polygenic Scores</h2>'
            + "".join(blocks) + "</div>")


def _qc_html(qc: dict) -> str:
    """Data-quality summary built from the scan's own QC counters."""
    if not qc:
        return ""
    bits: list[str] = []
    pairs = (
        ("flipped", "calls had alleles complemented to match the reference, "
                    "which is routine"),
        ("ambiguous", "calls are palindromic (A/T or C/G) and their strand "
                      "cannot be verified"),
        ("no_call", "positions returned no call and score 0"),
        ("conflicts", "positions where two pooled source files disagree, both "
                      "readings retained"),
    )
    for key, phrase in pairs:
        value = qc.get(key)
        if isinstance(value, int) and value > 0:
            bits.append(f"<li><strong>{value}</strong> {_esc(phrase)}</li>")
    if not bits:
        return ""
    return ('<div class="section"><h2>Data Quality</h2>'
            '<ul style="margin:0 0 10px 20px;padding:0;font-size:.88em;">'
            + "".join(bits) + "</ul>"
            "<p>A palindromic site is an A/T or C/G genotype: complementing both "
            "alleles returns the same pair, so the reporting strand cannot be "
            "recovered from the file. Those calls are capped and flagged rather "
            "than trusted.</p></div>")


def _ai_prompt(profile: dict, findings: list[dict]) -> str:
    """Build the pasteable AI analysis prompt.

    Extended for v2: magnitude, repute, CPIC level, review stars, confidence and
    carrier status travel with each finding, because those are the fields that
    let a model rank and qualify what it is given. Without carrier status in
    particular, a model will confidently discuss a variant the patient does not
    have.
    """
    selected = _prescription_relevant(findings)
    seen = {id(f) for f in selected}
    extra = [f for f in _positional(findings)
             if str(f.get("silo") or "") in ("pre_prescription", "actionable")
             and id(f) not in seen]
    payload_findings = []
    for finding in (selected + _by_magnitude(extra))[:30]:
        payload_findings.append({
            "rsid": finding.get("rsid"),
            "gene": finding.get("gene"),
            "genotype": finding.get("genotype"),
            "zygosity": finding.get("zygosity"),
            "silo": finding.get("silo"),
            "magnitude": finding.get("magnitude"),
            "magnitude_scale": "0 to 10, null means unscored",
            "repute": finding.get("repute") or "",
            "cpic_level": finding.get("cpic_level") or "",
            "review_stars": _star_count(finding),
            "confidence": finding.get("confidence") or "none",
            "carrier": finding.get("carrier"),
            "clinical_sig": finding.get("clinical_sig"),
            "strand_ambiguous": bool(finding.get("ambiguous")
                                     or finding.get("freq_ambiguous")),
            "population_frequency_percent": finding.get("freq"),
            "interpretation": str(finding.get("interpretation") or "")[:150],
        })
    return json.dumps({
        "patient": profile.get("name", ""),
        "sex": str(profile.get("sex") or ""),
        "source": "consumer DNA array, not clinical grade",
        "score_provenance": ("magnitude and repute are computed by DNAInsight "
                             "from CC0 and public-domain evidence, they are not "
                             "SNPedia values"),
        "findings": payload_findings,
    }, indent=2)


# ---------------------------------------------------------------------------
# The document
# ---------------------------------------------------------------------------

def generate_doctor_report(profile: dict, findings: list[dict],
                           extras: dict | None = None) -> str:
    """Generate the Doctor Discussion Report as one self-contained HTML string.

    Args:
        profile:  dict from database.get_profile()
        findings: finding dicts, either v1.2 or v2 shaped. A missing v2 key is
                  omitted rather than guessed at.
        extras:   optional scan payload. Only ``qc`` and ``population`` are read,
                  and both are optional, so every existing two-argument caller
                  keeps working unchanged.

    Returns: HTML string with no external reference of any kind.
    """
    extras = extras or {}
    rows = [f for f in (findings or []) if isinstance(f, dict)]

    name = _esc(profile.get("name") or "Unknown")
    dob = _esc(profile.get("dob") or "N/A")
    sex = _esc(str(profile.get("sex") or "N/A").title())
    provider = _esc(str(profile.get("provider") or "Unknown").title())
    generated = _esc(datetime.now(timezone.utc).strftime("%Y-%m-%d"))
    population = _esc(extras.get("population") or "") or "not recorded"

    positional = _positional(rows)
    pre_rx = [f for f in positional
              if str(f.get("silo") or "") == "pre_prescription"]
    actionable = [f for f in positional if str(f.get("silo") or "") == "actionable"]
    genes = sorted({str(f.get("gene") or "")
                    for f in pre_rx + actionable if f.get("gene")})
    genes_str = _esc(", ".join(genes)) if genes else "None"

    high = [f for f in positional
            if (_num(f.get("magnitude")) or 0.0) >= CONFIRM_THRESHOLD]
    non_carrier = [f for f in positional if f.get("carrier") is False]
    ambiguous = [f for f in positional
                 if f.get("ambiguous") or f.get("freq_ambiguous")]
    ai_json = _esc(_ai_prompt(profile, rows))

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Doctor Discussion Report, {name}</title>
<style>{_CSS}</style>
</head>
<body>
<div class="container">

  <div class="header">
    <h1>Doctor Discussion Report</h1>
    <p style="margin:0;opacity:0.85;font-size:0.9em;">
      Genetic summary for clinical review &nbsp;|&nbsp; {generated}
      &nbsp;|&nbsp; DNAInsight v{_esc(APP_VERSION)}
    </p>
  </div>

  <div class="disclaimer">
    <strong>FOR HEALTHCARE PROVIDER USE:</strong> This document summarizes
    consumer DNA array findings for discussion with a licensed provider.
    Consumer arrays are NOT clinical-grade tests. DNAInsight is not a medical
    device. Results require clinical validation before any prescribing or
    diagnostic decision. This document does not constitute a medical record.
  </div>

  <div class="section">
    <h2>Patient Information</h2>
    <div class="info-grid">
      <div class="info-cell"><label>Full Name</label><span class="value">{name}</span></div>
      <div class="info-cell"><label>Date of Birth</label><span class="value">{dob}</span></div>
      <div class="info-cell"><label>Biological Sex</label><span class="value">{sex}</span></div>
      <div class="info-cell"><label>DNA Source</label><span class="value">{provider}</span></div>
      <div class="info-cell"><label>Frequency Population</label><span class="value">{population}</span></div>
      <div class="info-cell"><label>Total Findings</label><span class="value">{len(rows)}</span></div>
      <div class="info-cell"><label>Prescription-Critical</label><span class="value" style="color:#c0392b;">{len(pre_rx)}</span></div>
      <div class="info-cell"><label>Actionable</label><span class="value" style="color:#e67e22;">{len(actionable)}</span></div>
      <div class="info-cell"><label>Magnitude {_trim(CONFIRM_THRESHOLD)} or Above</label><span class="value">{len(high)}</span></div>
      <div class="info-cell"><label>Non-Carrier Entries</label><span class="value">{len(non_carrier)}</span></div>
      <div class="info-cell"><label>Strand-Ambiguous</label><span class="value">{len(ambiguous)}</span></div>
      <div class="info-cell"><label>Genes of Interest</label><span class="value">{len(genes)}</span></div>
    </div>
    <p style="margin-top:14px;font-size:0.9em;">
      <strong>Genes with variants of clinical interest:</strong> {genes_str}
    </p>
    <p style="font-size:0.88em;">
      Magnitude below is a DNAInsight 0 to 10 interest score, not a severity
      grade. A finding with no scoreable evidence reads as unscored, never as 0.
      {len(non_carrier)} entries are reference genotypes where the patient does
      not carry the reported variant, and those classifications do not apply to
      them.
    </p>
  </div>

  <div class="section">
    <h2>Prescription-Critical Variants</h2>
    <div class="alert-box">
      Sorted by CPIC assignment level, then by magnitude. These may require dose
      adjustment, alternative drug selection, or enhanced monitoring before
      prescribing the affected drug classes. Rows marked NOT A CARRIER are shown
      for completeness and carry no prescribing implication.
    </div>
    {_prescription_table(rows)}
  </div>

  <div class="section">
    <h2>Drug Interaction Summary by Medicine</h2>
    <p style="font-size:0.86em;color:#666;">
      Grouped from the drug annotations carried on each finding. Highest CPIC
      level first within each medicine.
    </p>
    {_medicine_summary(rows)}
  </div>

  <div class="section">
    <h2>Affected Drug Classes by Gene</h2>
    <p style="font-size:0.86em;color:#666;">
      Gene level view, retained from the previous report format for the pathways
      that carry no per-variant drug annotation.
    </p>
    {_drug_class_table(rows)}
  </div>

  <div class="section">
    <h2>Recommended Laboratory Follow-Ups</h2>
    <p style="font-size:0.86em;color:#666;">
      Each line is generated by a finding in this scan and names it. Nothing here
      is a standing panel.
    </p>
    {_lab_section(rows)}
  </div>

  <div class="section">
    <h2>Confidence and Limitations</h2>
    {_limitations_html()}
  </div>

  {_qc_html(extras.get("qc") or {})}

  <div class="section">
    <h2>Full Finding Detail</h2>
    <p style="font-size:0.86em;color:#666;">
      All positional findings, in the three silos, magnitude descending within
      each. The coloured left border encodes repute.
    </p>
    {_all_silos(rows)}
  </div>

  {_non_snp_section(rows)}

  <div class="section">
    <h2>AI-Assisted Analysis Prompt (Grok / Claude / ChatGPT)</h2>
    <p style="font-size:0.85em;color:#555;margin-bottom:10px;">
      Copy the text below and paste it into Grok, Claude, or ChatGPT for a
      detailed AI-assisted pharmacogenomics interpretation. Do not share this
      with services you do not trust. The magnitude, repute, CPIC level, review
      stars and carrier status are included because they are what let a model
      rank and qualify what it is given.
    </p>
    <div class="code-block">You are a clinical pharmacogenomics specialist. Analyze the following patient DNA findings and provide:
1. A plain-language summary of each finding and its clinical significance.
2. Specific drug classes the patient should discuss with their prescriber.
3. Any lifestyle or supplement considerations based on the variants.
4. Recommended monitoring or follow-up labs.
5. Questions the patient should bring to their next appointment.

Rank your answer by cpic_level first and then by magnitude. Field notes:
- magnitude is a DNAInsight 0 to 10 interest score, null means unscored, and it is not a severity grade.
- repute is the direction of effect: Good, Bad, or empty for neutral, mixed or not applicable.
- review_stars is the ClinVar review level, 0 to 4.
- carrier false means the patient does NOT carry the reported variant, so the classification does not apply to them. Say so plainly and do not discuss risk for those entries.
- strand_ambiguous true means the genotype is A/T or C/G and the strand cannot be verified, so the call may be the complement.

IMPORTANT: Separate confirmed clinical findings from consumer-array limitations. Flag any finding where chip-based detection is unreliable. Consumer arrays do not call star alleles completely, so a normal result does not rule out a variant.

Patient findings (JSON):
{ai_json}
</div>
  </div>

  <div class="section">
    <h2>Notes for Prescriber</h2>
    <p style="font-size:0.9em;">
      The patient presents consumer genetic data from a <strong>{provider}</strong>
      DNA array. Consumer arrays genotype approximately 600,000 to 700,000 SNPs
      using a chip-based method. They <strong>cannot</strong> reliably detect:
      copy number variants (for example CYP2D6 duplications), star allele
      haplotypes requiring phasing, structural variants, or rare variants not on
      the array.
    </p>
    <p style="font-size:0.9em;">
      For clinical prescribing decisions, order a
      <strong>validated clinical PGx panel</strong> (for example GeneSight,
      Genomind, Invitae PGx, or equivalent) before making medication changes
      solely based on this report.
    </p>
  </div>

  <div class="footer">
    DNAInsight v{_esc(APP_VERSION)} &nbsp;|&nbsp; Generated {generated} from
    {len(rows)} findings &nbsp;|&nbsp; Data sourced from ClinVar, PharmGKB, CPIC
    and MyVariant.info<br>
    DNAInsight is not a medical device and this document is not medical advice.
    Not for clinical use without independent validation.
  </div>

</div>
</body>
</html>"""
