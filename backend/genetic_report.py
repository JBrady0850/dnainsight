"""
genetic_report.py -- the printable Genetic Health Report, v2 field set.

WHAT THIS IS
------------
The static report the person tested reads about their own data. One
self-contained HTML file, inline CSS, printable, and it makes no network request
of any kind. Plain language throughout: this is the artifact that gets read
without anyone standing there to explain it.

WHY IT LOOKS LIKE THIS
----------------------
interactive_report.py is the reference for how a v2 finding is presented, and
this file matches its standards on purpose:

  * the same repute colours, so green means the same thing in both artifacts;
  * the same word, "unscored", for a null magnitude, never the number 0;
  * the same refusal to hide a caveat behind a toggle;
  * the same escaping rule. Every value that reaches the document goes through
    _esc(). A raw DNA file is attacker controllable in principle, so an
    interpretation string is never trusted as markup.

The v1.2 shape is preserved deliberately: the same three silos, the same
disclaimers, the same single printable column. Everything v2 added is layered
into each finding block rather than parked in an appendix, because a caveat that
lives somewhere else is a caveat nobody reads.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from . import APP_VERSION

__all__ = ["generate_genetic_report"]


CATEGORY_LABELS: dict[str, str] = {
    "PHARM":  "Pharmacogenomics & Drug Metabolism",
    "METAB":  "Metabolic Health & Diabetes Risk",
    "INFLAM": "Inflammation & Immune Response",
    "NEURO":  "Neurological & Cognitive Markers",
    "DETOX":  "Detox & Toxin Clearance",
    "CARDIO": "Cardiovascular Risk",
    "OTHER":  "Other Notable Variants",
}

SILO_LABELS: dict[str, tuple[str, str, str]] = {
    "pre_prescription": ("Prescription-Critical", "#c0392b",
                         "Requires prescriber review before any medication in this pathway"),
    "actionable":       ("Actionable Health Finding", "#e67e22",
                         "Lifestyle, supplement, or monitoring action recommended"),
    "informational":    ("Informational", "#2980b9",
                         "Background information; no immediate action required"),
}

# Repute colours are fixed by docs/API_V2.md section 5 and are shared with the
# interactive report so the two artifacts cannot disagree about what green means.
REPUTE_GOOD = "#60B060"
REPUTE_BAD = "#FF9090"
REPUTE_UNSET = "#C0C0C0"

# Entity types that never get a repute colour. A trait is not good or bad and a
# polygenic score is a statistic, so colouring either would be editorialising
# about a person rather than reporting their data.
NEUTRAL_ENTITIES = ("trait", "prs")

# Magnitude at or above which a finding earns a false-positive warning. A strong
# claim from a consumer array is exactly the case worth double checking.
CONFIRM_THRESHOLD = 6.0

PRS_DISCLAIMER = (
    "A polygenic score is a statistical predictor, not a diagnostic test. It is "
    "computed from a small subset of known variants, it is calibrated on a "
    "European reference population and transfers poorly to others, and it "
    "ignores every environmental and lifestyle factor. A high score does not "
    "mean you will develop the condition, and a low score does not mean you "
    "will not."
)

WHAT_NEXT = (
    ("n1", "Prescription-critical",
     "Show these to a prescriber or pharmacist BEFORE any medication change. "
     "Do not start, stop or adjust a medicine on your own because of this "
     "report."),
    ("n2", "Actionable",
     "Diet, exercise, supplement or monitoring items. Worth raising at your "
     "next appointment. Nothing here is an emergency."),
    ("n3", "Informational",
     "Background only. Interesting to know, and not a reason to do anything."),
)

GLOSSARY = (
    ("zygosity",
     "whether your two copies of a position match. Homozygous means both copies "
     "are the same, heterozygous means they differ, hemizygous means there is "
     "only one copy to read, and a no-call means the test could not read the "
     "position at all."),
    ("carrier",
     "whether you actually have the variant that a classification describes. A "
     "classification describes an allele, not a position, so if you do not have "
     "the allele then the classification does not apply to you."),
    ("magnitude",
     "a 0 to 10 estimate of how interesting a finding is. DNAInsight computes it "
     "from public evidence. It is not a severity score and not a probability, "
     "and a finding nobody could score is shown as unscored rather than as 0."),
    ("repute",
     "the direction of a finding: Good, Bad, or left blank when it is neutral, "
     "mixed or simply unknown. Traits and polygenic scores are always left "
     "blank."),
    ("review stars",
     "how much expert agreement the ClinVar database records for a "
     "classification, from 0 to 4. Four stars is a practice guideline. Zero "
     "stars means a single submitter who stated no criteria."),
    ("palindromic site",
     "an A/T or C/G genotype. Swapping the two alleles for their complements "
     "gives the same pair back, so the strand cannot be verified from the data "
     "and the reading shown may be the complement."),
)


_CSS = """
:root{--good:#60B060;--bad:#FF9090;--unset:#C0C0C0;--blue:#1a3a6b;
--mid:#2980b9;--rx:#c0392b;--act:#e67e22;--info:#2980b9;--bg:#f4f6f9;
--text:#2c3e50;--muted:#7f8c8d;--line:#e3e9ef}
*{box-sizing:border-box}
body{font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:0;
background:var(--bg);color:var(--text);line-height:1.5}
.container{max-width:920px;margin:0 auto;padding:24px}
.header{background:linear-gradient(135deg,#1a3a6b,#2980b9);color:#fff;
padding:30px;border-radius:10px;margin-bottom:20px}
.header h1{margin:0 0 8px 0;font-size:1.8em}
.header .meta{font-size:.9em;opacity:.88;line-height:1.8}
.disclaimer{background:#fff3cd;border:1px solid #ffc107;border-radius:6px;
padding:14px 16px;margin-bottom:18px;font-size:.86em}
.card{background:#fff;border-radius:8px;padding:16px 18px;margin-bottom:18px;
box-shadow:0 1px 4px rgba(0,0,0,.1)}
.card h2{margin:0 0 10px 0;font-size:1.1em;color:var(--blue)}
.card p{margin:0 0 8px 0}
.summary-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;
margin-bottom:18px}
.summary-card{background:#fff;border-radius:8px;padding:16px;text-align:center;
box-shadow:0 1px 4px rgba(0,0,0,.1)}
.summary-card .num{font-size:2em;font-weight:bold}
.silo{margin-bottom:26px}
.silo-head{color:#fff;padding:10px 16px;border-radius:6px 6px 0 0}
.silo-head h2{margin:0;font-size:1.1em}
.silo-head .sub{font-size:.85em;opacity:.9}
.silo-body{border:1px solid var(--line);border-top:none;
border-radius:0 0 6px 6px;background:#fff;padding:12px}
.finding{border:1px solid var(--line);border-left:5px solid var(--unset);
background:#fff;border-radius:0 6px 6px 0;padding:12px 14px;margin-bottom:10px}
.finding.good{border-left-color:var(--good)}
.finding.bad{border-left-color:var(--bad)}
.finding.dub{border-left-style:dashed}
.fh{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:4px}
.rsid{font-family:Consolas,monospace;font-weight:700}
.tok{font-family:Consolas,monospace;color:var(--mid);font-size:.9em}
.gene{color:var(--muted)}
.mag{background:var(--blue);color:#fff;border-radius:11px;padding:1px 9px;
font-size:.76em;font-weight:700}
.mag.zero{background:#b9c1c9}.mag.low{background:#8fa6bd}
.mag.mid{background:var(--mid)}.mag.high{background:#b5341f}
.tag{font-size:.72em;padding:2px 7px;border-radius:9px;background:#ecf0f1;
color:#556}
.tag.cpic{background:#e8f0fe;color:#17408b}
.tag.star{background:#fef6d8;color:#7a5c00}
.tag.flag{background:#fff8e1;color:#7d5900}
.sum{font-weight:600}
.body{font-size:.9em;color:#445}
.crit{font-family:Consolas,monospace;font-size:.8em;background:#f7f9fb;
border:1px solid var(--line);border-radius:4px;padding:6px 8px;margin-top:6px;
white-space:pre-wrap;word-break:break-word}
.nocarry{background:#eef7ee;border:2px solid var(--good);border-radius:6px;
padding:10px 12px;margin:8px 0;font-size:.9em}
.warn{background:#fff8e1;border:1px solid #ffc107;border-radius:5px;
padding:8px 11px;margin-top:8px;font-size:.85em}
.quiet{color:var(--muted);font-size:.82em;margin-top:6px}
.cav{margin-top:8px;background:#f7f9fb;border-left:3px solid var(--muted);
padding:7px 10px;font-size:.85em}
.cav ul{margin:4px 0 0 16px;padding:0}
.cfl{margin-top:8px;background:#fff8e1;border:1px solid #ffc107;
border-radius:5px;padding:8px 11px;font-size:.85em}
.calls{display:flex;flex-wrap:wrap;gap:10px;margin-top:6px}
.callcell{border:1px solid var(--line);background:#fff;border-radius:5px;
padding:5px 9px;min-width:118px}
.calllabel{font-size:.78em;color:var(--muted)}
.callgt{font-family:Consolas,monospace;font-weight:700}
table.kv{width:100%;border-collapse:collapse;font-size:.82em;margin-top:8px}
table.kv td{padding:3px 6px;vertical-align:top;border-top:1px solid var(--line)}
table.kv td:first-child{color:var(--muted);width:32%}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:.86em;margin-bottom:6px}
.dot{display:inline-block;width:11px;height:11px;border-radius:3px;
margin-right:5px}
.next{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}
.next div{border-radius:6px;padding:11px 13px;font-size:.86em;
background:#f7f9fb;border-top:4px solid var(--muted)}
.next .n1{border-top-color:var(--rx)}
.next .n2{border-top-color:var(--act)}
.next .n3{border-top-color:var(--info)}
.next strong{display:block;margin-bottom:4px}
dl.gloss{margin:0;font-size:.87em}
dl.gloss dt{font-weight:700;margin-top:8px}
dl.gloss dd{margin:2px 0 0 0;color:#445}
.footer{text-align:center;font-size:.8em;color:var(--muted);margin-top:28px;
padding-top:14px;border-top:1px solid var(--line);line-height:1.6}
.none{color:var(--muted);font-style:italic;font-size:.9em}
@media(max-width:760px){.summary-grid{grid-template-columns:repeat(2,1fr)}
.next{grid-template-columns:1fr}}
@media print{body{background:#fff}.container{padding:0;max-width:100%}
.card,.summary-card{box-shadow:none;border:1px solid var(--line)}
.finding{break-inside:avoid}}
"""


# ---------------------------------------------------------------------------
# Escaping. The same rule as interactive_report.py, for the same reason.
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


def _magnitude_sort_value(finding: dict) -> float:
    """Sort key for magnitude.

    A null magnitude counts as 1, the convention documented in API_V2 section
    2.2: unscored means nobody assessed it, which belongs above "assessed and
    boring" and below "interesting".
    """
    value = _num(finding.get("magnitude"))
    return 1.0 if value is None else value


def _magnitude_text(finding: dict) -> str:
    """Magnitude with its scale spelled out, or the word unscored.

    Never renders a null as 0. Zero is a real score that a no-call earns; blank
    means nobody could score it, and the two must not read the same.
    """
    value = _num(finding.get("magnitude"))
    if value is None:
        return "unscored"
    return f"{_trim(value)} out of 10"


def _magnitude_class(finding: dict) -> str:
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
    """CSS class driving the coloured left border of a finding block."""
    if str(finding.get("entity_type") or "snp") in NEUTRAL_ENTITIES:
        return "unset"
    repute = str(finding.get("repute") or "")
    if repute == "Good":
        return "good"
    if repute == "Bad":
        return "bad"
    return "unset"


def _repute_text(finding: dict) -> str:
    """Repute in words, saying why it is blank when it is blank."""
    if str(finding.get("entity_type") or "snp") in NEUTRAL_ENTITIES:
        return "grey on purpose, a trait or a score is not good or bad"
    repute = str(finding.get("repute") or "")
    if repute in ("Good", "Bad"):
        return repute
    return "not set, this one is neutral, mixed or unknown"


def _stars_html(finding: dict) -> str:
    """Review stars as filled and hollow star characters, 0 to 4."""
    raw = finding.get("review_stars")
    count = raw if isinstance(raw, int) and 0 <= raw <= 4 else 0
    return "&#9733;" * count + "&#9734;" * (4 - count)


def _star_count(finding: dict) -> int:
    """Review stars as a plain integer, clamped to the documented 0 to 4."""
    raw = finding.get("review_stars")
    return raw if isinstance(raw, int) and 0 <= raw <= 4 else 0


def _confidence_text(finding: dict) -> str:
    """Confidence in words. 'none' is spelled out rather than left blank."""
    value = str(finding.get("confidence") or "").strip().lower()
    if value in ("high", "moderate", "low"):
        return value
    return "none, nothing here is well evidenced"


def _frequency_html(finding: dict) -> str:
    """Frequency: the number, the band, how it was obtained, and from where.

    A null frequency and a 0.0 frequency are different facts. Null means the
    panel holds no data for this genotype. Zero means the panel was checked and
    the genotype was never seen in it. Rendering both as "0%" would invent a
    certainty the data does not have, so they get different words.
    """
    value = _num(finding.get("freq"))
    if value is None:
        return '<span class="none">no data</span>'

    population = _esc(finding.get("freq_population") or "")
    where = f" in {population}" if population else ""
    method = str(finding.get("freq_method") or "")
    if finding.get("freq_derived") or method == "hardy_weinberg":
        how = "derived under Hardy-Weinberg, not counted directly"
    elif method == "observed":
        how = "observed directly in the panel"
    else:
        how = "source method not recorded"

    if value == 0.0:
        return (f"not observed in this panel{where} "
                f'<span class="none">({how})</span>')

    band = str(finding.get("freq_band") or "unknown").replace("_", " ")
    band_html = "" if band == "unknown" else f", {_esc(band)}"
    return (f"{_trim(value)}% of people{where}{band_html} "
            f'<span class="none">({how})</span>')


# ---------------------------------------------------------------------------
# Honesty blocks. Every one of these is inline and cannot be folded away.
# ---------------------------------------------------------------------------

def _carrier_html(finding: dict) -> str:
    """The non-carrier banner, or nothing.

    Presentation choice: a full-width bordered banner at the TOP of the finding,
    not a badge and not a footnote. A ClinVar classification describes an ALLELE,
    not a position, so printing "pathogenic" beside someone who holds two
    reference copies is simply wrong, and it is the most common way a report of
    this kind frightens people for no reason. Being loud is the whole point.
    """
    if finding.get("carrier") is not False:
        return ""
    return (
        '<div class="nocarry"><strong>You do not carry this variant.</strong> '
        "Your genotype does not carry the reported variant, so the "
        "classification described below does not apply to you. It is listed only "
        "so you can see that the position was checked.</div>"
    )


def _strand_html(finding: dict) -> str:
    """Strand notes: a loud one for ambiguity, a quiet one for a routine flip."""
    parts: list[str] = []
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        parts.append(
            '<div class="warn"><strong>The strand cannot be checked here.</strong> '
            "This genotype is A/T or C/G, so the strand cannot be verified from "
            "the data and the reading shown may be the complement of what you "
            "actually have.</div>")
    if finding.get("flipped") or finding.get("freq_flipped"):
        # Quiet by design. A flip is routine bookkeeping, and dressing it up as a
        # warning would train the reader to ignore the warnings that matter.
        parts.append(
            '<div class="quiet">The alleles here were complemented to match the '
            "reference strand. That is routine and is not a problem.</div>")
    return "".join(parts)


def _confirm_html(finding: dict) -> str:
    """False-positive warning for a high magnitude call."""
    value = _num(finding.get("magnitude"))
    if value is None or value < CONFIRM_THRESHOLD:
        return ""
    return (
        '<div class="warn"><strong>Confirm this one before acting on it.</strong> '
        f"A magnitude of {_trim(CONFIRM_THRESHOLD)} or more out of 10 is a strong "
        "claim, and a strong claim from a consumer array is sometimes a false "
        "positive. Ask for a clinically validated test before you or a "
        "prescriber act on it.</div>")


def _caveats_html(finding: dict) -> str:
    """Every caveat string, inline. Never behind a toggle."""
    items = [c for c in (finding.get("caveats") or []) if str(c).strip()]
    single = finding.get("caveat")
    if single and str(single).strip():
        items.append(single)
    if not items:
        return ""
    body = "".join(f"<li>{_esc(c)}</li>" for c in items)
    return f'<div class="cav"><strong>Read with care</strong><ul>{body}</ul></div>'


def _provenance_html(finding: dict) -> str:
    """How many of your files called this position, and any disagreement."""
    parts: list[str] = []
    count = finding.get("count")
    if isinstance(count, int) and count > 1:
        labels = [str(x) for x in (finding.get("labels") or []) if str(x).strip()]
        named = f' ({_esc(", ".join(labels))})' if labels else ""
        parts.append(f'<div class="quiet">{count} of your files called this '
                     f"position{named}.</div>")
    if finding.get("conflict"):
        cells = []
        for call in finding.get("calls") or []:
            genotype = call.get("genotype") or (
                str(call.get("allele1") or "") + str(call.get("allele2") or ""))
            cells.append('<div class="callcell"><div class="calllabel">'
                         + _esc(call.get("label"))
                         + '</div><div class="callgt">' + _esc(genotype)
                         + "</div></div>")
        parts.append('<div class="cfl"><strong>Your files disagree here.</strong> '
                     "Both calls are kept and neither was chosen, because a "
                     "disagreement between two arrays tells you something about "
                     'reliability that a silent merge would destroy.<div class="calls">'
                     + "".join(cells) + "</div></div>")
    return "".join(parts)


def _reliability_html(finding: dict) -> str:
    """Warn when a polygenic score did not have enough of its variants."""
    if str(finding.get("entity_type") or "") != "prs":
        return ""
    if finding.get("reliable") is not False:
        return ""
    return ('<div class="warn"><strong>This score is not reliable for you.</strong> '
            "Your array covered too few of the positions the model needs, so "
            "treat the number below as an illustration and not as a result.</div>")


# ---------------------------------------------------------------------------
# One finding
# ---------------------------------------------------------------------------

def _badges_html(finding: dict) -> str:
    """Small tags after the identifier. Worst news first, deliberately."""
    tags: list[tuple[str, str]] = []
    if finding.get("carrier") is False:
        tags.append(("flag", "you are not a carrier"))
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        tags.append(("flag", "strand ambiguous"))
    if finding.get("conflict"):
        tags.append(("flag", "your files disagree"))
    if str(finding.get("zygosity") or "") == "no_call":
        tags.append(("flag", "no call"))
    if finding.get("flipped"):
        tags.append(("", "strand flipped"))
    stars = _star_count(finding)
    if stars:
        tags.append(("star", f"{_stars_html(finding)} {stars} of 4"))
    level = str(finding.get("cpic_level") or "").strip()
    if level:
        tags.append(("cpic", f"CPIC {_esc(level)}"))
    entity = str(finding.get("entity_type") or "snp")
    if entity != "snp":
        tags.append(("", _esc(entity)))
    return "".join(f'<span class="tag {cls}">{text}</span>' for cls, text in tags)


def _detail_rows(finding: dict) -> str:
    """The per-finding detail table. Empty values are dropped, not left blank."""
    entity = str(finding.get("entity_type") or "snp")
    rows: list[tuple[str, str]] = [
        ("Magnitude", _esc(_magnitude_text(finding))),
        ("Repute", _esc(_repute_text(finding))),
        ("Confidence", _esc(_confidence_text(finding))),
        ("Review stars", f"{_stars_html(finding)} {_star_count(finding)} of 4"),
    ]

    def add(label: str, value: str) -> None:
        if value:
            rows.append((label, value))

    add("CPIC level", _esc(str(finding.get("cpic_level") or "").strip()))
    add("Evidence", _esc(finding.get("evidence")))
    if entity == "snp":
        rows.append(("How common", _frequency_html(finding)))
        gmaf = _num(finding.get("gmaf"))
        add("Global minor allele frequency", "" if gmaf is None else _trim(gmaf))
        if finding.get("chromosome"):
            add("Position", _esc(f"chromosome {finding.get('chromosome')} at "
                                 f"{finding.get('position')}"))
        pubs = finding.get("publications")
        add("Publications", str(pubs) if isinstance(pubs, int) and pubs else "")

    carrier = finding.get("carrier")
    if carrier is True:
        add("Carrier", "yes, you have the reported variant")
    elif carrier is False:
        add("Carrier", "no, you do not carry the reported variant")
    add("Zygosity", _esc(str(finding.get("zygosity") or "").replace("_", " ")))
    copies = finding.get("variant_copies")
    add("Copies of the variant", f"{copies} of 2" if isinstance(copies, int) else "")
    count = finding.get("count")
    add("Files that called it",
        str(count) if isinstance(count, int) and count > 1 else "")
    if finding.get("ambiguous") or finding.get("freq_ambiguous"):
        add("Strand", "cannot be verified, this is a palindromic site")
    elif finding.get("flipped"):
        add("Strand", "complemented to match the reference, which is routine")

    coverage = _num(finding.get("coverage"))
    add("Coverage of the model", "" if coverage is None
        else f"{_trim(coverage * 100)}% of the positions it needs")
    percentile = _num(finding.get("percentile"))
    add("Percentile", "" if percentile is None else _trim(percentile))
    add("Band", _esc(str(finding.get("band") or "").replace("_", " ")))
    if entity == "prs":
        add("Reliable", "yes" if finding.get("reliable") else
            "no, the coverage is below the 90% this model needs")
    matched = finding.get("matched_rsids") or []
    add("Positions matched", str(len(matched)) if matched else "")
    add("Conditions", _esc(finding.get("conditions")))
    for label, key in (("Topics", "topics"), ("Medicines", "medicines")):
        values = [str(v) for v in (finding.get(key) or []) if str(v).strip()]
        add(label, _esc(", ".join(values[:14])))

    body = "".join(f"<tr><td>{label}</td><td>{value}</td></tr>"
                   for label, value in rows)
    return f'<table class="kv"><tbody>{body}</tbody></table>'


def _finding_html(finding: dict) -> str:
    """Render one finding block: honesty banners first, detail table last."""
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
    header.append(f'<span class="mag {_magnitude_class(finding)}">'
                  f"{_esc(_magnitude_text(finding))}</span>")
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
# Sections
# ---------------------------------------------------------------------------

def _snps(findings: list[dict]) -> list[dict]:
    """Positional findings only."""
    return [f for f in findings
            if str(f.get("entity_type") or "snp") not in ("genoset", "trait", "prs")]


def _by_magnitude(items: list[dict]) -> list[dict]:
    """Highest magnitude first, with a null counting as 1."""
    return sorted(items, key=_magnitude_sort_value, reverse=True)


def _silo_section(findings: list[dict], silo: str) -> str:
    """One silo, sorted by magnitude descending.

    v1.2 grouped these by gene. That grouping is dropped here because magnitude
    order within the silo is what a reader needs, and the two orderings fight
    each other. The gene is printed on every card instead, so nothing is lost.
    """
    items = [f for f in _snps(findings) if str(f.get("silo") or "") == silo]
    if not items:
        return ""
    label, colour, description = SILO_LABELS.get(silo, ("Other", "#7f8c8d", ""))
    cards = "".join(_finding_html(f) for f in _by_magnitude(items))
    return (f'<div class="silo"><div class="silo-head" style="background:{colour};">'
            f"<h2>{_esc(label)} ({len(items)})</h2>"
            f'<div class="sub">{_esc(description)}</div></div>'
            f'<div class="silo-body">{cards}</div></div>')


def _non_snp_section(findings: list[dict]) -> str:
    """Genosets, traits and polygenic scores, grouped by kind.

    Kept out of the three silos on purpose. A genoset is a rule across several
    positions and a polygenic score is a statistic over many, so the
    per-position fields a silo card leans on do not apply, and mixing them in
    would imply a comparability that is not there.
    """
    groups = (
        ("genoset", "Rule-based findings",
         "Each of these is a rule over several positions. The rule text is "
         "printed so you can see exactly what had to be true.", ""),
        ("trait", "Traits",
         "A trait is a description, not a verdict. These are always grey, never "
         "green or red, because a trait is neither good nor bad.", ""),
        ("prs", "Polygenic scores",
         "These estimate where you sit in a distribution. Read the warning "
         "first.", PRS_DISCLAIMER),
    )
    blocks: list[str] = []
    for entity, title, intro, disclaimer in groups:
        items = [f for f in findings if str(f.get("entity_type") or "") == entity]
        if not items:
            continue
        # The polygenic disclaimer sits ABOVE the numbers deliberately. A caveat
        # printed after a percentile is a caveat read after the damage is done.
        warning = f'<div class="warn">{_esc(disclaimer)}</div>' if disclaimer else ""
        cards = "".join(_finding_html(f) for f in _by_magnitude(items))
        blocks.append(f'<div class="card"><h2>{_esc(title)} ({len(items)})</h2>'
                      f"<p>{_esc(intro)}</p>{warning}{cards}</div>")
    return "".join(blocks)


def _qc_html(qc: dict) -> str:
    """Data-quality summary, built only from what the scan actually reported."""
    if not qc:
        return ""
    bits: list[str] = []
    pairs = (
        ("flipped", "had their alleles complemented to match the reference, "
                    "which is routine"),
        ("ambiguous", "are A/T or C/G genotypes whose strand cannot be verified "
                      "from the data"),
        ("no_call", "could not be read at all and so score 0"),
        ("conflicts", "are positions where two of your own files disagree"),
    )
    for key, phrase in pairs:
        value = qc.get(key)
        if isinstance(value, int) and value > 0:
            bits.append(f"<li><strong>{value}</strong> {_esc(phrase)}</li>")
    if not bits:
        return ""
    return ('<div class="card"><h2>Data quality</h2>'
            '<ul style="margin:0 0 8px 18px;padding:0;font-size:.88em;">'
            + "".join(bits) + "</ul>"
            '<p style="font-size:.86em;">A palindromic site is an A/T or C/G '
            "genotype: complementing the two alleles gives the same pair back, "
            "so nothing in the file can tell you which strand was reported. "
            "Those calls are held back rather than trusted.</p></div>")


def _how_to_read_html() -> str:
    """The legend, stated in the same terms the interactive report uses."""
    return (
        '<div class="card"><h2>How to read this report</h2>'
        '<div class="legend">'
        f'<span><span class="dot" style="background:{REPUTE_GOOD};"></span>'
        "Green edge: generally favourable</span>"
        f'<span><span class="dot" style="background:{REPUTE_BAD};"></span>'
        "Red edge: generally unfavourable</span>"
        f'<span><span class="dot" style="background:{REPUTE_UNSET};"></span>'
        "Grey edge: neutral, mixed or unknown</span></div>"
        "<p>Traits and polygenic scores are ALWAYS grey, whatever they say, "
        "because a trait is not good or bad. Colouring one would be passing "
        "judgement on you rather than reporting your data.</p>"
        "<p>Magnitude is a 0 to 10 estimate of how interesting a finding is, "
        "computed here from public evidence. Higher means more worth reading, "
        "not more dangerous. A finding nobody could score is shown as unscored, "
        "never as 0, because 0 is a real score that a failed reading earns.</p>"
        "<p>Where a block says your genotype does not carry the reported "
        "variant, the classification attached to that variant does not apply to "
        "you at all.</p>"
        "<p>If two of your own files disagree at a position, no winner is "
        "picked. Both calls are kept and neither was chosen, and you can see "
        "them side by side.</p></div>")


def generate_genetic_report(profile: dict, findings: list[dict],
                            extras: dict | None = None) -> str:
    """Generate the full Genetic Health Report as one self-contained HTML string.

    Args:
        profile:  dict from database.get_profile()
        findings: finding dicts, either v1.2 or v2 shaped. A missing v2 key is
                  simply not rendered rather than guessed at.
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
    generated = _esc(datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"))
    population = _esc(extras.get("population") or "")
    population_note = (f" &nbsp;|&nbsp; <strong>Frequencies for:</strong> {population}"
                       if population else "")

    positional = _snps(rows)
    counts = {key: len([f for f in positional if str(f.get("silo") or "") == key])
              for key in ("pre_prescription", "actionable", "informational")}
    other = len(rows) - len(positional)
    other_note = (f'<p class="none">{other} further entries are rules, traits or '
                  "polygenic scores rather than single positions. They are "
                  "listed in their own section further down.</p>") if other else ""

    next_html = "".join(
        f'<div class="{cls}"><strong>{_esc(title)}</strong>{_esc(text)}</div>'
        for cls, title, text in WHAT_NEXT)
    glossary_html = "".join(f"<dt>{_esc(term)}</dt><dd>{_esc(text)}</dd>"
                            for term, text in GLOSSARY)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genetic Health Report, {name}</title>
<style>{_CSS}</style>
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
      <strong>Generated:</strong> {generated}{population_note}
    </div>
  </div>

  <div class="disclaimer">
    <strong>IMPORTANT DISCLAIMER:</strong> This report is for personal
    informational use only. It is not a substitute for professional medical
    advice, diagnosis, or treatment. Consumer DNA arrays do not perform
    clinical-grade sequencing. Always consult a licensed healthcare provider or
    clinical geneticist before making any medical decisions based on genetic
    data.
  </div>

  <div class="summary-grid">
    <div class="summary-card">
      <div class="num" style="color:#2c3e50;">{len(rows)}</div>
      <div>Total Findings</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#c0392b;">{counts["pre_prescription"]}</div>
      <div>Prescription-Critical</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#e67e22;">{counts["actionable"]}</div>
      <div>Actionable</div>
    </div>
    <div class="summary-card">
      <div class="num" style="color:#2980b9;">{counts["informational"]}</div>
      <div>Informational</div>
    </div>
  </div>

  {_qc_html(extras.get("qc") or {})}

  {_how_to_read_html()}

  <div class="card">
    <h2>What to do next</h2>
    <div class="next">{next_html}</div>
  </div>

  {other_note}

  {_silo_section(rows, "pre_prescription")}
  {_silo_section(rows, "actionable")}
  {_silo_section(rows, "informational")}

  {_non_snp_section(rows)}

  <div class="card">
    <h2>Words used in this report</h2>
    <dl class="gloss">{glossary_html}</dl>
  </div>

  <div class="footer">
    DNAInsight v{_esc(APP_VERSION)} &nbsp;|&nbsp; Open-source personal DNA
    analysis tool &nbsp;|&nbsp; Data sourced from ClinVar, PharmGKB, and
    MyVariant.info<br>
    This report does not constitute medical advice. Not for clinical use.
    DNAInsight is not a medical device.<br><br>
    Magnitude and repute in this report are computed by DNAInsight from CC0 and
    public-domain evidence, including CPIC levels, ClinVar review status,
    population frequency and publication counts. They are not SNPedia values of
    the same name, and they are not clinical severity scores.
  </div>

</div>
</body>
</html>"""
