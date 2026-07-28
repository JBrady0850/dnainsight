"""
vreports.py -- verify the two static reports end to end.

Boots Flask in-process with the test client, creates a profile from a real
upload, runs a v2 scan, generates BOTH static reports through the documented v1
endpoint, fetches each through the view endpoint, and validates the HTML:
structure, offline guarantees, the v2 honesty features, and escaping.

A second pass calls the two generators directly with a crafted finding set, so
the cases a sample upload may not contain (a non-carrier, a null magnitude, a
0.0 frequency, a pooled conflict, a palindromic call, a genoset, a trait and an
unreliable polygenic score) are always exercised, and so a hostile
interpretation string can be proven to come out escaped.

NOTE ON THE RETRY LOOP: importing backend.database probes its candidate paths by
creating and then UNLINKING a test database at the same location, so any other
python process that imports it deletes this run's dnainsight.db. That is
pre-existing behaviour and not what this script is testing, so the endpoint pass
re-creates the schema before each step and is retried if the file is pulled out
from under it.
"""

import io
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

import isolated_db
isolated_db.use_temp_db()
import app as app_module
from backend import database as db
from backend import genetic_report as gr
from backend import doctor_report as dr

flask_app = app_module.create_app()
flask_app.config["TESTING"] = True
client = flask_app.test_client()

fails: list[str] = []


def check(label: str, cond: bool, detail: str = "") -> bool:
    print(("  [ ok ] " if cond else "  [FAIL] ") + label
          + (("  " + detail) if detail else ""))
    if not cond:
        fails.append(label + (" " + detail if detail else ""))
    return bool(cond)


def ensure_schema() -> None:
    """Recreate the schema if something removed the database file underneath us."""
    db.init_db()


HOSTILE = "<script>alert('xss')</script> & \"quoted\" <b>bold</b>"

REQUIRED_STRINGS = (
    "unscored", "does not carry", "cannot be verified", "Both calls are kept",
    "not SNPedia", "not a medical device",
)

REPUTE_HEXES = ("#60B060", "#FF9090", "#C0C0C0")

FORBIDDEN_PATTERNS = (
    (r"https?://", "no absolute http(s) URL"),
    (r"<link[^>]+href", "no external stylesheet link"),
    (r"src\s*=\s*[\"']?https?:", "no remote src"),
    (r"@import", "no CSS @import"),
)


def structure(kind: str, html: str) -> None:
    """Structural, offline and required-wording checks common to both reports."""
    check(f"{kind}: one <html", html.count("<html") == 1, str(html.count("<html")))
    check(f"{kind}: one </html>", html.count("</html>") == 1)
    check(f"{kind}: one <body>", html.count("<body>") == 1)
    check(f"{kind}: one </body>", html.count("</body>") == 1)
    for tag in ("div", "table", "style", "dl", "ul"):
        opens = len(re.findall(rf"<{tag}[\s>]", html))
        closes = html.count(f"</{tag}>")
        check(f"{kind}: <{tag}> tags balanced", opens == closes,
              f"{opens} open, {closes} close")
    for pattern, label in FORBIDDEN_PATTERNS:
        hits = re.findall(pattern, html)
        check(f"{kind}: {label}", not hits, f"{len(hits)} hit(s)")
    check(f"{kind}: no script tag anywhere, so no finding text can carry one",
          "<script" not in html)
    for hexcode in REPUTE_HEXES:
        check(f"{kind}: repute colour {hexcode} present", hexcode in html)
    for needle in REQUIRED_STRINGS:
        check(f"{kind}: says '{needle}'", needle in html)


def endpoint_pass(sample: Path) -> dict | None:
    """One full endpoint run. Returns None when the database was removed."""
    ensure_schema()
    resp = client.post("/api/profiles", data={
        "name": "Static Report Verify", "dob": "1970-01-01", "sex": "female",
        "file": (io.BytesIO(sample.read_bytes()), sample.name),
    }, content_type="multipart/form-data")
    pid = (resp.get_json() or {}).get("profile_id")
    if resp.status_code not in (200, 201) or not pid:
        print(f"        profile creation failed, HTTP {resp.status_code}")
        return None

    ensure_schema()
    resp = client.post(f"/api/profiles/{pid}/scan/v2", json={"use_api": False})
    if resp.status_code != 200:
        print(f"        scan refused, HTTP {resp.status_code}")
        return None
    state: dict = {}
    for _ in range(120):
        state = client.get(f"/api/profiles/{pid}/scan/v2/status").get_json() or {}
        if state.get("done"):
            break
        time.sleep(0.25)
    if not state.get("done") or state.get("error"):
        print(f"        scan did not finish: {state}")
        return None

    result: dict = {"pid": pid, "findings": state.get("findings"),
                    "status": {}, "view": {}, "html": {}}
    for kind in ("genetic", "doctor"):
        ensure_schema()
        made = client.post(f"/api/profiles/{pid}/reports", json={"type": kind})
        result["status"][kind] = made.status_code
        body = made.get_json() or {}
        if made.status_code != 201:
            print(f"        {kind} generation failed: {body}")
            return None
        view = client.get(f"/api/reports/{body.get('report_id')}/view")
        result["view"][kind] = view.status_code
        if view.status_code != 200:
            print(f"        {kind} view failed, HTTP {view.status_code}")
            return None
        result["html"][kind] = view.data.decode("utf-8", errors="replace")

    ensure_schema()
    found = client.get(f"/api/profiles/{pid}/findings/v2?limit=0").get_json() or {}
    result["non_carriers"] = [f for f in found.get("findings", [])
                              if f.get("carrier") is False]
    return result


print("=" * 78)
print("STATIC REPORT VERIFICATION (genetic_report.py, doctor_report.py)")
print("=" * 78)

# Largest upload, because it carries the most positions and therefore the best
# chance of a real non-carrier and a real palindromic call. The shortest-name
# tie-break keeps the selection self-limiting; see tools/sample.py.
from sample import pick_richest_sample, ensure_fixture

ensure_fixture(ROOT)
sample = pick_richest_sample(ROOT)
if sample is None:
    print("no sample upload, cannot verify")
    sys.exit(1)
print(f"\nsample upload: {sample.name} ({sample.stat().st_size:,} bytes)")

run: dict | None = None
for attempt in (1, 2, 3):
    print(f"\n-- endpoint pass, attempt {attempt} --")
    run = endpoint_pass(sample)
    if run:
        break
if not run:
    print("\nFAIL: the endpoint pass could not complete, see messages above")
    sys.exit(1)

print(f"        profile {run['pid']}, {run['findings']} findings")
served = run["html"]
for kind in ("genetic", "doctor"):
    print(f"\n-- {kind} report --")
    check(f"{kind}: POST returns 201", run["status"][kind] == 201,
          str(run["status"][kind]))
    check(f"{kind}: GET view returns 200", run["view"][kind] == 200,
          str(run["view"][kind]))
    html = served[kind]
    print(f"        {len(html):,} bytes, {len(html.splitlines())} lines")
    structure(kind, html)

print("\n-- doctor report specifics --")
doc = served["doctor"]
for needle, label in (
    ("Prescription-Critical Variants", "prescription-critical section"),
    ("Phenotype implication", "prescription table phenotype column"),
    ("Affected drugs", "prescription table affected drugs column"),
    ("CPIC level", "prescription table CPIC column"),
    ("Review stars", "prescription table review stars column"),
    ("Zygosity", "prescription table zygosity column"),
    ("AI-Assisted Analysis Prompt", "AI prompt block"),
    ("Patient findings (JSON):", "AI prompt payload"),
    ("Confidence and Limitations", "confidence and limitations block"),
    ("do not call star alleles completely", "star allele limitation stated"),
    ("Recommended Laboratory Follow-Ups", "lab follow-up section"),
    ("Drug Interaction Summary by Medicine", "per-medicine drug grouping"),
):
    check(f"doctor: {label}", needle in doc)

# The prompt JSON is HTML-escaped like everything else, so the quotes around each
# key arrive as &quot;. The browser renders them back to quotes, so what the user
# copies is valid JSON, but the check has to look for the escaped form.
for key in ("magnitude", "repute", "cpic_level", "review_stars", "carrier",
            "confidence", "strand_ambiguous"):
    check(f"doctor: AI prompt carries {key}",
          f"&quot;{key}&quot;:" in doc or f'"{key}":' in doc)

print("\n-- genetic report specifics --")
gen = served["genetic"]
for needle, label in (
    ("What to do next", "what-to-do-next block"),
    ("Words used in this report", "glossary block"),
    ("palindromic site", "glossary defines palindromic site"),
    ("zygosity", "glossary defines zygosity"),
    ("review stars", "glossary defines review stars"),
    ("Prescription-Critical", "prescription silo retained"),
    ("Informational", "informational silo retained"),
    ("IMPORTANT DISCLAIMER", "v1.2 disclaimer retained"),
    ("out of 10", "magnitude scale made explicit"),
):
    check(f"genetic: {label}", needle in gen)

print("\n-- a real non-carrier finding from the scan --")
non_carriers = run.get("non_carriers") or []
print(f"        {len(non_carriers)} non-carrier finding(s) in this scan")
if non_carriers:
    rsid = str(non_carriers[0].get("rsid"))
    check("genetic: the non-carrier finding is rendered", rsid in gen, rsid)
    check("genetic: non-carrier wording present",
          "You do not carry this variant" in gen)
    check("doctor: non-carrier wording present",
          "NOT a carrier of this variant" in doc)
else:
    print("        sample has none, the crafted pass below covers this case")

print("\n-- crafted findings, generators called directly --")
crafted = [
    {"rsid": "rs_noncarrier", "entity_type": "snp", "gene": "CYP2C19",
     "genotype": "GG", "token": "(G;G)", "zygosity": "homozygous",
     "silo": "pre_prescription", "magnitude": 6.5, "repute": "Bad",
     "confidence": "high", "review_stars": 3, "cpic_level": "A",
     "evidence": "CPIC Level A", "carrier": False, "variant_copies": 0,
     "freq": None, "summary": HOSTILE, "interpretation": HOSTILE,
     "caveats": [HOSTILE], "medicines": ["clopidogrel"], "count": 1,
     "clinical_sig": "pathogenic"},
    {"rsid": "rs_palindromic", "entity_type": "snp", "gene": "MTHFR",
     "genotype": "AT", "zygosity": "heterozygous", "silo": "actionable",
     "magnitude": None, "repute": "", "confidence": "none", "review_stars": 0,
     "ambiguous": True, "freq": 0.0, "freq_population": "CEU",
     "freq_method": "observed", "carrier": True, "variant_copies": 1,
     "summary": "Palindromic case", "flipped": True},
    {"rsid": "rs_conflict", "entity_type": "snp", "gene": "VKORC1",
     "genotype": "AG", "silo": "informational", "magnitude": 5.0,
     "repute": "Good", "conflict": True, "count": 2,
     "labels": ["file A", "file B"], "review_stars": 1,
     "calls": [{"label": "file A", "genotype": "AG"},
               {"label": "file B", "genotype": "AA"}],
     "freq": 12.5, "freq_band": "common", "freq_population": "CEU",
     "freq_derived": True, "freq_method": "hardy_weinberg"},
    {"rsid": "rs_nullmag", "entity_type": "snp", "gene": "APOE",
     "genotype": "CT", "silo": "informational", "magnitude": None},
    {"rsid": "rs_zeromag", "entity_type": "snp", "gene": "NAT2",
     "genotype": "NN", "zygosity": "no_call", "silo": "informational",
     "magnitude": 0.0},
    {"rsid": "dgs_demo", "entity_type": "genoset", "silo": "informational",
     "magnitude": 3.0, "repute": "Bad", "criteria": HOSTILE,
     "summary": "Rule based demo", "matched_rsids": ["rs1", "rs2"]},
    {"rsid": "trait_demo", "entity_type": "trait", "silo": "informational",
     "magnitude": 1.0, "repute": "Bad", "summary": "Trait demo",
     "caveats": ["A trait caveat that must appear inline."]},
    {"rsid": "prs_demo", "entity_type": "prs", "silo": "informational",
     "magnitude": None, "repute": "Bad", "percentile": 88.0, "band": "high",
     "coverage": 0.42, "reliable": False, "summary": "Score demo",
     "caveats": ["A score caveat that must appear inline."]},
]
profile = {"name": "Crafted " + HOSTILE, "dob": "1980-02-03", "sex": "male",
           "provider": "generic"}
qc = {"flipped": 7, "ambiguous": 3, "no_call": 5, "conflicts": 2}
extras = {"qc": qc, "population": "CEU"}

direct = {
    "genetic-direct": gr.generate_genetic_report(profile, crafted, extras),
    "doctor-direct": dr.generate_doctor_report(profile, crafted, extras),
}
for kind, html in direct.items():
    print(f"\n   {kind}: {len(html):,} bytes")
    structure(kind, html)
    check(f"{kind}: hostile markup escaped", "&lt;script&gt;" in html)
    check(f"{kind}: hostile markup is not live", "<b>bold</b>" not in html)
    check(f"{kind}: null magnitude reads unscored", "unscored" in html)
    check(f"{kind}: null frequency reads no data", "no data" in html)
    check(f"{kind}: zero frequency reads not observed",
          "not observed in this panel" in html)
    check(f"{kind}: Hardy-Weinberg derivation stated", "Hardy-Weinberg" in html)
    check(f"{kind}: source population named", "CEU" in html)
    check(f"{kind}: non-carrier wording present", "does not carry" in html)
    check(f"{kind}: conflict keeps both calls", "Both calls are kept" in html)
    check(f"{kind}: both conflicting genotypes shown",
          "file A" in html and "file B" in html)
    check(f"{kind}: strand ambiguity warned", "cannot be verified" in html)
    check(f"{kind}: routine flip noted quietly", "complemented" in html)
    check(f"{kind}: high magnitude flagged for confirmation",
          "Confirm" in html and "false positive" in html)
    check(f"{kind}: filled and hollow stars rendered",
          "&#9733;" in html and "&#9734;" in html)
    check(f"{kind}: CPIC level shown", "CPIC A" in html)
    check(f"{kind}: evidence label shown", "CPIC Level A" in html)
    check(f"{kind}: confidence shown", "high" in html)
    check(f"{kind}: caveats inline, nothing folded away",
          "caveat that must appear inline" in html and "<details" not in html)
    check(f"{kind}: multi-file provenance stated",
          ("2 of your files" in html) or ("Called by 2 pooled" in html))
    check(f"{kind}: qc counters rendered",
          ">7<" in html and ">3<" in html and ">5<" in html)
    check(f"{kind}: palindromic sites explained", "palindromic" in html.lower())
    check(f"{kind}: genoset section carries the rule text", "dgs_demo" in html)
    check(f"{kind}: trait grouped separately", "trait_demo" in html)
    check(f"{kind}: prs percentile shown", "88" in html)
    check(f"{kind}: prs band shown", "high" in html)
    check(f"{kind}: prs coverage shown", "42" in html)
    check(f"{kind}: prs reliability stated", "reliable" in html.lower())
    disclaimer_at = html.find("statistical predictor")
    score_at = html.find("prs_demo")
    check(f"{kind}: prs disclaimer sits above the results",
          0 <= disclaimer_at < score_at, f"{disclaimer_at} vs {score_at}")
    # Card order, not first mention: the doctor report legitimately names rsIDs
    # in the lab follow-up list before the cards appear.
    order = [html.find(f'<span class="rsid">{x}</span>')
             for x in ("rs_conflict", "rs_nullmag", "rs_zeromag")]
    check(f"{kind}: silo sorted by magnitude with null as 1",
          order == sorted(order) and -1 not in order, str(order))

print("\n-- traits and polygenic scores are always grey --")
for name, module in (("genetic", gr), ("doctor", dr)):
    for entity in ("trait", "prs"):
        cls = module._repute_class({"entity_type": entity, "repute": "Bad"})
        check(f"{name}: {entity} with repute Bad still renders grey",
              cls == "unset", cls)
    check(f"{name}: a Bad SNP still renders red",
          module._repute_class({"entity_type": "snp", "repute": "Bad"}) == "bad")
    body = direct[f"{name}-direct"]
    check(f"{name}: says why traits and scores are grey",
          ("not good or bad" in body) or ("no direction of effect" in body))

print("\n-- backward compatibility of the signatures --")
check("genetic: two-argument call still works",
      gr.generate_genetic_report(profile, crafted).startswith("<!DOCTYPE"))
check("doctor: two-argument call still works",
      dr.generate_doctor_report(profile, crafted).startswith("<!DOCTYPE"))
check("genetic: an empty finding list does not raise",
      gr.generate_genetic_report(profile, []).count("<html") == 1)
check("doctor: an empty finding list does not raise",
      dr.generate_doctor_report(profile, []).count("<html") == 1)

print("\n-- no em or en dashes --")
for kind, html in list(served.items()) + list(direct.items()):
    check(f"{kind}: no em or en dash",
          ("—" not in html) and ("–" not in html))

print("\n-- cleanup --")
try:
    client.delete(f"/api/profiles/{run['pid']}")
except Exception as exc:
    print(f"        cleanup skipped: {exc}")

print("\n" + "=" * 78)
if fails:
    print(f"FAIL: {len(fails)} check(s) failed")
    for item in fails:
        print("  -", item)
else:
    print("PASS: every static report check passed")
print("=" * 78)
sys.exit(0 if not fails else 1)
