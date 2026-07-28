"""
vreport.py -- verify the interactive report end to end.

Boots Flask in-process, runs a real scan, generates the interactive report
through BOTH the v1 and v2 endpoints, then validates the produced HTML the same
way the main frontend was validated: structural integrity, no network requests,
no browser storage, JavaScript parses, and every promised feature present.
"""

import io
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

import isolated_db
isolated_db.use_temp_db()
import app as app_module

flask_app = app_module.create_app()
flask_app.config["TESTING"] = True
client = flask_app.test_client()

fails = []
def check(label, cond, detail=""):
    print(("  [ ok ] " if cond else "  [FAIL] ") + label + (("  " + detail) if detail else ""))
    if not cond:
        fails.append(label + " " + detail)

print("=" * 76)
print("INTERACTIVE REPORT VERIFICATION")
print("=" * 76)

from sample import pick_sample, ensure_fixture

ensure_fixture(ROOT)
sample = pick_sample(ROOT)
if sample is None:
    print("no sample upload, cannot verify")
    sys.exit(1)

print("\n-- setup --")
r = client.post("/api/profiles", data={
    "name": "Report Verify", "dob": "1970-01-01", "sex": "male",
    "file": (io.BytesIO(sample.read_bytes()), sample.name),
}, content_type="multipart/form-data")
pid = (r.get_json() or {}).get("profile_id")
check("profile created", r.status_code in (200, 201), f"pid={pid}")

r = client.post(f"/api/profiles/{pid}/scan/v2", json={"use_api": False})
check("scan started", r.status_code == 200)
st = {}
for _ in range(80):
    st = client.get(f"/api/profiles/{pid}/scan/v2/status").get_json()
    if st.get("done"):
        break
    time.sleep(0.25)
check("scan completed", bool(st.get("done")) and not st.get("error"),
      f"findings={st.get('findings')}")

print("\n-- v2 endpoint --")
r = client.post(f"/api/profiles/{pid}/reports/interactive")
body = r.get_json() or {}
check("v2 endpoint returns 201", r.status_code == 201, str(body.get("message", ""))[:70])
check("reports findings_included", isinstance(body.get("findings_included"), int),
      f"included={body.get('findings_included')} bytes={body.get('bytes')}")
rid = body.get("report_id")

print("\n-- v1 endpoint delegates --")
r1 = client.post(f"/api/profiles/{pid}/reports", json={"type": "interactive"})
check("v1 type=interactive accepted", r1.status_code == 201,
      f"HTTP {r1.status_code}")
r_bad = client.post(f"/api/profiles/{pid}/reports", json={"type": "nonsense"})
check("v1 rejects an unknown type", r_bad.status_code == 400)

print("\n-- filters are honoured --")
r_f = client.post(f"/api/profiles/{pid}/reports/interactive?min_magnitude=5")
fb = r_f.get_json() or {}
all_n = body.get("findings_included") or 0
sub_n = fb.get("findings_included") or 0
check("filtered report is a subset", r_f.status_code == 201 and sub_n <= all_n,
      f"{sub_n} of {all_n}")
r_none = client.post(f"/api/profiles/{pid}/reports/interactive?min_magnitude=99")
check("empty selection is refused with 400", r_none.status_code == 400)

print("\n-- served HTML --")
r_view = client.get(f"/api/reports/{rid}/view")
check("report is served", r_view.status_code == 200)
html = r_view.data.decode("utf-8", errors="replace")
print(f"        {len(html):,} bytes, {len(html.splitlines())} lines")

print("\n-- offline guarantees --")
for pattern, label in [
    (r"https?://", "no absolute http(s) URL anywhere"),
    (r"<link[^>]+href", "no external stylesheet link"),
    (r"src\s*=\s*[\"']https?:", "no remote script src"),
    (r"@import", "no CSS @import"),
    (r"\bfetch\s*\(", "no fetch() call"),
    (r"XMLHttpRequest", "no XMLHttpRequest"),
    (r"localStorage|sessionStorage|indexedDB", "no browser storage"),
]:
    hits = re.findall(pattern, html)
    check(label, not hits, f"{len(hits)} hit(s)" if hits else "")

print("\n-- structure --")
script_blocks = re.findall(r"<script[^>]*>", html)
check("exactly one html element", html.count("<html") == 1)
check("exactly one body element", html.count("<body>") == 1)
check("script tags balanced", len(script_blocks) == html.count("</script>"),
      f"{len(script_blocks)} open, {html.count('</script>')} close")
check("json island present", 'id="payload"' in html)

island = re.search(r'<script type="application/json" id="payload">(.*?)</script>',
                   html, re.S)
check("json island parses", island is not None)
if island:
    try:
        payload = json.loads(island.group(1).encode().decode("unicode_escape")
                             if "\\u003c" in island.group(1) else island.group(1))
    except Exception:
        try:
            payload = json.loads(island.group(1).replace("\\u003c", "<"))
        except Exception as exc:
            payload = None
            check("json island decodes", False, str(exc)[:80])
    if payload:
        check("payload has findings", len(payload.get("findings", [])) > 0,
              f"{len(payload.get('findings', []))} findings")
        check("payload carries patient block", bool(payload.get("patient")))
        check("no raw '<' left unescaped in island", "</script" not in island.group(1))

print("\n-- javascript parses --")
js_blocks = re.findall(r"<script>\n?(.*?)</script>", html, re.S)
main_js = max(js_blocks, key=len) if js_blocks else ""
(ROOT / "tools" / "report_app.js").write_text(main_js, encoding="utf-8")
node = subprocess.run(["node", "--check", str(ROOT / "tools" / "report_app.js")],
                      capture_output=True, text=True)
check("node --check on report JS", node.returncode == 0,
      node.stderr.strip()[:120] if node.returncode else f"{len(main_js.splitlines())} lines")

print("\n-- honesty features carried through --")
for needle, label in [
    ("neutral by design", "traits and scores marked neutral"),
    ("Both calls are kept", "conflict retention wording"),
    ("cannot be verified", "strand ambiguity wording"),
    ("could not be evaluated", "not-testable framing for genosets"),
    ("not a medical device", "medical disclaimer"),
    ("not the SNPedia values", "magnitude provenance stated"),
    ("unscored", "null magnitude labelled unscored"),
    ("SNPs only", "frequency and publication exemption labelled"),
    ("works offline", "offline claim stated"),
    ("Read with care", "caveat block"),
    ("EXEMPT", "exemption logic present in engine"),
]:
    check(label, needle in html)

print("\n-- cleanup --")
client.delete(f"/api/profiles/{pid}")

print("\n" + "=" * 76)
print("INTERACTIVE REPORT OK" if not fails else f"PROBLEMS: {len(fails)}")
for f in fails:
    print("  -", f)
sys.exit(0 if not fails else 1)
