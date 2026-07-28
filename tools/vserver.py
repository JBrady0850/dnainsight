"""vserver.py -- boot the Flask app in-process and exercise every v2 endpoint.

Uses Flask's test client, so nothing binds a port and nothing touches the
network. Verifies the app actually starts with both blueprints registered and
that the v2 contract paths respond.
"""
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

import isolated_db
isolated_db.use_temp_db()
import app as app_module

flask_app = app_module.create_app()
flask_app.config["TESTING"] = True
client = flask_app.test_client()

print("=" * 76)
print("DNAInsight v2 SERVER VERIFICATION (in-process test client)")
print("=" * 76)

rules = sorted({str(r.rule) for r in flask_app.url_map.iter_rules()})
print(f"\nregistered routes: {len(rules)}")
v2_only = [r for r in rules if "/v2" in r or any(
    k in r for k in ("capabilities", "populations", "sources", "facets",
                     "genosets", "traits", "prs", "pgx", "conflicts", "trio",
                     "qc", "snpedia"))]
print(f"v2 routes: {len(v2_only)}")
for r in v2_only:
    print("   ", r)

results = []


def check(label, method, path, expect, **kw):
    fn = getattr(client, method)
    resp = fn(path, **kw)
    ok = resp.status_code in (expect if isinstance(expect, tuple) else (expect,))
    results.append((label, resp.status_code, ok))
    body = ""
    try:
        data = resp.get_json()
        if isinstance(data, dict):
            keys = list(data)[:7]
            body = f"keys={keys}"
        elif isinstance(data, list):
            body = f"list of {len(data)}"
    except Exception:
        body = f"{len(resp.data)} bytes"
    flag = "ok  " if ok else "FAIL"
    print(f"  [{flag}] {resp.status_code} {method.upper():<6} {path:<52} {body}")
    return resp


print("\n-- system --")
check("status", "get", "/api/status", 200)
check("version", "get", "/api/version", 200)
check("capabilities", "get", "/api/capabilities", 200)
r = check("populations", "get", "/api/populations", 200)
pops = r.get_json()
print(f"        populations={len(pops['populations'])} "
      f"available={sum(1 for p in pops['populations'] if p['available'])} "
      f"default={pops['default']} modes={pops['aggregate_modes']}")

print("\n-- create a profile from a real upload file --")
from sample import pick_sample, ensure_fixture

ensure_fixture(ROOT)
sample = pick_sample(ROOT)
if sample is None:
    print("  no sample upload available, cannot continue")
    sys.exit(1)
data = {
    "name": "V2 Verify",
    "dob": "1970-01-01",
    "sex": "male",
    "file": (io.BytesIO(sample.read_bytes()), sample.name),
}
r = check("create profile", "post", "/api/profiles", (200, 201),
          data=data, content_type="multipart/form-data")
pid = (r.get_json() or {}).get("profile_id")
print(f"        profile_id={pid}")

print("\n-- sources --")
check("list sources", "get", f"/api/profiles/{pid}/sources", 200)
data2 = {
    "file": (io.BytesIO(sample.read_bytes()), "second_self.txt"),
    "role": "self",
    "label": "second self file",
}
check("add self source (pooling)", "post", f"/api/profiles/{pid}/sources",
      (200, 201), data=data2, content_type="multipart/form-data")
data3 = {
    "file": (io.BytesIO(sample.read_bytes()), "mother.txt"),
    "role": "mother",
}
check("add mother source", "post", f"/api/profiles/{pid}/sources",
      (200, 201), data=data3, content_type="multipart/form-data")
check("bad role rejected", "post", f"/api/profiles/{pid}/sources", 400,
      data={"file": (io.BytesIO(b"x"), "x.txt"), "role": "cousin"},
      content_type="multipart/form-data")
r = check("list sources again", "get", f"/api/profiles/{pid}/sources", 200)
print(f"        sources={len(r.get_json()['sources'])}")

print("\n-- v2 scan (synchronous wait) --")
check("start scan", "post", f"/api/profiles/{pid}/scan/v2", 200,
      json={"use_api": False, "population": "CEU"})
import time
for _ in range(60):
    st = client.get(f"/api/profiles/{pid}/scan/v2/status").get_json()
    if st.get("done"):
        break
    time.sleep(0.25)
print(f"        final status: {json.dumps(st)}")
if st.get("error"):
    print("        SCAN ERROR, aborting")
    sys.exit(1)

print("\n-- findings and filtering --")
r = check("findings", "get", f"/api/profiles/{pid}/findings/v2", 200)
body = r.get_json()
print(f"        total={body['total']} returned={body['returned']} "
      f"ranges={body['ranges']} population={body['population']}")
print(f"        summary={json.dumps(body['summary'])}")
print(f"        qc flipped={body['qc'].get('flipped')} "
      f"ambiguous={body['qc'].get('ambiguous')} conflicts={body['qc'].get('conflicts')}")
check("filter min_magnitude", "get",
      f"/api/profiles/{pid}/findings/v2?min_magnitude=2", 200)
check("filter repute", "get", f"/api/profiles/{pid}/findings/v2?repute=bad", 200)
check("filter entity", "get",
      f"/api/profiles/{pid}/findings/v2?entity_type=genoset,prs", 200)
check("sort asc", "get",
      f"/api/profiles/{pid}/findings/v2?sort=frequency&order=asc", 200)
check("region query", "get", f"/api/profiles/{pid}/findings/v2?q=chr1", 200)
check("operator query", "get", f"/api/profiles/{pid}/findings/v2?q=/MAG>=2", 200)
check("population switch", "get",
      f"/api/profiles/{pid}/findings/v2?population=JPT", 200)
check("facets", "get", f"/api/profiles/{pid}/facets", 200)

print("\n-- subsystem views --")
r = check("genosets", "get", f"/api/profiles/{pid}/genosets", 200)
g = r.get_json()
print(f"        matched={len(g['matched'])} unmatched={len(g['unmatched'])} "
      f"not_testable={len(g['incomplete'])}")
check("traits", "get", f"/api/profiles/{pid}/traits", 200)
check("prs", "get", f"/api/profiles/{pid}/prs", 200)
check("pgx", "get", f"/api/profiles/{pid}/pgx", 200)
check("conflicts", "get", f"/api/profiles/{pid}/conflicts", 200)
check("trio", "get", f"/api/profiles/{pid}/trio", 200)
check("qc", "get", f"/api/profiles/{pid}/qc", 200)

print("\n-- export honouring filters --")
for fmt in ("json", "csv", "tsv"):
    r = check(f"export {fmt}", "get",
              f"/api/profiles/{pid}/export/v2/{fmt}?min_magnitude=1", 200)
    print(f"        {len(r.data):,} bytes")
check("bad export format", "get", f"/api/profiles/{pid}/export/v2/xml", 400)

print("\n-- snpedia licence gate --")
check("snpedia status", "get", "/api/admin/snpedia/status", 200)
r = check("harvest without acceptance MUST be 403", "post",
          "/api/admin/snpedia/harvest", 403, json={})
gate = r.get_json() or {}
print(f"        license={gate.get('license')}")
print(f"        notice present: {bool(gate.get('notice'))}")

print("\n-- v1 endpoints still intact --")
check("v1 findings", "get", f"/api/profiles/{pid}/findings", 200)
check("v1 lookup", "get", f"/api/profiles/{pid}/lookup/rs1801133", 200)
check("v1 db-status", "get", "/api/admin/db-status", 200)
check("v1 export csv", "get", f"/api/profiles/{pid}/export/csv", 200)

print("\n-- reports --")
for kind in ("genetic", "doctor"):
    r = check(f"generate {kind} report", "post", f"/api/profiles/{pid}/reports",
              (200, 201), json={"type": kind})
    rid = (r.get_json() or {}).get("report_id")
    if rid:
        check(f"view {kind} report", "get", f"/api/reports/{rid}/view", 200)

print("\n-- cleanup --")
check("delete profile", "delete", f"/api/profiles/{pid}", 200)

failed = [r for r in results if not r[2]]
print("\n" + "=" * 76)
print(f"{len(results) - len(failed)}/{len(results)} endpoint checks passed")
if failed:
    print("FAILED:")
    for label, code, _ in failed:
        print(f"  {label}: HTTP {code}")
print("SERVER OK" if not failed else "SERVER PROBLEMS FOUND")
sys.exit(0 if not failed else 1)
