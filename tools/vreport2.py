"""vreport2.py -- prove report filenames no longer collide within one second."""
import io
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
        fails.append(label)

print("=" * 70)
print("REPORT FILENAME COLLISION TEST")
print("=" * 70)

from sample import pick_sample, ensure_fixture

ensure_fixture(ROOT)
sample = pick_sample(ROOT)
r = client.post("/api/profiles", data={
    "name": "Collision Test", "dob": "1980-01-01", "sex": "female",
    "file": (io.BytesIO(sample.read_bytes()), sample.name),
}, content_type="multipart/form-data")
pid = r.get_json()["profile_id"]

client.post(f"/api/profiles/{pid}/scan/v2", json={"use_api": False})
for _ in range(80):
    if client.get(f"/api/profiles/{pid}/scan/v2/status").get_json().get("done"):
        break
    time.sleep(0.25)

# Fire six reports as fast as possible, deliberately inside one second.
print("\n-- six reports generated back to back --")
made = []
for i in range(6):
    kind = ["genetic", "doctor", "interactive"][i % 3]
    resp = client.post(f"/api/profiles/{pid}/reports", json={"type": kind})
    body = resp.get_json() or {}
    made.append((body.get("report_id"), kind, body.get("filename")))
    print(f"     rid={body.get('report_id')} {kind:<12} {body.get('filename')}")

names = [m[2] for m in made if m[2]]
check("all six produced a filename", len(names) == 6, f"{len(names)} of 6")
check("every filename is distinct", len(set(names)) == len(names),
      f"{len(set(names))} distinct of {len(names)}")

# Each report id must serve its OWN content, not the last one written.
print("\n-- each report id serves its own content --")
bodies = {}
for rid, kind, _name in made:
    if rid is None:
        continue
    resp = client.get(f"/api/reports/{rid}/view")
    bodies[rid] = resp.data
    check(f"rid {rid} serves 200", resp.status_code == 200)

interactive_ids = [rid for rid, kind, _ in made if kind == "interactive"]
genetic_ids = [rid for rid, kind, _ in made if kind == "genetic"]
if interactive_ids and genetic_ids:
    a, b = bodies.get(interactive_ids[0], b""), bodies.get(genetic_ids[0], b"")
    check("an interactive report is not serving a genetic one",
          (b'id="payload"' in a) and (b'id="payload"' not in b),
          f"interactive {len(a):,}B, genetic {len(b):,}B")

check("no two ids serve byte-identical content across types",
      len({bodies[k] for k in bodies}) >= 3,
      f"{len({bodies[k] for k in bodies})} distinct bodies from {len(bodies)} ids")

print("\n-- cleanup --")
client.delete(f"/api/profiles/{pid}")

print("\n" + "=" * 70)
print("FILENAME COLLISION FIXED" if not fails else f"PROBLEMS: {len(fails)}")
for f in fails:
    print("  -", f)
sys.exit(0 if not fails else 1)
