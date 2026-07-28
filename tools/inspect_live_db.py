"""Read-only inspection of the REAL production database and uploads folder.

Answers: what profiles exist right now, which ones are verification-harness
artifacts, and what's left over in uploads/. Read-only connection throughout;
nothing here can write to dnainsight.db.
"""
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / "dnainsight.db"
UPLOADS = ROOT / "uploads"

HARNESS_NAMES = {
    "DURABILITY CANARY", "Report Verify", "Collision Test", "V2 Verify",
    "Static Report Verify", "Probe",
}
HARNESS_PREFIXES = ("Crafted",)

print("DATABASE: %s" % DB)
print("exists=%s size=%s bytes" % (DB.exists(), DB.stat().st_size if DB.exists() else "n/a"))

if not DB.exists():
    raise SystemExit("no database file, nothing to inspect")

conn = sqlite3.connect("file:%s?mode=ro" % DB.as_posix(), uri=True)
conn.row_factory = sqlite3.Row

profiles = conn.execute(
    "SELECT id, name, dob, sex, provider, created_at FROM profiles ORDER BY id"
).fetchall()
print("\n%d profile(s) total\n" % len(profiles))

harness_ids = []
real_ids = []
for r in profiles:
    name = r["name"] or ""
    is_harness = name in HARNESS_NAMES or name.startswith(HARNESS_PREFIXES)
    tag = "HARNESS" if is_harness else "user-created"
    (harness_ids if is_harness else real_ids).append(r["id"])
    findings_n = conn.execute(
        "SELECT COUNT(*) FROM findings WHERE profile_id=?", (r["id"],)
    ).fetchone()[0]
    uploads_n = conn.execute(
        "SELECT COUNT(*) FROM snp_uploads WHERE profile_id=?", (r["id"],)
    ).fetchone()[0]
    reports_n = conn.execute(
        "SELECT COUNT(*) FROM reports WHERE profile_id=?", (r["id"],)
    ).fetchone()[0]
    print("  id=%-4d %-9s name=%-28r dob=%-11s sex=%-7s provider=%-8s "
          "created=%s  findings=%d uploads=%d reports=%d"
          % (r["id"], tag, name, r["dob"], r["sex"], r["provider"],
             r["created_at"], findings_n, uploads_n, reports_n))

print("\nSUMMARY")
print("  harness/verification profiles: %d  ids=%s" % (len(harness_ids), harness_ids))
print("  everything else (real or unrecognized): %d  ids=%s" % (len(real_ids), real_ids))

conn.close()

print("\nUPLOADS DIRECTORY: %s" % UPLOADS)
if UPLOADS.exists():
    files = sorted(UPLOADS.iterdir(), key=lambda p: p.stat().st_mtime)
    print("  %d file(s)" % len(files))
    for p in files[-15:]:
        print("    %9d bytes  %s" % (p.stat().st_size, p.name[:90]))
else:
    print("  (missing)")
