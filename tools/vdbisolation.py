"""Prove the harnesses no longer touch the real dnainsight.db.

Snapshots the real DB's profile count and mtime, runs every harness that used
to write to it, then re-checks. Any change means isolation failed.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
DB = ROOT / "dnainsight.db"
UPLOADS = ROOT / "uploads"

sys.path.insert(0, str(ROOT))
from backend import database as db_module
from sample import ensure_fixture

# init_db() is idempotent (CREATE TABLE IF NOT EXISTS), so this is safe whether
# the real dnainsight.db already exists or this is a genuinely fresh clone
# where it does not exist yet.
db_module.init_db()

# Seed the one persistent, tracked-fixture-derived sample BEFORE snapshotting,
# so this proves no pollution beyond that intentional, idempotent seed rather
# than flagging the seed itself as a violation on a clean uploads/.
ensure_fixture(ROOT)

before_profiles = len(db_module.list_profiles())
before_mtime = DB.stat().st_mtime if DB.exists() else None
before_uploads = len(list(UPLOADS.iterdir())) if UPLOADS.exists() else 0

print("BEFORE: real db profiles=%d  uploads/ files=%d  db mtime=%s"
      % (before_profiles, before_uploads, before_mtime))

HARNESSES = ["vdeps.py", "vreport.py", "vreport2.py", "vreports.py",
            "vserver.py", "vdbloss.py"]

for name in HARNESSES:
    p = subprocess.run([PY, str(ROOT / "tools" / name)], cwd=str(ROOT),
                       capture_output=True, text=True, timeout=120)
    tail = [ln for ln in (p.stdout + p.stderr).strip().splitlines() if ln.strip()][-3:]
    print("  ran %-16s exit=%d  %s" % (name, p.returncode, tail[-1][:90] if tail else ""))

after_profiles = len(db_module.list_profiles())
after_mtime = DB.stat().st_mtime if DB.exists() else None
after_uploads = len(list(UPLOADS.iterdir())) if UPLOADS.exists() else 0

print("\nAFTER:  real db profiles=%d  uploads/ files=%d  db mtime=%s"
      % (after_profiles, after_uploads, after_mtime))

ok = (before_profiles == after_profiles == 0
     and before_mtime == after_mtime
     and before_uploads == after_uploads)
print("\n%s" % ("ISOLATION HOLDS: the real database and uploads/ were not touched"
                if ok else "ISOLATION FAILED: the real database or uploads/ changed"))
sys.exit(0 if ok else 1)
