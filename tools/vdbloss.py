"""
vdbloss.py -- prove or disprove that importing backend.database destroys data.

`backend/database.py` resolves DB_PATH at import time by calling _test_sqlite()
on the real database path, and _test_sqlite() ends with path.unlink(). If that
runs against an existing database, every profile and finding is gone.

This script does not reason about it. It writes real data, then imports the
module in a SEPARATE process, then checks whether the data survived.

ISOLATION: every subprocess is launched with DNAINSIGHT_DB_PATH set to a
throwaway file (see isolated_db.py), never the app's real dnainsight.db. This
test exercises the exact same _resolve_db_path()/_test_sqlite() code path
either way -- nothing in that logic special-cases "the app directory" -- so
redirecting it proves the same thing without any chance of writing a canary
profile into a real user's database. Before this, running this script left
"DURABILITY CANARY" rows in the live database, which surfaced in the app's own
profile list with no way to tell them apart from real data.

Run: python tools/vdbloss.py
Exit 0 means data survives. Exit 1 means data loss is real.
"""

import shutil
import subprocess
import sys
import sqlite3
import tempfile
import uuid
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
PY = sys.executable

TEST_DIR = Path(tempfile.gettempdir()) / "dnainsight_harness" / f"vdbloss_{uuid.uuid4().hex[:8]}"
TEST_DIR.mkdir(parents=True, exist_ok=True)
DB = TEST_DIR / "dnainsight.db"

ENV = {"DNAINSIGHT_DB_PATH": str(DB)}


def run(script: Path) -> subprocess.CompletedProcess:
    import os
    env = dict(os.environ)
    env.update(ENV)
    return subprocess.run([PY, str(script)], cwd=str(ROOT), capture_output=True,
                          text=True, env=env)


print("=" * 74)
print("DATABASE DURABILITY TEST")
print(f"isolated test database: {DB}")
print("=" * 74)

# Step 1. Create a profile through the real code path, against the ISOLATED path.
print("\n1. writing a profile through the application")
setup = ROOT / "tools" / "_dbloss_setup.py"
setup.write_text(
    "import sys\n"
    f"sys.path.insert(0, r'{ROOT}')\n"
    "from backend import database as db\n"
    "db.init_db()\n"
    "pid = db.create_profile('DURABILITY CANARY', '1990-01-01', 'other', 'test')\n"
    "db.upsert_finding(pid, None, {'rsid': 'rs_canary', 'gene': 'CANARY',\n"
    "    'genotype': 'AA', 'silo': 'informational'})\n"
    "print('created pid', pid)\n"
    "print('profiles now', len(db.list_profiles()))\n"
    "print('db path', db.DB_PATH)\n",
    encoding="utf-8")
p = run(setup)
print("   " + "\n   ".join(p.stdout.strip().splitlines()))
if p.returncode != 0:
    print("   setup failed:", p.stderr.strip()[:300])
    shutil.rmtree(TEST_DIR, ignore_errors=True)
    sys.exit(1)

size_before = DB.stat().st_size if DB.exists() else 0
print(f"   isolated db exists={DB.exists()} size={size_before:,} bytes")


def canary_count() -> int:
    if not DB.exists():
        return -1
    try:
        conn = sqlite3.connect(str(DB))
        n = conn.execute(
            "SELECT COUNT(*) FROM profiles WHERE name='DURABILITY CANARY'"
        ).fetchone()[0]
        conn.close()
        return n
    except sqlite3.Error as exc:
        print(f"   sqlite error reading canary: {exc}")
        return -1


before = canary_count()
print(f"   canary profiles on disk: {before}")
if before < 1:
    print("   could not establish a baseline, aborting")
    shutil.rmtree(TEST_DIR, ignore_errors=True)
    sys.exit(1)

# Step 2. THE TEST. Import backend.database in a fresh process, isolated path,
# and do nothing else. If DB_PATH resolution deletes the file, the canary
# vanishes.
print("\n2. importing backend.database in a separate process, doing nothing else")
probe = ROOT / "tools" / "_dbloss_probe.py"
probe.write_text(
    "import sys\n"
    f"sys.path.insert(0, r'{ROOT}')\n"
    "import backend.database as db\n"
    "print('imported, DB_PATH =', db.DB_PATH)\n",
    encoding="utf-8")
p = run(probe)
print("   " + "\n   ".join(p.stdout.strip().splitlines()))

after = canary_count()
size_after = DB.stat().st_size if DB.exists() else 0
print(f"\n   isolated db exists={DB.exists()} size={size_after:,} bytes")
print(f"   canary profiles after a bare import: {after}")

# Step 3. And through create_app(), which is what actually happens on launch.
print("\n3. calling create_app(), which is what a real launch does")
probe2 = ROOT / "tools" / "_dbloss_probe2.py"
probe2.write_text(
    "import sys\n"
    f"sys.path.insert(0, r'{ROOT}')\n"
    "import app as a\n"
    "a.create_app()\n"
    "from backend import database as db\n"
    "print('profiles after create_app:', len(db.list_profiles()))\n",
    encoding="utf-8")
p = run(probe2)
print("   " + "\n   ".join((p.stdout + p.stderr).strip().splitlines()[:8]))
after_app = canary_count()
print(f"   canary profiles after create_app(): {after_app}")

# Cleanup: the probe scripts AND the entire isolated test directory. Nothing
# this script writes should outlive it.
for f in (setup, probe, probe2):
    f.unlink(missing_ok=True)
shutil.rmtree(TEST_DIR, ignore_errors=True)
print(f"\n   isolated test directory removed: {not TEST_DIR.exists()}")

print("\n" + "=" * 74)
lost_on_import = before >= 1 and after < 1
lost_on_launch = before >= 1 and after_app < 1
if lost_on_import or lost_on_launch:
    print("DATA LOSS CONFIRMED")
    if lost_on_import:
        print(f"  importing backend.database destroyed the data ({before} -> {after})")
    if lost_on_launch:
        print(f"  create_app() destroyed the data ({before} -> {after_app})")
    print("\n  Every launch of DNAInsight would wipe every stored profile,")
    print("  finding and report. This is a data-destroying defect.")
    sys.exit(1)
print("DATA SURVIVES")
print(f"  canary intact through a bare import and through create_app() ({after_app})")
sys.exit(0)
