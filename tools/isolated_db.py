"""Point backend.database / routes / routes_v2 at throwaway files and folders
BEFORE any harness imports them.

ROOT CAUSE THIS FIXES
backend.database._resolve_db_path()'s first candidate is
<app dir>/dnainsight.db, and backend.routes/routes_v2's UPLOAD_DIR and
REPORTS_DIR are hardcoded to <app dir>/uploads and <app dir>/reports_output --
the exact files and folders the real, installed app uses. Any tools script
that imports backend.database or app (which imports both route modules)
resolves to those SAME live locations unless something redirects them first.
Six harnesses did this: vdeps.py, vreport.py, vreport2.py, vreports.py,
vserver.py, vdbloss.py. Running them left "DURABILITY CANARY" /
"Report Verify" / "Collision Test" / "Static Report Verify" / "V2 Verify"
profiles in the real database and 99+ files in the real uploads/ directory,
indistinguishable from data a real user created, because neither
list_profiles() nor the uploads folder carries any test/real marker.

THE FIX
- database._resolve_db_path() checks DNAINSIGHT_DB_PATH first.
- routes.py and routes_v2.py check DNAINSIGHT_UPLOAD_DIR / DNAINSIGHT_REPORTS_DIR
  first.
All three are resolved once, at import time, so they must be set BEFORE the
harness's own `import app` / `from backend import database` line.

USAGE
    import isolated_db; isolated_db.use_temp_db()
    import app as app_module          # only now, so it sees the redirected paths
"""
import os
import tempfile
import uuid


def use_temp_db() -> str:
    """Redirect DB_PATH, UPLOAD_DIR and REPORTS_DIR to a fresh directory this
    process owns, and return the database path.

    A uuid-suffixed name, not just the pid, because pids get reused and two
    harnesses launched close together must never collide.
    """
    root = os.path.join(tempfile.gettempdir(), "dnainsight_harness",
                        f"verify_{os.getpid()}_{uuid.uuid4().hex[:8]}")
    db_path = os.path.join(root, "dnainsight.db")
    upload_dir = os.path.join(root, "uploads")
    reports_dir = os.path.join(root, "reports_output")
    os.makedirs(upload_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.environ["DNAINSIGHT_DB_PATH"] = db_path
    os.environ["DNAINSIGHT_UPLOAD_DIR"] = upload_dir
    os.environ["DNAINSIGHT_REPORTS_DIR"] = reports_dir
    return db_path
