"""Prove the new rate-limit tests are not vacuous.

Copies the tree to a temp dir, reverts the limiter threading to exactly the code
that shipped before this fix, and runs tests/test_snpedia_ratelimit.py there. A
test that cannot fail against the bug is not a regression guard.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable

REVERT = [
    ("def fetch_subject(subject: str, session: Any = None,\n"
     "                  limiter: _RateLimiter | None = None) -> dict:\n",
     "def fetch_subject(subject: str, session: Any = None) -> dict:\n"),
    ("def fetch_wikitext(title: str, session: Any = None,\n"
     "                   limiter: _RateLimiter | None = None) -> str:\n",
     "def fetch_wikitext(title: str, session: Any = None) -> str:\n"),
    ("        session=session,\n        limiter=limiter,\n    )\n",
     "        session=session,\n    )\n"),
    ("fetch_subject(title, session=session, limiter=limiter)",
     "fetch_subject(title, session=session)"),
    ('fetch_subject(f"{title}{token}", session=session, limiter=limiter)',
     'fetch_subject(f"{title}{token}", session=session)'),
    ("fetch_wikitext(title, session=session, limiter=limiter)",
     "fetch_wikitext(title, session=session)"),
    ("fetch_subject(name, session=session, limiter=limiter)",
     "fetch_subject(name, session=session)"),
    ('fetch_wikitext(f"{name}/criteria", session=session, limiter=limiter)',
     'fetch_wikitext(f"{name}/criteria", session=session)'),
]

tmp = Path(tempfile.mkdtemp(prefix="dnai_mut_"))
for sub in ("backend", "tests", "data"):
    shutil.copytree(ROOT / sub, tmp / sub,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
for f in ("app.py", "requirements.txt", ".flake8"):
    if (ROOT / f).exists():
        shutil.copy2(ROOT / f, tmp / f)

target = tmp / "backend" / "snpedia.py"
text = target.read_text(encoding="utf-8")
reverted = 0
for new, old in REVERT:
    n = text.count(new)
    if n:
        text = text.replace(new, old)
        reverted += n
target.write_text(text, encoding="utf-8", newline="\n")
print("reverted %d occurrences of the limiter threading in the temp copy" % reverted)
if reverted == 0:
    print("ABORT: nothing to revert, the mutation would be a no-op")
    shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(1)

p = subprocess.run([PY, "-m", "pytest", "tests/test_snpedia_ratelimit.py",
                    "-q", "--no-header", "-p", "no:cacheprovider"],
                   cwd=str(tmp), capture_output=True, text=True)
out = (p.stdout + p.stderr).strip().splitlines()
print("\npytest against the reverted code, exit %d:" % p.returncode)
for ln in out[-14:]:
    print("  " + ln[:150])

shutil.rmtree(tmp, ignore_errors=True)
print()
if p.returncode == 0:
    print("VACUOUS: the tests pass even with the bug restored. They guard nothing.")
    sys.exit(1)
print("NOT VACUOUS: restoring the bug fails the new tests, as a guard must.")
