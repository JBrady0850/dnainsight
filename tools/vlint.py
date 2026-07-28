"""Run the CI lint job's exact flake8 invocation and fail on any finding.

The go-live gate ran the test suite but never ran flake8, so a lint failure that
would have turned CI red was invisible locally. This closes that hole.
"""
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGETS = ["backend/", "data/", "app.py"]

print("LINT (the CI lint job's own command)")
print("  flake8 " + " ".join(TARGETS))

p = subprocess.run([sys.executable, "-m", "flake8"] + TARGETS,
                   cwd=str(ROOT), capture_output=True, text=True)
lines = [ln for ln in (p.stdout + p.stderr).splitlines() if ln.strip()]

if p.returncode == 127 or "No module named flake8" in "".join(lines):
    print("  [FAIL] flake8 is not installed, so the CI lint job cannot be verified")
    sys.exit(1)

codes = Counter()
byfile = defaultdict(list)
for ln in lines:
    parts = ln.split(":", 3)
    if len(parts) < 4:
        continue
    f, row, msg = parts[0], parts[1], parts[3].strip()
    codes[msg.split(" ", 1)[0]] += 1
    byfile[f].append((int(row), msg))

if not lines:
    print("  [ ok ] 0 findings")
    print("  LINT OK")
    sys.exit(0)

print("  [FAIL] %d finding(s), exit %d" % (len(lines), p.returncode))
for code, n in codes.most_common():
    print("      %-6s x%d" % (code, n))
for f in sorted(byfile):
    print("      %s (%d)" % (f, len(byfile[f])))
    for row, msg in sorted(byfile[f])[:8]:
        print("          %5d  %s" % (row, msg))
print("  LINT WOULD FAIL CI")
sys.exit(1)
