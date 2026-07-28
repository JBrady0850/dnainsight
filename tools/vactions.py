"""Fail if any workflow pins a GitHub Action major version that bundles Node 20.

GitHub deprecated the Node 20 action runtime on 2025-09-19 and now force-runs
those actions on Node 24 while printing a warning on every job. The fix is to
pin a major version that declares Node 24 natively.

Sources for the version boundaries, checked 2026-07-27:
  actions/checkout      v5.0.0 "Update actions checkout to use node 24"
  actions/setup-python  v6.0.0 upgraded to Node 24
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WF = ROOT / ".github" / "workflows"

# action name -> first major version that declares Node 24
MIN_MAJOR = {
    "actions/checkout": 5,
    "actions/setup-python": 6,
    "actions/setup-node": 5,
    "actions/cache": 4,
    "actions/upload-artifact": 5,
    "actions/download-artifact": 5,
}

USES = re.compile(r"uses:\s*([A-Za-z0-9_.\-]+/[A-Za-z0-9_.\-]+)@v(\d+)")

print("GITHUB ACTIONS RUNTIME CHECK")
if not WF.exists():
    print("  no .github/workflows directory, nothing to check")
    sys.exit(0)

files = sorted(list(WF.glob("*.yml")) + list(WF.glob("*.yaml")))
if not files:
    print("  no workflow files found")
    sys.exit(0)

stale = []
seen = 0
for f in files:
    for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
        m = USES.search(line)
        if not m:
            continue
        action, major = m.group(1), int(m.group(2))
        seen += 1
        floor = MIN_MAJOR.get(action)
        if floor is None:
            print("  [info] %s:%d  %s@v%d  (no Node floor recorded)"
                  % (f.name, i, action, major))
        elif major < floor:
            print("  [FAIL] %s:%d  %s@v%d bundles Node 20, needs v%d or later"
                  % (f.name, i, action, major, floor))
            stale.append("%s:%d %s@v%d" % (f.name, i, action, major))
        else:
            print("  [ ok ] %s:%d  %s@v%d  (Node 24)"
                  % (f.name, i, action, major))

print("  %d action pin(s) across %d workflow file(s)" % (seen, len(files)))
if stale:
    print("  NODE 20 ACTION PINS: %d" % len(stale))
    sys.exit(1)
print("  ACTIONS RUNTIME OK, no Node 20 pins")
