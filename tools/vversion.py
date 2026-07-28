"""vversion.py -- one version string across every place that declares one."""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
found: dict = {}


def grab(label: str, rel: str, pattern: str) -> None:
    p = ROOT / rel
    if not p.exists():
        return
    m = re.search(pattern, p.read_text(encoding="utf-8", errors="replace"), re.M)
    if m:
        found[label] = m.group(1)


def grab_json(label: str, rel: str) -> None:
    p = ROOT / rel
    if not p.exists():
        return
    try:
        meta = json.loads(p.read_text(encoding="utf-8")).get("_meta", {})
    except ValueError:
        return
    if meta.get("version"):
        found[label] = meta["version"]


print("=" * 74)
print("VERSION CONSISTENCY")
print("=" * 74)

grab("backend/__init__.py", "backend/__init__.py", r'__version__\s*=\s*"([^"]+)"')
grab("data/build_reference.py", "data/build_reference.py",
     r'REFERENCE_VERSION\s*=\s*"([^"]+)"')
grab("CHANGELOG.md newest", "CHANGELOG.md", r"^##\s*\[([0-9][^\]]*)\]")
grab("README.md heading", "README.md", r"^DNAInsight v([0-9][0-9A-Za-z.\-]*)")
grab_json("snp_reference.json", "data/snp_reference.json")
grab_json("genosets.json", "data/genosets.json")
grab_json("prs_models.json", "data/prs_models.json")
grab_json("frequencies.json", "data/frequencies.json")

for where, value in found.items():
    print(f"  {where:<30} {value}")

distinct = sorted(set(found.values()))
print()
if len(distinct) == 1:
    print(f"  ok: all {len(found)} declarations agree on {distinct[0]}")
    sys.exit(0)
print(f"  MISMATCH across {len(found)} declarations: {distinct}")
print("  Not a blocker for a commit, but pick one before tagging a release.")
sys.exit(1)
