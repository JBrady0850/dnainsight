"""vversion.py -- version declarations, checked within their own group.

WHY THIS IS TWO GROUPS AND NOT ONE
----------------------------------
This script used to compare every version string in the repository against
every other one and report a mismatch whenever they differed. They always
differ, because two unrelated things declare versions here:

  APPLICATION   backend/__init__.py, CHANGELOG.md, the README heading.
                What release of DNAInsight this is.

  DATA ARTEFACT REFERENCE_VERSION and the `_meta.version` field inside each
                built JSON. What build of the bundled reference data this is.

Those are deliberately independent. `backend/provenance.py` says so directly:
an artefact semver says nothing about which ClinVar release went into it, and
bumping the app for an installer fix must not imply the reference data was
rebuilt. Forcing them to agree would make the artefact version a lie.

So the old check could never pass, and a warning that is always on is a warning
nobody reads. It fired on every run of `tools/golive.py` from v2.0 onward and
was carried as a permanent "safe to ship" note, which is exactly how a real
version skew would have been missed.

Each group is now checked against itself. Cross-group differences are expected
and are reported as information, not as a failure.
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()

APPLICATION: dict = {}
ARTEFACT: dict = {}


def grab(into: dict, label: str, rel: str, pattern: str) -> None:
    p = ROOT / rel
    if not p.exists():
        return
    m = re.search(pattern, p.read_text(encoding="utf-8", errors="replace"), re.M)
    if m:
        into[label] = m.group(1)


def grab_json(into: dict, label: str, rel: str) -> None:
    p = ROOT / rel
    if not p.exists():
        return
    try:
        meta = json.loads(p.read_text(encoding="utf-8")).get("_meta", {})
    except ValueError:
        return
    if meta.get("version"):
        into[label] = meta["version"]


def collect() -> tuple[dict, dict]:
    """Fill both groups. Importable so the test suite checks the real thing."""
    APPLICATION.clear()
    ARTEFACT.clear()

    grab(APPLICATION, "backend/__init__.py", "backend/__init__.py",
         r'__version__\s*=\s*"([^"]+)"')
    grab(APPLICATION, "CHANGELOG.md newest", "CHANGELOG.md",
         r"^##\s*\[([0-9][^\]]*)\]")
    grab(APPLICATION, "README.md heading", "README.md",
         r"^DNAInsight v([0-9][0-9A-Za-z.\-]*)")

    grab(ARTEFACT, "data/build_reference.py", "data/build_reference.py",
         r'REFERENCE_VERSION\s*=\s*"([^"]+)"')
    grab_json(ARTEFACT, "snp_reference.json", "data/snp_reference.json")
    grab_json(ARTEFACT, "genosets.json", "data/genosets.json")
    grab_json(ARTEFACT, "prs_models.json", "data/prs_models.json")
    grab_json(ARTEFACT, "frequencies.json", "data/frequencies.json")

    return APPLICATION, ARTEFACT


def report(name: str, group: dict) -> bool:
    """Print one group and return True when it is internally consistent."""
    print(f"  {name}")
    for where, value in group.items():
        print(f"    {where:<30} {value}")
    if not group:
        print("    (nothing declared)")
        return True
    distinct = sorted(set(group.values()))
    if len(distinct) == 1:
        print(f"    ok: all {len(group)} agree on {distinct[0]}")
        return True
    print(f"    MISMATCH across {len(group)} declarations: {distinct}")
    return False


def main() -> int:
    print("=" * 74)
    print("VERSION CONSISTENCY")
    print("=" * 74)

    application, artefact = collect()
    app_ok = report("APPLICATION", application)
    print()
    art_ok = report("DATA ARTEFACT", artefact)
    print()

    app_values = sorted(set(application.values()))
    art_values = sorted(set(artefact.values()))
    if app_ok and art_ok and app_values != art_values:
        print("  info: application and data artefact versions differ, which is "
              "expected.")
        print("        They track different things and are not required to "
              "match.")

    if app_ok and art_ok:
        print("  ok: every version group is internally consistent")
        return 0

    print("  Not a blocker for a commit, but fix it before tagging a release.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
