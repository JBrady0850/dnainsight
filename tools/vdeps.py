"""
vdeps.py -- verify requirements.txt is satisfiable by what is actually installed.

Does NOT install anything. Parses the requirement specifiers and compares them
against the installed distributions, which is what actually matters: the point of
the pin change was to stop the file claiming a version nobody had verified.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()

try:
    from importlib.metadata import version as installed_version, PackageNotFoundError
except ImportError:
    print("importlib.metadata unavailable")
    sys.exit(1)


def parse(line: str):
    """Return (name, [(op, version), ...]) or None for a non-requirement line."""
    line = line.split("#")[0].strip()
    if not line or line.startswith("-"):
        return None
    m = re.match(r"^([A-Za-z0-9_.\-]+)\s*(.*)$", line)
    if not m:
        return None
    name, rest = m.group(1), m.group(2)
    specs = re.findall(r"(>=|<=|==|!=|~=|<|>)\s*([0-9][0-9A-Za-z.\-]*)", rest)
    return name, specs


def tup(v: str):
    """Version to comparable tuple, ignoring non-numeric suffixes."""
    parts = []
    for p in re.split(r"[.\-]", v):
        m = re.match(r"^(\d+)", p)
        parts.append(int(m.group(1)) if m else 0)
    return tuple(parts + [0] * (4 - len(parts)))[:4]


def satisfies(have: str, op: str, want: str) -> bool:
    a, b = tup(have), tup(want)
    return {">=": a >= b, "<=": a <= b, "==": a == b, "!=": a != b,
            "<": a < b, ">": a > b, "~=": a >= b and a[0] == b[0]}[op]


fails = []
print("=" * 70)
print("DEPENDENCY VERIFICATION")
print("=" * 70)

for rel in ("requirements.txt", "requirements-dev.txt"):
    path = ROOT / rel
    if not path.exists():
        print(f"\n{rel}: MISSING")
        continue
    print(f"\n{rel}")
    for raw in path.read_text(encoding="utf-8").splitlines():
        parsed = parse(raw)
        if not parsed:
            continue
        name, specs = parsed
        try:
            have = installed_version(name)
        except PackageNotFoundError:
            print(f"  [FAIL] {name:<12} NOT INSTALLED (required {specs})")
            fails.append(f"{name} not installed")
            continue
        bad = [(op, want) for op, want in specs if not satisfies(have, op, want)]
        spec_text = ",".join(f"{op}{want}" for op, want in specs) or "any"
        if bad:
            print(f"  [FAIL] {name:<12} installed {have}, violates {bad}")
            fails.append(f"{name} {have} violates {bad}")
        else:
            print(f"  [ ok ] {name:<12} installed {have}, satisfies {spec_text}")

# The runtime must import with only what requirements.txt asks for.
print("\nruntime import check")
for mod in ("flask", "requests"):
    try:
        __import__(mod)
        print(f"  [ ok ] import {mod}")
    except Exception as exc:
        print(f"  [FAIL] import {mod}: {exc}")
        fails.append(f"import {mod}")

# The app itself must construct.
print("\napplication construction")
sys.path.insert(0, str(ROOT))
try:
    import isolated_db
    isolated_db.use_temp_db()
    import app as app_module
    a = app_module.create_app()
    rules = {str(r.rule) for r in a.url_map.iter_rules()}
    print(f"  [ ok ] create_app() built {len(rules)} routes")
    for needed in ("/api/status", "/api/capabilities",
                   "/api/profiles/<int:pid>/findings/v2",
                   "/api/profiles/<int:pid>/reports/interactive"):
        if needed in rules:
            print(f"  [ ok ] route present: {needed}")
        else:
            print(f"  [FAIL] route MISSING: {needed}")
            fails.append(f"route {needed}")
except Exception as exc:
    print(f"  [FAIL] create_app(): {type(exc).__name__}: {exc}")
    fails.append("create_app")

print("\n" + "=" * 70)
print("DEPENDENCIES OK" if not fails else f"PROBLEMS: {len(fails)}")
for f in fails:
    print("  -", f)
sys.exit(0 if not fails else 1)
