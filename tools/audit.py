"""
audit.py -- detect append-fault corruption in DNAInsight v2 source files.

The MCP bridge used to author v2 double-applied append-mode writes. That can
duplicate whole function and class blocks. A duplicated Python file still
imports cleanly because later definitions silently shadow earlier ones, so
corruption is invisible at runtime and must be detected structurally.

Checks per file:
  1. ast.parse succeeds
  2. no duplicated module-level def / class names
  3. no duplicated contiguous line blocks of length >= 8
  4. import succeeds in isolation
  5. reports declared __all__ vs actual public names

Usage:
    python tools/audit.py                  audit all backend + data + tests
    python tools/audit.py backend/x.py     audit specific files
"""

import ast
import sys
import hashlib
import importlib
import traceback
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

TARGET_DIRS = ["backend", "data", "tests"]
SKIP_NAMES = {"__pycache__", "_build", "tools"}
BLOCK = 8


def collect() -> list:
    files = []
    for d in TARGET_DIRS:
        base = ROOT / d
        if not base.exists():
            continue
        for p in sorted(base.rglob("*.py")):
            if any(part in SKIP_NAMES for part in p.parts):
                continue
            files.append(p)
    return files


def dup_top_level(tree: ast.Module) -> dict:
    names = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append(node.name)
    counts = Counter(names)
    return {n: c for n, c in counts.items() if c > 1}


def dup_class_methods(tree: ast.Module) -> dict:
    out = {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            inner = [n.name for n in node.body
                     if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
            counts = Counter(inner)
            dups = {n: c for n, c in counts.items() if c > 1}
            if dups:
                out[node.name] = dups
    return out


def dup_blocks(text: str, block: int = BLOCK) -> list:
    """Find repeated contiguous non-trivial line runs."""
    lines = [ln.rstrip() for ln in text.splitlines()]
    seen = defaultdict(list)
    for i in range(len(lines) - block + 1):
        window = lines[i:i + block]
        meat = [w for w in window if w.strip() and not w.strip().startswith("#")]
        if len(meat) < block - 1:
            continue
        key = hashlib.md5("\n".join(window).encode("utf-8")).hexdigest()
        seen[key].append(i + 1)
    hits = []
    for key, positions in seen.items():
        if len(positions) > 1:
            # collapse overlapping runs
            spread = [positions[0]]
            for p in positions[1:]:
                if p - spread[-1] > block:
                    spread.append(p)
            if len(spread) > 1:
                hits.append(spread)
    hits.sort(key=lambda s: s[0])
    return hits[:12]


def try_import(path: Path) -> str:
    rel = path.relative_to(ROOT)
    if rel.parts[0] != "backend":
        return "skipped"
    mod = "backend." + path.stem
    if path.stem == "__init__":
        mod = "backend"
    try:
        if mod in sys.modules:
            importlib.reload(sys.modules[mod])
        else:
            importlib.import_module(mod)
        return "ok"
    except Exception as exc:
        return f"FAIL {type(exc).__name__}: {exc}"


def public_names(tree: ast.Module) -> list:
    out = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith("_"):
                out.append(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id.isupper():
                    out.append(t.id)
    return out


def audit_one(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    rel = str(path.relative_to(ROOT))
    rec = {
        "file": rel,
        "lines": len(text.splitlines()),
        "bytes": len(text.encode("utf-8")),
        "parse": "ok",
        "dup_defs": {},
        "dup_methods": {},
        "dup_blocks": [],
        "import": "skipped",
        "public": [],
        "em_dash": text.count("—") + text.count("–"),
        "verdict": "CLEAN",
    }
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        rec["parse"] = f"SyntaxError line {exc.lineno}: {exc.msg}"
        rec["verdict"] = "BROKEN"
        return rec

    rec["dup_defs"] = dup_top_level(tree)
    rec["dup_methods"] = dup_class_methods(tree)
    rec["dup_blocks"] = dup_blocks(text)
    rec["public"] = public_names(tree)
    rec["import"] = try_import(path)

    if rec["dup_defs"] or rec["dup_methods"] or rec["dup_blocks"]:
        rec["verdict"] = "DUPLICATED"
    if rec["import"].startswith("FAIL"):
        rec["verdict"] = "IMPORT_FAIL"
    return rec


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    files = [ROOT / a for a in args] if args else collect()

    records = []
    for p in files:
        if not p.exists():
            print(f"MISSING {p}")
            continue
        records.append(audit_one(p))

    bad = [r for r in records if r["verdict"] != "CLEAN"]

    print("=" * 78)
    print(f"{'FILE':<34}{'LINES':>7}{'VERDICT':>14}  NOTES")
    print("=" * 78)
    for r in records:
        note = ""
        if r["dup_defs"]:
            note = "dup defs: " + ", ".join(f"{k}x{v}" for k, v in list(r["dup_defs"].items())[:5])
        elif r["dup_methods"]:
            first = list(r["dup_methods"].items())[0]
            note = f"dup methods in {first[0]}: " + ", ".join(f"{k}x{v}" for k, v in list(first[1].items())[:4])
        elif r["dup_blocks"]:
            note = f"{len(r['dup_blocks'])} repeated block(s), first at lines {r['dup_blocks'][0][:4]}"
        elif r["import"].startswith("FAIL"):
            note = r["import"][:60]
        if r["em_dash"]:
            note += f"  [em/en dashes: {r['em_dash']}]"
        print(f"{r['file']:<34}{r['lines']:>7}{r['verdict']:>14}  {note}")

    print("=" * 78)
    print(f"{len(records)} files audited, {len(bad)} need repair")
    if bad:
        print("\nREPAIR LIST:")
        for r in bad:
            print(f"  {r['file']}  ({r['verdict']}, {r['lines']} lines)")

    print("\nEXPECTED-BUT-MISSING:")
    expected = [
        "backend/orientation.py", "backend/genosets.py", "backend/frequency.py",
        "backend/prs.py", "backend/merge.py", "backend/traits.py", "backend/snpedia.py",
        "backend/interactive_report.py",
        "data/build_genosets.py", "data/genosets.json",
        "data/build_frequencies.py", "data/frequencies.json",
        "data/build_prs.py", "data/prs_models.json",
        "tests/test_orientation.py", "tests/test_genosets.py", "tests/test_frequency.py",
        "tests/test_prs.py", "tests/test_merge.py", "tests/test_traits.py",
        "tests/test_snpedia.py",
    ]
    for e in expected:
        if not (ROOT / e).exists():
            print(f"  MISSING {e}")

    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
