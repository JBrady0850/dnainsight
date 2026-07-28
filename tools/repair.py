"""
repair.py -- undo append-fault duplication in DNAInsight v2 source files.

FAULT MODEL
The MCP file bridge applied every append-mode write twice. So a file authored
as rewrite(A), append(B), append(C) landed on disk as:

    A + B + B + C + C

Duplicated regions are therefore always ADJACENT and exactly equal. This is
recoverable losslessly: find every maximal adjacent repeated run and drop the
second copy.

ALGORITHM
Walk the line list left to right. At each index i, find the largest k such that
lines[i:i+k] == lines[i+k:i+2k] and k >= MIN_RUN. Prefer the largest k, because
a doubled 200-line chunk also technically contains many smaller doubled runs,
and collapsing the largest first is what actually reverses the fault. Delete the
second copy and re-test at the same index so nested duplication is caught.

SAFETY
  - Writes a .bak of every file it touches before modifying.
  - Refuses to act if the result fails ast.parse.
  - Refuses to act if the result loses any unique top-level def or class name.
  - Reports every collapse with line numbers so the change is reviewable.
  - --dry-run prints the plan and writes nothing.

Usage:
    python tools/repair.py --dry-run
    python tools/repair.py
    python tools/repair.py backend/traits.py
"""

import ast
import sys
import shutil
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).parent.parent.resolve()

MIN_RUN = 4          # do not collapse runs shorter than this
MIN_MEAT = 2         # a run must contain at least this many code-bearing lines

DEFAULT_TARGETS = [
    "backend/frequency.py",
    "backend/merge.py",
    "backend/prs.py",
    "backend/snpedia.py",
    "backend/traits.py",
    "data/build_frequencies.py",
    "data/build_genosets.py",
]


def top_level_names(text: str) -> Counter:
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return Counter()
    names = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append(node.name)
    return Counter(names)


def is_meaty(run: list) -> bool:
    meat = [ln for ln in run if ln.strip() and not ln.strip().startswith("#")]
    return len(meat) >= MIN_MEAT


def collapse(lines: list) -> tuple:
    """Return (new_lines, collapse_log)."""
    out = list(lines)
    log = []
    i = 0
    guard = 0
    while i < len(out):
        guard += 1
        if guard > 500000:
            log.append(("ABORT", i, 0, "guard tripped"))
            break
        n = len(out)
        best_k = 0
        # largest possible duplicate run starting at i
        max_k = (n - i) // 2
        k = max_k
        while k >= MIN_RUN:
            if out[i:i + k] == out[i + k:i + 2 * k] and is_meaty(out[i:i + k]):
                best_k = k
                break
            k -= 1
        if best_k:
            log.append(("COLLAPSE", i + 1, best_k, out[i].strip()[:60]))
            del out[i + best_k:i + 2 * best_k]
            continue          # retest same index for nested duplication
        i += 1
    return out, log


def repair_file(rel: str, dry_run: bool = False) -> dict:
    path = ROOT / rel
    rec = {"file": rel, "status": "skipped", "before": 0, "after": 0, "collapses": 0, "note": ""}
    if not path.exists():
        rec["note"] = "missing"
        return rec

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    rec["before"] = len(lines)

    before_names = top_level_names(text)
    if not before_names:
        rec["status"] = "PARSE_FAIL_BEFORE"
        rec["note"] = "file does not parse; not touching it"
        return rec

    new_lines, log = collapse(lines)
    rec["collapses"] = sum(1 for e in log if e[0] == "COLLAPSE")
    rec["after"] = len(new_lines)

    if rec["after"] == rec["before"]:
        rec["status"] = "already clean"
        return rec

    new_text = "\n".join(new_lines)
    if not new_text.endswith("\n"):
        new_text += "\n"

    try:
        ast.parse(new_text)
    except SyntaxError as exc:
        rec["status"] = "REJECTED"
        rec["note"] = f"result would not parse: line {exc.lineno} {exc.msg}"
        return rec

    after_names = top_level_names(new_text)
    lost = set(before_names) - set(after_names)
    if lost:
        rec["status"] = "REJECTED"
        rec["note"] = "would lose definitions: " + ", ".join(sorted(lost))
        return rec

    still_dup = {k: v for k, v in after_names.items() if v > 1}
    rec["note"] = ("residual duplicate defs: " + ", ".join(f"{k}x{v}" for k, v in still_dup.items())) if still_dup else "unique defs restored"

    if dry_run:
        rec["status"] = "would repair"
        for kind, start, k, preview in log[:20]:
            print(f"    {kind} at line {start}: {k} lines  | {preview}")
        return rec

    shutil.copy2(path, path.with_suffix(path.suffix + ".bak"))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(new_text)
    rec["status"] = "repaired"
    return rec


def main() -> int:
    dry = "--dry-run" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    targets = args or DEFAULT_TARGETS

    print("=" * 78)
    print("APPEND-FAULT REPAIR" + ("  (DRY RUN)" if dry else ""))
    print("=" * 78)

    records = []
    for rel in targets:
        print(f"\n{rel}")
        rec = repair_file(rel.replace("\\", "/"), dry_run=dry)
        records.append(rec)
        print(f"  {rec['status']}: {rec['before']} -> {rec['after']} lines, "
              f"{rec['collapses']} collapse(s). {rec['note']}")

    print("\n" + "=" * 78)
    ok = [r for r in records if r["status"] in ("repaired", "would repair", "already clean")]
    print(f"{len(ok)}/{len(records)} handled")
    rejected = [r for r in records if r["status"] == "REJECTED"]
    if rejected:
        print("REJECTED (needs manual rewrite):")
        for r in rejected:
            print(f"  {r['file']}: {r['note']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
