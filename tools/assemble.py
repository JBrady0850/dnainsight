"""
assemble.py -- build-time helper for DNAInsight v2 development.

The MCP file bridge used to author this release double-applies append-mode
writes, so large source files are authored as numbered part files and then
concatenated here. Parts are written with rewrite mode only, which is
idempotent under that fault.

Usage:
    python tools/assemble.py <target_relative_path> [--keep]

Looks for tools/<basename>.partNN files, concatenates them in sorted order,
writes the target as UTF-8 with LF line endings, removes the parts unless
--keep is given, and prints the resulting line count.
"""

import sys
import glob
from pathlib import Path

BUILD_DIR = Path(__file__).parent.resolve()
ROOT = BUILD_DIR.parent


def assemble(target_rel: str, keep: bool = False) -> int:
    target = (ROOT / target_rel).resolve()
    stem = target.name
    pattern = str(BUILD_DIR / f"{stem}.part*")
    parts = sorted(glob.glob(pattern))
    if not parts:
        print(f"ERROR: no parts found matching {pattern}")
        return 1

    chunks = []
    for p in parts:
        text = Path(p).read_text(encoding="utf-8")
        chunks.append(text)
        print(f"  + {Path(p).name}  ({len(text.splitlines())} lines)")

    body = "".join(chunks)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(body)

    if not keep:
        for p in parts:
            Path(p).unlink(missing_ok=True)

    lines = len(body.splitlines())
    print(f"WROTE {target}  ({lines} lines, {len(body)} bytes, {len(parts)} parts)")
    return 0


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    keep = "--keep" in sys.argv[1:]
    if not args:
        print(__doc__)
        return 2
    rc = 0
    for target_rel in args:
        rc |= assemble(target_rel, keep=keep)
    return rc


if __name__ == "__main__":
    sys.exit(main())
