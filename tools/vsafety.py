"""
vsafety.py -- licence and payload safety, run before any push.

Three failures this catches, each of which would be hard to undo once pushed:

  1. Committing anything SNPedia-derived. SNPedia is CC-BY-NC-SA-3.0-US, so a
     harvested cache in the tree would relicense the whole repository and
     foreclose commercial use permanently.
  2. Committing personal genetic data. An uploads/ file or a .db is somebody's
     genome. Once pushed to a public remote it cannot be recalled.
  3. Committing a file GitHub will reject. The hard limit is 100 MB and it warns
     well before that.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
BIG_MB = 45

NC_MARKERS = ("snpedia_cache", "snpedia.db", "snpedia_export", "snpedia_harvest.db")
DNA_MARKERS = ("uploads/", "uploads\\")
DB_SUFFIXES = (".db", ".db-wal", ".db-shm", ".sqlite", ".sqlite3")

fails, notes = [], []


def git(*args):
    p = subprocess.run(["git"] + list(args), cwd=str(ROOT),
                       capture_output=True, text=True, errors="replace")
    return p.returncode, (p.stdout or "")


print("=" * 74)
print("LICENCE AND PAYLOAD SAFETY")
print("=" * 74)

code, out = git("ls-files", "--cached", "--others", "--exclude-standard")
if code != 0:
    print("  git unavailable, cannot verify. Treating as a BLOCKER.")
    sys.exit(1)

candidates = [c.strip() for c in out.splitlines() if c.strip()]
print(f"\n  {len(candidates)} path(s) are tracked or would be added by 'git add -A'")

# 1. SNPedia
nc = [c for c in candidates if any(m in c.lower() for m in NC_MARKERS)]
if nc:
    fails.append(f"SNPedia-derived file(s) would be committed: {nc}")
    print(f"  [FAIL] SNPedia-derived: {nc}")
else:
    print("  [ ok ] nothing SNPedia-derived would be committed")

# Confirm the real cache exists outside the repo, which is where it belongs.
home_cache = Path.home() / ".dnainsight"
print(f"  [info] harvest cache location: {home_cache} "
      f"({'exists' if home_cache.exists() else 'not created yet'}), "
      f"{'OUTSIDE' if str(ROOT).lower() not in str(home_cache).lower() else 'INSIDE'} the repo")

# 2. Personal data and databases
dna = [c for c in candidates
       if any(m in c.replace("\\", "/") for m in DNA_MARKERS)
       or c.lower().endswith(DB_SUFFIXES)]
if dna:
    fails.append(f"personal data or a database would be committed: {dna[:8]}")
    print(f"  [FAIL] personal data or database: {dna[:8]}")
else:
    print("  [ ok ] no uploads/ content and no database file would be committed")

# Prove the ignore rules actually fire, rather than trusting the file's text.
print("\n  ignore rules, verified with git check-ignore:")
for probe in ("uploads/whatever.txt", "dnainsight.db", "data/reference.db",
              "data/reference.db-wal", "_build/scratch.py",
              "backend/anything.bak", "snpedia_cache.db",
              "data/variant_summary.txt.gz"):
    rc, why = git("check-ignore", "-v", probe)
    if rc == 0:
        print(f"    [ ok ] {probe:<34} {why.strip().split(':')[0]}:{why.strip().split(':')[1]}")
    else:
        fails.append(f"{probe} is NOT ignored")
        print(f"    [FAIL] {probe:<34} NOT IGNORED")

# 3. Size
print("\n  file sizes:")
big = []
for rel in candidates:
    p = ROOT / rel
    try:
        size = p.stat().st_size
    except OSError:
        continue
    if size > BIG_MB * 1024 * 1024:
        big.append((rel, size))
if big:
    for rel, size in big:
        fails.append(f"{rel} is {size/1048576:.1f} MB")
        print(f"    [FAIL] {rel} {size/1048576:.1f} MB, above the {BIG_MB} MB threshold")
else:
    largest = sorted(
        ((ROOT / c).stat().st_size, c) for c in candidates
        if (ROOT / c).exists() and (ROOT / c).is_file())[-3:]
    print(f"    [ ok ] nothing above {BIG_MB} MB. Largest committable files:")
    for size, rel in reversed(largest):
        print(f"           {size/1024:>10,.0f} KB  {rel}")

# 4. The repo must actually claim MIT, since that is what the data licences allow.
lic = ROOT / "LICENSE"
if lic.exists() and "MIT" in lic.read_text(encoding="utf-8", errors="replace")[:400]:
    print("\n  [ ok ] LICENSE is MIT, consistent with CC0 and public domain data only")
else:
    notes.append("LICENSE does not obviously state MIT")
    print("\n  [warn] LICENSE does not obviously state MIT")

ds = ROOT / "data" / "DATA_SOURCES.md"
if ds.exists():
    text = ds.read_text(encoding="utf-8", errors="replace")
    required = ["CC0", "public domain", "CC-BY-NC-SA", "SNPedia", "CPIC", "ClinVar"]
    missing = [r for r in required if r not in text]
    if missing:
        notes.append(f"DATA_SOURCES.md does not mention: {missing}")
        print(f"  [warn] DATA_SOURCES.md does not mention: {missing}")
    else:
        print(f"  [ ok ] DATA_SOURCES.md present, {len(text.splitlines())} lines, "
              "covers every licence class")
else:
    fails.append("data/DATA_SOURCES.md is missing but the changelog references it")
    print("  [FAIL] data/DATA_SOURCES.md missing")

print("\n" + "=" * 74)
print("SAFETY OK" if not fails else f"SAFETY BLOCKERS: {len(fails)}")
for f in fails:
    print("  -", f)
sys.exit(0 if not fails else 1)
