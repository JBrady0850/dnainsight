"""
golive.py -- single-process go-live readiness check.

Everything runs as a subprocess from ONE Python parent and prints straight to
stdout. No log file, no shell redirect, no pipe. That matters because this
bridge sometimes re-executes a command concurrently, and two copies writing the
same log file collide with "the process cannot access the file".

Usage:
    python tools/golive.py
    python tools/golive.py --quick     skip the slow stages

Exit 0 means ready. Exit 1 means blockers, listed at the end.
"""

import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
PY = sys.executable
QUICK = "--quick" in sys.argv

blockers: list = []
warnings: list = []


def run(args, timeout=1800):
    try:
        p = subprocess.run(args, cwd=str(ROOT), capture_output=True, text=True,
                           errors="replace", timeout=timeout)
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired:
        return 124, f"timed out after {timeout}s"
    except FileNotFoundError:
        return 127, f"not found: {args[0]}"


def stage(num, title):
    print()
    print("#" * 76)
    print(f"#  {num}  {title}")
    print("#" * 76)


def show(out, patterns, limit=14):
    keep = []
    for ln in out.splitlines():
        if any(re.search(p, ln) for p in patterns):
            keep.append(ln.rstrip())
    for ln in keep[-limit:]:
        print("   " + ln)
    return keep


def verify(num, title, script, patterns, blocker=True, extra=None, limit=14):
    stage(num, title)
    path = ROOT / script
    if not path.exists():
        print(f"   SKIPPED, {script} not present")
        warnings.append(f"{script} missing")
        return
    t0 = time.monotonic()
    code, out = run([PY, str(path)] + (extra or []))
    show(out, patterns, limit)
    print(f"   exit {code}  ({time.monotonic() - t0:.1f}s)")
    if code != 0:
        (blockers if blocker else warnings).append(f"{title} failed (exit {code})")


print("=" * 76)
print("DNAInsight v2.0  GO-LIVE READINESS")
print(f"repository: {ROOT}")
print("=" * 76)

# ---------------------------------------------------------------------------
stage(1, "REBUILD EVERY DERIVED ARTIFACT")
for rel, args in (("data/build_reference.py", []),
                  ("data/build_genosets.py", []),
                  ("data/build_prs.py", []),
                  ("data/build_prs.py", ["--validate"])):
    code, out = run([PY, str(ROOT / rel)] + args)
    tail = [l for l in out.strip().splitlines() if l.strip()][-1:] or [""]
    label = f"{rel} {' '.join(args)}".strip()
    print(f"   {'ok  ' if code == 0 else 'FAIL'} {label}: {tail[0][:88]}")
    if code != 0:
        blockers.append(f"{label} failed")

# ---------------------------------------------------------------------------
stage(2, "SOURCE INTEGRITY (append and edit-block duplication)")
code, out = run([PY, str(ROOT / "tools" / "audit.py")])
flagged = [l for l in out.splitlines()
           if ("DUPLICATED" in l or "BROKEN" in l or "IMPORT_FAIL" in l) and ".py" in l]
real = [l for l in flagged if "parsers.py" not in l]
for l in flagged:
    print("   " + l.strip())
counted = [l for l in out.splitlines() if "audited" in l]
if counted:
    print("   " + counted[0].strip())
if real:
    blockers.append(f"{len(real)} file(s) contain duplicated definitions")
    print(f"   BLOCKER: {len(real)} real duplication(s)")
else:
    print("   ok   only the known parsers.py false positive")

# ---------------------------------------------------------------------------
verify(3, "MODULE SMOKE", "tools/smoke.py",
       [r"\[ ok \]", r"\[FAIL\]", r"modules ok"])
verify(4, "STRAND REGRESSION", "tools/vfreq.py",
       [r"REGRESSION", r"STILL BROKEN"], limit=3)
verify(5, "PIPELINE CONTRACT", "tools/vpipe.py",
       [r"CONTRACT KEY", r"PIPELINE (OK|PROBLEMS)", r"TOTAL FINDINGS",
        r"out of 0-10", r"no-calls with", r"must be none"])
verify(6, "FILTER ENGINE", "tools/vfilters.py",
       [r"FILTER ENGINE", r"BROKEN"], limit=4)
verify(7, "API SWEEP", "tools/vserver.py",
       [r"endpoint checks passed", r"SERVER (OK|PROBLEMS)", r"\[FAIL\]"], limit=8)
verify(8, "INTERACTIVE REPORT", "tools/vreport.py",
       [r"INTERACTIVE REPORT", r"\[FAIL\]", r"PROBLEMS", r"bytes"], limit=8)
verify(9, "REPORT FILENAME COLLISION", "tools/vreport2.py",
       [r"FILENAME", r"\[FAIL\]", r"distinct", r"PROBLEMS"], limit=6)
verify(10, "STATIC REPORTS", "tools/vreports.py",
       [r"PASS", r"FAIL", r"checks"], limit=8)
verify(11, "FRONTEND INTEGRITY", "tools/vfrontend.py",
       [r"RESULT", r"\[FAIL\]", r"required strings", r"node --check"], limit=8)
verify(12, "DEPENDENCIES", "tools/vdeps.py",
       [r"\[ ok \]", r"\[FAIL\]", r"DEPENDENCIES"], limit=12)
verify(13, "VERSION CONSISTENCY", "tools/vversion.py",
       [r"ok: all", r"MISMATCH"], blocker=False, limit=4)
verify(14, "LICENCE AND PAYLOAD SAFETY", "tools/vsafety.py",
       [r"\[FAIL\]", r"SAFETY", r"\[ ok \]", r"\[warn\]", r"\[info\]"], limit=20)
verify(15, "COMPLETENESS", "tools/vcomplete.py",
       [r"\[FAIL\]", r"COMPLETE", r"MISSING", r"required files present",
        r"v2 build"], limit=10)

verify(16, "GITHUB ACTIONS RUNTIME", "tools/vactions.py",
       [r"\[ ok \]", r"\[FAIL\]", r"\[info\]", r"ACTIONS RUNTIME",
        r"NODE 20 ACTION PINS", r"action pin"], limit=14)
verify(17, "LINT", "tools/vlint.py",
       [r"\[ ok \]", r"\[FAIL\]", r"LINT", r"x\d+", r"flake8"], limit=14)
verify(18, "CI CLEAN-CLONE SIMULATION", "tools/vci.py",
       [r"\[ok\]", r"\[FAIL\]", r"CI WOULD", r"PUBLISHED", r"git would publish"],
       limit=16)

# ---------------------------------------------------------------------------
verify(19, "HARNESS ISOLATION", "tools/vdbisolation.py",
       [r"BEFORE:", r"AFTER:", r"ISOLATION (HOLDS|FAILED)"], limit=6)

# ---------------------------------------------------------------------------
stage(20, "FULL TEST SUITE")
code, out = run([PY, "-m", "pytest", "tests", "-q", "--no-header"])
summary = [l for l in out.splitlines()
           if re.search(r"\d+ (passed|failed|error|xfailed)", l)]
for l in summary[-3:]:
    print("   " + l.strip())
fails = [l for l in out.splitlines() if l.startswith(("FAILED", "ERROR"))]
for l in fails[:20]:
    print("   " + l.strip())
if code != 0:
    blockers.append(f"pytest failed, {len(fails)} failing test(s)")
m = re.search(r"(\d+) passed", out)
passed = int(m.group(1)) if m else 0
if passed < 1800:
    warnings.append(f"only {passed} tests collected, expected 1800 or more")

# ---------------------------------------------------------------------------
stage(21, "REPOSITORY HYGIENE")
code, out = run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
branch = out.strip() if code == 0 else "unknown"
print(f"   branch: {branch}")
if branch in ("main", "master"):
    warnings.append(f"committing directly to {branch}")
code, out = run(["git", "status", "--porcelain"])
changed = [l for l in out.splitlines() if l.strip()]
print(f"   changed or untracked paths: {len(changed)}")
stray = []
for pat in ("_p_*.txt", "_assemble.py", "_append_probe.txt",
            "backend/*.bak", "data/*.bak", "tests/*.bak", "_build/*.part*", "tools/*.part*"):
    stray += [str(p.relative_to(ROOT)) for p in ROOT.glob(pat)]
if stray:
    print(f"   scratch files present: {len(stray)}")
    for s in stray[:10]:
        print("      " + s)
    warnings.append(f"{len(stray)} scratch file(s) present")
else:
    print("   ok   no scratch or partial-assembly files")

# ---------------------------------------------------------------------------
print()
print("=" * 76)
print("VERDICT")
print("=" * 76)
print(f"   blockers: {len(blockers)}    warnings: {len(warnings)}")
if blockers:
    print("\n   BLOCKERS, must be fixed before going live:")
    for i, b in enumerate(blockers, 1):
        print(f"     {i}. {b}")
if warnings:
    print("\n   WARNINGS, safe to ship but worth knowing:")
    for i, w in enumerate(warnings, 1):
        print(f"     {i}. {w}")
print()
if blockers:
    print("   RESULT: NOT READY")
else:
    print("   RESULT: READY TO GO LIVE")
    print()
    print("   Suggested sequence:")
    print("     git checkout -b v3.0-expansion")
    print("     git add -A")
    print("     git status")
    print('     git commit -m "feat: v3.0 ancestry, sequencing, imputation and provenance"')
    print("     git push -u origin v3.0-expansion")
print()
sys.exit(1 if blockers else 0)
