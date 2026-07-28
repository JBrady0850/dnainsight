"""Simulate the CI workflow against ONLY the files git would actually publish.

A green local run proves nothing about CI: the working tree holds gitignored
artifacts that a fresh clone will not have. This copies exactly the set
`git ls-files --cached --others --exclude-standard` returns into a temp tree,
then runs the workflow's own commands there, in order.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
fails = []


def run(cmd, cwd, label):
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, shell=False)
    ok = p.returncode == 0
    print("  [%s] %s  (exit %d)" % ("ok" if ok else "FAIL", label, p.returncode))
    out = (p.stdout or "") + (p.stderr or "")
    tail = [ln for ln in out.strip().splitlines() if ln.strip()][-6:]
    for ln in tail:
        print("        " + ln[:150])
    if not ok:
        fails.append(label)
    return p


print("=" * 76)
print("CI SIMULATION AGAINST THE PUBLISHABLE FILE SET")
print("=" * 76)

listed = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                        cwd=str(ROOT), capture_output=True, text=True)
files = [ln for ln in listed.stdout.splitlines() if ln.strip()]
print("git would publish %d files" % len(files))

need = ["data/snp_reference.json", "data/genosets.json", "data/prs_models.json",
        "data/frequencies.json", "requirements.txt", ".flake8",
        ".github/workflows/ci.yml"]
for n in need:
    print("  %-34s %s" % (n, "PUBLISHED" if n in files else "NOT PUBLISHED"))

tmp = Path(tempfile.mkdtemp(prefix="dnai_ci_"))
print("\nclean tree: %s" % tmp)
for rel in files:
    src = ROOT / rel
    if not src.exists():
        continue
    dst = tmp / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
print("copied %d files\n" % len(list(tmp.rglob("*.py"))))

print("JOB: test")
run([PY, "data/build_reference.py"], tmp, "python data/build_reference.py")
run([PY, "-m", "pytest", "tests/", "-q", "--no-header", "-p", "no:cacheprovider"],
    tmp, "python -m pytest tests/")

print("\nJOB: lint")
run([PY, "-m", "flake8", "backend/", "data/", "app.py"], tmp,
    "flake8 backend/ data/ app.py")

print("\n" + "=" * 76)
if fails:
    print("CI WOULD FAIL: " + ", ".join(fails))
else:
    print("CI WOULD PASS: every workflow step green on the publishable file set")
print("=" * 76)
shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if fails else 0)
