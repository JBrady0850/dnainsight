"""capture_screenshot.py -- regenerate the README dashboard screenshot.

WHY THIS IS A SCRIPT AND NOT A MANUAL TASK
------------------------------------------
The previous screenshot showed a "v3.1.0" banner and shipped unchanged through
v3.2, v3.3 and v3.4. It was the first visual on the README and it was wrong for
three releases, because re-capturing it was a manual step nobody owned and
nothing could detect: no test reads pixels.

So the capture is scripted, its version is recorded in docs/SCREENSHOT.md, and
tests/test_readme_currency.py fails when that recorded version falls behind
backend.__version__. The staleness is now impossible to miss even though the
image itself is still opaque to the suite.

ISOLATION, WHICH IS NOT OPTIONAL
--------------------------------
This boots the real application and creates a profile. Run against the default
paths it would write into the user's own dnainsight.db and uploads/ directory,
which is the exact fault tools/isolated_db.py exists to prevent. Everything here
runs in a temporary directory that is removed afterwards, and the fixture it
uploads is the synthetic one from tests/, never a real export.

WHAT IS IN THE FRAME
--------------------
The dashboard for a synthetic profile: the finding counts, the next-step
actions, and the capability table. That table is the honest part and the reason
this is the chosen shot. It shows three states rather than two, so a reader sees
at a glance which capabilities ship, which need a builder run, and which need a
third-party tool they install themselves.

Usage:
    python tools/capture_screenshot.py
    python tools/capture_screenshot.py --out /tmp/shot.png --keep
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
FIXTURE = ROOT / "tests" / "fixtures" / "sample_23andme.txt"
DEFAULT_OUT = ROOT / "DNAInsight.png"
PORT = 5077
VIEWPORT = (1366, 900)
SCALE = 2  # a 2x capture stays legible when GitHub scales it down

# tests/test_readme_assets.py caps the committed file at 400 KB, because GitHub
# serves it to every visitor. A full-colour 2x capture lands well above that, so
# the palettise step below is required rather than cosmetic.
MAX_KB = 400


def _log(msg: str) -> None:
    print(f"  {msg}", flush=True)


def _start_app(workdir: Path):
    """Boot Flask in a thread, pointed entirely at throwaway paths."""
    os.environ["DNAINSIGHT_DB_PATH"] = str(workdir / "shot.db")
    os.environ["DNAINSIGHT_UPLOAD_DIR"] = str(workdir / "uploads")
    os.environ["DNAINSIGHT_REPORTS_DIR"] = str(workdir / "reports")
    for key in ("DNAINSIGHT_UPLOAD_DIR", "DNAINSIGHT_REPORTS_DIR"):
        Path(os.environ[key]).mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(ROOT))
    import app as app_module  # imported AFTER the env vars, never before

    # create_app() also initialises the schema, so the throwaway database is
    # built here rather than on the first request.
    flask_app = app_module.create_app()
    thread = threading.Thread(
        target=lambda: flask_app.run(host="127.0.0.1", port=PORT,
                                     debug=False, use_reloader=False),
        daemon=True,
    )
    thread.start()
    return flask_app


def _wait_for_server(timeout: float = 30.0) -> None:
    import urllib.error
    import urllib.request
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{PORT}/api/version", timeout=2) as r:
                if r.status == 200:
                    return
        except (urllib.error.URLError, OSError):
            time.sleep(0.25)
    raise RuntimeError("the application did not start")


def _seed_profile() -> int:
    """Create the synthetic profile the screenshot shows, then scan it.

    POST /api/profiles takes the raw file in the SAME multipart request as the
    details; there is no create-then-upload pair. The scan afterwards is what
    puts real counts on the dashboard, and a screenshot of an unscanned profile
    would show four zeros and advertise nothing.
    """
    import json
    import urllib.request
    import uuid

    boundary = uuid.uuid4().hex

    def part(name: str, value: str) -> bytes:
        return (f"--{boundary}\r\n"
                f'Content-Disposition: form-data; name="{name}"\r\n\r\n'
                f"{value}\r\n").encode()

    # A synthetic date of birth on purpose. tests/test_readme_assets.py checks
    # the alt text for personal identifiers; this keeps the pixels clean too.
    body = (
        part("name", "Sample Profile")
        + part("dob", "1980-01-01")
        + part("sex", "M")
        + part("provider", "23andme")
        + (f"--{boundary}\r\n"
           f'Content-Disposition: form-data; name="file"; filename="sample.txt"\r\n'
           f"Content-Type: text/plain\r\n\r\n").encode()
        + FIXTURE.read_bytes()
        + f"\r\n--{boundary}--\r\n".encode()
    )
    req = urllib.request.Request(
        f"http://127.0.0.1:{PORT}/api/profiles",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        created = json.loads(r.read())
    pid = created.get("id") or created.get("profile_id") or (created.get("profile") or {}).get("id")
    if not pid:
        raise RuntimeError(f"no profile id in {created!r}")
    _log(f"profile {pid} created, {created.get('snp_count', '?')} variants read")

    # Offline subsystems only. The API pass batches every rsID and would make
    # this take an hour for no visual difference.
    scan = urllib.request.Request(
        f"http://127.0.0.1:{PORT}/api/profiles/{pid}/scan",
        data=json.dumps({"use_api": False}).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(scan, timeout=600) as r:
        result = json.loads(r.read())
    _log(f"scan complete: {result.get('total_findings', result.get('count', '?'))} findings")
    return pid


def _capture(pid: int, out: Path) -> None:
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(
            viewport={"width": VIEWPORT[0], "height": VIEWPORT[1]},
            device_scale_factor=SCALE,
        )
        page.goto(f"http://127.0.0.1:{PORT}/", wait_until="networkidle")
        # The dashboard renders after /api/version and /api/capabilities
        # resolve, so waiting on the capability table is waiting on both.
        page.wait_for_selector("#caps-list table tr", timeout=30_000)
        page.wait_for_timeout(600)
        page.screenshot(path=str(out), full_page=True)
        browser.close()


def _palettise(path: Path) -> None:
    """Quantise to a palette. The UI is flat colour, so this is lossless enough."""
    from PIL import Image
    before = path.stat().st_size / 1024
    img = Image.open(path).convert("RGB")
    for colours in (256, 128, 64):
        img.convert("P", palette=Image.ADAPTIVE, colors=colours).save(
            path, optimize=True)
        after = path.stat().st_size / 1024
        if after < MAX_KB:
            _log(f"palettised at {colours} colours: {before:.0f} KB -> {after:.0f} KB")
            return
        img = Image.open(path).convert("RGB")
    raise RuntimeError(f"could not get below {MAX_KB} KB, now {after:.0f} KB")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--keep", action="store_true",
                        help="leave the throwaway workdir in place for inspection")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    workdir = Path(tempfile.mkdtemp(prefix="dnainsight-shot-"))
    _log(f"workdir {workdir}")
    try:
        _start_app(workdir)
        _wait_for_server()
        _log("application up")
        pid = _seed_profile()
        _log(f"profile {pid} seeded")
        _capture(pid, out)
        _log(f"captured {out}")
        _palettise(out)
    finally:
        if not args.keep:
            shutil.rmtree(workdir, ignore_errors=True)

    import backend
    _log(f"captured at v{backend.__version__}. "
         f"Record it in docs/SCREENSHOT.md or the suite will fail.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
