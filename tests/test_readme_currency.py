"""The README is part of the build, so it is tested like part of the build.

WHY THIS FILE EXISTS
--------------------
`tests/test_readme_assets.py` already proves the README's images will RENDER.
Nothing proved the README was TRUE. Those are different failures and the second
one is worse, because a README that renders perfectly while describing the
previous release is believed.

Every number in a README rots on a schedule nobody sets. Before v3.4.0 this one
carried three different test counts in three places, two of them stale and none
of them equal, plus an endpoint count that was one short from the release that
added the concordance route. Each was correct when written. None was wrong in a
way a human reader would notice.

So every count the README states is checked against the repository, and every
module and document it lists is checked for existence in both directions. The
screenshot is covered too, through a recorded capture version, because an image
cannot be parsed but it can lie just as loudly: the v3.1.0 banner in the
dashboard screenshot survived three releases.

WHAT THIS DELIBERATELY DOES NOT DO
----------------------------------
It does not check prose. A claim like "nothing leaves your machine" is enforced
by `tests/test_uploads.py` and the privacy tests, not here. This file only
covers facts the repository can answer for itself, because a test that guesses
at intent produces false failures and gets deleted.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import backend  # noqa: E402

README = (ROOT / "README.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------

def _collected_test_count() -> int:
    """Ask pytest, rather than counting `def test_` and being wrong.

    An AST count misses every parametrised case, and this suite parametrises
    heavily, so the only honest number is the one the collector produces.
    `--collect-only` does not execute anything, so this is not recursive.
    """
    result = subprocess.run(
        [sys.executable, "-m", "pytest", str(ROOT / "tests"), "--collect-only", "-q",
         "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=str(ROOT),
    )
    match = re.search(r"(\d+) tests? collected", result.stdout)
    if not match:
        pytest.skip(f"could not read a collection count from pytest:\n{result.stdout[-500:]}")
    return int(match.group(1))


def _readme_test_counts() -> list[tuple[str, int]]:
    """Every place the README states how many tests there are."""
    found = []
    badge = re.search(r"Tests-(\d+)_passing", README)
    if badge:
        found.append(("badge", int(badge.group(1))))
    table = re.search(r"^\| Tests \| \d+ \| (\d+) \|", README, re.M)
    if table:
        found.append(("comparison table", int(table.group(1))))
    tree = re.search(r"^├── tests/\s+(\d+) tests", README, re.M)
    if tree:
        found.append(("project layout", int(tree.group(1))))
    return found


def test_the_readme_states_a_test_count_somewhere():
    # Without this, deleting every count would make the checks below vacuous.
    assert _readme_test_counts(), "README states no test count at all"


def test_every_stated_test_count_agrees_with_every_other():
    """Three places used to hold three different numbers."""
    counts = _readme_test_counts()
    distinct = {n for _, n in counts}
    assert len(distinct) == 1, (
        "README disagrees with itself about the test count: "
        + ", ".join(f"{where}={n}" for where, n in counts)
    )


def test_the_stated_test_count_is_the_count_pytest_collects():
    stated = _readme_test_counts()[0][1]
    actual = _collected_test_count()
    assert stated == actual, (
        f"README says {stated} tests, pytest collects {actual}. "
        f"Update the badge, the comparison table and the project layout together."
    )


def test_the_stated_v3_endpoint_count_matches_the_routes_that_exist():
    stated = re.search(r"plus (\d+) v3 paths", README)
    assert stated, "README no longer states a v3 path count"
    source = (ROOT / "backend" / "routes_v3.py").read_text(encoding="utf-8")
    actual = source.count("@api_v3.route(")
    assert int(stated.group(1)) == actual, (
        f"README says {stated.group(1)} v3 paths, routes_v3.py defines {actual}"
    )


# ---------------------------------------------------------------------------
# The project layout tree, checked in both directions
# ---------------------------------------------------------------------------

# Listed for a reason rather than by accident. A private helper module is not
# worth a line in a public README; anything here that stops existing should
# fail rather than mislead.
_LAYOUT_EXEMPT_MODULES = {"__init__.py"}


def _layout_block() -> str:
    match = re.search(r"```\ndnainsight/\n(.*?)```", README, re.S)
    assert match, "the project layout block is gone from the README"
    return match.group(1)


def test_every_backend_module_appears_in_the_project_layout():
    listed = _layout_block()
    missing = sorted(
        p.name for p in (ROOT / "backend").glob("*.py")
        if p.name not in _LAYOUT_EXEMPT_MODULES and p.name not in listed
    )
    assert not missing, (
        "backend modules exist that the README's project layout does not list: "
        + ", ".join(missing)
    )


def test_every_document_appears_in_the_project_layout():
    listed = _layout_block()
    missing = sorted(
        p.name for p in (ROOT / "docs").glob("*.md") if p.name not in listed
    )
    assert not missing, (
        "documents exist that the README's project layout does not list: "
        + ", ".join(missing)
    )


def test_the_layout_lists_nothing_that_stopped_existing():
    """The other direction. A removed module leaves a line behind silently."""
    listed = _layout_block()
    ghosts = []
    for name in re.findall(r"([A-Za-z0-9_.-]+\.(?:py|md|html|txt))", listed):
        if not list(ROOT.rglob(name)):
            ghosts.append(name)
    assert not ghosts, (
        "the README's project layout names files that do not exist: "
        + ", ".join(sorted(set(ghosts)))
    )


# ---------------------------------------------------------------------------
# The screenshot
# ---------------------------------------------------------------------------

SCREENSHOT_DOC = ROOT / "docs" / "SCREENSHOT.md"


def test_the_screenshot_records_the_version_it_was_captured_at():
    """An image cannot be parsed, so the capture records its own version.

    This is not bureaucracy. The dashboard screenshot carried a "v3.1.0" banner
    through v3.2, v3.3 and v3.4, and every reader of the README saw a version
    number three releases old as the first visual on the page. Nothing in the
    suite could see it, because nothing in the suite reads pixels.
    """
    assert SCREENSHOT_DOC.is_file(), "docs/SCREENSHOT.md is missing"
    text = SCREENSHOT_DOC.read_text(encoding="utf-8")
    match = re.search(r"^- \*\*Captured at version:\*\* (\S+)", text, re.M)
    assert match, "docs/SCREENSHOT.md does not record a capture version"
    assert match.group(1) == backend.__version__, (
        f"the screenshot was captured at v{match.group(1)} and the application "
        f"is v{backend.__version__}. Re-capture it: the version banner is the "
        f"first thing in the image."
    )


def test_the_screenshot_doc_records_how_to_regenerate_it():
    # A recorded version with no recorded method just moves the staleness.
    text = SCREENSHOT_DOC.read_text(encoding="utf-8")
    assert "tools/capture_screenshot.py" in text
    assert (ROOT / "tools" / "capture_screenshot.py").is_file()
