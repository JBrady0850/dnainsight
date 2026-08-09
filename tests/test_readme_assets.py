"""Tests that every image the documentation references will actually render on GitHub.

WHY THIS FILE EXISTS
--------------------
A broken README image is invisible to the person who broke it. It renders fine
on their machine and 404s for everybody else, and nothing in the test suite or
the release gate notices, because the image is not code.

Three ways it breaks, all of them silent on Windows:

  1. The file gets gitignored. `.gitignore` grew binary patterns in v3.0 for
     panels and alignments, and one careless `*.png` would remove the only
     screenshot from the repository while leaving it on disk.

  2. The case stops matching. `git config core.ignorecase` is `true` on Windows,
     so `DNAInsight.png` on disk satisfies a README that says `dnainsight.png`
     locally, and GitHub, which is case-sensitive, shows a broken image.

  3. The file stops being the format its extension claims, usually after a
     conversion or an optimisation pass writes the wrong container.

None of this needs the network or a GitHub account. It is all checkable from
the repository itself, which is the only kind of check that survives.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

# Markdown inline images: ![alt](path). Also catches HTML <img src="path">,
# because README files pick up raw HTML over time.
_MD_IMAGE = re.compile(r"!\[[^\]]*\]\(([^)\s]+)")
_HTML_IMAGE = re.compile(r"<img[^>]+src=[\"']([^\"']+)[\"']", re.IGNORECASE)

# First bytes of each container this project is likely to reference.
_SIGNATURES = {
    ".png": [b"\x89PNG\r\n\x1a\n"],
    ".jpg": [b"\xff\xd8\xff"],
    ".jpeg": [b"\xff\xd8\xff"],
    ".gif": [b"GIF87a", b"GIF89a"],
    ".webp": [b"RIFF"],
}

DOC_FILES = sorted(
    p for p in ROOT.rglob("*.md")
    if ".git" not in p.parts and "node_modules" not in p.parts
)


def _local_image_refs(path: Path) -> list[str]:
    """Every image reference in one document that points inside the repository."""
    text = path.read_text(encoding="utf-8", errors="replace")
    refs = _MD_IMAGE.findall(text) + _HTML_IMAGE.findall(text)
    out = []
    for ref in refs:
        # Shields.io badges and any other absolute URL are somebody else's
        # problem, and a network check here would make the suite flaky.
        if ref.startswith(("http://", "https://", "data:", "//", "#")):
            continue
        out.append(ref.split("#")[0].split("?")[0])
    return out


def _all_refs() -> list[tuple[Path, str]]:
    pairs = []
    for doc in DOC_FILES:
        for ref in _local_image_refs(doc):
            pairs.append((doc, ref))
    return pairs


ALL_REFS = _all_refs()


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True)


def _git_available() -> bool:
    return _git("rev-parse", "--git-dir").returncode == 0


# ---------------------------------------------------------------------------
# The reference set itself
# ---------------------------------------------------------------------------

def test_the_readme_references_at_least_one_local_image():
    # Guards against this whole file quietly becoming a no-op if the screenshot
    # is ever removed from the README without anyone deciding to remove it.
    readme = ROOT / "README.md"
    assert _local_image_refs(readme), "README.md references no local image"


def test_every_referenced_image_exists():
    missing = [f"{doc.relative_to(ROOT)} -> {ref}"
               for doc, ref in ALL_REFS
               if not (doc.parent / ref).is_file()]
    assert not missing, "referenced image does not exist:\n  " + "\n  ".join(missing)


def test_every_referenced_image_matches_on_disk_case_exactly():
    """The Windows trap, stated precisely.

    Path.is_file() returns True for the wrong case on a case-insensitive
    filesystem, so existence alone proves nothing. Comparing the reference
    against the real directory listing is the only check that catches a
    mismatch before GitHub does.
    """
    wrong = []
    for doc, ref in ALL_REFS:
        target = (doc.parent / ref)
        if not target.is_file():
            continue
        actual = {p.name for p in target.parent.iterdir()}
        if target.name not in actual:
            match = next((n for n in actual if n.lower() == target.name.lower()), "?")
            wrong.append(
                f"{doc.relative_to(ROOT)} references {target.name!r} "
                f"but the file on disk is {match!r}"
            )
    assert not wrong, (
        "case mismatch, works on Windows and 404s on GitHub:\n  " + "\n  ".join(wrong)
    )


def test_no_referenced_image_uses_a_backslash():
    # A Windows path separator in Markdown is not a path, it is an escape.
    bad = [f"{doc.relative_to(ROOT)} -> {ref}" for doc, ref in ALL_REFS if "\\" in ref]
    assert not bad, "backslash in an image path:\n  " + "\n  ".join(bad)


def test_no_referenced_image_escapes_the_repository():
    outside = []
    for doc, ref in ALL_REFS:
        resolved = (doc.parent / ref).resolve()
        if not str(resolved).startswith(str(ROOT.resolve())):
            outside.append(f"{doc.relative_to(ROOT)} -> {ref}")
    assert not outside, "image path leaves the repository:\n  " + "\n  ".join(outside)


def test_no_referenced_image_is_an_absolute_path():
    # An absolute path renders on the author's machine and nowhere else.
    bad = [f"{doc.relative_to(ROOT)} -> {ref}" for doc, ref in ALL_REFS
           if ref.startswith("/") or re.match(r"^[A-Za-z]:", ref)]
    assert not bad, "absolute image path:\n  " + "\n  ".join(bad)


# ---------------------------------------------------------------------------
# File contents
# ---------------------------------------------------------------------------

def test_every_referenced_image_has_the_signature_its_extension_claims():
    wrong = []
    for doc, ref in ALL_REFS:
        target = doc.parent / ref
        if not target.is_file():
            continue
        expected = _SIGNATURES.get(target.suffix.lower())
        if not expected:
            continue
        head = target.read_bytes()[:8]
        if not any(head.startswith(sig) for sig in expected):
            wrong.append(f"{target.relative_to(ROOT)} starts with {head!r}")
    assert not wrong, "file is not the format its extension claims:\n  " + "\n  ".join(wrong)


def test_no_referenced_image_is_empty():
    empty = [str((doc.parent / ref).relative_to(ROOT)) for doc, ref in ALL_REFS
             if (doc.parent / ref).is_file() and (doc.parent / ref).stat().st_size == 0]
    assert not empty, "zero-byte image:\n  " + "\n  ".join(empty)


def test_the_readme_screenshot_stays_a_reasonable_size():
    """GitHub serves the raw file to every visitor, so this is a real cost.

    Not a correctness failure, but a ceiling worth having: the screenshot was
    425 KB before it was palettised and there is no reason for it to drift back.
    """
    shot = ROOT / "DNAInsight.png"
    if not shot.is_file():
        pytest.skip("no screenshot in this checkout")
    size_kb = shot.stat().st_size / 1024
    assert size_kb < 400, f"DNAInsight.png is {size_kb:.0f} KB, palettise it"


def test_the_readme_screenshot_carries_no_personal_identifiers():
    """The capture before v3.1 showed a real name and date of birth.

    That image ships in a public repository, so this is a privacy check rather
    than a style one. The alt text is what a reader and a screen reader both
    receive, so it is checked directly; the pixels were verified at capture
    time by asserting the rendered page contained none of these strings.
    """
    readme = (ROOT / "README.md").read_text(encoding="utf-8", errors="replace")
    alt_texts = re.findall(r"!\[([^\]]*)\]\([^)]*\)", readme)
    forbidden = ("brady", "1970-06-03", "bradyj")
    for alt in alt_texts:
        lowered = alt.lower()
        for token in forbidden:
            assert token not in lowered, f"personal identifier {token!r} in image alt text"


# ---------------------------------------------------------------------------
# Git will actually publish it
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _git_available(), reason="not a git checkout")
def test_no_referenced_image_is_gitignored():
    """An ignored image is present locally and absent from the clone GitHub renders."""
    ignored = []
    for doc, ref in ALL_REFS:
        target = doc.parent / ref
        if not target.is_file():
            continue
        result = _git("check-ignore", "-q", str(target))
        if result.returncode == 0:
            ignored.append(str(target.relative_to(ROOT)))
    assert not ignored, (
        "image is gitignored and will not exist in the pushed repository:\n  "
        + "\n  ".join(ignored)
    )


@pytest.mark.skipif(not _git_available(), reason="not a git checkout")
def test_the_readme_screenshot_is_tracked_under_the_referenced_name():
    """Compare against git's index, which is what actually gets pushed.

    The working tree can disagree with the index on case while
    core.ignorecase is true, and the index is the version that wins.
    """
    shot = ROOT / "DNAInsight.png"
    if not shot.is_file():
        pytest.skip("no screenshot in this checkout")
    tracked = _git("ls-files", "--", "*.png").stdout.split()
    assert "DNAInsight.png" in tracked, (
        f"README references DNAInsight.png but git tracks {tracked!r}"
    )


@pytest.mark.skipif(not _git_available(), reason="not a git checkout")
def test_no_referenced_image_needs_git_lfs():
    """LFS pointers render as broken images wherever LFS is not configured."""
    needs_lfs = []
    for doc, ref in ALL_REFS:
        target = doc.parent / ref
        if not target.is_file():
            continue
        result = _git("check-attr", "filter", "--", str(target))
        if "filter: lfs" in result.stdout:
            needs_lfs.append(str(target.relative_to(ROOT)))
    assert not needs_lfs, "image is routed through git-lfs:\n  " + "\n  ".join(needs_lfs)


# ---------------------------------------------------------------------------
# The support link
#
# It is the one element in this README that exists to be clicked, so a typo in
# it fails silently: the badge still renders and the link still looks right.
# ---------------------------------------------------------------------------

SUPPORT_URL = "https://buymeacoffee.com/jbrady2852"


def test_the_readme_carries_the_support_link():
    readme = (ROOT / "README.md").read_text(encoding="utf-8", errors="replace")
    assert SUPPORT_URL in readme


def test_the_support_link_sits_in_the_badge_block_at_the_top():
    # Above the screenshot, so it is visible without scrolling on any screen.
    lines = (ROOT / "README.md").read_text(encoding="utf-8", errors="replace").splitlines()
    link_line = next(i for i, l in enumerate(lines) if SUPPORT_URL in l)
    shot_line = next(i for i, l in enumerate(lines) if "DNAInsight.png" in l)
    assert link_line < shot_line, "support link fell below the screenshot"
    assert link_line < 15, f"support link drifted to line {link_line + 1}, expected the badge block"


def test_the_support_badge_is_a_clickable_link_not_a_bare_image():
    # ![alt](badge) renders the badge and goes nowhere. The wrapping [ ]( ) is
    # the entire point, and it is easy to lose when reordering badges.
    readme = (ROOT / "README.md").read_text(encoding="utf-8", errors="replace")
    line = next(l for l in readme.splitlines() if SUPPORT_URL in l)
    assert line.startswith("[!["), f"support badge is not wrapped in a link: {line[:40]}"
    assert line.rstrip().endswith(f"({SUPPORT_URL})"), "link target is not the support URL"
