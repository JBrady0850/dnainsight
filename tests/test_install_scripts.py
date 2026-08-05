"""Tests that lock install.bat and install.sh to reality.

WHY THIS FILE EXISTS
--------------------
Install scripts are the one part of a project nobody runs during development,
so they rot silently and the first person to notice is a new user on a clean
machine. Every assertion here encodes a way these two scripts have already gone
stale or could:

  - referencing a builder or requirements file that was renamed or deleted
  - the two scripts drifting apart, so Windows and macOS install different things
  - step counters saying "[3/4]" after a fifth step was added
  - claiming success without ever importing the application
  - pointing at a document that does not exist

None of this needs the scripts to be executed. It is all statically checkable,
which means it can run in CI on any platform, which is the only way a Windows
batch file ever gets tested from a Linux runner.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

BAT = ROOT / "install.bat"
SH = ROOT / "install.sh"

EM_DASH = "—"


@pytest.fixture(scope="module")
def bat_text() -> str:
    return BAT.read_text(encoding="utf-8", errors="replace")


@pytest.fixture(scope="module")
def sh_text() -> str:
    return SH.read_text(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Existence
# ---------------------------------------------------------------------------

def test_windows_installer_exists():
    assert BAT.is_file()


def test_unix_installer_exists():
    assert SH.is_file()


def test_unix_installer_has_a_shebang(sh_text):
    assert sh_text.startswith("#!")


def test_unix_installer_shebang_is_portable(sh_text):
    # /usr/bin/env bash rather than /bin/bash, because macOS ships an ancient
    # bash at /bin/bash and Homebrew installs a modern one elsewhere.
    assert sh_text.splitlines()[0] == "#!/usr/bin/env bash"


def test_unix_installer_sets_errexit(sh_text):
    assert re.search(r"^set -e\s*$", sh_text, re.MULTILINE)


# ---------------------------------------------------------------------------
# Referenced files must exist
# ---------------------------------------------------------------------------

def _referenced_paths(text: str) -> set[str]:
    """Every repository-relative path either script hands to Python or pip."""
    found = set()
    for match in re.finditer(r"(data[\\/][a-z_]+\.py)", text):
        found.add(match.group(1).replace("\\", "/"))
    for match in re.finditer(r"(requirements[a-z-]*\.txt)", text):
        found.add(match.group(1))
    for match in re.finditer(r"(docs[\\/][A-Z_]+\.md)", text):
        found.add(match.group(1).replace("\\", "/"))
    for match in re.finditer(r"open\('([^']+)'\)", text):
        found.add(match.group(1))
    return found


def test_windows_installer_references_only_files_that_exist(bat_text):
    for rel in sorted(_referenced_paths(bat_text)):
        assert (ROOT / rel).exists(), f"install.bat references missing {rel}"


def test_unix_installer_references_only_files_that_exist(sh_text):
    for rel in sorted(_referenced_paths(sh_text)):
        assert (ROOT / rel).exists(), f"install.sh references missing {rel}"


def test_both_installers_reference_the_same_builders(bat_text, sh_text):
    # Windows and macOS must install the same thing. A builder added to one
    # script and not the other produces two different products with one name.
    builders_bat = {p for p in _referenced_paths(bat_text) if p.startswith("data/")}
    builders_sh = {p for p in _referenced_paths(sh_text) if p.startswith("data/")}
    assert builders_bat == builders_sh


def test_both_installers_reference_the_same_requirements(bat_text, sh_text):
    reqs_bat = {p for p in _referenced_paths(bat_text) if p.startswith("requirements")}
    reqs_sh = {p for p in _referenced_paths(sh_text) if p.startswith("requirements")}
    assert reqs_bat == reqs_sh


def test_the_bundled_reference_builder_is_run(bat_text, sh_text):
    assert "data/build_reference.py" in _referenced_paths(bat_text)
    assert "data/build_reference.py" in _referenced_paths(sh_text)


# ---------------------------------------------------------------------------
# Step numbering
# ---------------------------------------------------------------------------

def _steps(text: str) -> list[tuple[int, int]]:
    return [(int(a), int(b)) for a, b in re.findall(r"\[(\d+)/(\d+)\]", text)]


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_step_counters_are_contiguous_and_agree_on_the_total(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    steps = _steps(text)
    assert steps, "no step counters found"
    totals = {total for _, total in steps}
    assert len(totals) == 1, f"step counters disagree on the total: {totals}"
    total = totals.pop()
    numbers = [n for n, _ in steps]
    assert numbers == sorted(numbers), "step numbers are out of order"
    assert numbers == list(range(1, total + 1)), (
        f"expected steps 1..{total}, got {numbers}. Renumber after adding a step."
    )


def test_both_installers_have_the_same_number_of_steps(bat_text, sh_text):
    assert _steps(bat_text)[-1][1] == _steps(sh_text)[-1][1]


# ---------------------------------------------------------------------------
# Verification step
#
# The whole point of the v3.0 change. An installer that prints a success banner
# without importing the application has verified nothing.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["bat", "sh"])
def test_installer_verifies_before_declaring_success(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    assert "Verifying install" in text


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_verification_builds_the_flask_app(name, bat_text, sh_text):
    # Importing backend alone proves very little. create_app() is where a
    # missing data file, a broken blueprint or a bad schema migration actually
    # surfaces.
    text = bat_text if name == "bat" else sh_text
    assert "create_app()" in text


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_verification_reports_the_application_version(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    assert "APP_VERSION" in text


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_verification_failure_is_fatal(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    assert "Verification failed" in text
    if name == "bat":
        assert "exit /b 1" in text
    else:
        assert "exit 1" in text


def test_the_verification_command_actually_works():
    # Run the same three checks the installers run, in this process. If this
    # fails, both installers would fail on a clean machine, and finding that
    # out here beats finding it out from a user.
    import json

    import backend
    from app import create_app

    create_app()
    with open(ROOT / "data" / "snp_reference.json", encoding="utf-8") as fh:
        payload = json.load(fh)
    assert payload
    assert backend.APP_VERSION


# ---------------------------------------------------------------------------
# Launcher creation
# ---------------------------------------------------------------------------

def test_windows_installer_creates_a_launcher(bat_text):
    assert "launch.bat" in bat_text


def test_unix_installer_creates_a_launcher(sh_text):
    assert "launch.sh" in sh_text
    assert "chmod +x launch.sh" in sh_text


def test_unix_launcher_changes_directory_before_running(sh_text):
    # launch.sh is double-clickable from a file manager, which starts it with an
    # arbitrary working directory. Without the cd, app.py cannot find frontend/.
    assert 'cd "$(dirname "$0")"' in sh_text or 'cd "\\$(dirname "\\$0")"' in sh_text


def test_windows_launcher_changes_directory_before_running(bat_text):
    assert "cd /d" in bat_text


# ---------------------------------------------------------------------------
# The v3.0 external-tool disclosure
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["bat", "sh"])
def test_installer_discloses_that_external_tools_are_not_installed(name, bat_text, sh_text):
    # A user who installs DNAInsight and then cannot find ancestry has hit a
    # licence boundary, not a bug. Saying so at install time is cheaper than
    # answering it later.
    text = bat_text if name == "bat" else sh_text
    assert "EXTERNAL_TOOLS.md" in text
    assert "NOT installed" in text


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_installer_names_the_out_of_tree_tools_directory(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    assert ".dnainsight" in text


def test_the_external_tools_document_exists():
    assert (ROOT / "docs" / "EXTERNAL_TOOLS.md").is_file()


@pytest.mark.parametrize("name", ["bat", "sh"])
def test_installer_promises_offline_operation(name, bat_text, sh_text):
    text = bat_text if name == "bat" else sh_text
    assert "offline" in text.lower()


# ---------------------------------------------------------------------------
# Portability traps
# ---------------------------------------------------------------------------

def test_windows_installer_invokes_pip_through_python(bat_text):
    # A machine can have python.exe on PATH without pip.exe on PATH. Calling
    # bare "pip" then fails with "not recognized" two steps after we already
    # proved python works, which is a confusing place to fail.
    for match in re.finditer(r"^\s*(?!.*python -m )pip\s+install", bat_text, re.MULTILINE):
        pytest.fail(f"bare pip call in install.bat: {match.group(0)!r}")


def test_unix_installer_invokes_pip_through_python(sh_text):
    for match in re.finditer(r"^\s*(?!.*python3 -m )(?!.*\.venv/bin/)pip\s+install",
                             sh_text, re.MULTILINE):
        pytest.fail(f"bare pip call in install.sh: {match.group(0)!r}")


def test_unix_installer_handles_pep668_managed_environments(sh_text):
    # Ubuntu 23.04+, Debian 12+ and Homebrew Python refuse a plain pip install.
    # Both escape hatches must be present or the installer dies on a modern Mac.
    assert "--break-system-packages" in sh_text
    assert "venv" in sh_text


def test_unix_installer_offers_a_package_manager_for_macos(sh_text):
    assert "brew" in sh_text


def test_windows_installer_enables_delayed_expansion_if_it_uses_it(bat_text):
    if "!" in bat_text and re.search(r"![A-Z_]+!", bat_text):
        assert "EnableDelayedExpansion" in bat_text


def test_windows_installer_does_not_assume_a_refreshed_path(bat_text):
    # Installing Python does not update PATH in the already-running shell. The
    # script has to either find the interpreter itself or tell the user to
    # reopen the window.
    assert "LOCALAPPDATA" in bat_text or "reopen" in bat_text.lower() or "again" in bat_text.lower()


# ---------------------------------------------------------------------------
# House style
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", [BAT, SH])
def test_no_em_dashes(path):
    assert EM_DASH not in path.read_text(encoding="utf-8", errors="replace")


@pytest.mark.parametrize("path", [BAT, SH])
def test_scripts_are_ascii(path):
    raw = path.read_bytes()
    non_ascii = [b for b in raw if b > 127]
    assert not non_ascii, (
        f"{path.name} contains non-ASCII bytes. Windows console code pages "
        f"mangle them and the installer output becomes unreadable."
    )


def test_unix_installer_parses_under_bash():
    """Syntax-check install.sh with the real shell that will run it.

    Skipped on Windows, and the reason is not squeamishness. shutil.which finds
    C:\\Windows\\System32\\bash.EXE, which is the WSL launcher. That bash lives in
    a different filesystem namespace, so handing it a Windows path produces
    "D:dnainsightinstall.sh: No such file or directory" and the test fails on a
    file that is perfectly valid. Translating the path to /mnt/d/... would then
    be checking whether WSL is installed and configured, which is a fact about
    the developer's laptop and not about this repository.

    The check still runs everywhere it means something: any Linux or macOS
    machine, and the CI job, which is where install.sh is actually executed.
    """
    import shutil
    import subprocess

    if sys.platform.startswith("win"):
        pytest.skip("bash on Windows is the WSL launcher and cannot read a Windows path")
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash not available on this platform")
    result = subprocess.run([bash, "-n", str(SH)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# ---------------------------------------------------------------------------
# Batch file quoting traps
#
# Both of these shipped in v1 and v2 and neither shows up in a code review,
# because the line looks like prose. They only surface when somebody actually
# runs the installer and reads the output, which is why these are tests now.
# ---------------------------------------------------------------------------

_LAUNCHER_REDIRECTS = ("launch.bat",)


def test_no_echo_line_accidentally_redirects_to_a_file(bat_text):
    """An unescaped > in an echo is a redirect, not a greater-than sign.

    The shipped v2 installer contained:

        echo TIP: Open DNAInsight and use Settings > Database to update

    which printed a truncated tip and silently created a file literally named
    "Database" in every Windows install directory. Escape it as ^> when it is
    meant as punctuation.
    """
    offenders = []
    for number, line in enumerate(bat_text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped.lower().startswith("echo"):
            continue
        # Writing the launcher is a deliberate redirect and is allowed.
        if any(target in stripped for target in _LAUNCHER_REDIRECTS):
            continue
        for match in re.finditer(r"(.)?([<>])", stripped):
            if match.group(1) != "^":
                offenders.append(f"line {number}: {stripped}")
                break
    assert not offenders, (
        "unescaped redirect in an echo, use ^> or ^< :\n  " + "\n  ".join(offenders)
    )


def test_no_echo_line_loses_a_bang_to_delayed_expansion(bat_text):
    """With EnableDelayedExpansion on, a bare ! is consumed by the parser.

    "echo Installation Complete!" printed "Installation Complete". Escape it as
    ^! when the exclamation mark is meant to reach the screen.
    """
    if "EnableDelayedExpansion" not in bat_text:
        pytest.skip("delayed expansion is not enabled, so a bare ! is literal")
    offenders = []
    for number, line in enumerate(bat_text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped.lower().startswith("echo"):
            continue
        for match in re.finditer(r"(.)?!", stripped):
            if match.group(1) != "^":
                offenders.append(f"line {number}: {stripped}")
                break
    assert not offenders, (
        "bare ! in an echo under delayed expansion, use ^! :\n  " + "\n  ".join(offenders)
    )


def test_the_installer_creates_no_stray_files_named_after_a_word(bat_text):
    # Direct regression guard for the exact file the old bug produced.
    assert "Settings > Database" not in bat_text
