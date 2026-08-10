# The README screenshot

- **File:** `DNAInsight.png`
- **Captured at version:** 3.4.1
- **Captured on:** 2026-08-10
- **Viewport:** 1366 x 900 at device scale factor 2, full page
- **Regenerate with:** `python tools/capture_screenshot.py`

## Why this file exists

An image cannot be parsed by the test suite, so it cannot be checked for truth
the way a number in the README can. The previous capture proved the cost of
that: it carried a **v3.1.0** banner in the page header and shipped unchanged
through v3.2, v3.3 and v3.4. The first visual any reader saw was a version three
releases old, and nothing in the repository could notice, because nothing in the
repository reads pixels.

This file is the machine-readable part of the image. `Captured at version` above
is checked against `backend.__version__` by
`tests/test_readme_currency.py::test_the_screenshot_records_the_version_it_was_captured_at`,
so the suite fails on the release that makes the screenshot stale rather than
three releases later.

Recording the version without recording the method would only move the problem,
so the capture is a script rather than a manual step, and a second test asserts
that script still exists.

## What is in the frame, and why this shot

The dashboard for a synthetic profile:

- the four finding counts, split by category
- the next-step actions
- the **What this build can do** table

That table is the reason this is the chosen shot. It reports three states rather
than two. `available` means DNAInsight ships the capability and it has data.
`not built` means the code is here and a builder has not been run. `separate
tool` means a third-party program is required, which the user installs
themselves under its own licence. Collapsing the last two would tell somebody to
run a builder that does not exist for that feature, and the distinction is the
same one the genoset engine draws between a rule that was checked and found
absent and a rule that was never testable.

## What the capture deliberately shows

**The profile is synthetic and reads `generic`, not `23andme`.** The fixture at
`tests/fixtures/sample_23andme.txt` is a fabricated file with placeholder
positions, so provider detection correctly declines to claim it came from
23andMe. Faking that field for a nicer screenshot would put a false detection
result in the most visible image in the project.

**The counts are small** (97 findings, not the hundreds a real export produces)
because the fixture is roughly a hundred rows and omits every fifth rsID on
purpose, so several multi-SNP genosets are genuinely incomplete.

**The date of birth is 1980-01-01 and the name is "Sample Profile".** The
capture before v3.1 showed a real name and a real date of birth in a public
repository. `tests/test_readme_assets.py` now checks the alt text for personal
identifiers, and the capture script hardcodes synthetic values so the pixels
cannot drift back.

## Isolation

`tools/capture_screenshot.py` boots the real application and creates a profile,
which against default paths would write into the user's own `dnainsight.db` and
`uploads/` directory. That is the exact fault `tools/isolated_db.py` documents.
The script sets `DNAINSIGHT_DB_PATH`, `DNAINSIGHT_UPLOAD_DIR` and
`DNAINSIGHT_REPORTS_DIR` to a temporary directory **before importing `app`**,
and removes that directory afterwards.

## Size

GitHub serves the raw file to every visitor, so
`tests/test_readme_assets.py::test_the_readme_screenshot_stays_a_reasonable_size`
caps it at 400 KB. A full-colour 2x capture is around 450 KB, so the script
palettises, stepping 256, 128 then 64 colours until it fits. The interface is
flat colour, so the loss is not visible. The current file is 166 KB.
