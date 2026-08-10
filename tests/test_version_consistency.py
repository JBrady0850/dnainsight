"""The version gate must be able to pass, and must still catch real skew.

WHY THIS FILE EXISTS
--------------------
`tools/vversion.py` compared every version string in the repository against
every other one. The repository declares two unrelated kinds of version, the
application release and the built data artefact build, and those are
deliberately independent, so the check could never pass. It failed on every run
of the release gate from v2.0 onward and was carried as a standing "safe to
ship" warning.

That is worse than not checking. A warning that is always on is a warning
nobody reads, so the one run where the application version genuinely had
skewed would have looked exactly like every other run.

These tests hold both halves of the fix: the gate passes when the repository is
correct, and it still fails when a group is genuinely inconsistent.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

import vversion  # noqa: E402


class TestTheRepositoryIsConsistent:
    def test_the_gate_passes_on_this_checkout(self):
        assert vversion.main() == 0

    def test_the_application_group_agrees_with_itself(self):
        application, _ = vversion.collect()
        assert application, "no application version was found at all"
        assert len(set(application.values())) == 1, (
            f"application version skew: {application}"
        )

    def test_the_artefact_group_agrees_with_itself(self):
        _, artefact = vversion.collect()
        assert artefact, "no artefact version was found at all"
        assert len(set(artefact.values())) == 1, (
            f"data artefact version skew: {artefact}"
        )

    def test_the_backend_version_is_the_one_the_readme_advertises(self):
        # The specific skew a user would actually notice: the badge and heading
        # claiming a release the code does not report.
        application, _ = vversion.collect()
        assert application["backend/__init__.py"] == application["README.md heading"]

    def test_the_changelog_documents_the_version_the_code_reports(self):
        application, _ = vversion.collect()
        assert application["CHANGELOG.md newest"] == application["backend/__init__.py"]


class TestTheGateStillCatchesSkew:
    def test_a_group_with_two_values_is_reported_as_a_mismatch(self):
        assert vversion.report("TEST", {"a": "1.0.0", "b": "2.0.0"}) is False

    def test_a_group_with_one_value_is_reported_as_consistent(self):
        assert vversion.report("TEST", {"a": "1.0.0", "b": "1.0.0"}) is True

    def test_an_empty_group_does_not_fail_the_gate(self):
        # A file being absent is not the same as a version being wrong, and the
        # completeness gate is what checks for missing files.
        assert vversion.report("TEST", {}) is True

    def test_the_two_groups_are_allowed_to_differ_from_each_other(self):
        """The whole point of the split, asserted directly.

        If someone later reinstates a cross-group comparison, this fails and the
        module docstring explains why it must not.
        """
        application, artefact = vversion.collect()
        assert set(application.values()) != set(artefact.values()), (
            "this checkout happens to have matching groups, so the test cannot "
            "prove they are permitted to differ; adjust it rather than deleting it"
        )
        assert vversion.main() == 0
