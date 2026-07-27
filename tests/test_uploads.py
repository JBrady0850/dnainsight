"""Regression guard for the upload destination path.

Two faults were observed in this repository and both are user-facing:

1. UNBOUNDED COMPONENT. The destination was f"{profile_name}_{filename}" with no
   length bound. Re-uploading a file this app had itself named prepended the
   profile name again, so the component grew every cycle. At 246 characters the
   next write crossed the filesystem's 255-byte limit for a single component and
   Windows raised OSError EINVAL, which surfaced to the caller as HTTP 500.
2. SILENT OVERWRITE. Two profiles sharing a name and a filename resolved to one
   path, so the second upload replaced the first person's raw export with no
   error and no warning. Raw DNA cannot be recovered once replaced, which makes a
   silent overwrite strictly worse than a refused upload.

Nothing here touches the real uploads directory; every case runs in tmp_path.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.routes import MAX_UPLOAD_STEM, _bounded_upload_path

# The limit every mainstream filesystem enforces on ONE path component.
COMPONENT_LIMIT = 255
NULL = chr(0)
DNA_EMOJI = chr(0x1F9EC)


class TestTheComponentIsAlwaysWithinTheFilesystemLimit:
    def test_a_short_name_is_left_alone(self, tmp_path):
        dest = _bounded_upload_path(tmp_path, "Jane_Doe_ancestry", ".txt")
        assert dest.name == "Jane_Doe_ancestry.txt"
        assert dest.parent == tmp_path

    def test_a_very_long_stem_is_truncated(self, tmp_path):
        dest = _bounded_upload_path(tmp_path, "x" * 4000, ".txt")
        assert len(dest.name) <= MAX_UPLOAD_STEM + len(".txt")
        assert len(dest.name) < COMPONENT_LIMIT

    def test_the_truncated_name_is_actually_writable(self, tmp_path):
        dest = _bounded_upload_path(tmp_path, "y" * 4000, ".txt")
        dest.write_bytes(b"rsid\tchromosome\tposition\tgenotype\n")
        assert dest.exists()

    @pytest.mark.parametrize("length", [1, 50, 99, 100, 101, 254, 255, 256, 1000])
    def test_no_length_produces_an_oversized_component(self, tmp_path, length):
        dest = _bounded_upload_path(tmp_path, "n" * length, ".txt")
        assert len(dest.name.encode("utf-8")) < COMPONENT_LIMIT

    def test_a_multibyte_name_is_bounded_in_bytes_not_just_characters(self, tmp_path):
        """A UTF-8 character can cost four bytes, so a character cap is not enough."""
        dest = _bounded_upload_path(tmp_path, DNA_EMOJI * 500, ".txt")
        assert len(dest.name.encode("utf-8")) < COMPONENT_LIMIT

    def test_the_extension_survives_truncation(self, tmp_path):
        for ext in (".txt", ".csv", ".tsv"):
            dest = _bounded_upload_path(tmp_path, "z" * 500, ext)
            assert dest.suffix == ext


class TestAnExistingUploadIsNeverOverwritten:
    def test_a_second_identical_request_gets_its_own_path(self, tmp_path):
        first = _bounded_upload_path(tmp_path, "Jane_Doe_ancestry", ".txt")
        first.write_bytes(b"first person")
        second = _bounded_upload_path(tmp_path, "Jane_Doe_ancestry", ".txt")
        assert second != first
        second.write_bytes(b"second person")
        assert first.read_bytes() == b"first person"

    def test_many_collisions_all_resolve_distinctly(self, tmp_path):
        made = []
        for _ in range(12):
            p = _bounded_upload_path(tmp_path, "Same_Name", ".txt")
            p.write_bytes(b"x")
            made.append(p)
        assert len({p.name for p in made}) == 12

    def test_collisions_stay_within_the_limit_even_when_truncated(self, tmp_path):
        made = []
        for _ in range(12):
            p = _bounded_upload_path(tmp_path, "q" * 4000, ".txt")
            p.write_bytes(b"x")
            made.append(p)
        assert len({p.name for p in made}) == 12
        for p in made:
            assert len(p.name.encode("utf-8")) < COMPONENT_LIMIT


FWD = chr(47)
BACK = chr(92)

HOSTILE_STEMS = [
    ".." + FWD + "escape",
    ".." + BACK + "escape",
    ".." + FWD + ".." + FWD + "etc" + FWD + "passwd",
    "sub" + FWD + "dir" + FWD + "name",
    "sub" + BACK + "dir" + BACK + "name",
    "a:b",
    "name" + NULL + "null",
]


class TestThePathCannotEscapeTheUploadDirectory:
    @pytest.mark.parametrize("stem", HOSTILE_STEMS)
    def test_no_separator_or_traversal_survives(self, tmp_path, stem):
        dest = _bounded_upload_path(tmp_path, stem, ".txt")
        assert dest.parent == tmp_path
        assert FWD not in dest.name
        assert BACK not in dest.name
        assert ".." not in dest.name
        assert dest.resolve().is_relative_to(tmp_path.resolve())

    def test_an_empty_stem_still_yields_a_usable_name(self, tmp_path):
        for stem in ("", "   ", "...", "___", FWD * 3):
            dest = _bounded_upload_path(tmp_path, stem, ".txt")
            assert dest.stem, "an empty stem produced a nameless file"
            assert dest.parent == tmp_path

    def test_a_missing_extension_defaults_rather_than_producing_none(self, tmp_path):
        dest = _bounded_upload_path(tmp_path, "profile", "")
        assert dest.suffix == ".txt"


class TestTheGrowthLoopThatCausedTheOutage:
    """Feed the output back in as input: the exact cycle that reached 246 chars."""

    def test_repeated_reupload_converges_instead_of_growing(self, tmp_path):
        prefix = "Collision_Test"
        name = "p1_3_mother"
        lengths = []
        for _ in range(40):
            dest = _bounded_upload_path(tmp_path, prefix + "_" + name, ".txt")
            dest.write_bytes(b"x")
            name = dest.stem
            lengths.append(len(dest.name))
        assert max(lengths) < COMPONENT_LIMIT
        assert max(lengths) <= MAX_UPLOAD_STEM + len(".txt")

    def test_forty_cycles_never_reuse_a_path(self, tmp_path):
        prefix = "Report_Verify"
        name = "sample"
        seen = set()
        for _ in range(40):
            dest = _bounded_upload_path(tmp_path, prefix + "_" + name, ".txt")
            assert dest not in seen, "a cycle reused a path and would overwrite raw DNA"
            dest.write_bytes(b"x")
            seen.add(dest)
            name = dest.stem
        assert len(seen) == 40


class TestBothRouteModulesUseTheOneImplementation:
    BACKEND = Path(__file__).parent.parent / "backend"

    def test_routes_v2_imports_it_rather_than_copying_it(self):
        src = (self.BACKEND / "routes_v2.py").read_text(encoding="utf-8")
        assert "from .routes import _bounded_upload_path" in src
        assert "def _bounded_upload_path" not in src, (
            "a second copy of the rule will drift, exactly as the report filename "
            "helper did before it was consolidated")

    def test_neither_module_still_builds_a_raw_destination(self):
        v1 = (self.BACKEND / "routes.py").read_text(encoding="utf-8")
        v2 = (self.BACKEND / "routes_v2.py").read_text(encoding="utf-8")
        assert 'UPLOAD_DIR / f"{safe_name}_{orig_name}"' not in v1
        assert 'UPLOAD_DIR / f"p{pid}_' not in v2

    def test_the_cap_is_low_enough_to_leave_room_for_a_counter(self):
        assert MAX_UPLOAD_STEM + len("_9999.txt") < COMPONENT_LIMIT
