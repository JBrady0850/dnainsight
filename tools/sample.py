"""Choose a raw-DNA sample from uploads/ that cannot be a harness's own output.

THE FAULT THIS REPLACES
Every harness used ``sorted((ROOT / "uploads").glob("*.txt"))[0]``, the
alphabetically first file. But each harness POSTs that file back to
/api/profiles, and the app stores the upload as
``<ProfileName>_<original name>``, inside the same uploads/ directory. So after
one run the alphabetically first file was one of the harness's own artifacts,
and every subsequent run prepended another prefix:

    p1_3_mother.txt
    Collision_Test_p1_3_mother.txt
    Collision_Test_Collision_Test_p1_3_mother.txt          ... and so on

The component grew by about fifteen characters per gate run. At 246 characters
the next run crossed the filesystem's 255-byte limit for a single component and
Windows raised OSError EINVAL, which surfaced as three simultaneous stage
failures. The gate had been green only because it had not been run enough times
yet: its result depended on its own history, which means it was not a gate.

THE INVARIANT THAT FIXES IT
Selection by SHORTEST NAME is self-limiting. A derived artifact is strictly
longer than the file it was derived from, because derivation only ever prepends
a prefix. So the shortest name can never be output derived from itself, no
matter how many times the gate runs. No prefix blacklist is needed, and none is
used here, because a blacklist has to be updated every time a harness is added
and silently stops working when someone forgets.
"""

import shutil
from pathlib import Path


def _candidates(root: Path) -> list:
    return [p for p in (Path(root) / "uploads").glob("*.txt") if p.is_file()]


_FIXTURE = Path(__file__).parent.parent / "tests" / "fixtures" / "sample_23andme.txt"


def ensure_fixture(root: Path) -> Path:
    """Guarantee at least one raw-DNA sample exists in uploads/ before a harness runs.

    uploads/ is gitignored and rightly stays empty on a fresh clone, so
    end-to-end harnesses that need a real upload to POST would otherwise fail
    the moment someone clears uploads/ or clones fresh. This seeds a synthetic,
    tracked fixture (tests/fixtures/sample_23andme.txt, built from the bundled
    reference's rsIDs, containing no real person's genetic data) into uploads/
    the first time it is needed, then returns the path selection normally would.
    """
    updir = Path(root) / "uploads"
    updir.mkdir(parents=True, exist_ok=True)
    if not _candidates(root):
        shutil.copyfile(_FIXTURE, updir / _FIXTURE.name)
    return pick_sample(root)


def pick_sample(root: Path):
    """Shortest upload name, ties broken alphabetically. None if there are none."""
    files = _candidates(root)
    if not files:
        return None
    return min(files, key=lambda p: (len(p.name), p.name))


def pick_richest_sample(root: Path):
    """Largest upload, then shortest name among equals.

    Some checks want the sample carrying the most positions, for the best chance
    of a real non-carrier call and a real palindromic call. Size decides first;
    the shortest-name tie-break keeps the self-limiting property, which matters
    because every fixture in uploads/ is currently the same size.
    """
    files = _candidates(root)
    if not files:
        return None
    biggest = max(p.stat().st_size for p in files)
    return min((p for p in files if p.stat().st_size == biggest),
               key=lambda p: (len(p.name), p.name))
