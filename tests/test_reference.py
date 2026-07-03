"""Unit tests for data/build_reference.py — reference table integrity."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from data.build_reference import REFERENCE, build_reference, find_duplicates


class TestReferenceIntegrity:
    def test_no_duplicate_rsids(self):
        """Every rsID in the curated table must be unique; duplicates silently
        drop annotations at build time."""
        dupes = find_duplicates()
        assert dupes == [], f"Duplicate rsIDs in REFERENCE: {sorted(set(dupes))}"

    def test_build_count_matches_table(self):
        """With no duplicates, the built dict size equals the table length."""
        ref = build_reference()
        assert len(ref) == len(REFERENCE)

    def test_all_rows_have_five_fields(self):
        for row in REFERENCE:
            assert len(row) == 5, f"Malformed row: {row!r}"

    def test_all_rsids_well_formed(self):
        for rsid, *_ in REFERENCE:
            assert rsid.startswith("rs"), f"Bad rsID: {rsid}"

    def test_category_values_valid(self):
        valid = {"PHARM", "METAB", "INFLAM", "NEURO", "DETOX", "CARDIO"}
        for rsid, gene, category, *_ in REFERENCE:
            assert category in valid, f"{rsid}: unknown category {category}"
