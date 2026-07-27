"""Tests for backend.filters, the findings filter, sort and facet engine.

Covers docs/API_V2.md sections 3.4, 3.5 and 3.6. The two rules that break most
easily are asserted first and hardest:

  1. Frequency and publication filters must never drop a genoset, a trait or a
     polygenic score. Those entities have no position, so they have no
     frequency and no citation count, and dropping them makes a user think the
     slider is broken.
  2. A null magnitude filters and sorts as 1, not as 0.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.filters import (
    EXEMPT_ENTITIES,
    SORT_KEYS,
    UNSCORED_MAGNITUDE,
    apply_filters,
    build_facets,
    filter_and_sort,
    matches_query,
    paginate,
    parse_bool,
    parse_csv,
    parse_float,
    parse_int,
    parse_query,
    sort_findings,
    summarise,
)

# ---------------------------------------------------------------------------
# A synthetic finding set covering all four entity types. Read only: nothing in
# filters.py mutates a finding, so every test shares this list.
# ---------------------------------------------------------------------------

MTHFR = {
    "rsid": "rs1801133", "entity_type": "snp", "gene": "MTHFR",
    "chromosome": "1", "position": 11856378, "genotype": "CT",
    "token": "(C;T)", "zygosity": "heterozygous", "magnitude": 4.5,
    "magnitude_source": "computed", "repute": "Bad", "confidence": "moderate",
    "summary": "MTHFR C677T: reduced enzyme activity.",
    "interpretation": "Impairs folate conversion and homocysteine clearance.",
    "clinical_sig": "risk factor", "clinvar_sig_code": 255,
    "review_status": "criteria provided, multiple submitters, no conflicts",
    "review_stars": 2, "cpic_level": "", "pgx_level": "",
    "evidence": "ClinVar risk factor, 2 stars", "publications": 720,
    "conditions": "Homocysteinemia", "conditions_list": ["Homocysteinemia"],
    "sources": ["bundled_reference"], "freq": 44.2, "freq_band": "common",
    "gmaf": 0.24, "silo": "actionable", "category": "NEURO",
    "topics": ["MTHFR", "homocysteine", "folate"],
    "medicines": ["methotrexate", "folic acid"],
    "variant_copies": 1, "carrier": True, "count": 2,
    "labels": ["self_a.txt", "self_b.txt"], "conflict": True,
    "flipped": False, "ambiguous": False, "dubious": False,
    "coverage": None, "discovered_at": "2026-07-01T09:00:00",
}

CYP2C19 = {
    "rsid": "rs4244285", "entity_type": "snp", "gene": "CYP2C19",
    "chromosome": "10", "position": 96541616, "genotype": "AA",
    "token": "(A;A)", "zygosity": "homozygous", "magnitude": 7.8,
    "magnitude_source": "computed", "repute": "Bad", "confidence": "high",
    "summary": "CYP2C19 *2: no enzyme function.",
    "interpretation": "Clopidogrel is not activated. Poor metabolizer.",
    "clinical_sig": "drug response", "clinvar_sig_code": 6,
    "review_status": "reviewed by expert panel", "review_stars": 3,
    "cpic_level": "A", "pgx_level": "Testing Recommended",
    "evidence": "CPIC Level A", "publications": 410,
    "conditions": "Clopidogrel response", "conditions_list": ["Clopidogrel response"],
    "sources": ["bundled_reference", "cpic"], "freq": 15.0,
    "freq_band": "uncommon", "gmaf": 0.15, "silo": "pre_prescription",
    "category": "PHARM", "topics": ["clopidogrel", "antiplatelet"],
    "medicines": ["clopidogrel", "prasugrel"],
    "variant_copies": 2, "carrier": True, "count": 1,
    "labels": ["self_a.txt"], "conflict": False,
    "flipped": True, "ambiguous": False, "dubious": False,
    "coverage": None, "discovered_at": "2026-07-02T09:00:00",
}

UNSCORED = {
    "rsid": "rs9999", "entity_type": "snp", "gene": "", "chromosome": "X",
    "position": 1000, "genotype": "AT", "token": "(A;T)",
    "zygosity": "heterozygous", "magnitude": None, "magnitude_source": "",
    "repute": "", "confidence": "none", "summary": "",
    "interpretation": "", "clinical_sig": "", "clinvar_sig_code": None,
    "review_status": "", "review_stars": 0, "cpic_level": "", "pgx_level": "",
    "evidence": "", "publications": 0, "conditions": "",
    "conditions_list": [], "sources": [], "freq": None,
    "freq_band": "unknown", "gmaf": None, "silo": "informational",
    "category": "", "topics": [], "medicines": [],
    "variant_copies": 0, "carrier": False, "count": 1, "labels": ["self_a.txt"],
    "conflict": False, "flipped": False, "ambiguous": True, "dubious": True,
    "coverage": None, "discovered_at": "2026-06-01T09:00:00",
}

FACTOR_V = {
    "rsid": "rs6025", "entity_type": "snp", "gene": "F5", "chromosome": "1",
    "position": 169519049, "genotype": "CT", "token": "(C;T)",
    "zygosity": "heterozygous", "magnitude": 2.0,
    "magnitude_source": "computed", "repute": "Bad", "confidence": "high",
    "summary": "Factor V Leiden carrier.",
    "interpretation": "Increased risk of venous thrombosis.",
    "clinical_sig": "pathogenic", "clinvar_sig_code": 5,
    "review_status": "practice guideline", "review_stars": 4,
    "cpic_level": "", "pgx_level": "", "evidence": "ClinVar pathogenic, 4 stars",
    "publications": 900, "conditions": "Thrombophilia",
    "conditions_list": ["Thrombophilia"], "sources": ["bundled_reference"],
    "freq": 0.0, "freq_band": "very_rare", "gmaf": 0.02, "silo": "actionable",
    "category": "CARDIO", "topics": ["thrombosis"], "medicines": ["warfarin"],
    "variant_copies": 1, "carrier": True, "count": 2,
    "labels": ["self_a.txt", "self_b.txt"], "conflict": False,
    "flipped": False, "ambiguous": False, "dubious": False,
    "coverage": None, "discovered_at": "2026-07-03T09:00:00",
}

NO_CALL = {
    "rsid": "rs4680", "entity_type": "snp", "gene": "COMT",
    "chromosome": "22", "position": 19951271, "genotype": "NN",
    "token": "(-;-)", "zygosity": "no_call", "magnitude": 0.0,
    "magnitude_source": "computed", "repute": "", "confidence": "none",
    "summary": "COMT Val158Met.", "interpretation": "Probe failed.",
    "clinical_sig": "uncertain significance", "clinvar_sig_code": 1,
    "review_status": "criteria provided, single submitter", "review_stars": 1,
    "cpic_level": "", "pgx_level": "", "evidence": "ClinVar, 1 star",
    "publications": 12, "conditions": "", "conditions_list": [],
    "sources": ["bundled_reference"], "freq": None, "freq_band": "unknown",
    "gmaf": None, "silo": "informational", "category": "NEURO",
    "topics": ["dopamine"], "medicines": [],
    "variant_copies": None, "carrier": None, "count": 0, "labels": [],
    "conflict": False, "flipped": False, "ambiguous": False, "dubious": True,
    "coverage": None, "discovered_at": "2026-05-01T09:00:00",
}

PROTECTIVE = {
    "rsid": "rs73598374", "entity_type": "snp", "gene": "ADA",
    "chromosome": "20", "position": 43255049, "genotype": "GG",
    "token": "(G;G)", "zygosity": "homozygous", "magnitude": 1.5,
    "magnitude_source": "computed", "repute": "Good", "confidence": "low",
    "summary": "Reference genotype, normal adenosine deaminase.",
    "interpretation": "Associated with reduced risk of poor sleep quality.",
    "clinical_sig": "likely benign", "clinvar_sig_code": 3,
    "review_status": "criteria provided, single submitter", "review_stars": 1,
    "cpic_level": "C", "pgx_level": "Informative PGx",
    "evidence": "ClinVar likely benign, 1 star", "publications": 8,
    "conditions": "Sleep quality", "conditions_list": ["Sleep quality"],
    "sources": ["bundled_reference"], "freq": 3.0, "freq_band": "rare",
    "gmaf": 0.03, "silo": "informational", "category": "NEURO",
    "topics": ["sleep"], "medicines": [],
    "variant_copies": 2, "carrier": True, "count": 1, "labels": ["self_a.txt"],
    "conflict": False, "flipped": False, "ambiguous": False, "dubious": False,
    "coverage": None, "discovered_at": "2026-04-01T09:00:00",
}

GENOSET = {
    "rsid": "dgs001", "entity_type": "genoset", "gene": "", "chromosome": "",
    "position": 0, "genotype": "", "token": "", "zygosity": "",
    "magnitude": 4.0, "magnitude_source": "computed", "repute": "",
    "confidence": "low", "summary": "You carry two copies of APOE e2.",
    "interpretation": "Lowest genetic Alzheimer risk of any APOE pair.",
    "clinical_sig": "", "clinvar_sig_code": None, "review_status": "",
    "review_stars": 0, "cpic_level": "", "pgx_level": "",
    "evidence": "genoset rule", "publications": 0, "conditions": "",
    "conditions_list": [], "sources": ["genoset"], "freq": None,
    "freq_band": "", "gmaf": None, "silo": "informational",
    "category": "LIPID", "topics": ["cholesterol"], "medicines": [],
    "variant_copies": None, "carrier": None, "count": 0, "labels": [],
    "conflict": False, "flipped": False, "ambiguous": False, "dubious": False,
    "criteria": "and(rs429358(T;T), rs7412(T;T))",
    "matched_rsids": ["rs429358", "rs7412"], "aka": "APOE e2/e2",
    "coverage": 1.0, "discovered_at": "2026-07-04T09:00:00",
}

TRAIT = {
    "rsid": "lactase_persistence", "entity_type": "trait", "gene": "LCT",
    "chromosome": "", "position": 0, "genotype": "", "token": "",
    "zygosity": "", "magnitude": 1.0, "magnitude_source": "computed",
    "repute": "", "confidence": "low",
    "summary": "Lactase persistent.", "interpretation": "Milk is digestible.",
    "clinical_sig": "", "clinvar_sig_code": None, "review_status": "",
    "review_stars": 0, "cpic_level": "", "pgx_level": "", "evidence": "trait",
    "publications": 0, "conditions": "", "conditions_list": [],
    "sources": ["trait"], "freq": None, "freq_band": "", "gmaf": None,
    "silo": "informational", "category": "Diet", "topics": ["diet"],
    "medicines": [], "variant_copies": None, "carrier": None, "count": 0,
    "labels": [], "conflict": False, "flipped": False, "ambiguous": False,
    "dubious": False, "coverage": 1.0, "caveats": [],
    "discovered_at": "2026-07-05T09:00:00",
}

PRS = {
    "rsid": "t2d", "entity_type": "prs", "gene": "", "chromosome": "",
    "position": 0, "genotype": "", "token": "", "zygosity": "",
    "magnitude": 3.0, "magnitude_source": "computed", "repute": "",
    "confidence": "low", "summary": "Type 2 diabetes polygenic score.",
    "interpretation": "62nd percentile against the reference panel.",
    "clinical_sig": "", "clinvar_sig_code": None, "review_status": "",
    "review_stars": 0, "cpic_level": "", "pgx_level": "", "evidence": "prs",
    "publications": 0, "conditions": "Type 2 diabetes",
    "conditions_list": ["Type 2 diabetes"], "sources": ["prs"], "freq": None,
    "freq_band": "", "gmaf": None, "silo": "informational",
    "category": "METAB", "topics": ["diabetes"], "medicines": ["metformin"],
    "variant_copies": None, "carrier": None, "count": 0, "labels": [],
    "conflict": False, "flipped": False, "ambiguous": False, "dubious": False,
    "coverage": 0.95, "percentile": 62.0, "band": "moderate",
    "reliable": True, "discovered_at": "2026-07-06T09:00:00",
}

SAMPLE = [MTHFR, CYP2C19, UNSCORED, FACTOR_V, NO_CALL, PROTECTIVE,
          GENOSET, TRAIT, PRS]

SNP_IDS = {"rs1801133", "rs4244285", "rs9999", "rs6025", "rs4680", "rs73598374"}
EXEMPT_IDS = {"dgs001", "lactase_persistence", "t2d"}

# A null magnitude, an explicit zero and a real score, for the sort convention.
MAG_PROBE = [
    {"rsid": "mag_null", "entity_type": "snp", "magnitude": None},
    {"rsid": "mag_zero", "entity_type": "snp", "magnitude": 0.0},
    {"rsid": "mag_two", "entity_type": "snp", "magnitude": 2.0},
]


def ids(rows):
    """rsIDs of a result set, in order."""
    return [row["rsid"] for row in rows]


def id_set(rows):
    """rsIDs of a result set, unordered."""
    return {row["rsid"] for row in rows}


class TestExemptEntities:
    """Rule 1. The frequency and publication filters must not touch these."""

    def test_exempt_entities_is_exactly_the_three_documented_types(self):
        assert EXEMPT_ENTITIES == {"genoset", "trait", "prs"}

    def test_a_snp_is_not_exempt(self):
        assert "snp" not in EXEMPT_ENTITIES

    def test_the_fixture_covers_all_four_entity_types(self):
        assert {f["entity_type"] for f in SAMPLE} == {"snp", "genoset",
                                                      "trait", "prs"}

    def test_min_publications_keeps_every_exempt_entity(self):
        for floor in (1, 5, 50, 500, 100000):
            result = apply_filters(SAMPLE, {"min_publications": floor})
            assert EXEMPT_IDS <= id_set(result)

    def test_min_publications_still_filters_snps(self):
        result = apply_filters(SAMPLE, {"min_publications": 1})
        assert "rs9999" not in id_set(result)

    def test_an_impossible_publication_floor_leaves_only_exempt_entities(self):
        result = apply_filters(SAMPLE, {"min_publications": 1000})
        assert id_set(result) == EXEMPT_IDS

    def test_max_publications_keeps_every_exempt_entity(self):
        for ceiling in (0, 1, 10, 500):
            result = apply_filters(SAMPLE, {"max_publications": ceiling})
            assert EXEMPT_IDS <= id_set(result)

    def test_max_publications_still_filters_snps(self):
        result = apply_filters(SAMPLE, {"max_publications": 0})
        assert id_set(result) == EXEMPT_IDS | {"rs9999"}

    def test_min_freq_keeps_every_exempt_entity(self):
        for floor in (0.1, 1, 10, 50, 99.9):
            result = apply_filters(SAMPLE, {"min_freq": floor})
            assert EXEMPT_IDS <= id_set(result)

    def test_min_freq_still_filters_snps(self):
        result = apply_filters(SAMPLE, {"min_freq": 10})
        assert id_set(result) == EXEMPT_IDS | {"rs1801133", "rs4244285"}

    def test_an_impossible_frequency_floor_leaves_only_exempt_entities(self):
        result = apply_filters(SAMPLE, {"min_freq": 99.9})
        assert id_set(result) == EXEMPT_IDS

    def test_max_freq_keeps_every_exempt_entity(self):
        for ceiling in (0.0, 1, 10, 50):
            result = apply_filters(SAMPLE, {"max_freq": ceiling})
            assert EXEMPT_IDS <= id_set(result)

    def test_max_freq_still_filters_snps(self):
        result = apply_filters(SAMPLE, {"max_freq": 10})
        assert "rs1801133" not in id_set(result)
        assert "rs4244285" not in id_set(result)

    def test_require_frequency_keeps_every_exempt_entity(self):
        result = apply_filters(SAMPLE, {"require_frequency": True})
        assert EXEMPT_IDS <= id_set(result)

    def test_both_filters_together_keep_every_exempt_entity(self):
        result = apply_filters(SAMPLE, {"min_freq": 20, "min_publications": 50,
                                        "require_frequency": True})
        assert EXEMPT_IDS <= id_set(result)

    def test_nudging_the_sliders_never_empties_the_exempt_entities(self):
        for floor in (0, 1, 2, 5, 25, 75, 100):
            params = {"min_freq": floor, "min_publications": floor}
            result = apply_filters(SAMPLE, params)
            assert len(id_set(result) & EXEMPT_IDS) == 3

    def test_an_exempt_entity_with_no_frequency_key_is_still_kept(self):
        rows = [{"rsid": "dgs002", "entity_type": "genoset"},
                {"rsid": "rs1", "entity_type": "snp"}]
        result = apply_filters(rows, {"min_freq": 5, "require_frequency": True})
        assert ids(result) == ["dgs002"]


class TestNullMagnitudeSortsAsOne:
    """Rule 2. An unscored variant is 1, not 0."""

    def test_the_unscored_constant_is_one(self):
        assert UNSCORED_MAGNITUDE == 1.0

    def test_min_magnitude_one_keeps_a_null_magnitude(self):
        result = apply_filters(MAG_PROBE, {"min_magnitude": 1})
        assert "mag_null" in id_set(result)

    def test_min_magnitude_one_point_one_drops_a_null_magnitude(self):
        result = apply_filters(MAG_PROBE, {"min_magnitude": 1.1})
        assert ids(result) == ["mag_two"]

    def test_min_magnitude_one_drops_an_explicit_zero(self):
        result = apply_filters(MAG_PROBE, {"min_magnitude": 1})
        assert "mag_zero" not in id_set(result)

    def test_a_null_magnitude_is_not_treated_as_zero(self):
        kept = id_set(apply_filters(MAG_PROBE, {"min_magnitude": 0.5}))
        assert "mag_null" in kept
        assert "mag_zero" not in kept

    def test_a_null_magnitude_is_not_treated_as_ten(self):
        result = apply_filters(MAG_PROBE, {"max_magnitude": 1.0})
        assert "mag_null" in id_set(result)

    def test_max_magnitude_below_one_drops_a_null_magnitude(self):
        result = apply_filters(MAG_PROBE, {"max_magnitude": 0.9})
        assert ids(result) == ["mag_zero"]

    def test_a_null_magnitude_sits_between_two_and_zero_descending(self):
        assert ids(sort_findings(MAG_PROBE, "magnitude", "desc")) == [
            "mag_two", "mag_null", "mag_zero"]

    def test_a_null_magnitude_sits_between_zero_and_two_ascending(self):
        assert ids(sort_findings(MAG_PROBE, "magnitude", "asc")) == [
            "mag_zero", "mag_null", "mag_two"]

    def test_a_missing_magnitude_key_behaves_like_a_null_one(self):
        rows = [{"rsid": "absent", "entity_type": "snp"}] + MAG_PROBE
        ordered = ids(sort_findings(rows, "magnitude", "desc"))
        assert ordered.index("absent") < ordered.index("mag_zero")

    def test_junk_magnitude_is_treated_as_unscored(self):
        rows = [{"rsid": "junk", "magnitude": "high"}]
        assert ids(apply_filters(rows, {"min_magnitude": 1})) == ["junk"]

    def test_a_boolean_magnitude_is_treated_as_unscored(self):
        rows = [{"rsid": "flagged", "magnitude": True}]
        assert ids(apply_filters(rows, {"min_magnitude": 1})) == ["flagged"]

    def test_the_fixture_unscored_snp_survives_a_floor_of_one(self):
        result = apply_filters(SAMPLE, {"min_magnitude": 1})
        assert "rs9999" in id_set(result)

    def test_the_fixture_no_call_is_dropped_by_a_floor_of_one(self):
        result = apply_filters(SAMPLE, {"min_magnitude": 1})
        assert "rs4680" not in id_set(result)


class TestParseCsv:
    def test_a_plain_list_is_split(self):
        assert parse_csv("a,b,c") == ["a", "b", "c"]

    def test_whitespace_is_stripped(self):
        assert parse_csv(" a , b ") == ["a", "b"]

    def test_empty_items_are_dropped(self):
        assert parse_csv("a,,b,") == ["a", "b"]

    def test_none_gives_an_empty_list(self):
        assert parse_csv(None) == []

    def test_an_empty_string_gives_an_empty_list(self):
        assert parse_csv("") == []

    def test_a_list_input_is_accepted(self):
        assert parse_csv(["a", "b"]) == ["a", "b"]

    def test_a_tuple_input_is_accepted(self):
        assert parse_csv(("a", "b")) == ["a", "b"]

    def test_numbers_are_stringified(self):
        assert parse_csv([5, 4]) == ["5", "4"]


class TestParseBool:
    def test_every_truthy_spelling(self):
        for text in ("1", "true", "TRUE", "yes", "Yes", "on", "y", "t"):
            assert parse_bool(text) is True

    def test_falsy_spellings(self):
        for text in ("0", "false", "no", "off", "n", "f", "maybe"):
            assert parse_bool(text) is False

    def test_none_returns_the_default(self):
        assert parse_bool(None) is False
        assert parse_bool(None, True) is True

    def test_an_empty_string_returns_the_default(self):
        assert parse_bool("") is False
        assert parse_bool("", True) is True

    def test_a_real_bool_passes_through(self):
        assert parse_bool(True) is True
        assert parse_bool(False) is False

    def test_surrounding_whitespace_is_tolerated(self):
        assert parse_bool("  yes  ") is True


class TestParseFloat:
    def test_a_number_is_parsed(self):
        assert parse_float("3.5") == 3.5

    def test_an_integer_string_is_parsed(self):
        assert parse_float("4") == 4.0

    def test_zero_is_preserved(self):
        assert parse_float("0") == 0.0

    def test_junk_returns_the_default(self):
        assert parse_float("banana") is None
        assert parse_float("banana", 1.0) == 1.0

    def test_none_returns_the_default(self):
        assert parse_float(None) is None
        assert parse_float(None, 2.0) == 2.0

    def test_an_empty_string_returns_the_default(self):
        assert parse_float("", 7.0) == 7.0

    def test_a_negative_number_is_parsed(self):
        assert parse_float("-2.5") == -2.5


class TestParseInt:
    def test_an_integer_is_parsed(self):
        assert parse_int("42") == 42

    def test_a_float_string_is_truncated(self):
        assert parse_int("3.7") == 3

    def test_zero_is_preserved(self):
        assert parse_int("0") == 0

    def test_junk_returns_the_default(self):
        assert parse_int("banana") is None
        assert parse_int("banana", 5) == 5

    def test_none_returns_the_default(self):
        assert parse_int(None) is None
        assert parse_int(None, 200) == 200

    def test_an_empty_string_returns_the_default(self):
        assert parse_int("", 200) == 200

    def test_a_negative_number_is_parsed(self):
        assert parse_int("-3") == -3


class TestParseQueryText:
    def test_the_result_carries_the_five_documented_keys(self):
        parsed = parse_query("MTHFR")
        assert set(parsed) == {"text", "regex", "region", "ops", "flags"}

    def test_bare_text_is_kept(self):
        assert parse_query("MTHFR")["text"] == "MTHFR"

    def test_bare_text_is_compiled_to_a_regex(self):
        parsed = parse_query("MTH.R")
        assert parsed["regex"] is not None
        assert parsed["regex"].search("mthfr")

    def test_the_regex_is_case_insensitive(self):
        assert parse_query("mthfr")["regex"].search("MTHFR")

    def test_an_invalid_regex_degrades_to_a_literal_match(self):
        parsed = parse_query("[unbalanced(")
        assert parsed["regex"] is not None
        assert parsed["regex"].search("a [unbalanced( thing")

    def test_an_invalid_regex_does_not_raise(self):
        for bad in ("[unbalanced(", "*", "(", "a{2,1}", "[z-a]"):
            assert parse_query(bad)["regex"] is not None

    def test_an_invalid_regex_keeps_the_original_text(self):
        assert parse_query("[unbalanced(")["text"] == "[unbalanced("

    def test_empty_input_gives_the_empty_shape(self):
        parsed = parse_query("")
        assert parsed["text"] == ""
        assert parsed["regex"] is None
        assert parsed["region"] is None
        assert parsed["ops"] == {}
        assert parsed["flags"] == set()

    def test_none_input_gives_the_empty_shape(self):
        assert parse_query(None)["regex"] is None

    def test_whitespace_only_input_gives_the_empty_shape(self):
        assert parse_query("   ")["text"] == ""


class TestParseQueryRegion:
    def test_a_bare_chromosome(self):
        assert parse_query("chr7")["region"] == {
            "chromosome": "7", "start": None, "end": None}

    def test_an_exact_position(self):
        assert parse_query("chr7:1234")["region"] == {
            "chromosome": "7", "start": 1234, "end": 1234}

    def test_an_inclusive_range(self):
        assert parse_query("chr7:1000-2000")["region"] == {
            "chromosome": "7", "start": 1000, "end": 2000}

    def test_a_position_plus_an_offset(self):
        assert parse_query("chr7:1000+500")["region"] == {
            "chromosome": "7", "start": 1000, "end": 1500}

    def test_the_chro_spelling(self):
        assert parse_query("chro7")["region"]["chromosome"] == "7"

    def test_the_chrom_spelling(self):
        assert parse_query("chrom7")["region"]["chromosome"] == "7"

    def test_the_chromosome_spelling(self):
        assert parse_query("chromosome7")["region"]["chromosome"] == "7"

    def test_chromosome_x(self):
        assert parse_query("chrX")["region"]["chromosome"] == "X"

    def test_chromosome_y(self):
        assert parse_query("chrY")["region"]["chromosome"] == "Y"

    def test_chromosome_m(self):
        assert parse_query("chrM")["region"]["chromosome"] == "M"

    def test_chromosome_mt(self):
        assert parse_query("chrMT")["region"]["chromosome"] == "MT"

    def test_a_lowercase_chromosome_is_upper_cased(self):
        assert parse_query("chrx")["region"]["chromosome"] == "X"

    def test_a_two_digit_chromosome(self):
        assert parse_query("chr22:1-99")["region"]["chromosome"] == "22"

    def test_a_region_clears_the_free_text(self):
        assert parse_query("chr7:1000-2000")["text"] == ""

    def test_a_region_with_trailing_text_is_not_a_region(self):
        parsed = parse_query("chr7:1000-2000 extra")
        assert parsed["region"] is None
        assert parsed["text"] == "chr7:1000-2000 extra"

    def test_a_whitespace_separated_chromosome_is_accepted(self):
        assert parse_query("chr 7")["region"]["chromosome"] == "7"


class TestParseQueryOperators:
    def test_a_clinvar_significance_whitelist(self):
        assert parse_query("/CLNSIG=5,4")["ops"] == {"clnsig": ("=", ["5", "4"])}

    def test_a_single_clinvar_code(self):
        assert parse_query("/CLNSIG=6")["ops"] == {"clnsig": ("=", ["6"])}

    def test_a_minimum_star_count(self):
        assert parse_query("/STARS>=2")["ops"] == {"stars": (">=", ["2"])}

    def test_a_minimum_magnitude(self):
        assert parse_query("/MAG>=3")["ops"] == {"mag": (">=", ["3"])}

    def test_a_minimum_pooled_call_count(self):
        assert parse_query("/COUNT>=2")["ops"] == {"count": (">=", ["2"])}

    def test_a_frequency_operator(self):
        assert parse_query("/FREQ<=10")["ops"] == {"freq": ("<=", ["10"])}

    def test_operators_are_case_insensitive(self):
        assert parse_query("/mag>=3")["ops"] == {"mag": (">=", ["3"])}

    def test_every_comparison_operator_parses(self):
        for operator in (">=", "<=", ">", "<", "="):
            parsed = parse_query("/MAG%s3" % operator)
            assert parsed["ops"]["mag"][0] == operator

    def test_whitespace_around_the_operator_is_tolerated(self):
        assert parse_query("/MAG >= 3")["ops"]["mag"] == (">=", ["3"])

    def test_a_decimal_value_is_kept(self):
        assert parse_query("/MAG>=3.5")["ops"]["mag"] == (">=", ["3.5"])

    def test_operators_are_stripped_from_the_text(self):
        assert parse_query("/MAG>=3")["text"] == ""

    def test_several_operators_combine(self):
        parsed = parse_query("/STARS>=2 /MAG>=3")
        assert parsed["ops"]["stars"] == (">=", ["2"])
        assert parsed["ops"]["mag"] == (">=", ["3"])

    def test_an_operator_leaves_the_surrounding_text_searchable(self):
        parsed = parse_query("MTHFR /MAG>=3")
        assert parsed["text"] == "MTHFR"
        assert parsed["ops"]["mag"] == (">=", ["3"])
        assert parsed["regex"].search("mthfr")

    def test_text_after_an_operator_is_kept(self):
        parsed = parse_query("/STARS>=2 warfarin")
        assert parsed["text"] == "warfarin"


class TestParseQueryFlags:
    def test_the_dubious_flag(self):
        assert parse_query("/dubious")["flags"] == {"dubious"}

    def test_the_flipped_flag(self):
        assert parse_query("/flipped")["flags"] == {"flipped"}

    def test_the_ambiguous_flag(self):
        assert parse_query("/ambiguous")["flags"] == {"ambiguous"}

    def test_the_conflict_flag(self):
        assert parse_query("/conflict")["flags"] == {"conflict"}

    def test_the_carrier_flag(self):
        assert parse_query("/carrier")["flags"] == {"carrier"}

    def test_the_nocall_flag(self):
        assert parse_query("/nocall")["flags"] == {"nocall"}

    def test_flags_are_case_insensitive(self):
        assert parse_query("/DUBIOUS")["flags"] == {"dubious"}

    def test_flags_are_stripped_from_the_text(self):
        assert parse_query("/dubious")["text"] == ""

    def test_several_flags_combine(self):
        assert parse_query("/dubious /flipped")["flags"] == {"dubious", "flipped"}

    def test_flags_operators_and_text_combine(self):
        parsed = parse_query("warfarin /STARS>=2 /nocall")
        assert parsed["text"] == "warfarin"
        assert parsed["flags"] == {"nocall"}
        assert parsed["ops"]["stars"] == (">=", ["2"])

    def test_an_unknown_flag_is_left_as_text(self):
        parsed = parse_query("/banana")
        assert parsed["flags"] == set()
        assert parsed["text"] == "/banana"


def matching(query, rows=None):
    """rsIDs of the rows that satisfy a free-text query."""
    parsed = parse_query(query)
    return [r["rsid"] for r in (rows if rows is not None else SAMPLE)
            if matches_query(r, parsed)]


class TestMatchesQueryFlags:
    def test_an_empty_query_matches_everything(self):
        assert len(matching("")) == len(SAMPLE)

    def test_the_dubious_flag_keeps_only_suspect_calls(self):
        assert set(matching("/dubious")) == {"rs9999", "rs4680"}

    def test_the_flipped_flag_keeps_only_flipped_calls(self):
        assert matching("/flipped") == ["rs4244285"]

    def test_the_ambiguous_flag_keeps_only_palindromic_calls(self):
        assert matching("/ambiguous") == ["rs9999"]

    def test_the_ambiguous_flag_also_accepts_freq_ambiguous(self):
        rows = [{"rsid": "rs1", "freq_ambiguous": True},
                {"rsid": "rs2", "freq_ambiguous": False}]
        assert matching("/ambiguous", rows) == ["rs1"]

    def test_the_conflict_flag_keeps_only_pooled_disagreements(self):
        assert matching("/conflict") == ["rs1801133"]

    def test_the_carrier_flag_keeps_only_carriers(self):
        assert set(matching("/carrier")) == {"rs1801133", "rs4244285",
                                             "rs6025", "rs73598374"}

    def test_the_carrier_flag_rejects_an_unknown_carrier_state(self):
        rows = [{"rsid": "rs1", "carrier": None}, {"rsid": "rs2", "carrier": 1}]
        assert matching("/carrier", rows) == []

    def test_the_nocall_flag_keeps_only_failed_probes(self):
        assert matching("/nocall") == ["rs4680"]

    def test_two_flags_are_combined_with_and(self):
        assert matching("/dubious /nocall") == ["rs4680"]


class TestMatchesQueryOperators:
    def test_a_clinvar_whitelist(self):
        assert set(matching("/CLNSIG=5,4")) == {"rs6025"}

    def test_a_single_clinvar_code(self):
        assert matching("/CLNSIG=6") == ["rs4244285"]

    def test_a_clinvar_whitelist_drops_findings_with_no_record(self):
        assert "rs9999" not in matching("/CLNSIG=1,5")

    def test_a_minimum_star_count(self):
        assert set(matching("/STARS>=3")) == {"rs4244285", "rs6025"}

    def test_a_maximum_star_count(self):
        assert "rs6025" not in matching("/STARS<=2")

    def test_a_minimum_magnitude(self):
        assert set(matching("/MAG>=4")) == {"rs1801133", "rs4244285", "dgs001"}

    def test_a_null_magnitude_is_compared_as_one(self):
        assert "mag_null" in matching("/MAG>=1", MAG_PROBE)

    def test_a_null_magnitude_fails_a_higher_floor(self):
        assert "mag_null" not in matching("/MAG>1", MAG_PROBE)

    def test_a_minimum_pooled_count(self):
        assert set(matching("/COUNT>=2")) == {"rs1801133", "rs6025"}

    def test_the_greater_than_operator_is_strict(self):
        rows = [{"rsid": "rs1", "review_stars": 2},
                {"rsid": "rs2", "review_stars": 3}]
        assert matching("/STARS>2", rows) == ["rs2"]

    def test_the_less_than_operator_is_strict(self):
        rows = [{"rsid": "rs1", "review_stars": 2},
                {"rsid": "rs2", "review_stars": 1}]
        assert matching("/STARS<2", rows) == ["rs2"]

    def test_the_equals_operator_is_exact(self):
        rows = [{"rsid": "rs1", "review_stars": 2},
                {"rsid": "rs2", "review_stars": 3}]
        assert matching("/STARS=2", rows) == ["rs1"]

    def test_the_equals_operator_tolerates_float_noise(self):
        rows = [{"rsid": "rs1", "magnitude": 3.0000000001}]
        assert matching("/MAG=3", rows) == ["rs1"]

    def test_a_missing_field_fails_the_comparison(self):
        rows = [{"rsid": "rs1"}]
        assert matching("/COUNT>=1", rows) == []


class TestMatchesQueryRegion:
    def test_a_bare_chromosome_keeps_that_chromosome(self):
        assert set(matching("chr1")) == {"rs1801133", "rs6025"}

    def test_a_bare_chromosome_drops_the_others(self):
        assert "rs4244285" not in matching("chr1")

    def test_chromosome_x_is_matched(self):
        assert matching("chrX") == ["rs9999"]

    def test_an_exact_position_is_matched(self):
        assert matching("chr1:11856378") == ["rs1801133"]

    def test_a_position_one_base_away_is_not_matched(self):
        assert matching("chr1:11856379") == []

    def test_an_inclusive_range_is_matched(self):
        assert matching("chr1:11856378-169519049") == ["rs1801133", "rs6025"]

    def test_a_range_boundary_is_inclusive(self):
        rows = [{"rsid": "rs1", "chromosome": "7", "position": 1000},
                {"rsid": "rs2", "chromosome": "7", "position": 2000},
                {"rsid": "rs3", "chromosome": "7", "position": 2001}]
        assert matching("chr7:1000-2000", rows) == ["rs1", "rs2"]

    def test_a_position_plus_offset_range_is_matched(self):
        rows = [{"rsid": "rs1", "chromosome": "7", "position": 1400},
                {"rsid": "rs2", "chromosome": "7", "position": 1600}]
        assert matching("chr7:1000+500", rows) == ["rs1"]

    def test_a_chr_prefix_on_the_finding_is_tolerated(self):
        rows = [{"rsid": "rs1", "chromosome": "chr7", "position": 10}]
        assert matching("chr7", rows) == ["rs1"]

    def test_a_lowercase_chromosome_on_the_finding_is_tolerated(self):
        rows = [{"rsid": "rs1", "chromosome": "x", "position": 10}]
        assert matching("chrX", rows) == ["rs1"]

    def test_a_finding_with_no_position_fails_a_positional_query(self):
        rows = [{"rsid": "rs1", "chromosome": "7"}]
        assert matching("chr7:100", rows) == []

    def test_an_entity_with_no_chromosome_is_dropped_by_a_region(self):
        assert "dgs001" not in matching("chr1")


class TestMatchesQueryText:
    def test_a_gene_symbol_is_searched(self):
        assert matching("MTHFR") == ["rs1801133"]

    def test_the_search_is_case_insensitive(self):
        assert matching("mthfr") == ["rs1801133"]

    def test_an_rsid_is_searched(self):
        assert matching("rs6025") == ["rs6025"]

    def test_the_summary_is_searched(self):
        assert matching("Factor V Leiden") == ["rs6025"]

    def test_the_interpretation_is_searched(self):
        assert matching("venous thrombosis") == ["rs6025"]

    def test_the_conditions_string_is_searched(self):
        assert matching("Homocysteinemia") == ["rs1801133"]

    def test_the_genotype_token_is_searched(self):
        assert matching("A;A") == ["rs4244285"]

    def test_the_clinical_significance_is_searched(self):
        assert set(matching("drug response")) == {"rs4244285"}

    def test_the_evidence_label_is_searched(self):
        assert matching("CPIC Level A") == ["rs4244285"]

    def test_the_topics_list_is_searched(self):
        assert matching("folate") == ["rs1801133"]

    def test_the_medicines_list_is_searched(self):
        assert matching("clopidogrel") == ["rs4244285"]

    def test_the_conditions_list_is_searched(self):
        assert matching("Type 2 diabetes") == ["t2d"]

    def test_a_genoset_criteria_string_is_searched(self):
        assert matching("rs7412") == ["dgs001"]

    def test_a_genoset_alias_is_searched(self):
        assert matching("APOE e2") == ["dgs001"]

    def test_a_regex_alternation_is_honoured(self):
        assert set(matching("MTHFR|CYP2C19")) == {"rs1801133", "rs4244285"}

    def test_an_unmatched_search_returns_nothing(self):
        assert matching("zzzznotpresent") == []

    def test_an_invalid_regex_matches_literally(self):
        rows = [{"rsid": "rs1", "gene": "A[B"}, {"rsid": "rs2", "gene": "AB"}]
        assert matching("A[B", rows) == ["rs1"]


class TestApplyFiltersNoParams:
    def test_no_parameters_filters_nothing(self):
        assert len(apply_filters(SAMPLE, {})) == len(SAMPLE)

    def test_an_empty_string_parameter_filters_nothing(self):
        params = {"silo": "", "gene": "", "entity_type": "", "repute": "",
                  "zygosity": "", "topic": "", "medicine": "", "condition": ""}
        assert len(apply_filters(SAMPLE, params)) == len(SAMPLE)

    def test_a_none_parameter_filters_nothing(self):
        params = {"min_magnitude": None, "max_magnitude": None,
                  "min_publications": None, "min_freq": None,
                  "max_freq": None, "min_stars": None, "q": None}
        assert len(apply_filters(SAMPLE, params)) == len(SAMPLE)

    def test_an_unrecognised_parameter_is_ignored(self):
        assert len(apply_filters(SAMPLE, {"banana": "yes"})) == len(SAMPLE)

    def test_the_population_parameter_is_not_a_filter(self):
        assert len(apply_filters(SAMPLE, {"population": "YRI"})) == len(SAMPLE)

    def test_an_empty_input_list_returns_an_empty_list(self):
        assert apply_filters([], {"min_magnitude": 1}) == []

    def test_a_none_input_list_returns_an_empty_list(self):
        assert apply_filters(None, {}) == []

    def test_the_returned_objects_are_the_same_dicts(self):
        result = apply_filters(SAMPLE, {"gene": "MTHFR"})
        assert result[0] is MTHFR


class TestApplyFiltersDocumentedParams:
    def test_silo(self):
        result = apply_filters(SAMPLE, {"silo": "actionable"})
        assert set(ids(result)) == {"rs1801133", "rs6025"}

    def test_an_unknown_silo_matches_nothing(self):
        assert apply_filters(SAMPLE, {"silo": "nowhere"}) == []

    def test_entity_type_single(self):
        assert ids(apply_filters(SAMPLE, {"entity_type": "genoset"})) == ["dgs001"]

    def test_entity_type_csv(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp,prs"})
        assert set(ids(result)) == SNP_IDS | {"t2d"}

    def test_min_magnitude(self):
        result = apply_filters(SAMPLE, {"min_magnitude": 4})
        assert set(ids(result)) == {"rs1801133", "rs4244285", "dgs001"}

    def test_max_magnitude(self):
        result = apply_filters(SAMPLE, {"max_magnitude": 1.0})
        assert set(ids(result)) == {"rs9999", "rs4680", "lactase_persistence"}

    def test_a_magnitude_window(self):
        result = apply_filters(SAMPLE, {"min_magnitude": 2,
                                       "max_magnitude": 4.5})
        assert set(ids(result)) == {"rs1801133", "rs6025", "dgs001", "t2d"}

    def test_repute_bad(self):
        result = apply_filters(SAMPLE, {"repute": "bad"})
        assert set(ids(result)) == {"rs1801133", "rs4244285", "rs6025"}

    def test_repute_good(self):
        assert ids(apply_filters(SAMPLE, {"repute": "good"})) == ["rs73598374"]

    def test_repute_unset(self):
        result = apply_filters(SAMPLE, {"repute": "unset"})
        assert set(ids(result)) == {"rs9999", "rs4680"} | EXEMPT_IDS

    def test_repute_csv_of_two_states(self):
        result = apply_filters(SAMPLE, {"repute": "good,bad"})
        assert len(result) == 4

    def test_repute_is_case_insensitive(self):
        assert len(apply_filters(SAMPLE, {"repute": "BAD"})) == 3

    def test_min_publications_on_snps(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp",
                                        "min_publications": 500})
        assert set(ids(result)) == {"rs1801133", "rs6025"}

    def test_max_publications_on_snps(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp",
                                        "max_publications": 20})
        assert set(ids(result)) == {"rs9999", "rs4680", "rs73598374"}

    def test_min_freq_on_snps(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "min_freq": 3})
        assert set(ids(result)) == {"rs1801133", "rs4244285", "rs73598374"}

    def test_max_freq_on_snps(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "max_freq": 20})
        assert "rs1801133" not in set(ids(result))
        assert "rs4244285" in set(ids(result))

    def test_a_frequency_window(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "min_freq": 1,
                                        "max_freq": 20})
        assert set(ids(result)) == {"rs4244285", "rs73598374"}

    def test_require_frequency_drops_a_null_frequency(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp",
                                        "require_frequency": True})
        assert "rs9999" not in set(ids(result))
        assert "rs4680" not in set(ids(result))

    def test_require_frequency_keeps_an_observed_zero(self):
        # 0.0 means "looked for and not seen"; null means "never looked".
        # They are different facts and must filter differently.
        result = apply_filters(SAMPLE, {"entity_type": "snp",
                                        "require_frequency": True})
        assert "rs6025" in set(ids(result))

    def test_require_frequency_defaults_to_off(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp"})
        assert "rs9999" in set(ids(result))

    def test_a_null_frequency_survives_a_floor_of_zero(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "min_freq": 0})
        assert "rs9999" in set(ids(result))

    def test_a_null_frequency_is_dropped_by_a_positive_floor(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "min_freq": 0.1})
        assert "rs9999" not in set(ids(result))

    def test_an_observed_zero_survives_a_floor_of_zero(self):
        result = apply_filters(SAMPLE, {"entity_type": "snp", "min_freq": 0})
        assert "rs6025" in set(ids(result))

    def test_clinvar_only_defaults_to_pathogenic_and_likely_pathogenic(self):
        assert ids(apply_filters(SAMPLE, {"clinvar_only": True})) == ["rs6025"]

    def test_clinvar_only_with_an_explicit_list_overrides_the_default(self):
        result = apply_filters(SAMPLE, {"clinvar_only": True,
                                        "clinvar_sig": "6"})
        assert ids(result) == ["rs4244285"]

    def test_clinvar_sig_without_clinvar_only(self):
        result = apply_filters(SAMPLE, {"clinvar_sig": "255,6"})
        assert set(ids(result)) == {"rs1801133", "rs4244285"}

    def test_clinvar_sig_drops_findings_with_no_record(self):
        result = apply_filters(SAMPLE, {"clinvar_sig": "1"})
        assert ids(result) == ["rs4680"]

    def test_clinvar_only_drops_every_finding_with_no_record(self):
        result = apply_filters(SAMPLE, {"clinvar_only": True,
                                        "clinvar_sig": "1,3,5,6,255"})
        assert not (EXEMPT_IDS & set(ids(result)))

    def test_min_stars(self):
        result = apply_filters(SAMPLE, {"min_stars": 3})
        assert set(ids(result)) == {"rs4244285", "rs6025"}

    def test_min_stars_of_zero_keeps_everything(self):
        assert len(apply_filters(SAMPLE, {"min_stars": 0})) == len(SAMPLE)

    def test_min_stars_of_four_keeps_only_a_practice_guideline(self):
        assert ids(apply_filters(SAMPLE, {"min_stars": 4})) == ["rs6025"]

    def test_gene(self):
        assert ids(apply_filters(SAMPLE, {"gene": "MTHFR"})) == ["rs1801133"]

    def test_gene_is_case_insensitive(self):
        assert ids(apply_filters(SAMPLE, {"gene": "mthfr"})) == ["rs1801133"]

    def test_gene_csv(self):
        result = apply_filters(SAMPLE, {"gene": "MTHFR,F5"})
        assert set(ids(result)) == {"rs1801133", "rs6025"}

    def test_topic(self):
        assert ids(apply_filters(SAMPLE, {"topic": "folate"})) == ["rs1801133"]

    def test_topic_is_case_insensitive(self):
        assert ids(apply_filters(SAMPLE, {"topic": "FOLATE"})) == ["rs1801133"]

    def test_topic_csv_matches_any(self):
        result = apply_filters(SAMPLE, {"topic": "folate,sleep"})
        assert set(ids(result)) == {"rs1801133", "rs73598374"}

    def test_a_finding_with_no_topics_is_dropped_by_a_topic_filter(self):
        assert "rs9999" not in ids(apply_filters(SAMPLE, {"topic": "folate"}))

    def test_medicine(self):
        assert ids(apply_filters(SAMPLE, {"medicine": "warfarin"})) == ["rs6025"]

    def test_medicine_csv_matches_any(self):
        result = apply_filters(SAMPLE, {"medicine": "warfarin,metformin"})
        assert set(ids(result)) == {"rs6025", "t2d"}

    def test_condition_matches_the_conditions_list(self):
        result = apply_filters(SAMPLE, {"condition": "thrombophilia"})
        assert ids(result) == ["rs6025"]

    def test_condition_matches_the_joined_conditions_string(self):
        rows = [{"rsid": "rs1", "conditions": "Long QT syndrome"}]
        result = apply_filters(rows, {"condition": "long qt syndrome"})
        assert ids(result) == ["rs1"]

    def test_condition_csv_matches_any(self):
        result = apply_filters(SAMPLE, {"condition": "thrombophilia,type 2 diabetes"})
        assert set(ids(result)) == {"rs6025", "t2d"}

    def test_zygosity_heterozygous(self):
        result = apply_filters(SAMPLE, {"zygosity": "heterozygous"})
        assert set(ids(result)) == {"rs1801133", "rs9999", "rs6025"}

    def test_zygosity_homozygous(self):
        result = apply_filters(SAMPLE, {"zygosity": "homozygous"})
        assert set(ids(result)) == {"rs4244285", "rs73598374"}

    def test_zygosity_no_call(self):
        assert ids(apply_filters(SAMPLE, {"zygosity": "no_call"})) == ["rs4680"]

    def test_zygosity_csv(self):
        result = apply_filters(SAMPLE, {"zygosity": "homozygous,no_call"})
        assert len(result) == 3

    def test_carrier_only(self):
        result = apply_filters(SAMPLE, {"carrier_only": True})
        assert set(ids(result)) == {"rs1801133", "rs4244285", "rs6025",
                                    "rs73598374"}

    def test_carrier_only_drops_a_confirmed_non_carrier(self):
        result = apply_filters(SAMPLE, {"carrier_only": True})
        assert "rs9999" not in set(ids(result))

    def test_carrier_only_drops_an_unknown_carrier_state(self):
        result = apply_filters(SAMPLE, {"carrier_only": True})
        assert "rs4680" not in set(ids(result))

    def test_conflicts_only(self):
        assert ids(apply_filters(SAMPLE, {"conflicts_only": True})) == ["rs1801133"]

    def test_ambiguous_only_includes_palindromic_and_flipped_calls(self):
        result = apply_filters(SAMPLE, {"ambiguous_only": True})
        assert set(ids(result)) == {"rs9999", "rs4244285"}

    def test_ambiguous_only_also_accepts_freq_ambiguous(self):
        rows = [{"rsid": "rs1", "freq_ambiguous": True}, {"rsid": "rs2"}]
        assert ids(apply_filters(rows, {"ambiguous_only": True})) == ["rs1"]

    def test_the_free_text_parameter(self):
        assert ids(apply_filters(SAMPLE, {"q": "MTHFR"})) == ["rs1801133"]

    def test_a_free_text_region(self):
        result = apply_filters(SAMPLE, {"q": "chr1"})
        assert set(ids(result)) == {"rs1801133", "rs6025"}

    def test_a_free_text_operator(self):
        result = apply_filters(SAMPLE, {"q": "/STARS>=3"})
        assert set(ids(result)) == {"rs4244285", "rs6025"}

    def test_filters_and_free_text_combine_with_and(self):
        result = apply_filters(SAMPLE, {"q": "/STARS>=3", "silo": "actionable"})
        assert ids(result) == ["rs6025"]

    def test_several_filters_combine_with_and(self):
        params = {"entity_type": "snp", "repute": "bad", "min_stars": 3,
                  "carrier_only": True}
        assert ids(apply_filters(SAMPLE, params)) == ["rs4244285", "rs6025"]

    def test_a_contradictory_combination_returns_nothing(self):
        params = {"min_magnitude": 9, "max_magnitude": 1}
        assert apply_filters(SAMPLE, params) == []


CHROM_PROBE = [{"rsid": "c_%s" % c, "chromosome": c, "position": 1}
               for c in ("MT", "X", "2", "10", "1", "Y", "22")]

RSID_PROBE = [{"rsid": r} for r in ("rs100", "rs99", "rs9", "i5000", "rs1000")]


class TestSortFindings:
    def test_sort_keys_has_the_ten_supported_keys(self):
        assert set(SORT_KEYS) == {"magnitude", "frequency", "publications",
                                  "location", "modified", "gmaf", "stars",
                                  "gene", "rsid", "coverage"}

    def test_sort_keys_covers_every_key_documented_in_3_4(self):
        for key in ("magnitude", "frequency", "publications", "location",
                    "modified", "gmaf", "stars", "gene", "rsid"):
            assert key in SORT_KEYS

    def test_every_key_sorts_ascending(self):
        for name, key in SORT_KEYS.items():
            ordered = sort_findings(SAMPLE, name, "asc")
            values = [key(row) for row in ordered]
            assert values == sorted(values)

    def test_every_key_sorts_descending(self):
        for name, key in SORT_KEYS.items():
            ordered = sort_findings(SAMPLE, name, "desc")
            values = [key(row) for row in ordered]
            assert values == sorted(values, reverse=True)

    def test_every_key_returns_the_whole_set(self):
        for name in SORT_KEYS:
            assert len(sort_findings(SAMPLE, name, "asc")) == len(SAMPLE)

    def test_both_directions_differ_for_a_multi_valued_column(self):
        for name in ("magnitude", "publications", "gene", "rsid"):
            first_desc = ids(sort_findings(SAMPLE, name, "desc"))[0]
            first_asc = ids(sort_findings(SAMPLE, name, "asc"))[0]
            assert first_desc != first_asc

    def test_magnitude_descending_is_the_default(self):
        assert ids(sort_findings(SAMPLE))[0] == "rs4244285"

    def test_magnitude_ascending(self):
        assert ids(sort_findings(SAMPLE, "magnitude", "asc"))[0] == "rs4680"

    def test_frequency_descending(self):
        assert ids(sort_findings(SAMPLE, "frequency", "desc"))[0] == "rs1801133"

    def test_publications_descending(self):
        assert ids(sort_findings(SAMPLE, "publications", "desc"))[0] == "rs6025"

    def test_stars_descending(self):
        assert ids(sort_findings(SAMPLE, "stars", "desc"))[0] == "rs6025"

    def test_gmaf_descending(self):
        assert ids(sort_findings(SAMPLE, "gmaf", "desc"))[0] == "rs1801133"

    def test_gene_ascending_starts_at_the_first_symbol(self):
        assert ids(sort_findings(SAMPLE, "gene", "asc"))[0] == "rs73598374"

    def test_gene_ascending_pushes_the_unnamed_entities_last(self):
        assert sort_findings(SAMPLE, "gene", "asc")[-1]["gene"] == ""

    def test_modified_descending_uses_the_discovery_timestamp(self):
        assert ids(sort_findings(SAMPLE, "modified", "desc"))[0] == "t2d"

    def test_coverage_descending(self):
        assert ids(sort_findings(SAMPLE, "coverage", "desc"))[0] == "dgs001"

    def test_coverage_ascending_puts_the_uncovered_entities_first(self):
        assert sort_findings(SAMPLE, "coverage", "asc")[0]["coverage"] is None

    def test_an_unknown_sort_key_falls_back_to_magnitude(self):
        assert ids(sort_findings(SAMPLE, "banana", "desc")) == \
            ids(sort_findings(SAMPLE, "magnitude", "desc"))

    def test_an_unknown_sort_key_does_not_raise(self):
        assert len(sort_findings(SAMPLE, "stale bookmark")) == len(SAMPLE)

    def test_an_empty_sort_key_falls_back_to_magnitude(self):
        assert ids(sort_findings(SAMPLE, "")) == \
            ids(sort_findings(SAMPLE, "magnitude"))

    def test_a_none_sort_key_falls_back_to_magnitude(self):
        assert ids(sort_findings(SAMPLE, None)) == \
            ids(sort_findings(SAMPLE, "magnitude"))

    def test_the_sort_key_is_case_insensitive(self):
        assert ids(sort_findings(SAMPLE, "MAGNITUDE")) == \
            ids(sort_findings(SAMPLE, "magnitude"))

    def test_an_unknown_order_is_treated_as_descending(self):
        assert ids(sort_findings(SAMPLE, "magnitude", "sideways")) == \
            ids(sort_findings(SAMPLE, "magnitude", "desc"))

    def test_a_none_order_is_treated_as_descending(self):
        assert ids(sort_findings(SAMPLE, "magnitude", None))[0] == "rs4244285"

    def test_sorting_does_not_mutate_the_input_list(self):
        before = ids(SAMPLE)
        sort_findings(SAMPLE, "rsid", "asc")
        assert ids(SAMPLE) == before

    def test_sorting_an_empty_list_is_safe(self):
        assert sort_findings([], "magnitude") == []

    def test_sorting_findings_with_missing_fields_is_safe(self):
        rows = [{"rsid": "rs1"}, {"rsid": "rs2"}]
        for name in SORT_KEYS:
            assert len(sort_findings(rows, name)) == 2

    def test_sorting_findings_with_junk_field_types_is_safe(self):
        rows = [{"rsid": "rs1", "magnitude": "high", "freq": "n/a",
                 "position": "unknown", "review_stars": None, "gmaf": [],
                 "publications": {}, "coverage": "full", "gene": None,
                 "chromosome": None, "discovered_at": None},
                {"rsid": "rs2", "magnitude": 1.0}]
        for name in SORT_KEYS:
            assert len(sort_findings(rows, name)) == 2


class TestChromosomeOrder:
    def test_numeric_chromosomes_sort_numerically(self):
        ordered = [f["chromosome"] for f in sort_findings(CHROM_PROBE,
                                                          "location", "asc")]
        assert ordered == ["1", "2", "10", "22", "X", "Y", "MT"]

    def test_ten_does_not_sort_before_two(self):
        ordered = [f["chromosome"] for f in sort_findings(CHROM_PROBE,
                                                          "location", "asc")]
        assert ordered.index("2") < ordered.index("10")

    def test_the_sex_chromosomes_follow_the_autosomes(self):
        ordered = [f["chromosome"] for f in sort_findings(CHROM_PROBE,
                                                          "location", "asc")]
        assert ordered.index("22") < ordered.index("X")

    def test_mitochondrial_dna_sorts_last(self):
        ordered = [f["chromosome"] for f in sort_findings(CHROM_PROBE,
                                                          "location", "asc")]
        assert ordered[-1] == "MT"

    def test_a_chr_prefix_is_ignored(self):
        rows = [{"rsid": "a", "chromosome": "chr10", "position": 1},
                {"rsid": "b", "chromosome": "chr2", "position": 1}]
        assert ids(sort_findings(rows, "location", "asc")) == ["b", "a"]

    def test_position_breaks_a_chromosome_tie(self):
        rows = [{"rsid": "a", "chromosome": "1", "position": 500},
                {"rsid": "b", "chromosome": "1", "position": 100}]
        assert ids(sort_findings(rows, "location", "asc")) == ["b", "a"]

    def test_an_unknown_chromosome_sorts_after_the_known_ones(self):
        rows = CHROM_PROBE + [{"rsid": "weird", "chromosome": "GL000209",
                               "position": 1}]
        assert ids(sort_findings(rows, "location", "asc"))[-1] == "weird"


class TestRsidOrder:
    def test_rsids_sort_numerically(self):
        assert ids(sort_findings(RSID_PROBE, "rsid", "asc")) == [
            "i5000", "rs9", "rs99", "rs100", "rs1000"]

    def test_rs99_precedes_rs100(self):
        ordered = ids(sort_findings(RSID_PROBE, "rsid", "asc"))
        assert ordered.index("rs99") < ordered.index("rs100")

    def test_the_i_prefix_is_grouped_separately(self):
        ordered = ids(sort_findings(RSID_PROBE, "rsid", "asc"))
        assert ordered[0] == "i5000"

    def test_a_non_standard_id_sorts_after_the_numbered_ones(self):
        rows = RSID_PROBE + [{"rsid": "lactase_persistence"}]
        assert ids(sort_findings(rows, "rsid", "asc"))[-1] == "lactase_persistence"

    def test_rsid_ordering_is_case_insensitive(self):
        rows = [{"rsid": "RS100"}, {"rsid": "rs9"}]
        assert ids(sort_findings(rows, "rsid", "asc")) == ["rs9", "RS100"]


class TestPaginate:
    def test_the_default_limit_returns_a_small_set_whole(self):
        assert len(paginate(SAMPLE)) == len(SAMPLE)

    def test_a_limit_of_zero_means_everything(self):
        assert len(paginate(SAMPLE, limit=0)) == len(SAMPLE)

    def test_a_limit_of_none_means_everything(self):
        assert len(paginate(SAMPLE, limit=None)) == len(SAMPLE)

    def test_a_limit_slices_the_head(self):
        assert ids(paginate(SAMPLE, limit=3)) == ids(SAMPLE)[:3]

    def test_a_limit_larger_than_the_set_returns_everything(self):
        assert len(paginate(SAMPLE, limit=500)) == len(SAMPLE)

    def test_an_offset_skips_the_head(self):
        assert ids(paginate(SAMPLE, limit=0, offset=5)) == ids(SAMPLE)[5:]

    def test_a_limit_and_an_offset_together(self):
        assert ids(paginate(SAMPLE, limit=2, offset=2)) == ids(SAMPLE)[2:4]

    def test_an_offset_past_the_end_returns_nothing(self):
        assert paginate(SAMPLE, limit=5, offset=999) == []

    def test_a_negative_offset_is_clamped_to_zero(self):
        assert ids(paginate(SAMPLE, limit=2, offset=-4)) == ids(SAMPLE)[:2]

    def test_a_none_offset_is_treated_as_zero(self):
        assert ids(paginate(SAMPLE, limit=2, offset=None)) == ids(SAMPLE)[:2]

    def test_an_offset_of_zero_with_no_limit_returns_everything(self):
        assert len(paginate(SAMPLE, limit=0, offset=0)) == len(SAMPLE)

    def test_paginating_an_empty_list_is_safe(self):
        assert paginate([], limit=10) == []


FACET_PROBE = [
    {"rsid": "rs1", "gene": "ZED", "entity_type": "snp", "silo": "actionable",
     "zygosity": "homozygous", "category": "CARDIO", "repute": "Bad",
     "clinical_sig": "pathogenic", "cpic_level": "A", "freq_band": "rare",
     "confidence": "high", "review_stars": 3, "topics": ["heart"],
     "medicines": ["warfarin"], "conditions_list": ["Thrombophilia"]},
    {"rsid": "rs2", "gene": "ABLE", "entity_type": "snp", "silo": "actionable",
     "zygosity": "heterozygous", "category": "CARDIO", "repute": "",
     "clinical_sig": "benign", "cpic_level": "", "freq_band": "common",
     "confidence": "low", "review_stars": 1, "topics": ["heart", "lipids"],
     "medicines": [], "conditions_list": []},
    {"rsid": "rs3", "gene": "ZED", "entity_type": "trait",
     "silo": "informational", "zygosity": "", "category": "Diet",
     "repute": "", "clinical_sig": "", "cpic_level": "", "freq_band": "",
     "confidence": "low", "review_stars": None, "topics": [],
     "medicines": [], "conditions_list": []},
    {"rsid": "rs4", "gene": "MID", "entity_type": "snp",
     "silo": "informational", "zygosity": "homozygous", "category": "",
     "repute": "Good", "clinical_sig": "benign", "cpic_level": "B",
     "freq_band": "rare", "confidence": "moderate", "review_stars": 1,
     "topics": ["lipids"], "medicines": ["statins"], "conditions_list": []},
]

DOCUMENTED_FACET_BUCKETS = ("genes", "topics", "medicines", "conditions",
                            "silos", "entity_types", "zygosity", "categories")


def facet_values(facets, bucket):
    return [row["value"] for row in facets[bucket]]


def facet_count(facets, bucket, value):
    for row in facets[bucket]:
        if row["value"] == value:
            return row["count"]
    return 0


class TestBuildFacets:
    def test_every_documented_bucket_is_rendered(self):
        facets = build_facets(SAMPLE)
        for bucket in DOCUMENTED_FACET_BUCKETS:
            assert bucket in facets

    def test_the_extra_scoring_buckets_are_rendered(self):
        facets = build_facets(SAMPLE)
        for bucket in ("reputes", "clinvar_significance", "review_stars",
                       "cpic_levels", "freq_bands", "confidence"):
            assert bucket in facets

    def test_every_bucket_is_a_list_of_value_and_count_rows(self):
        facets = build_facets(SAMPLE)
        for rows in facets.values():
            for row in rows:
                assert set(row) == {"value", "count"}

    def test_counts_are_correct(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "genes", "ZED") == 2
        assert facet_count(facets, "genes", "ABLE") == 1

    def test_buckets_sort_by_count_descending(self):
        facets = build_facets(FACET_PROBE)
        counts = [row["count"] for row in facets["genes"]]
        assert counts == sorted(counts, reverse=True)

    def test_ties_sort_by_value_ascending(self):
        facets = build_facets(FACET_PROBE)
        assert facet_values(facets, "genes") == ["ZED", "ABLE", "MID"]

    def test_blank_values_are_not_counted(self):
        facets = build_facets(FACET_PROBE)
        assert "" not in facet_values(facets, "categories")

    def test_a_missing_key_is_not_counted(self):
        facets = build_facets([{"rsid": "rs1"}])
        assert facets["genes"] == []

    def test_topics_are_exploded_from_the_list(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "topics", "heart") == 2
        assert facet_count(facets, "topics", "lipids") == 2

    def test_medicines_are_exploded_from_the_list(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "medicines", "warfarin") == 1

    def test_conditions_are_exploded_from_the_list(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "conditions", "Thrombophilia") == 1

    def test_silos_are_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "silos", "actionable") == 2

    def test_entity_types_are_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "entity_types", "snp") == 3
        assert facet_count(facets, "entity_types", "trait") == 1

    def test_zygosity_is_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "zygosity", "homozygous") == 2

    def test_categories_are_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "categories", "CARDIO") == 2

    def test_a_blank_repute_is_labelled_not_set(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "reputes", "Not Set") == 2

    def test_review_stars_are_counted_as_strings(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "review_stars", "1") == 2
        assert facet_count(facets, "review_stars", "3") == 1

    def test_a_null_star_count_is_not_counted(self):
        facets = build_facets(FACET_PROBE)
        total = sum(row["count"] for row in facets["review_stars"])
        assert total == 3

    def test_cpic_levels_are_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "cpic_levels", "A") == 1

    def test_freq_bands_are_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "freq_bands", "rare") == 2

    def test_confidence_is_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "confidence", "low") == 2

    def test_clinvar_significance_is_counted(self):
        facets = build_facets(FACET_PROBE)
        assert facet_count(facets, "clinvar_significance", "benign") == 2

    def test_an_empty_input_gives_every_bucket_as_an_empty_list(self):
        facets = build_facets([])
        assert all(rows == [] for rows in facets.values())

    def test_a_none_input_is_safe(self):
        assert build_facets(None)["genes"] == []

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   docs/API_V2.md section 3.6 lists clinvar_diseases in the /facets response
    #   and routes_v2.py returns build_facets verbatim, so the documented bucket
    #   is never rendered.
    def test_the_documented_clinvar_diseases_bucket_is_rendered(self):
        assert "clinvar_diseases" in build_facets(SAMPLE)


class TestSummarise:
    def test_the_documented_shape(self):
        assert set(summarise(SAMPLE)) == {"silos", "entity_types", "reputes"}

    def test_silo_counts(self):
        assert summarise(SAMPLE)["silos"]["actionable"] == 2

    def test_entity_type_counts(self):
        entities = summarise(SAMPLE)["entity_types"]
        assert entities["snp"] == 6
        assert entities["genoset"] == 1
        assert entities["trait"] == 1
        assert entities["prs"] == 1

    def test_repute_counts(self):
        reputes = summarise(SAMPLE)["reputes"]
        assert reputes == {"Good": 1, "Bad": 3, "unset": 5}

    def test_the_three_repute_keys_are_always_present(self):
        assert set(summarise([])["reputes"]) == {"Good", "Bad", "unset"}

    def test_a_missing_silo_counts_as_informational(self):
        assert summarise([{"rsid": "rs1"}])["silos"] == {"informational": 1}

    def test_a_missing_entity_type_counts_as_snp(self):
        assert summarise([{"rsid": "rs1"}])["entity_types"] == {"snp": 1}

    def test_an_unexpected_repute_counts_as_unset(self):
        reputes = summarise([{"rsid": "rs1", "repute": "Sideways"}])["reputes"]
        assert reputes["unset"] == 1

    def test_an_empty_input_gives_zero_counts(self):
        summary = summarise([])
        assert summary["silos"] == {}
        assert summary["entity_types"] == {}

    def test_a_none_input_is_safe(self):
        assert summarise(None)["reputes"]["unset"] == 0


class TestFilterAndSort:
    def test_the_documented_response_shape(self):
        out = filter_and_sort(SAMPLE, {})
        assert set(out) == {"findings", "total", "returned", "offset",
                            "filtered_summary", "facets"}

    def test_total_counts_the_filtered_set_before_pagination(self):
        out = filter_and_sort(SAMPLE, {"limit": 2})
        assert out["total"] == len(SAMPLE)

    def test_returned_counts_the_page(self):
        out = filter_and_sort(SAMPLE, {"limit": 2})
        assert out["returned"] == 2
        assert len(out["findings"]) == 2

    def test_the_offset_is_echoed(self):
        assert filter_and_sort(SAMPLE, {"offset": 3})["offset"] == 3

    def test_a_missing_offset_is_zero(self):
        assert filter_and_sort(SAMPLE, {})["offset"] == 0

    def test_the_default_limit_is_two_hundred(self):
        out = filter_and_sort(SAMPLE, {})
        assert out["returned"] == len(SAMPLE)

    def test_a_limit_of_zero_returns_everything(self):
        out = filter_and_sort(SAMPLE, {"limit": 0})
        assert out["returned"] == out["total"]

    def test_the_offset_moves_the_window(self):
        first = filter_and_sort(SAMPLE, {"limit": 1, "offset": 0})
        second = filter_and_sort(SAMPLE, {"limit": 1, "offset": 1})
        assert first["findings"][0]["rsid"] != second["findings"][0]["rsid"]

    def test_filters_are_applied_before_sorting(self):
        out = filter_and_sort(SAMPLE, {"entity_type": "snp",
                                       "sort": "magnitude", "order": "desc"})
        assert ids(out["findings"])[0] == "rs4244285"
        assert out["total"] == 6

    def test_the_sort_key_is_honoured(self):
        out = filter_and_sort(SAMPLE, {"sort": "publications", "order": "desc"})
        assert ids(out["findings"])[0] == "rs6025"

    def test_the_order_is_honoured(self):
        out = filter_and_sort(SAMPLE, {"sort": "magnitude", "order": "asc"})
        assert ids(out["findings"])[0] == "rs4680"

    def test_the_filtered_summary_describes_the_filtered_set(self):
        out = filter_and_sort(SAMPLE, {"entity_type": "snp"})
        assert out["filtered_summary"]["entity_types"] == {"snp": 6}

    def test_the_facets_describe_the_filtered_set(self):
        out = filter_and_sort(SAMPLE, {"gene": "MTHFR"})
        assert facet_values(out["facets"], "genes") == ["MTHFR"]

    def test_the_facets_are_built_before_pagination(self):
        out = filter_and_sort(SAMPLE, {"limit": 1})
        assert facet_count(out["facets"], "entity_types", "snp") == 6

    def test_an_empty_finding_set_is_safe(self):
        out = filter_and_sort([], {})
        assert out["total"] == 0
        assert out["returned"] == 0
        assert out["findings"] == []

    def test_a_free_text_query_flows_through(self):
        out = filter_and_sort(SAMPLE, {"q": "/nocall"})
        assert ids(out["findings"]) == ["rs4680"]
