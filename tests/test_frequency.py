"""Tests for backend.frequency population frequencies and strand handling."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

from backend import frequency
from backend.frequency import (
    AGGREGATE_MODES,
    ANNOTATION_KEYS,
    DEFAULT_POPULATION,
    POPULATION_CODES,
    RARITY_RAMP,
    UNKNOWN_COLOR,
    aggregate_frequency,
    allele_frequency,
    annotate,
    annotate_all,
    available_populations,
    genotype_frequency,
    genotype_frequency_detail,
    gmaf,
    is_palindromic,
    load_frequencies,
    minor_allele,
    observed_alleles,
    population_series,
    rarity_band,
    rarity_color,
    resolve_strand,
    sample_size,
)

FREQ_PATH = Path(__file__).parent.parent / "data" / "frequencies.json"

# MTHFR C677T: quoted C/T on a consumer array, stored A/G by dbSNP. The single
# most useful test case in the file, because a strand bug here reads as
# "no data" rather than as an error.
MTHFR = "rs1801133"

# 9p21, a C/G site, and FTO, an A/T site. Both palindromic: complementing
# either genotype yields the other observed allele, so no metadata can settle
# which strand the array reported.
PALINDROMIC_CG = "rs1333049"
PALINDROMIC_AT = "rs9939609"

# LPA. The PJL panel stores GG at exactly 0.0, which is a measurement and not
# an absence of data.
ZERO_RSID = "rs10455872"
ZERO_POPULATION = "PJL"

UNKNOWN_RSID = "rs00000000"


def require_data() -> dict:
    """Return the loaded frequency table, skipping when it is not built."""
    if not FREQ_PATH.exists():
        pytest.skip("data/frequencies.json has not been built")
    data = load_frequencies()
    if not data:
        pytest.skip("data/frequencies.json carries no entries")
    return data


def require_rsid(rsid: str) -> dict:
    """Return one entry, skipping when that rsID is not bundled."""
    data = require_data()
    if rsid not in data:
        pytest.skip(f"{rsid} is not in the bundled frequency data")
    return data[rsid]


class TestRarityBand:
    def test_none_is_unknown(self):
        assert rarity_band(None) == "unknown"

    def test_zero_is_very_rare_not_unknown(self):
        assert rarity_band(0.0) == "very_rare"

    def test_just_under_one_is_very_rare(self):
        assert rarity_band(0.99) == "very_rare"

    def test_exactly_one_is_rare(self):
        assert rarity_band(1.0) == "rare"

    def test_just_under_five_is_rare(self):
        assert rarity_band(4.99) == "rare"

    def test_exactly_five_is_uncommon(self):
        assert rarity_band(5.0) == "uncommon"

    def test_just_under_twenty_is_uncommon(self):
        assert rarity_band(19.99) == "uncommon"

    def test_exactly_twenty_is_common(self):
        assert rarity_band(20.0) == "common"

    def test_exactly_fifty_is_common(self):
        assert rarity_band(50.0) == "common"

    def test_just_over_fifty_is_majority(self):
        assert rarity_band(50.01) == "majority"

    def test_one_hundred_is_majority(self):
        assert rarity_band(100.0) == "majority"

    def test_unparsable_value_is_unknown(self):
        assert rarity_band("not a number") == "unknown"

    def test_band_is_monotonic_across_the_boundaries(self):
        order = ["very_rare", "rare", "uncommon", "common", "majority"]
        seen = [rarity_band(f) for f in (0.0, 1.0, 5.0, 20.0, 100.0)]
        assert seen == order


class TestRarityColor:
    def test_none_gets_the_unknown_color(self):
        assert rarity_color(None) == UNKNOWN_COLOR

    def test_unknown_color_is_a_light_neutral(self):
        # Absent data must never be rendered as an alarming finding.
        assert UNKNOWN_COLOR == "#EFEFEF"

    def test_every_color_is_seven_characters_starting_with_hash(self):
        for value in (None, 0.0, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0):
            color = rarity_color(value)
            assert isinstance(color, str)
            assert len(color) == 7
            assert color.startswith("#")

    def test_majority_genotype_is_white(self):
        assert rarity_color(50.0) == "#FFFFFF"
        assert rarity_color(99.0) == "#FFFFFF"

    def test_twenty_percent_is_the_blush_stop(self):
        assert rarity_color(20.0) == "#FFF3EE"

    def test_ten_percent_is_the_pale_salmon_stop(self):
        assert rarity_color(10.0) == "#FFE0D4"

    def test_five_percent_is_the_salmon_stop(self):
        assert rarity_color(5.0) == "#FFC5B2"

    def test_one_percent_is_the_coral_stop(self):
        assert rarity_color(1.0) == "#FA9A80"

    def test_a_tenth_of_a_percent_is_the_orange_red_stop(self):
        assert rarity_color(0.1) == "#EF6A4C"

    def test_zero_is_the_deep_red_stop_not_the_unknown_grey(self):
        assert rarity_color(0.0) == "#D8351B"
        assert rarity_color(0.0) != UNKNOWN_COLOR

    def test_a_negative_value_falls_through_to_the_last_stop(self):
        assert rarity_color(-1.0) == RARITY_RAMP[-1][1]

    def test_unparsable_value_gets_the_unknown_color(self):
        assert rarity_color("not a number") == UNKNOWN_COLOR

    def test_the_ramp_runs_from_common_to_rare(self):
        bounds = [lower for lower, _ in RARITY_RAMP]
        assert bounds == sorted(bounds, reverse=True)


class TestZeroVersusNone:
    def test_a_stored_zero_comes_back_as_zero(self):
        entry = require_rsid(ZERO_RSID)
        table = (entry.get("genotypes") or {}).get(ZERO_POPULATION) or {}
        if table.get("GG") != 0.0:
            pytest.skip("the bundled data no longer stores a zero here")
        assert genotype_frequency(ZERO_RSID, "G", "G", ZERO_POPULATION) == 0.0

    def test_a_stored_zero_is_not_none(self):
        entry = require_rsid(ZERO_RSID)
        table = (entry.get("genotypes") or {}).get(ZERO_POPULATION) or {}
        if table.get("GG") != 0.0:
            pytest.skip("the bundled data no longer stores a zero here")
        assert genotype_frequency(ZERO_RSID, "G", "G", ZERO_POPULATION) is not None

    def test_a_no_call_genotype_gives_none_for_a_known_rsid(self):
        require_rsid(MTHFR)
        assert genotype_frequency(MTHFR, "-", "-") is None

    def test_zero_and_none_land_in_different_bands(self):
        assert rarity_band(0.0) != rarity_band(None)

    def test_zero_and_none_get_different_colors(self):
        assert rarity_color(0.0) != rarity_color(None)

    def test_annotate_preserves_a_zero_rather_than_blanking_it(self):
        entry = require_rsid(ZERO_RSID)
        table = (entry.get("genotypes") or {}).get(ZERO_POPULATION) or {}
        if table.get("GG") != 0.0:
            pytest.skip("the bundled data no longer stores a zero here")
        finding = {"rsid": ZERO_RSID, "allele1": "G", "allele2": "G"}
        annotate(finding, ZERO_POPULATION)
        assert finding["freq"] == 0.0
        assert finding["freq_band"] == "very_rare"
        assert finding["freq_color"] != UNKNOWN_COLOR


class TestObservedAlleles:
    def test_mthfr_is_stored_in_dbsnp_orientation(self):
        require_rsid(MTHFR)
        assert observed_alleles(MTHFR) == {"A", "G"}

    def test_the_array_alleles_are_not_the_stored_alleles(self):
        # This is the whole reason the strand layer exists.
        require_rsid(MTHFR)
        assert "C" not in observed_alleles(MTHFR)
        assert "T" not in observed_alleles(MTHFR)

    def test_a_palindromic_site_stores_c_and_g(self):
        require_rsid(PALINDROMIC_CG)
        assert observed_alleles(PALINDROMIC_CG) == {"C", "G"}

    def test_a_palindromic_site_stores_a_and_t(self):
        require_rsid(PALINDROMIC_AT)
        assert observed_alleles(PALINDROMIC_AT) == {"A", "T"}

    def test_only_real_bases_are_returned(self):
        require_data()
        for rsid in (MTHFR, PALINDROMIC_CG, PALINDROMIC_AT):
            assert observed_alleles(rsid) <= {"A", "C", "G", "T"}

    def test_a_single_population_is_a_subset_of_the_union(self):
        require_rsid(MTHFR)
        assert observed_alleles(MTHFR, "CEU") <= observed_alleles(MTHFR)

    def test_an_unknown_rsid_gives_an_empty_set(self):
        assert observed_alleles(UNKNOWN_RSID) == set()


class TestIsPalindromic:
    def test_a_over_t_is_palindromic(self):
        assert is_palindromic("A", "T") is True

    def test_t_over_a_is_palindromic(self):
        assert is_palindromic("T", "A") is True

    def test_c_over_g_is_palindromic(self):
        assert is_palindromic("C", "G") is True

    def test_g_over_c_is_palindromic(self):
        assert is_palindromic("G", "C") is True

    def test_a_homozygote_is_not_palindromic(self):
        assert is_palindromic("A", "A") is False

    def test_a_c_g_and_t_homozygotes_are_not_palindromic(self):
        for base in ("A", "C", "G", "T"):
            assert is_palindromic(base, base) is False

    def test_a_over_g_is_not_palindromic(self):
        assert is_palindromic("A", "G") is False

    def test_c_over_t_is_not_palindromic(self):
        assert is_palindromic("C", "T") is False

    def test_a_no_call_is_not_palindromic(self):
        assert is_palindromic("A", "-") is False
        assert is_palindromic("", "T") is False

    def test_case_does_not_matter(self):
        assert is_palindromic("a", "t") is True
        assert is_palindromic("c", "g") is True


class TestResolveStrand:
    def test_mthfr_array_alleles_flip_to_the_stored_orientation(self):
        require_rsid(MTHFR)
        out = resolve_strand(MTHFR, "C", "T")
        assert out["flipped"] is True
        assert out["resolved"] is True
        assert out["allele1"] == "G"
        assert out["allele2"] == "A"

    def test_a_flipped_site_is_not_reported_as_ambiguous(self):
        require_rsid(MTHFR)
        assert resolve_strand(MTHFR, "C", "T")["ambiguous"] is False

    def test_the_observed_alleles_are_reported_back_sorted(self):
        require_rsid(MTHFR)
        assert resolve_strand(MTHFR, "C", "T")["observed"] == ["A", "G"]

    def test_already_oriented_alleles_are_left_alone(self):
        require_rsid(MTHFR)
        out = resolve_strand(MTHFR, "G", "A")
        assert out["flipped"] is False
        assert out["resolved"] is True
        assert (out["allele1"], out["allele2"]) == ("G", "A")

    def test_a_flipped_homozygote_resolves_too(self):
        require_rsid(MTHFR)
        out = resolve_strand(MTHFR, "C", "C")
        assert out["flipped"] is True
        assert (out["allele1"], out["allele2"]) == ("G", "G")

    def test_allele_order_is_preserved_through_a_flip(self):
        require_rsid(MTHFR)
        out = resolve_strand(MTHFR, "T", "C")
        assert (out["allele1"], out["allele2"]) == ("A", "G")

    def test_a_palindromic_cg_site_is_ambiguous_and_not_flipped(self):
        # The caller must be able to warn instead of quietly guessing.
        require_rsid(PALINDROMIC_CG)
        out = resolve_strand(PALINDROMIC_CG, "C", "G")
        assert out["ambiguous"] is True
        assert out["flipped"] is False
        assert out["resolved"] is True

    def test_a_palindromic_at_site_is_ambiguous_and_not_flipped(self):
        require_rsid(PALINDROMIC_AT)
        out = resolve_strand(PALINDROMIC_AT, "A", "T")
        assert out["ambiguous"] is True
        assert out["flipped"] is False
        assert out["resolved"] is True

    def test_a_palindromic_site_keeps_the_unflipped_reading(self):
        require_rsid(PALINDROMIC_CG)
        out = resolve_strand(PALINDROMIC_CG, "C", "G")
        assert (out["allele1"], out["allele2"]) == ("C", "G")

    def test_an_unknown_rsid_falls_through_unresolved(self):
        out = resolve_strand(UNKNOWN_RSID, "C", "T")
        assert out["resolved"] is False
        assert out["flipped"] is False
        assert (out["allele1"], out["allele2"]) == ("C", "T")
        assert out["observed"] == []

    def test_a_no_call_returns_a_well_formed_dict(self):
        out = resolve_strand(MTHFR, "-", "T")
        assert out["resolved"] is False
        assert out["allele1"] == ""

    def test_every_documented_key_is_always_present(self):
        for args in ((MTHFR, "C", "T"), (UNKNOWN_RSID, "A", "G"),
                     (MTHFR, "-", "-"), (PALINDROMIC_AT, "A", "T")):
            out = resolve_strand(*args)
            for key in ("allele1", "allele2", "flipped", "ambiguous",
                        "resolved", "observed"):
                assert key in out

    def test_lowercase_input_is_normalised_to_uppercase(self):
        require_rsid(MTHFR)
        out = resolve_strand(MTHFR, "c", "t")
        assert (out["allele1"], out["allele2"]) == ("G", "A")


class TestAlleleFrequency:
    def test_the_array_allele_resolves_by_default(self):
        require_rsid(MTHFR)
        value = allele_frequency(MTHFR, "C", "CEU")
        assert value is not None
        assert isinstance(value, float)

    def test_the_array_allele_is_absent_without_strand_tolerance(self):
        # The literal table has no C at all. Turning tolerance off must say so
        # rather than inventing a number.
        require_rsid(MTHFR)
        assert allele_frequency(MTHFR, "C", "CEU", strand_tolerant=False) is None

    def test_the_stored_allele_resolves_without_strand_tolerance(self):
        require_rsid(MTHFR)
        assert allele_frequency(MTHFR, "G", "CEU",
                                strand_tolerant=False) is not None

    def test_the_tolerant_lookup_returns_the_complements_frequency(self):
        require_rsid(MTHFR)
        tolerant = allele_frequency(MTHFR, "C", "CEU")
        literal = allele_frequency(MTHFR, "G", "CEU", strand_tolerant=False)
        assert tolerant == literal

    def test_the_other_array_allele_maps_to_the_other_stored_allele(self):
        require_rsid(MTHFR)
        assert allele_frequency(MTHFR, "T", "CEU") == \
            allele_frequency(MTHFR, "A", "CEU", strand_tolerant=False)

    def test_allele_frequencies_are_fractions_not_percentages(self):
        require_rsid(MTHFR)
        for allele in ("C", "T"):
            value = allele_frequency(MTHFR, allele, "CEU")
            assert 0.0 <= value <= 1.0

    def test_the_two_alleles_sum_to_one(self):
        require_rsid(MTHFR)
        total = (allele_frequency(MTHFR, "C", "CEU")
                 + allele_frequency(MTHFR, "T", "CEU"))
        assert abs(total - 1.0) < 0.01

    def test_an_unknown_rsid_gives_none(self):
        assert allele_frequency(UNKNOWN_RSID, "A", "CEU") is None

    def test_an_unknown_population_gives_none(self):
        require_rsid(MTHFR)
        assert allele_frequency(MTHFR, "C", "NOT_A_POPULATION") is None

    def test_a_no_call_allele_gives_none(self):
        require_rsid(MTHFR)
        for token in ("", "-", "N", "0"):
            assert allele_frequency(MTHFR, token, "CEU") is None


class TestGenotypeFrequencyDetail:
    def test_the_flipped_genotype_resolves_to_a_real_frequency(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "CEU")
        assert detail["frequency"] is not None
        assert detail["flipped"] is True
        assert detail["method"] == "hardy_weinberg"

    def test_a_derived_table_is_never_promoted_to_observed(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "CEU")
        assert detail["derived"] is True
        assert detail["method"] != "observed"

    def test_the_queried_genotype_is_reported_in_stored_orientation(self):
        require_rsid(MTHFR)
        assert genotype_frequency_detail(MTHFR, "C", "T", "CEU")["queried"] == "GA"

    def test_the_population_is_echoed_back(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "TSI")
        assert detail["population"] == "TSI"

    def test_the_sample_size_is_carried_alongside(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "CEU")
        assert detail["n"] == sample_size(MTHFR, "CEU")

    def test_a_heterozygote_is_order_insensitive(self):
        require_rsid(MTHFR)
        forward = genotype_frequency(MTHFR, "C", "T", "CEU")
        reverse = genotype_frequency(MTHFR, "T", "C", "CEU")
        assert forward == reverse

    def test_a_palindromic_site_is_flagged_ambiguous_in_the_detail(self):
        require_rsid(PALINDROMIC_CG)
        detail = genotype_frequency_detail(PALINDROMIC_CG, "C", "G", "CEU")
        assert detail["ambiguous"] is True
        assert detail["frequency"] is not None

    def test_an_unknown_rsid_is_reported_unavailable(self):
        detail = genotype_frequency_detail(UNKNOWN_RSID, "A", "G", "CEU")
        assert detail["frequency"] is None
        assert detail["method"] == "unavailable"
        assert detail["derived"] is False

    def test_a_no_call_is_reported_unavailable(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "-", "-", "CEU")
        assert detail["frequency"] is None
        assert detail["method"] == "unavailable"

    def test_the_convenience_wrapper_returns_the_detail_frequency(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "CEU")
        assert genotype_frequency(MTHFR, "C", "T", "CEU") == detail["frequency"]

    def test_every_documented_detail_key_is_present(self):
        require_rsid(MTHFR)
        detail = genotype_frequency_detail(MTHFR, "C", "T", "CEU")
        for key in ("frequency", "population", "derived", "method", "n",
                    "flipped", "ambiguous", "queried"):
            assert key in detail


class TestHardyWeinberg:
    def test_the_heterozygote_is_two_p_q(self):
        require_rsid(MTHFR)
        p = allele_frequency(MTHFR, "C", "CEU")
        q = allele_frequency(MTHFR, "T", "CEU")
        het = genotype_frequency(MTHFR, "C", "T", "CEU")
        assert abs(het - 2.0 * p * q * 100.0) < 0.05

    def test_the_major_homozygote_is_p_squared(self):
        require_rsid(MTHFR)
        p = allele_frequency(MTHFR, "C", "CEU")
        hom = genotype_frequency(MTHFR, "C", "C", "CEU")
        assert abs(hom - p * p * 100.0) < 0.05

    def test_the_minor_homozygote_is_q_squared(self):
        require_rsid(MTHFR)
        q = allele_frequency(MTHFR, "T", "CEU")
        hom = genotype_frequency(MTHFR, "T", "T", "CEU")
        assert abs(hom - q * q * 100.0) < 0.05

    def test_the_three_genotypes_sum_to_about_one_hundred(self):
        require_rsid(MTHFR)
        total = sum(genotype_frequency(MTHFR, a, b, "CEU")
                    for a, b in (("C", "C"), ("C", "T"), ("T", "T")))
        assert abs(total - 100.0) < 0.5

    def test_the_three_genotypes_sum_to_one_hundred_in_every_population(self):
        require_rsid(MTHFR)
        for population in POPULATION_CODES:
            values = [genotype_frequency(MTHFR, a, b, population)
                      for a, b in (("C", "C"), ("C", "T"), ("T", "T"))]
            if any(v is None for v in values):
                continue
            assert abs(sum(values) - 100.0) < 0.5, population

    def test_derivation_fires_when_no_genotype_table_exists(self, tmp_path,
                                                            monkeypatch):
        # Exercises the pure Hardy-Weinberg fallback: allele data only, so
        # 2pq and p squared have to be computed rather than looked up.
        # FREQUENCY_FILE has to be redirected as well as load_frequencies
        # called, because every internal accessor calls load_frequencies with
        # no argument and would otherwise pull the real file straight back in.
        real_path = frequency.FREQUENCY_FILE
        payload = (
            '{"_meta": {"version": "test"}, "frequencies": '
            '{"rs9999999": {"alleles": {"CEU": {"A": 0.2, "G": 0.8}}, '
            '"n": {"CEU": 100}}}}'
        )
        path = tmp_path / "frequencies.json"
        path.write_text(payload, encoding="utf-8")
        try:
            monkeypatch.setattr(frequency, "FREQUENCY_FILE", path)
            load_frequencies(path)
            het = genotype_frequency_detail("rs9999999", "A", "G", "CEU")
            assert het["method"] == "hardy_weinberg"
            assert het["derived"] is True
            assert abs(het["frequency"] - 32.0) < 0.01
            hom = genotype_frequency_detail("rs9999999", "G", "G", "CEU")
            assert abs(hom["frequency"] - 64.0) < 0.01
        finally:
            monkeypatch.setattr(frequency, "FREQUENCY_FILE", real_path)
            load_frequencies(real_path)

    def test_an_explicit_path_sticks_across_accessors(self, tmp_path):
        # An explicit path is remembered for the rest of the process. Every
        # accessor calls load_frequencies() with no argument, so if the override
        # were forgotten a caller who pointed the module at a fixture would
        # silently be reading the bundled file again on the very next lookup.
        require_data()
        path = tmp_path / "frequencies.json"
        path.write_text('{"frequencies": {"rs9999999": {}}}', encoding="utf-8")
        try:
            assert "rs9999999" in load_frequencies(path)
            # The fixture has no MTHFR, so the override is genuinely in force.
            assert observed_alleles(MTHFR) == set()
            assert gmaf(MTHFR) is None
        finally:
            frequency.reset_source()

    def test_reset_source_restores_the_bundled_file(self, tmp_path):
        require_rsid(MTHFR)
        path = tmp_path / "frequencies.json"
        path.write_text('{"frequencies": {"rs9999999": {}}}', encoding="utf-8")
        load_frequencies(path)
        assert observed_alleles(MTHFR) == set()
        frequency.reset_source()
        assert observed_alleles(MTHFR) == {"A", "G"}

    def test_reset_source_is_safe_to_call_twice(self):
        require_data()
        frequency.reset_source()
        frequency.reset_source()
        assert isinstance(load_frequencies(), dict)


class TestPopulationSeries:
    def test_exactly_one_entry_is_marked_yours(self):
        require_rsid(MTHFR)
        series = population_series(MTHFR, "C", "T", "CEU")
        assert sum(1 for entry in series if entry["yours"]) == 1

    def test_the_marked_entry_is_the_requested_population(self):
        require_rsid(MTHFR)
        series = population_series(MTHFR, "C", "T", "TSI")
        yours = [entry for entry in series if entry["yours"]]
        assert yours[0]["code"] == "TSI"

    def test_one_entry_per_available_population(self):
        require_rsid(MTHFR)
        series = population_series(MTHFR, "C", "T", "CEU")
        assert len(series) == len(available_populations())

    def test_every_entry_carries_the_documented_keys(self):
        require_rsid(MTHFR)
        for entry in population_series(MTHFR, "C", "T", "CEU"):
            for key in ("code", "label", "brief", "frequency", "yours"):
                assert key in entry

    def test_an_unavailable_population_still_anchors_exactly_one_entry(self):
        require_rsid(MTHFR)
        series = population_series(MTHFR, "C", "T", "NOT_A_POPULATION")
        assert sum(1 for entry in series if entry["yours"]) == 1

    def test_an_unknown_rsid_still_returns_a_well_formed_series(self):
        require_data()
        series = population_series(UNKNOWN_RSID, "A", "G", "CEU")
        assert len(series) == len(available_populations())
        assert all(entry["frequency"] is None for entry in series)
        assert sum(1 for entry in series if entry["yours"]) == 1


class TestAggregateFrequency:
    def test_max_is_at_least_min(self):
        require_rsid(MTHFR)
        assert aggregate_frequency(MTHFR, "C", "T", "MAX") >= \
            aggregate_frequency(MTHFR, "C", "T", "MIN")

    def test_avg_sits_between_min_and_max(self):
        require_rsid(MTHFR)
        low = aggregate_frequency(MTHFR, "C", "T", "MIN")
        mid = aggregate_frequency(MTHFR, "C", "T", "AVG")
        high = aggregate_frequency(MTHFR, "C", "T", "MAX")
        assert low <= mid <= high

    def test_max_equals_the_largest_value_in_the_series(self):
        require_rsid(MTHFR)
        values = [entry["frequency"]
                  for entry in population_series(MTHFR, "C", "T", "CEU")
                  if entry["frequency"] is not None]
        assert aggregate_frequency(MTHFR, "C", "T", "MAX") == round(max(values), 2)

    def test_min_equals_the_smallest_value_in_the_series(self):
        require_rsid(MTHFR)
        values = [entry["frequency"]
                  for entry in population_series(MTHFR, "C", "T", "CEU")
                  if entry["frequency"] is not None]
        assert aggregate_frequency(MTHFR, "C", "T", "MIN") == round(min(values), 2)

    def test_the_mode_is_case_insensitive(self):
        require_rsid(MTHFR)
        assert aggregate_frequency(MTHFR, "C", "T", "max") == \
            aggregate_frequency(MTHFR, "C", "T", "MAX")

    def test_an_unknown_mode_returns_none(self):
        require_rsid(MTHFR)
        assert aggregate_frequency(MTHFR, "C", "T", "MEDIAN") is None

    def test_the_four_documented_modes_are_declared(self):
        assert set(AGGREGATE_MODES) == {"MAX", "AVG", "MIN", "GLOBAL"}

    def test_global_derives_from_gmaf_in_stored_orientation(self):
        require_rsid(MTHFR)
        value = aggregate_frequency(MTHFR, "G", "A", "GLOBAL")
        assert value is not None
        assert 0.0 <= value <= 100.0

    def test_global_is_strand_tolerant_like_the_population_modes(self):
        # Regression guard. GLOBAL used to compare the queried alleles against
        # the stored minor allele directly, so an array-orientation query
        # returned None for every minus-strand variant while MAX returned a
        # number. That inconsistency reads as missing data rather than a bug.
        require_rsid(MTHFR)
        assert aggregate_frequency(MTHFR, "C", "T", "GLOBAL") is not None
        assert aggregate_frequency(MTHFR, "C", "T", "MAX") is not None

    def test_global_agrees_across_both_orientations(self):
        # C;T (array plus strand) and G;A (dbSNP orientation) are the same
        # genotype, so GLOBAL must return the same number for both.
        require_rsid(MTHFR)
        array_strand = aggregate_frequency(MTHFR, "C", "T", "GLOBAL")
        dbsnp_strand = aggregate_frequency(MTHFR, "G", "A", "GLOBAL")
        assert array_strand == dbsnp_strand

    def test_an_unknown_rsid_returns_none_for_every_mode(self):
        for mode in AGGREGATE_MODES:
            assert aggregate_frequency(UNKNOWN_RSID, "A", "G", mode) is None


class TestAnnotate:
    def test_a_bare_finding_gets_every_documented_key(self):
        require_rsid(MTHFR)
        finding = {"rsid": MTHFR, "allele1": "C", "allele2": "T"}
        annotate(finding)
        for key in ANNOTATION_KEYS:
            assert key in finding

    def test_annotate_returns_the_same_dict(self):
        require_rsid(MTHFR)
        finding = {"rsid": MTHFR, "allele1": "C", "allele2": "T"}
        assert annotate(finding) is finding

    def test_an_empty_dict_is_annotated_without_raising(self):
        require_data()
        finding: dict = {}
        annotate(finding)
        for key in ANNOTATION_KEYS:
            assert key in finding
        assert finding["freq"] is None
        assert finding["freq_band"] == "unknown"
        assert finding["freq_color"] == UNKNOWN_COLOR

    def test_junk_values_do_not_raise(self):
        require_data()
        finding = {"rsid": 7, "allele1": None, "allele2": ["X"]}
        annotate(finding)
        assert finding["freq"] is None

    def test_a_non_dict_is_returned_unchanged(self):
        assert annotate("not a finding") == "not a finding"

    def test_the_band_and_color_agree_with_the_frequency(self):
        require_rsid(MTHFR)
        finding = {"rsid": MTHFR, "allele1": "C", "allele2": "T"}
        annotate(finding)
        assert finding["freq_band"] == rarity_band(finding["freq"])
        assert finding["freq_color"] == rarity_color(finding["freq"])

    def test_the_requested_population_is_recorded(self):
        require_rsid(MTHFR)
        finding = {"rsid": MTHFR, "allele1": "C", "allele2": "T"}
        annotate(finding, "FIN")
        assert finding["freq_population"] == "FIN"

    def test_the_flipped_lookup_still_produces_a_frequency(self):
        require_rsid(MTHFR)
        finding = {"rsid": MTHFR, "allele1": "C", "allele2": "T"}
        annotate(finding)
        assert finding["freq"] is not None
        assert finding["freq_method"] == "hardy_weinberg"

    def test_annotate_all_annotates_every_finding(self):
        require_rsid(MTHFR)
        findings = [{"rsid": MTHFR, "allele1": "C", "allele2": "T"},
                    {"rsid": UNKNOWN_RSID, "allele1": "A", "allele2": "G"}]
        annotate_all(findings)
        assert all("freq_band" in f for f in findings)

    def test_annotate_all_tolerates_none_and_empty(self):
        assert annotate_all([]) == []
        assert annotate_all(None) is None


class TestUnknownRsid:
    def test_gmaf_is_none(self):
        assert gmaf(UNKNOWN_RSID) is None

    def test_minor_allele_is_the_empty_string(self):
        assert minor_allele(UNKNOWN_RSID) == ""

    def test_allele_frequency_is_none(self):
        assert allele_frequency(UNKNOWN_RSID, "A") is None

    def test_genotype_frequency_is_none(self):
        assert genotype_frequency(UNKNOWN_RSID, "A", "G") is None

    def test_sample_size_is_none(self):
        assert sample_size(UNKNOWN_RSID) is None

    def test_observed_alleles_is_empty(self):
        assert observed_alleles(UNKNOWN_RSID) == set()

    def test_annotation_degrades_to_unknown_rather_than_raising(self):
        require_data()
        finding = {"rsid": UNKNOWN_RSID, "allele1": "A", "allele2": "G"}
        annotate(finding)
        assert finding["freq"] is None
        assert finding["freq_band"] == "unknown"
        assert finding["freq_color"] == UNKNOWN_COLOR
        assert finding["gmaf"] is None
        assert finding["minor_allele"] == ""

    def test_an_empty_rsid_is_handled_the_same_way(self):
        assert gmaf("") is None
        assert genotype_frequency("", "A", "G") is None
        assert observed_alleles("") == set()


class TestLoadingAndMetadata:
    def test_load_frequencies_returns_a_dict(self):
        assert isinstance(load_frequencies(), dict)

    def test_available_populations_are_declared_panels(self):
        require_data()
        for population in available_populations():
            assert population["code"] in POPULATION_CODES

    def test_the_default_population_is_available(self):
        require_data()
        codes = [p["code"] for p in available_populations()]
        assert DEFAULT_POPULATION in codes

    def test_the_metadata_declares_a_version_and_a_source(self):
        require_data()
        meta = frequency.get_metadata()
        assert meta.get("version")
        assert meta.get("source")

    def test_a_missing_data_file_returns_an_empty_table(self):
        absent = Path(__file__).parent / "no_such_frequencies.json"
        try:
            assert load_frequencies(absent) == {}
        finally:
            # Restore the module cache for the tests that follow.
            load_frequencies(frequency.FREQUENCY_FILE)

    def test_the_coverage_report_counts_what_is_bundled(self):
        data = require_data()
        report = frequency.coverage_report()
        assert report["rsids"] == len(data)
        assert report["populations"] == len(available_populations())
        assert report["with_gmaf"] >= 0
