"""Tests for backend.snpedia: licence gating, cache location and offline lookups.

Every test in this file must pass with the machine offline. An autouse fixture
replaces the module's HTTP entry points with functions that fail loudly, so any
accidental network access shows up as a test failure rather than a hang, and
every test that touches the filesystem is redirected into tmp_path.
"""

import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import snpedia
from backend.snpedia import (
    ANNOTATION_KEYS,
    API_BASE,
    CACHE_DIR,
    CACHE_PATH,
    LICENSE_NAME,
    LICENSE_URL,
    NOTICE,
    SCOPES,
    annotate,
    cache_path,
    cache_status,
    export_cache,
    harvest,
    harvest_genosets,
    init_cache,
    is_available,
    lookup,
    lookup_genoset,
    lookup_genotype,
    normalize_token,
    page_title,
    parse_population_diversity,
    purge_cache,
)

REPO_ROOT = Path(snpedia.__file__).resolve().parent.parent
REAL_CACHE_FILE = Path.home() / ".dnainsight" / "snpedia_cache.db"
REAL_CACHE_EXISTED = REAL_CACHE_FILE.exists()

DOCUMENTED_TABLES = ("meta", "snps", "genotypes", "genosets", "harvest_log")
DOCUMENTED_INDICES = ("idx_genotypes_rsid", "idx_snps_gene")
STATUS_KEYS = {
    "available", "path", "snps", "genotypes", "genosets",
    "last_harvest", "scope", "license", "notice",
}

POPULATION_TEMPLATE = """{{Rsnum
|rsid=1815739
|gene=ACTN3
|Chromosome=11
}}

{{ population diversity
|geno1=(C;C)
|geno2=(C;T)
|geno3=(T;T)
|CEU|46.9|44.2|8.8
|HCB|100.0|0.0|0.0
|CHB|100.0|0.0|0.0
|JPT|97.7|2.3|0.0
|YRI|55.0|38.3|6.7
|ASW|58.5|33.9|7.6
|CHD|93.5|6.5|0.0
|GIH|64.7|30.6|4.7
|LWK|53.5|40.4|6.1
|MEX|68.0|28.0|4.0
|MKK|60.9|33.1|6.0
|TSI|48.9|41.1|10.0
|HapMapRevision=28
}}

[[Category:Is a snp]]
"""

MALFORMED_TEMPLATE = """{{ population diversity
|geno1=(A;A)
|geno2=(A;G)
|CEU|46.9|44.2
|YRI|55.0|38.3|6.7
|HapMapRevision=28
}}
"""

PERCENT_TEMPLATE = """{{ population diversity
|geno1=(A;A)
|CEU|46.9%|44.2%|8.9%
|HapMapRevision=27
}}
"""


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Make any network call an immediate, obvious failure."""

    def forbidden(*args, **kwargs):
        raise AssertionError("this test tried to reach the network")

    monkeypatch.setattr(snpedia, "_request", forbidden)
    monkeypatch.setattr(snpedia, "_http_get", forbidden)
    monkeypatch.setattr(snpedia, "make_session", lambda: None)


@pytest.fixture
def cache(monkeypatch, tmp_path):
    """Redirect the module's cache constant into tmp_path."""
    path = tmp_path / "cache" / "snpedia_cache.db"
    monkeypatch.setattr(snpedia, "CACHE_PATH", path)
    monkeypatch.setattr(snpedia, "CACHE_DIR", path.parent)
    return path


class TestHarvestIsLicenceGated:
    def test_harvest_refuses_by_default(self):
        with pytest.raises(PermissionError):
            harvest()

    def test_harvest_refuses_when_the_licence_is_declined(self):
        with pytest.raises(PermissionError):
            harvest(scope="restricted", rsids=["rs53576"], accept_license=False)

    def test_harvest_refuses_a_truthy_non_true_value(self):
        for value in (1, "yes", "true", [1]):
            with pytest.raises(PermissionError):
                harvest(rsids=["rs53576"], accept_license=value)

    def test_harvest_refuses_none(self):
        with pytest.raises(PermissionError):
            harvest(rsids=["rs53576"], accept_license=None)

    def test_the_harvest_refusal_carries_the_notice(self):
        with pytest.raises(PermissionError) as info:
            harvest()
        assert str(info.value) == NOTICE

    def test_the_harvest_refusal_says_noncommercial(self):
        with pytest.raises(PermissionError) as info:
            harvest()
        assert "Noncommercial" in str(info.value)

    def test_the_harvest_refusal_carries_the_licence_url(self):
        with pytest.raises(PermissionError) as info:
            harvest()
        assert LICENSE_URL in str(info.value)

    def test_harvest_genosets_refuses_by_default(self):
        with pytest.raises(PermissionError):
            harvest_genosets()

    def test_harvest_genosets_refuses_when_the_licence_is_declined(self):
        with pytest.raises(PermissionError):
            harvest_genosets(accept_license=False)

    def test_harvest_genosets_refuses_a_truthy_non_true_value(self):
        for value in (1, "yes", 0.0):
            with pytest.raises(PermissionError):
                harvest_genosets(accept_license=value)

    def test_the_genoset_refusal_says_noncommercial_and_links_the_licence(self):
        with pytest.raises(PermissionError) as info:
            harvest_genosets()
        message = str(info.value)
        assert "Noncommercial" in message
        assert LICENSE_URL in message

    def test_the_licence_gate_runs_before_anything_else(self, cache):
        with pytest.raises(PermissionError):
            harvest(scope="all", accept_license=False)
        assert not cache.exists()

    def test_an_unknown_scope_is_rejected_without_touching_the_network(self, cache):
        with pytest.raises(ValueError):
            harvest(scope="everything", accept_license=True)

    def test_the_restricted_scope_demands_an_explicit_rsid_list(self, cache):
        with pytest.raises(ValueError):
            harvest(scope="restricted", rsids=None, accept_license=True)

    def test_the_documented_scopes(self):
        assert SCOPES == ("restricted", "chip_23andme_v5", "chip_ancestry_v2", "all")


class TestNotice:
    def test_the_notice_names_snpedia(self):
        assert "snpedia.com" in NOTICE

    def test_the_notice_links_the_licence(self):
        assert LICENSE_URL in NOTICE
        assert "by-nc-sa/3.0/us" in NOTICE

    def test_the_notice_names_the_licence(self):
        assert "Attribution-Noncommercial-Share" in NOTICE
        assert "3.0 United States" in NOTICE

    def test_the_notice_says_personal(self):
        assert "personal" in NOTICE

    def test_the_notice_says_non_commercial(self):
        assert "non-commercial" in NOTICE

    def test_the_notice_forbids_redistribution(self):
        assert "not redistribute" in NOTICE

    def test_the_notice_says_the_cache_lives_outside_the_repository(self):
        assert "outside the DNAInsight repository" in NOTICE
        assert "home" in NOTICE

    def test_the_notice_says_dnainsight_ships_no_snpedia_data(self):
        assert "ships no SNPedia data" in NOTICE

    def test_the_notice_requires_share_alike_on_any_derived_database(self):
        assert "must itself carry this same" in NOTICE

    def test_the_notice_forbids_probing_rs_numbers(self):
        assert "probe arbitrary rs numbers" in NOTICE

    def test_the_licence_name_constant(self):
        assert "Noncommercial" in LICENSE_NAME
        assert "3.0 United States" in LICENSE_NAME

    def test_the_harvester_targets_the_bots_host(self):
        assert API_BASE.startswith("https://bots.snpedia.com")


class TestCacheLivesOutsideTheRepository:
    def test_the_cache_path_is_under_the_home_directory(self):
        assert cache_path().is_relative_to(Path.home())

    def test_the_repository_path_never_appears_in_the_cache_path(self):
        assert str(REPO_ROOT) not in str(cache_path())

    def test_the_cache_path_is_not_inside_the_repository(self):
        assert not cache_path().is_relative_to(REPO_ROOT)

    def test_the_cache_directory_is_the_documented_dot_folder(self):
        assert CACHE_DIR == Path.home() / ".dnainsight"
        assert CACHE_PATH == CACHE_DIR / "snpedia_cache.db"

    def test_cache_path_returns_the_module_constant(self):
        assert cache_path() == Path(CACHE_PATH)

    def test_cache_path_returns_a_path_object(self):
        assert isinstance(cache_path(), Path)

    def test_no_cache_file_is_committed_anywhere_in_the_repository(self):
        assert list(REPO_ROOT.rglob("snpedia_cache.db")) == []

    def test_the_repository_ships_no_snpedia_licence_file(self):
        assert list(REPO_ROOT.rglob("LICENSE-SNPEDIA.txt")) == []

    def test_the_redirected_cache_is_still_outside_the_repository(self, cache):
        assert str(REPO_ROOT) not in str(cache_path())


class TestInitCache:
    def test_the_cache_file_is_created(self, cache):
        init_cache()
        assert cache.exists()

    def test_the_parent_directory_is_created(self, cache):
        assert not cache.parent.exists()
        init_cache()
        assert cache.parent.is_dir()

    def test_every_documented_table_exists(self, cache):
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            names = {row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'")}
        finally:
            conn.close()
        for table in DOCUMENTED_TABLES:
            assert table in names, table

    def test_the_frequencies_table_exists_for_the_diversity_template(self, cache):
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            names = {row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'")}
        finally:
            conn.close()
        assert "frequencies" in names

    def test_the_documented_indices_exist(self, cache):
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            names = {row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type = 'index'")}
        finally:
            conn.close()
        for index in DOCUMENTED_INDICES:
            assert index in names, index

    def test_the_notice_travels_with_the_data(self, cache):
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
        finally:
            conn.close()
        assert meta["notice"] == NOTICE
        assert meta["license"] == LICENSE_NAME
        assert meta["license_url"] == LICENSE_URL
        assert "snpedia.com" in meta["source"]

    def test_the_meta_table_records_the_redistribution_rule(self, cache):
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
        finally:
            conn.close()
        assert "Do not redistribute" in meta["redistribution"]
        assert meta["schema_version"]

    def test_init_cache_is_idempotent(self, cache):
        init_cache()
        init_cache()
        conn = sqlite3.connect(str(cache))
        try:
            rows = conn.execute("SELECT COUNT(*) FROM meta").fetchone()[0]
        finally:
            conn.close()
        assert rows == 7

    def test_init_cache_writes_nothing_to_the_real_home_directory(self, cache):
        init_cache()
        assert cache.exists()
        assert REAL_CACHE_FILE.exists() == REAL_CACHE_EXISTED


class TestCacheStatus:
    def test_a_missing_cache_reports_unavailable(self, cache):
        status = cache_status()
        assert status["available"] is False

    def test_a_missing_cache_still_returns_every_key(self, cache):
        assert set(cache_status()) == STATUS_KEYS

    def test_a_missing_cache_reports_zero_counts(self, cache):
        status = cache_status()
        assert status["snps"] == 0
        assert status["genotypes"] == 0
        assert status["genosets"] == 0

    def test_a_missing_cache_reports_the_path_it_looked_for(self, cache):
        assert cache_status()["path"] == str(cache)

    def test_a_fresh_empty_cache_is_still_unavailable(self, cache):
        init_cache()
        status = cache_status()
        assert status["available"] is False
        assert status["snps"] == 0

    def test_a_fresh_empty_cache_returns_every_key(self, cache):
        init_cache()
        assert set(cache_status()) == STATUS_KEYS

    def test_a_fresh_empty_cache_has_no_harvest_history(self, cache):
        init_cache()
        status = cache_status()
        assert status["last_harvest"] is None
        assert status["scope"] is None

    def test_the_status_always_carries_the_licence(self, cache):
        init_cache()
        status = cache_status()
        assert status["license"] == LICENSE_NAME
        assert status["notice"] == NOTICE

    def test_is_available_is_false_without_a_cache(self, cache):
        assert is_available() is False

    def test_is_available_is_false_for_an_empty_cache(self, cache):
        init_cache()
        assert is_available() is False


class TestLookupsWithoutACache:
    def test_lookup_returns_none_when_there_is_no_cache_file(self, cache):
        assert lookup("rs53576") is None

    def test_lookup_returns_none_against_an_empty_cache(self, cache):
        init_cache()
        assert lookup("rs53576") is None

    def test_lookup_tolerates_an_empty_rsid(self, cache):
        assert lookup("") is None
        assert lookup(None) is None

    def test_lookup_genotype_returns_none_when_there_is_no_cache_file(self, cache):
        assert lookup_genotype("rs53576", "A", "G") is None

    def test_lookup_genotype_returns_none_against_an_empty_cache(self, cache):
        init_cache()
        assert lookup_genotype("rs53576", "A", "G") is None

    def test_lookup_genotype_tolerates_an_empty_rsid(self, cache):
        assert lookup_genotype("", "A", "G") is None

    def test_lookup_genotype_tolerates_no_call_alleles(self, cache):
        init_cache()
        assert lookup_genotype("rs53576", "-", "-") is None

    def test_lookup_genoset_returns_none_when_there_is_no_cache_file(self, cache):
        assert lookup_genoset("gs123") is None

    def test_lookup_genoset_returns_none_against_an_empty_cache(self, cache):
        init_cache()
        assert lookup_genoset("gs123") is None

    def test_lookup_genoset_tolerates_an_empty_name(self, cache):
        assert lookup_genoset("") is None
        assert lookup_genoset(None) is None

    def test_no_lookup_creates_a_cache_file(self, cache):
        lookup("rs53576")
        lookup_genotype("rs53576", "A", "G")
        lookup_genoset("gs123")
        assert not cache.exists()


class TestAnnotate:
    def _finding(self):
        return {"rsid": "rs53576", "allele1": "A", "allele2": "G"}

    def test_annotate_returns_the_same_dict(self, cache):
        finding = self._finding()
        assert annotate(finding) is finding

    def test_annotate_never_raises_without_a_cache(self, cache):
        assert annotate(self._finding())["rsid"] == "rs53576"

    def test_annotate_sets_every_documented_key(self, cache):
        finding = annotate(self._finding())
        for key in ANNOTATION_KEYS:
            assert key in finding, key

    def test_annotate_leaves_the_original_call_alone(self, cache):
        finding = annotate(self._finding())
        assert (finding["allele1"], finding["allele2"]) == ("A", "G")

    def test_annotate_fills_nothing_in_when_the_cache_is_absent(self, cache):
        finding = annotate(self._finding())
        assert finding["magnitude"] is None
        assert finding["max_magnitude"] is None
        assert finding["publications"] is None
        assert finding["gmaf"] is None
        assert finding["orientation"] is None
        assert finding["stabilized_orientation"] is None

    def test_annotate_uses_empty_strings_and_lists_not_missing_keys(self, cache):
        finding = annotate(self._finding())
        assert finding["repute"] == ""
        assert finding["summary"] == ""
        assert finding["snpedia_topics"] == []
        assert finding["snpedia_medicines"] == []
        assert finding["snpedia_conditions"] == []

    def test_annotate_defaults_the_strand_flags_to_false(self, cache):
        finding = annotate(self._finding())
        assert finding["flipped"] is False
        assert finding["ambiguous"] is False

    def test_annotate_still_links_to_the_snpedia_page(self, cache):
        finding = annotate(self._finding())
        assert "snpedia.com" in finding["snpedia_url"]
        assert finding["snpedia_url"].endswith("Rs53576")

    def test_annotate_does_not_link_a_non_rs_identifier(self, cache):
        finding = annotate({"rsid": "blood_type", "allele1": "N", "allele2": "N"})
        assert finding["snpedia_url"] is None

    def test_annotate_tolerates_an_empty_finding(self, cache):
        finding = annotate({})
        assert set(ANNOTATION_KEYS).issubset(finding)
        assert finding["snpedia_url"] is None

    def test_annotate_tolerates_a_missing_rsid(self, cache):
        finding = annotate({"allele1": "A", "allele2": "G"})
        assert finding["magnitude"] is None

    def test_annotate_tolerates_a_none_rsid(self, cache):
        assert annotate({"rsid": None})["repute"] == ""

    def test_annotate_of_a_non_dict_returns_all_keys_as_none(self, cache):
        result = annotate("rs53576")
        assert set(result) == set(ANNOTATION_KEYS)
        assert all(value is None for value in result.values())

    def test_annotate_preserves_an_existing_summary(self, cache):
        finding = annotate({"rsid": "rs53576", "summary": "already known"})
        assert finding["summary"] == "already known"

    def test_annotate_against_an_empty_cache_behaves_the_same(self, cache):
        init_cache()
        finding = annotate(self._finding())
        assert finding["magnitude"] is None
        assert finding["repute"] == ""
        assert finding["flipped"] is False

    def test_annotate_creates_no_cache_file(self, cache):
        annotate(self._finding())
        assert not cache.exists()


class TestSmallHelpers:
    def test_page_title_capitalises_the_first_letter_only(self):
        assert page_title("rs53576") == "Rs53576"
        assert page_title("i12345") == "I12345"

    def test_page_title_tolerates_empty_input(self):
        assert page_title("") == ""
        assert page_title(None) == ""

    def test_page_title_leaves_an_already_capitalised_title(self):
        assert page_title("Rs53576") == "Rs53576"

    def test_normalize_token_uppercases_and_removes_spaces(self):
        assert normalize_token("(a; g)") == "(A;G)"

    def test_normalize_token_closes_an_unclosed_token(self):
        assert normalize_token("(A;G") == "(A;G)"

    def test_normalize_token_tolerates_empty_input(self):
        assert normalize_token(None) == ""
        assert normalize_token("") == ""

    def test_the_annotation_keys_constant_is_a_tuple(self):
        assert isinstance(ANNOTATION_KEYS, tuple)
        assert "magnitude" in ANNOTATION_KEYS
        assert "snpedia_url" in ANNOTATION_KEYS


class TestParsePopulationDiversity:
    def test_the_result_carries_the_documented_keys(self):
        parsed = parse_population_diversity(POPULATION_TEMPLATE)
        assert set(parsed) == {"genotypes", "populations", "revision"}

    def test_the_three_genotypes_are_parsed_in_order(self):
        parsed = parse_population_diversity(POPULATION_TEMPLATE)
        assert parsed["genotypes"] == ["(C;C)", "(C;T)", "(T;T)"]

    def test_the_hapmap_revision_is_parsed(self):
        assert parse_population_diversity(POPULATION_TEMPLATE)["revision"] == "28"

    def test_every_population_row_is_parsed(self):
        populations = parse_population_diversity(POPULATION_TEMPLATE)["populations"]
        assert len(populations) == 12
        assert set(populations) == {
            "CEU", "HCB", "CHB", "JPT", "YRI", "ASW",
            "CHD", "GIH", "LWK", "MEX", "MKK", "TSI",
        }

    def test_the_percentages_are_floats_in_page_order(self):
        populations = parse_population_diversity(POPULATION_TEMPLATE)["populations"]
        assert populations["CEU"] == [46.9, 44.2, 8.8]
        assert populations["YRI"] == [55.0, 38.3, 6.7]

    def test_every_row_has_exactly_three_percentages(self):
        populations = parse_population_diversity(POPULATION_TEMPLATE)["populations"]
        for code, values in populations.items():
            assert len(values) == 3, code
            for value in values:
                assert isinstance(value, float), code

    def test_the_legacy_hcb_code_is_kept_alongside_chb(self):
        populations = parse_population_diversity(POPULATION_TEMPLATE)["populations"]
        assert "HCB" in populations
        assert "CHB" in populations
        assert populations["HCB"] == populations["CHB"]
        assert populations["HCB"] == [100.0, 0.0, 0.0]

    def test_named_parameters_never_leak_into_the_populations(self):
        populations = parse_population_diversity(POPULATION_TEMPLATE)["populations"]
        for leaked in ("geno1", "geno2", "geno3", "HapMapRevision", "(C;C)"):
            assert leaked not in populations

    def test_the_template_name_is_matched_case_insensitively(self):
        text = "{{ Population Diversity\n|geno1=(A;A)\n|CEU|10.0|20.0|70.0\n}}\n"
        parsed = parse_population_diversity(text)
        assert parsed["genotypes"] == ["(A;A)"]
        assert parsed["populations"] == {"CEU": [10.0, 20.0, 70.0]}

    def test_a_missing_revision_is_an_empty_string(self):
        text = "{{ population diversity\n|geno1=(A;A)\n|CEU|10.0|20.0|70.0\n}}\n"
        assert parse_population_diversity(text)["revision"] == ""

    def test_percent_signs_are_tolerated(self):
        parsed = parse_population_diversity(PERCENT_TEMPLATE)
        assert parsed["populations"] == {"CEU": [46.9, 44.2, 8.9]}
        assert parsed["revision"] == "27"

    def test_a_page_without_the_template_yields_empty_containers(self):
        parsed = parse_population_diversity("{{Rsnum\n|rsid=1815739\n}}\n")
        assert parsed == {"genotypes": [], "populations": {}, "revision": ""}

    def test_empty_input_yields_empty_containers(self):
        assert parse_population_diversity("") == {
            "genotypes": [], "populations": {}, "revision": "",
        }

    def test_none_input_yields_empty_containers(self):
        assert parse_population_diversity(None)["populations"] == {}

    def test_a_malformed_row_does_not_crash_the_parser(self):
        parsed = parse_population_diversity(MALFORMED_TEMPLATE)
        assert isinstance(parsed["populations"], dict)

    def test_a_malformed_row_is_skipped_rather_than_guessed(self):
        parsed = parse_population_diversity(MALFORMED_TEMPLATE)
        assert "CEU" not in parsed["populations"]

    def test_a_malformed_row_does_not_stop_later_rows(self):
        parsed = parse_population_diversity(MALFORMED_TEMPLATE)
        assert parsed["populations"]["YRI"] == [55.0, 38.3, 6.7]
        assert parsed["revision"] == "28"

    def test_a_truncated_template_does_not_crash(self):
        text = "{{ population diversity\n|geno1=(A;A)\n|CEU|10.0|20.0"
        parsed = parse_population_diversity(text)
        assert parsed["populations"] == {}
        assert parsed["genotypes"] == ["(A;A)"]

    def test_parsing_touches_no_files(self, cache):
        parse_population_diversity(POPULATION_TEMPLATE)
        assert not cache.exists()


class TestExportCache:
    def test_export_refuses_when_there_is_no_cache(self, cache, tmp_path):
        with pytest.raises(FileNotFoundError):
            export_cache(tmp_path / "out" / "copy.db")

    def test_export_copies_the_database(self, cache, tmp_path):
        init_cache()
        target = tmp_path / "out" / "copy.db"
        result = export_cache(target)
        assert result["exported"] is True
        assert Path(result["destination"]) == target
        assert target.exists()
        assert result["bytes"] > 0

    def test_export_writes_a_sibling_licence_file(self, cache, tmp_path):
        init_cache()
        result = export_cache(tmp_path / "out" / "copy.db")
        licence = Path(result["license_file"])
        assert licence.name == "LICENSE-SNPEDIA.txt"
        assert licence.parent == tmp_path / "out"
        assert licence.exists()

    def test_the_licence_file_contains_the_notice(self, cache, tmp_path):
        init_cache()
        result = export_cache(tmp_path / "out" / "copy.db")
        text = Path(result["license_file"]).read_text(encoding="utf-8")
        assert text == NOTICE
        assert "Noncommercial" in text
        assert LICENSE_URL in text

    def test_export_reports_the_licence_terms(self, cache, tmp_path):
        init_cache()
        result = export_cache(tmp_path / "out" / "copy.db")
        assert result["license"] == LICENSE_NAME
        assert result["license_url"] == LICENSE_URL
        assert result["source"] == str(cache)

    def test_export_to_a_directory_keeps_the_cache_filename(self, cache, tmp_path):
        init_cache()
        destination = tmp_path / "outdir"
        destination.mkdir()
        result = export_cache(destination)
        assert Path(result["destination"]) == destination / cache.name
        assert (destination / "LICENSE-SNPEDIA.txt").exists()

    def test_export_creates_the_destination_directory(self, cache, tmp_path):
        init_cache()
        result = export_cache(tmp_path / "deep" / "nested" / "copy.db")
        assert Path(result["destination"]).exists()

    def test_the_exported_copy_is_a_usable_database(self, cache, tmp_path):
        init_cache()
        result = export_cache(tmp_path / "out" / "copy.db")
        conn = sqlite3.connect(result["destination"])
        try:
            value = conn.execute(
                "SELECT value FROM meta WHERE key = 'notice'").fetchone()[0]
        finally:
            conn.close()
        assert value == NOTICE

    def test_export_writes_nothing_into_the_repository(self, cache, tmp_path):
        init_cache()
        export_cache(tmp_path / "out" / "copy.db")
        assert list(REPO_ROOT.rglob("LICENSE-SNPEDIA.txt")) == []


class TestPurgeCache:
    def test_purge_removes_the_cache_file(self, cache):
        init_cache()
        assert cache.exists()
        result = purge_cache()
        assert result["removed"] is True
        assert not cache.exists()

    def test_purge_reports_the_documented_keys(self, cache):
        init_cache()
        result = purge_cache()
        assert set(result) == {"removed", "path", "snps", "genotypes", "genosets"}

    def test_purge_reports_what_was_in_the_cache(self, cache):
        init_cache()
        result = purge_cache()
        assert result["path"] == str(cache)
        assert result["snps"] == 0
        assert result["genotypes"] == 0
        assert result["genosets"] == 0

    def test_purging_twice_reports_nothing_removed(self, cache):
        init_cache()
        purge_cache()
        result = purge_cache()
        assert result["removed"] is False

    def test_purge_without_a_cache_is_harmless(self, cache):
        result = purge_cache()
        assert result["removed"] is False
        assert result["path"] == str(cache)

    def test_purge_leaves_the_status_reporting_unavailable(self, cache):
        init_cache()
        purge_cache()
        assert cache_status()["available"] is False
        assert is_available() is False


class TestTheRealHomeDirectoryIsUntouched:
    def test_the_redirected_cache_really_is_used(self, cache, tmp_path):
        init_cache()
        assert str(tmp_path) in cache_status()["path"]

    def test_no_test_created_the_real_cache_file(self):
        assert REAL_CACHE_FILE.exists() == REAL_CACHE_EXISTED

    def test_the_default_cache_path_is_still_under_the_home_directory(self):
        assert cache_path() == REAL_CACHE_FILE
        assert cache_path().is_relative_to(Path.home())

    def test_no_snpedia_artefact_was_left_in_the_repository(self):
        assert list(REPO_ROOT.rglob("snpedia_cache.db")) == []
        assert list(REPO_ROOT.rglob("LICENSE-SNPEDIA.txt")) == []
