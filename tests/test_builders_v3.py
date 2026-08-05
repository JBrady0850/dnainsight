"""Tests for the v3.0 data builders and the external tool manifest.

Three things are being protected here and they are worth naming, because a
reader who does not know why these tests exist will delete half of them:

  1. NO NETWORK, EVER, WITHOUT CONSENT. Both builders must be importable and
     --help-able offline, and a run without --accept-terms must be a dry run
     that fetches nothing and writes nothing. Every test that could otherwise
     reach the network monkeypatches requests so that any call raises, which
     means a regression shows up as a failure rather than as a slow test.

  2. THE LICENCE GATES. The HGDP double opt-in, the permanent refusal of the
     restricted SGDP tier and of the Allen Ancient DNA Resource, and the CPIC
     forbidden-column assertion. These are the parts of the builders that exist
     for legal rather than technical reasons, which makes them exactly the
     parts somebody will "simplify" later.

  3. MANIFEST DRIFT. data/tools_manifest.json is what a licence auditor reads
     first. The set equality tests below fail in BOTH directions, so adding a
     tool to backend/external.py without regenerating the manifest is a test
     failure and not a quiet inconsistency.

No test in this file touches the network, writes into the repository, or needs
an external tool installed.
"""

import importlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

import requests

from backend import ancestry, external
from backend.diplotype import unverified_entries
from data import build_panel, build_pgx_alleles


ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = ROOT / "data" / "tools_manifest.json"
NEW_FILES = (
    ROOT / "data" / "build_panel.py",
    ROOT / "data" / "build_pgx_alleles.py",
    ROOT / "data" / "tools_manifest.json",
    ROOT / "tests" / "test_builders_v3.py",
)


class NetworkTouched(AssertionError):
    """Raised by the stubbed transport when a test tries to reach the network."""


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """Point ~/.dnainsight at a fresh temp directory for every test.

    Same fixture the ancestry and imputation suites use. Without it a builder
    test would write into the developer's real panel directory.
    """
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    external.reset_cache()
    yield
    external.reset_cache()


@pytest.fixture()
def no_network(monkeypatch):
    """Make every HTTP verb raise, so an accidental fetch fails loudly."""

    def boom(*args, **kwargs):
        raise NetworkTouched(f"network call attempted: {args[:2]}")

    for verb in ("get", "post", "put", "head", "request"):
        monkeypatch.setattr(requests, verb, boom, raising=False)
    monkeypatch.setattr(requests.Session, "request", boom, raising=False)
    return boom


def repo_snapshot() -> set:
    """Every tracked-looking file in the repository, with its size.

    Used to prove a dry run wrote nothing into the tree. __pycache__ and the
    pytest cache are excluded because importing a module legitimately creates
    them and they are gitignored anyway.
    """
    skip = {".git", "__pycache__", ".pytest_cache", "uploads",
            "reports_output", ".venv", "venv"}
    out = set()
    for path in ROOT.rglob("*"):
        if any(part in skip for part in path.parts):
            continue
        if path.is_file():
            out.add((str(path.relative_to(ROOT)), path.stat().st_size))
    return out


# ===========================================================================
# 1. Import and --help work offline
# ===========================================================================

def test_build_panel_imports_without_network(no_network):
    module = importlib.reload(importlib.import_module("data.build_panel"))
    assert module.PANEL_ID == "onekg_sgdp"


def test_build_pgx_imports_without_network(no_network):
    module = importlib.reload(importlib.import_module("data.build_pgx_alleles"))
    assert module.CPIC_SPDX == "CC0-1.0"


def test_build_panel_help_exits_zero():
    result = subprocess.run(
        [sys.executable, "data/build_panel.py", "--help"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert "--accept-terms" in result.stdout


def test_build_pgx_help_exits_zero():
    result = subprocess.run(
        [sys.executable, "data/build_pgx_alleles.py", "--help"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert "--accept-terms" in result.stdout


def test_panel_parser_carries_every_documented_flag():
    parser = build_panel.build_parser()
    flags = {action.option_strings[0] for action in parser._actions
             if action.option_strings}
    for expected in ("--accept-terms", "--dry-run", "--array-file",
                     "--chromosomes", "--statistic", "--min-maf",
                     "--include-hgdp", "--accept-consent-caveat",
                     "--skip-sgdp", "--sgdp-dir", "--sgdp-vcf-template",
                     "--sgdp-metadata", "--onekg-source", "--limit",
                     "--markers-per-population", "--panel-id"):
        assert expected in flags


def test_pgx_parser_carries_every_documented_flag():
    parser = build_pgx_alleles.build_parser()
    flags = {action.option_strings[0] for action in parser._actions
             if action.option_strings}
    for expected in ("--accept-terms", "--dry-run", "--genes", "--out",
                     "--report", "--fail-on-conflict"):
        assert expected in flags


# ===========================================================================
# 2. The consent gate: no --accept-terms means dry run and no fetch
# ===========================================================================

def test_panel_without_accept_terms_fetches_nothing(no_network, capsys):
    assert build_panel.main([]) == 0
    assert "DRY RUN" in capsys.readouterr().out


def test_pgx_without_accept_terms_fetches_nothing(no_network, capsys):
    assert build_pgx_alleles.main([]) == 0
    assert "DRY RUN" in capsys.readouterr().out


def test_panel_dry_run_prints_every_url(no_network, capsys):
    build_panel.main([])
    out = capsys.readouterr().out
    assert build_panel.ONEKG_PANEL_URL in out
    assert build_panel.GENETIC_MAP_URL in out
    assert build_panel.SGDP_METADATA_URL in out


def test_panel_dry_run_prints_the_licences(no_network, capsys):
    build_panel.main([])
    out = capsys.readouterr().out
    assert "EMBL-EBI terms" in out
    assert "Simons Genome Diversity Project" in out


def test_panel_dry_run_prints_an_estimated_size(no_network, capsys):
    build_panel.main(["--chromosomes", "22"])
    out = capsys.readouterr().out
    assert "Total of the measurable parts" in out
    assert "MiB" in out or "GiB" in out


def test_pgx_dry_run_prints_the_endpoints_and_licence(no_network, capsys):
    build_pgx_alleles.main([])
    out = capsys.readouterr().out
    assert build_pgx_alleles.PAIR_URL in out
    assert "CC0-1.0" in out


def test_panel_dry_run_writes_nothing_into_the_repository(no_network):
    before = repo_snapshot()
    build_panel.main(["--chromosomes", "21,22"])
    assert repo_snapshot() == before


def test_pgx_dry_run_writes_nothing_into_the_repository(no_network):
    before = repo_snapshot()
    build_pgx_alleles.main([])
    assert repo_snapshot() == before


def test_panel_dry_run_writes_nothing_into_the_panel_root(no_network):
    build_panel.main([])
    assert not Path(external.panel_root()).exists()


def test_pgx_dry_run_does_not_create_its_output(no_network, tmp_path):
    target = tmp_path / "pgx_alleles.json"
    build_pgx_alleles.main(["--out", str(target)])
    assert not target.exists()


def test_dry_run_flag_forces_a_dry_run_even_with_accept_terms(no_network, capsys):
    assert build_panel.main(["--accept-terms", "--dry-run"]) == 0
    assert "DRY RUN" in capsys.readouterr().out


def test_pgx_dry_run_flag_forces_a_dry_run_even_with_accept_terms(no_network,
                                                                 capsys):
    assert build_pgx_alleles.main(["--accept-terms", "--dry-run"]) == 0
    assert "DRY RUN" in capsys.readouterr().out


def test_accept_terms_actually_reaches_the_transport(no_network):
    # Proves the gate is the only thing holding the builder back: with consent
    # given, the very next thing it does is a network call.
    with pytest.raises(NetworkTouched):
        build_panel.main(["--accept-terms", "--chromosomes", "22",
                          "--skip-sgdp"])


# ===========================================================================
# 3. The HGDP double gate
# ===========================================================================

def test_hgdp_refused_when_neither_flag_is_passed():
    allowed, reason = build_panel.hgdp_gate(False, False)
    assert allowed is False
    assert "REFUSED BY DEFAULT" in reason


def test_hgdp_refused_with_include_flag_alone():
    allowed, reason = build_panel.hgdp_gate(True, False)
    assert allowed is False
    assert "--accept-consent-caveat" in reason


def test_hgdp_refused_with_consent_flag_alone():
    allowed, reason = build_panel.hgdp_gate(False, True)
    assert allowed is False
    assert "--include-hgdp" in reason


def test_hgdp_allowed_only_with_both_flags():
    allowed, reason = build_panel.hgdp_gate(True, True)
    assert allowed is True
    assert "double opt-in" in reason


def test_hgdp_refusal_cites_nature_genetics_and_primed():
    text = build_panel.REFUSAL_HGDP
    assert "Nature Genetics" in text
    assert "24 November 2025" in text
    assert "PRIMED" in text
    assert "21 August 2024" in text


def test_hgdp_refusal_says_it_is_not_a_licence_objection():
    assert "NOT a licence objection" in build_panel.REFUSAL_HGDP


def test_hgdp_half_gate_exits_two_and_never_fetches(no_network, capsys):
    assert build_panel.main(["--accept-terms", "--include-hgdp"]) == 2
    assert "half-answered consent gate" in capsys.readouterr().out


def test_hgdp_consent_flag_alone_exits_two(no_network, capsys):
    assert build_panel.main(["--accept-terms", "--accept-consent-caveat"]) == 2
    assert "half-answered consent gate" in capsys.readouterr().out


def test_hgdp_absent_from_the_plan_by_default():
    args = build_panel.build_parser().parse_args(["--chromosomes", "22"])
    plan = build_panel.download_plan(args)
    assert plan["hgdp_included"] is False
    assert "hgdp" in plan["refusals"]


def test_hgdp_present_in_the_plan_with_both_flags():
    args = build_panel.build_parser().parse_args(
        ["--chromosomes", "22", "--include-hgdp", "--accept-consent-caveat"])
    plan = build_panel.download_plan(args)
    assert plan["hgdp_included"] is True
    assert "hgdp" not in plan["refusals"]
    assert any("HGDP" in item["source"] for item in plan["items"])


# ===========================================================================
# 4. The permanent refusals
# ===========================================================================

def test_restricted_sgdp_refusal_quotes_the_commercial_term():
    text = build_panel.refuse_restricted_sgdp()
    assert "will not use the data for any commercial purposes" in text
    assert "21-sample" in text


def test_restricted_sgdp_has_no_enabling_flag():
    parser = build_panel.build_parser()
    flags = " ".join(f for action in parser._actions
                     for f in action.option_strings)
    assert "restricted" not in flags


def test_restricted_sgdp_is_refused_in_every_plan():
    for extra in ([], ["--skip-sgdp"], ["--include-hgdp",
                                        "--accept-consent-caveat"]):
        args = build_panel.build_parser().parse_args(["--chromosomes", "22",
                                                      *extra])
        plan = build_panel.download_plan(args)
        assert "commercial purposes" in plan["refusals"]["sgdp_restricted"]


def test_aadr_refused_entirely():
    text = build_panel.refuse_aadr()
    assert "Allen Ancient DNA Resource" in text
    assert "Unread is not permissive" in text
    assert build_panel.AADR_URL in text


def test_aadr_has_no_enabling_flag():
    parser = build_panel.build_parser()
    flags = " ".join(f for action in parser._actions
                     for f in action.option_strings)
    assert "aadr" not in flags.lower()


def test_aadr_is_refused_in_every_plan():
    args = build_panel.build_parser().parse_args([])
    assert "aadr" in build_panel.download_plan(args)["refusals"]


def test_all_three_refusals_are_printed_in_a_dry_run(no_network, capsys):
    build_panel.main([])
    out = capsys.readouterr().out
    assert "WILL NOT DOWNLOAD" in out
    assert "Allen Ancient DNA Resource" in out
    assert "restricted tier" in out
    assert "Human Genome Diversity Project" in out


def test_panel_refusal_text_agrees_with_external_registry():
    # The registry and the builder must not tell two different stories about
    # why the same thing was excluded.
    excluded = external.PANELS["onekg_sgdp"]["excluded"]
    assert set(excluded) == set(build_panel.REFUSALS)


# ===========================================================================
# 5. The CPIC forbidden-column assertion
# ===========================================================================

def test_forbidden_columns_are_the_two_pharmgkb_ones():
    assert build_pgx_alleles.FORBIDDEN_COLUMNS == ("clinpgxlevel", "pgxtesting")


def test_pair_select_names_only_the_three_permitted_columns():
    assert build_pgx_alleles.PAIR_SELECT == "genesymbol,drugname,cpiclevel"


def test_no_select_list_mentions_a_forbidden_column():
    selects = (build_pgx_alleles.PAIR_SELECT,
               build_pgx_alleles.ALLELE_SELECT,
               build_pgx_alleles.ALLELE_DEFINITION_SELECT,
               build_pgx_alleles.SEQUENCE_LOCATION_SELECT,
               build_pgx_alleles.ALLELE_LOCATION_VALUE_SELECT)
    for select in selects:
        for column in build_pgx_alleles.FORBIDDEN_COLUMNS:
            assert column not in select


def test_clinpgxlevel_in_a_row_raises():
    row = {"genesymbol": "CYP2C19", "drugname": "clopidogrel",
           "cpiclevel": "A", "clinpgxlevel": "1A"}
    with pytest.raises(build_pgx_alleles.ForbiddenColumn) as exc:
        build_pgx_alleles.assert_no_forbidden_columns(row, "pair_view")
    assert exc.value.column == "clinpgxlevel"


def test_pgxtesting_in_a_row_raises():
    with pytest.raises(build_pgx_alleles.ForbiddenColumn):
        build_pgx_alleles.assert_no_forbidden_columns({"pgxtesting": None})


def test_forbidden_column_check_is_case_insensitive():
    with pytest.raises(build_pgx_alleles.ForbiddenColumn):
        build_pgx_alleles.assert_no_forbidden_columns({"ClinPGxLevel": "1A"})


def test_forbidden_column_check_ignores_surrounding_whitespace():
    with pytest.raises(build_pgx_alleles.ForbiddenColumn):
        build_pgx_alleles.assert_no_forbidden_columns({" pgxtesting ": "x"})


def test_forbidden_column_error_explains_the_licence():
    with pytest.raises(build_pgx_alleles.ForbiddenColumn) as exc:
        build_pgx_alleles.assert_no_forbidden_columns({"clinpgxlevel": "1A"})
    assert "PharmGKB" in str(exc.value)
    assert "prohibiting sale" in str(exc.value)


def test_a_clean_row_passes_and_is_returned_unchanged():
    row = {"genesymbol": "TPMT", "drugname": "azathioprine", "cpiclevel": "A"}
    assert build_pgx_alleles.assert_no_forbidden_columns(row) is row


def test_check_rows_raises_on_the_one_bad_row_in_a_batch():
    rows = [{"genesymbol": "TPMT"}, {"genesymbol": "DPYD"},
            {"genesymbol": "CYP2D6", "pgxtesting": "Actionable"}]
    with pytest.raises(build_pgx_alleles.ForbiddenColumn):
        build_pgx_alleles.check_rows(rows, "allele")


def test_check_rows_returns_every_row_when_all_are_clean():
    rows = [{"genesymbol": "TPMT"}, {"genesymbol": "DPYD"}]
    assert build_pgx_alleles.check_rows(rows) == rows


def test_a_request_without_a_select_list_is_refused(no_network):
    with pytest.raises(ValueError):
        build_pgx_alleles._get(build_pgx_alleles.PAIR_URL, "", "pair_view")


# ===========================================================================
# 6. The marker statistics, against hand-computed examples
# ===========================================================================

def test_fst_is_one_for_a_fixed_difference():
    # Two populations fixed for opposite alleles: Hs = 0, Ht = 0.5, Fst = 1.
    assert build_panel.wright_fst({"A": (20, 20), "B": (0, 20)}) == pytest.approx(1.0)


def test_fst_is_zero_when_frequencies_are_identical():
    # p = 0.5 in both: Hs = Ht = 0.5, so (Ht - Hs) / Ht = 0.
    assert build_panel.wright_fst({"A": (10, 20), "B": (10, 20)}) == pytest.approx(0.0)


def test_fst_hand_computed_equal_sample_sizes():
    # p1 = 0.8, p2 = 0.2, equal n. pbar = 0.5 so Ht = 0.5.
    # Hs = (2*0.8*0.2 + 2*0.2*0.8) / 2 = 0.32. Fst = (0.5 - 0.32) / 0.5 = 0.36.
    assert build_panel.wright_fst(
        {"A": (16, 20), "B": (4, 20)}) == pytest.approx(0.36)


def test_fst_hand_computed_unequal_sample_sizes():
    # p1 = 0.8 over 20 copies, p2 = 0.2 over 10 copies.
    # pbar = 18/30 = 0.6, Ht = 2*0.6*0.4 = 0.48.
    # Hs = (20*0.32 + 10*0.32) / 30 = 0.32.
    # Fst = (0.48 - 0.32) / 0.48 = 1/3.
    assert build_panel.wright_fst(
        {"A": (16, 20), "B": (2, 10)}) == pytest.approx(1.0 / 3.0)


def test_fst_is_zero_for_an_invariant_site():
    assert build_panel.wright_fst({"A": (0, 20), "B": (0, 20)}) == 0.0


def test_fst_is_zero_with_fewer_than_two_populations():
    assert build_panel.wright_fst({"A": (10, 20)}) == 0.0
    assert build_panel.wright_fst({}) == 0.0


def test_fst_ignores_populations_with_no_called_copies():
    with_empty = build_panel.wright_fst(
        {"A": (20, 20), "B": (0, 20), "C": (0, 0)})
    without = build_panel.wright_fst({"A": (20, 20), "B": (0, 20)})
    assert with_empty == pytest.approx(without)


def test_population_fst_is_one_versus_rest():
    counts = {"A": (20, 20), "B": (0, 20), "C": (0, 20)}
    # A against the pooled rest, which is 0 alt out of 40 copies.
    expected = build_panel.wright_fst({"A": (20, 20), "REST": (0, 40)})
    assert build_panel.population_fst(counts, "A") == pytest.approx(expected)


def test_population_fst_is_zero_for_an_absent_population():
    assert build_panel.population_fst({"A": (10, 20)}, "ZZZ") == 0.0


def test_informativeness_is_ln_two_for_a_fixed_difference():
    import math
    value = build_panel.informativeness({"A": (20, 20), "B": (0, 20)})
    assert value == pytest.approx(math.log(2.0))


def test_informativeness_is_zero_when_frequencies_are_identical():
    assert build_panel.informativeness(
        {"A": (10, 20), "B": (10, 20)}) == pytest.approx(0.0)


def test_informativeness_hand_computed():
    # p1 = 0.8, p2 = 0.2, pbar = 0.5 for both alleles by symmetry.
    # Per allele: -0.5 ln 0.5 + 0.5(0.8 ln 0.8) + 0.5(0.2 ln 0.2)
    #           = 0.346574 - 0.089257 - 0.160944 = 0.096373
    # Two alleles: In = 0.192745 nats.
    assert build_panel.informativeness(
        {"A": (16, 20), "B": (4, 20)}) == pytest.approx(0.192745, abs=1e-6)


def test_informativeness_uses_an_unweighted_mean_unlike_fst():
    # With unequal sample sizes the two statistics must disagree, because Fst
    # weights by called copies and In does not. If they ever agree here, one of
    # them has quietly changed its definition.
    counts = {"A": (16, 20), "B": (2, 10)}
    assert build_panel.informativeness(counts) != pytest.approx(
        build_panel.wright_fst(counts))


def test_marker_statistic_dispatches_on_name():
    counts = {"A": (16, 20), "B": (4, 20)}
    assert build_panel.marker_statistic(counts, "fst") == pytest.approx(0.36)
    assert build_panel.marker_statistic(counts, "informativeness") == \
        pytest.approx(build_panel.informativeness(counts))


def test_marker_statistic_rejects_an_unknown_name():
    with pytest.raises(ValueError):
        build_panel.marker_statistic({"A": (1, 2)}, "vibes")


def test_allele_frequency_handles_a_zero_denominator():
    assert build_panel.allele_frequency((0, 0)) == 0.0
    assert build_panel.allele_frequency((3, 12)) == pytest.approx(0.25)


# ===========================================================================
# 7. Bounded selection and small helpers
# ===========================================================================

def test_topn_keeps_only_the_highest_scores():
    top = build_panel.TopN(2)
    for score, name in ((0.1, "a"), (0.9, "b"), (0.5, "c"), (0.2, "d")):
        top.add(score, name)
    assert [name for _score, name in top.items()] == ["b", "c"]


def test_topn_breaks_ties_by_insertion_order():
    top = build_panel.TopN(2)
    for name in ("first", "second", "third"):
        top.add(0.5, name)
    assert [name for _score, name in top.items()] == ["first", "second"]


def test_topn_of_size_zero_keeps_nothing():
    top = build_panel.TopN(0)
    top.add(1.0, "x")
    assert top.items() == []


def test_interpolate_cm_is_linear_between_two_map_points():
    assert build_panel.interpolate_cm(150, [100, 200], [1.0, 2.0]) == \
        pytest.approx(1.5)


def test_interpolate_cm_clamps_outside_the_mapped_range():
    assert build_panel.interpolate_cm(10, [100, 200], [1.0, 2.0]) == 1.0
    assert build_panel.interpolate_cm(900, [100, 200], [1.0, 2.0]) == 2.0


def test_interpolate_cm_with_no_map_returns_zero():
    assert build_panel.interpolate_cm(150, [], []) == 0.0


def test_sha256_of_a_known_file(tmp_path):
    import hashlib
    target = tmp_path / "x.txt"
    target.write_bytes(b"dnainsight")
    assert build_panel.sha256_of(target) == \
        hashlib.sha256(b"dnainsight").hexdigest()


# ===========================================================================
# 8. Source parsing
# ===========================================================================

PANEL_TEXT = (
    "sample\tpop\tsuper_pop\tgender\n"
    "HG00096\tGBR\tEUR\tmale\n"
    "HG00097\tGBR\tEUR\tfemale\n"
    "NA18525\tCHB\tEAS\tfemale\n"
)


def test_parse_panel_file_reads_population_and_superpopulation():
    parsed = build_panel.parse_panel_file(PANEL_TEXT.splitlines())
    assert parsed["HG00096"] == ("GBR", "EUR")
    assert parsed["NA18525"] == ("CHB", "EAS")
    assert len(parsed) == 3


def test_parse_panel_file_resolves_columns_by_name_not_position():
    reordered = ("gender\tsuper_pop\tpop\tsample\n"
                 "male\tEUR\tGBR\tHG00096\n")
    parsed = build_panel.parse_panel_file(reordered.splitlines())
    assert parsed["HG00096"] == ("GBR", "EUR")


def test_parse_sgdp_metadata_reads_sample_and_population():
    csv = ("Sample ID (Illumina),Population ID,Region\n"
           "LP6005441-DNA_A01,Yoruba,Africa\n"
           "LP6005441-DNA_B01,French,West Eurasia\n")
    parsed = build_panel.parse_sgdp_metadata(csv.splitlines())
    assert parsed["LP6005441-DNA_A01"] == ("YORUBA", "AFRICA")
    assert parsed["LP6005441-DNA_B01"] == ("FRENCH", "WEST_EURASIA")


def test_parse_sgdp_metadata_skips_a_row_with_no_population():
    csv = ("Sample ID (Illumina),Population ID,Region\n"
           "LP1,,Africa\n"
           "LP2,French,West Eurasia\n")
    parsed = build_panel.parse_sgdp_metadata(csv.splitlines())
    assert "LP1" not in parsed
    assert "LP2" in parsed


def test_vcf_sample_names_reads_everything_after_format():
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                        "INFO", "FORMAT", "HG00096", "HG00097"])
    assert build_panel.vcf_sample_names(header) == ["HG00096", "HG00097"]


def test_vcf_sample_names_returns_nothing_for_a_sites_only_header():
    header = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                        "INFO"])
    assert build_panel.vcf_sample_names(header) == []


def _vcf_line(rsid="rs1", ref="A", alt="G", calls=("0|0", "0|1")):
    return "\t".join(["22", "16050115", rsid, ref, alt, ".", "PASS", ".",
                      "GT", *calls])


def test_parse_vcf_record_reads_a_biallelic_snv():
    record = build_panel.parse_vcf_record(_vcf_line())
    assert record["chrom"] == "22"
    assert record["pos"] == 16050115
    assert record["rsid"] == "rs1"
    assert record["ref"] == "A"
    assert record["alt"] == "G"


def test_parse_vcf_record_skips_headers():
    assert build_panel.parse_vcf_record("##fileformat=VCFv4.1") is None
    assert build_panel.parse_vcf_record("#CHROM\tPOS") is None


def test_parse_vcf_record_skips_multiallelic_and_indels():
    assert build_panel.parse_vcf_record(_vcf_line(alt="G,T")) is None
    assert build_panel.parse_vcf_record(_vcf_line(ref="A", alt="AT")) is None


def test_parse_vcf_record_blanks_a_missing_identifier():
    assert build_panel.parse_vcf_record(_vcf_line(rsid="."))["rsid"] == ""


def test_allele_counts_counts_diploid_calls_per_population():
    samples = ["s1", "s2", "s3"]
    pops = {"s1": ("CEU", "EUR"), "s2": ("CEU", "EUR"), "s3": ("YRI", "AFR")}
    counts = build_panel.allele_counts(["0|0", "1|1", "0|1"], samples, pops)
    assert counts["CEU"] == (2, 4)
    assert counts["YRI"] == (1, 2)


def test_allele_counts_ignores_missing_calls():
    samples = ["s1", "s2"]
    pops = {"s1": ("CEU", "EUR"), "s2": ("CEU", "EUR")}
    assert build_panel.allele_counts(["./.", "0|1"], samples, pops)["CEU"] == (1, 2)


def test_allele_counts_handles_haploid_and_unphased_calls():
    samples = ["s1", "s2"]
    pops = {"s1": ("CEU", "EUR"), "s2": ("CEU", "EUR")}
    assert build_panel.allele_counts(["1", "0/1"], samples, pops)["CEU"] == (2, 3)


def test_allele_counts_reads_only_the_gt_subfield():
    samples = ["s1"]
    pops = {"s1": ("CEU", "EUR")}
    assert build_panel.allele_counts(["0|1:0.98:1,0,0"], samples,
                                     pops)["CEU"] == (1, 2)


def test_allele_counts_ignores_unlabelled_samples():
    samples = ["s1", "stranger"]
    pops = {"s1": ("CEU", "EUR")}
    counts = build_panel.allele_counts(["1|1", "1|1"], samples, pops)
    assert counts == {"CEU": (2, 2)}


# ===========================================================================
# 9. Output shapes, checked against the readers in backend/ancestry.py
# ===========================================================================

def test_panel_outputs_cover_what_external_declares():
    declared = set(external.PANELS["onekg_sgdp"]["files"])
    assert declared <= set(build_panel.PANEL_OUTPUTS)


def test_panel_outputs_include_the_v3_additions():
    for name in ("informative_markers.tsv", "q_columns.tsv", "BUILD.txt"):
        assert name in build_panel.PANEL_OUTPUTS


def test_populations_tsv_round_trips_through_the_ancestry_reader(tmp_path):
    target = tmp_path / "populations.tsv"
    build_panel.write_populations_tsv(
        target,
        {"HG1": ("CEU", "EUR"), "HG2": ("YRI", "AFR")},
        ["HG1", "HG2"])
    parsed = ancestry.parse_population_map(target)
    assert parsed["order"] == ["CEU", "YRI"]
    assert parsed["samples"]["HG2"] == "YRI"


def test_q_columns_round_trips_through_the_ancestry_reader(tmp_path):
    build_panel.write_q_columns(tmp_path / "q_columns.tsv",
                                ["CEU", "YRI", "CHB"],
                                {"CEU": "EUR", "YRI": "AFR", "CHB": "EAS"})
    labels, labelled = ancestry._column_labels(tmp_path, 3)
    assert labelled is True
    assert labels == ["CEU", "YRI", "CHB"]


def test_informative_markers_round_trips_through_the_ancestry_reader(tmp_path):
    target = tmp_path / "informative_markers.tsv"
    build_panel.write_informative_markers(
        target, {"CEU": [(0.4, "rs1"), (0.3, "rs2")], "YRI": [(0.5, "rs3")]},
        "fst")
    parsed = ancestry._read_informative_markers(target)
    assert parsed["CEU"] == ["rs1", "rs2"]
    assert parsed["YRI"] == ["rs3"]


def test_build_txt_records_provenance_licences_and_hashes(tmp_path):
    target = tmp_path / "BUILD.txt"
    build_panel.write_build_txt(target, {
        "panel_id": "onekg_sgdp", "build": "GRCh37",
        "builder_version": "3.0.0", "retrieved_at": "2026-08-04T00:00:00Z",
        "chromosomes": ["22"], "sample_count": 2504, "population_count": 26,
        "marker_count": 12345, "statistic": "fst", "min_maf": 0.01,
        "array_filter": "", "source_versions": {"1000_genomes": "phase 3"},
        "licences": {"onekg": build_panel.LICENCES["onekg"]},
        "refusals": {"aadr": build_panel.REFUSAL_AADR},
        "hashes": {"panel.map": "abc123"},
    })
    text = target.read_text(encoding="utf-8")
    assert "genome_build:    GRCh37" in text
    assert "2026-08-04T00:00:00Z" in text
    assert "LICENCES, VERBATIM" in text
    assert "abc123  panel.map" in text
    assert "Allen Ancient DNA Resource" in text


def test_panel_dir_is_outside_the_repository():
    target = build_panel.panel_dir()
    assert ROOT not in target.parents
    assert target.parent == Path(external.panel_root())


def test_onekg_source_default_is_the_distribution_that_has_rsids():
    args = build_panel.build_parser().parse_args([])
    assert args.onekg_source == "beagle"
    assert "1kg.phase3.v5a" in build_panel._onekg_vcf_url("22", "beagle")


def test_igsr_source_is_still_reachable_and_uses_the_v1c_name_for_x():
    assert "v5b" in build_panel._onekg_vcf_url("22", "igsr")
    assert "v1c" in build_panel._onekg_vcf_url("X", "igsr")


def test_sgdp_column_aligns_one_sample_to_the_panel_marker_index(tmp_path):
    import gzip
    source = tmp_path / "SGDP_A.vcf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.1\n")
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT", "SGDP_A"]) + "\n")
        for pos, gt in ((100, "1|1"), (200, "0|1"), (300, "./.")):
            fh.write("\t".join(["22", str(pos), f"rs{pos}", "A", "G", ".",
                                "PASS", ".", "GT", gt]) + "\n")
    index = {("22", 100, "A", "G"): 0, ("22", 200, "A", "G"): 1,
             ("22", 300, "A", "G"): 2, ("22", 400, "A", "G"): 3}
    column = build_panel._sgdp_column(str(source), index, 4, local=True)
    # 2 copies, 1 copy, no call, and a marker the sample's VCF never mentions.
    assert list(column) == [2, 1, 255, 255]


def test_sgdp_column_will_not_merge_a_different_alternate_allele(tmp_path):
    import gzip
    source = tmp_path / "SGDP_B.vcf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT", "SGDP_B"]) + "\n")
        fh.write("\t".join(["22", "100", "rs100", "A", "T", ".", "PASS", ".",
                            "GT", "1|1"]) + "\n")
    # The panel has A>G at this position; the sample carries A>T. Matching on
    # position alone would record a homozygous alternate call for the wrong
    # allele, which inverts the genotype.
    column = build_panel._sgdp_column(str(source), {("22", 100, "A", "G"): 0},
                                      1, local=True)
    assert list(column) == [255]


def test_merge_columns_appends_samples_without_disturbing_the_records(tmp_path):
    import gzip
    panel = tmp_path / "panel.vcf.gz"
    with gzip.open(panel, "wt", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.1\n")
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT", "HG1"]) + "\n")
        for pos in (100, 200):
            fh.write("\t".join(["22", str(pos), f"rs{pos}", "A", "G", ".",
                                "PASS", ".", "GT", "0|0"]) + "\n")
    column = tmp_path / "SGDP_A.gt"
    column.write_bytes(bytes([2, 255]))
    merged = tmp_path / "merged.vcf.gz"
    assert build_panel._merge_columns(panel, merged, [("SGDP_A", column)]) == 2
    with gzip.open(merged, "rt", encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if not line.startswith("##")]
    assert lines[0].endswith("\tHG1\tSGDP_A")
    assert lines[1].endswith("\t0|0\t1|1")
    assert lines[2].endswith("\t0|0\t.|.")


def test_merge_columns_is_a_no_op_with_no_extra_samples(tmp_path):
    assert build_panel._merge_columns(tmp_path / "a.gz",
                                      tmp_path / "b.gz", []) == 0


def test_sgdp_sources_prefers_a_local_directory(tmp_path):
    (tmp_path / "SGDP_A.variants.vcf.gz").write_bytes(b"")
    args = build_panel.build_parser().parse_args(["--sgdp-dir", str(tmp_path)])
    found = build_panel._sgdp_sources(args, {"SGDP_A": ("YORUBA", "AFRICA")})
    assert len(found) == 1
    assert found[0][0] == "SGDP_A" and found[0][2] is True


def test_sgdp_sources_falls_back_to_the_url_template():
    args = build_panel.build_parser().parse_args(
        ["--sgdp-vcf-template", "https://example.invalid/{sample}.vcf.gz"])
    found = build_panel._sgdp_sources(args, {"SGDP_A": ("YORUBA", "AFRICA")})
    assert found == [("SGDP_A", "https://example.invalid/SGDP_A.vcf.gz", False)]


def test_sgdp_metadata_override_appears_in_the_plan():
    args = build_panel.build_parser().parse_args(
        ["--sgdp-metadata", "/local/sgdp.csv"])
    plan = build_panel.download_plan(args)
    sgdp = [i for i in plan["items"] if "SGDP" in i["source"]][0]
    assert "/local/sgdp.csv" in sgdp["urls"]


def test_sgdp_plan_states_that_its_url_shape_was_unverified():
    args = build_panel.build_parser().parse_args([])
    plan = build_panel.download_plan(args)
    sgdp = [i for i in plan["items"] if "SGDP" in i["source"]][0]
    assert "could NOT be verified" in sgdp["caveat"]


def test_skip_sgdp_removes_it_from_the_plan_entirely():
    args = build_panel.build_parser().parse_args(["--skip-sgdp"])
    plan = build_panel.download_plan(args)
    assert not any("SGDP" in item["source"] for item in plan["items"])


def test_igsr_source_carries_the_empty_id_column_caveat():
    args = build_panel.build_parser().parse_args(["--onekg-source", "igsr"])
    plan = build_panel.download_plan(args)
    onekg = [i for i in plan["items"] if "1000 Genomes" in i["source"]][0]
    assert "ID column" in onekg["caveat"]


# ===========================================================================
# 10. CPIC allele handling and reconciliation
# ===========================================================================

def test_is_concrete_base_accepts_only_single_unambiguous_bases():
    for base in ("A", "C", "G", "T", "a"):
        assert build_pgx_alleles.is_concrete_base(base)
    for junk in ("R", "Y", "S", "N", "", "AT", None, "-"):
        assert not build_pgx_alleles.is_concrete_base(junk)


def test_match_cpic_allele_matches_exactly():
    assert build_pgx_alleles.match_cpic_allele(["*1", "*2", "*17"], "*2") == "*2"


def test_match_cpic_allele_matches_a_parenthetical_common_name():
    names = ["Reference", "c.1129-5923C>G, c.1236G>A (HapB3)"]
    assert build_pgx_alleles.match_cpic_allele(names, "HapB3") == names[1]


def test_match_cpic_allele_does_not_substring_match_star_numbers():
    # "*1" must not match "*10", "*11" or "*100".
    assert build_pgx_alleles.match_cpic_allele(["*10", "*11", "*100"], "*1") is None


def test_match_cpic_allele_returns_none_for_a_retired_allele():
    assert build_pgx_alleles.match_cpic_allele(["*1", "*5"], "*1B") is None


def _fake_cpic(definitions, alleles=None, structural=()):
    return {"definitions": definitions, "alleles": alleles or {},
            "structural": set(structural)}


def test_reconcile_confirms_a_matching_allele():
    local = [{"gene": "CYP2C19", "allele": "*6",
              "variants": {"rs72552267": "A"}, "note": ""}]
    cpic = {"CYP2C19": _fake_cpic({"*6": {"rs72552267": "A"}})}
    row = build_pgx_alleles.reconcile(local, cpic)[0]
    assert row["status"] == "confirmed"


def test_reconcile_flags_a_conflicting_plus_strand_base():
    local = [{"gene": "TPMT", "allele": "*2",
              "variants": {"rs1800462": "C"}, "note": ""}]
    cpic = {"TPMT": _fake_cpic({"*2": {"rs1800462": "G"}})}
    row = build_pgx_alleles.reconcile(local, cpic)[0]
    assert row["status"] == "conflict"
    assert row["mismatches"] == [{"rsid": "rs1800462", "local": "C",
                                  "cpic": "G"}]


def test_reconcile_flags_an_allele_cpic_no_longer_defines():
    local = [{"gene": "SLCO1B1", "allele": "*1B",
              "variants": {"rs2306283": "G"}, "note": ""}]
    cpic = {"SLCO1B1": _fake_cpic({"*1": {}, "*5": {"rs4149056": "C"}})}
    row = build_pgx_alleles.reconcile(local, cpic)[0]
    assert row["status"] == "allele_not_in_cpic"


def test_reconcile_flags_a_missing_rsid_inside_a_known_allele():
    local = [{"gene": "DPYD", "allele": "HapB3",
              "variants": {"rs56038477": "T"}, "note": ""}]
    cpic = {"DPYD": _fake_cpic({"HapB3": {"rs67376798": "A"}})}
    row = build_pgx_alleles.reconcile(local, cpic)[0]
    assert row["status"] == "rsid_not_in_cpic"


def test_reconcile_refuses_to_confirm_an_iupac_ambiguity_code():
    local = [{"gene": "CYP2D6", "allele": "*2",
              "variants": {"rs1135840": "C"}, "note": ""}]
    cpic = {"CYP2D6": _fake_cpic({"*2": {"rs1135840": "S"}})}
    row = build_pgx_alleles.reconcile(local, cpic)[0]
    assert row["status"] == "ambiguous"


def test_reconcile_reports_a_gene_that_was_not_fetched():
    local = [{"gene": "NUDT15", "allele": "*5", "variants": {"rs1": "A"},
              "note": ""}]
    row = build_pgx_alleles.reconcile(local, {})[0]
    assert row["status"] == "gene_not_fetched"


def test_reconcile_covers_every_currently_unverified_allele():
    entries = unverified_entries()
    assert len(entries) == 17, (
        "backend/diplotype.py changed its unverified list; the reconciliation "
        "report and this test both need re-reading rather than re-numbering."
    )
    report = build_pgx_alleles.reconcile(entries, {})
    assert len(report) == len(entries)
    assert {(r["gene"], r["allele"]) for r in report} == \
        {(e["gene"], e["allele"]) for e in entries}


def test_reconciliation_markdown_names_every_status_it_found():
    local = [
        {"gene": "CYP2C19", "allele": "*6", "variants": {"rs1": "A"}, "note": ""},
        {"gene": "TPMT", "allele": "*2", "variants": {"rs2": "C"}, "note": ""},
    ]
    cpic = {"CYP2C19": _fake_cpic({"*6": {"rs1": "A"}}),
            "TPMT": _fake_cpic({"*2": {"rs2": "G"}})}
    text = build_pgx_alleles.reconciliation_markdown(
        build_pgx_alleles.reconcile(local, cpic), when="2026-08-04T00:00:00Z")
    assert "## confirmed (1)" in text
    assert "## conflict (1)" in text
    assert "was modified" in text and "Nothing in" in text


def test_build_records_marks_an_ambiguous_allele_unverified():
    cpic = {"CYP2D6": _fake_cpic(
        {"*2": {"rs1135840": "S"}, "*4": {"rs3892097": "T"}},
        {"*4": {"function": "no function", "activity": 0.0, "strength": ""}})}
    records = {r["allele"]: r
               for r in build_pgx_alleles.build_records(
                   cpic, when="2026-08-04T00:00:00Z")}
    assert records["*2"]["verified"] is False
    assert records["*2"]["ambiguous_positions"] == ["rs1135840"]
    assert records["*4"]["verified"] is True
    assert records["*4"]["variants"] == {"rs3892097": "T"}


def test_build_records_carries_the_verification_source_and_date():
    cpic = {"TPMT": _fake_cpic({"*3C": {"rs1142345": "C"}})}
    record = build_pgx_alleles.build_records(
        cpic, when="2026-08-04T12:00:00Z")[0]
    assert record["verified_date"] == "2026-08-04"
    assert "cpicpgx.org" in record["verified_source"]


def test_write_output_embeds_the_reconciliation_and_the_licence(tmp_path):
    target = tmp_path / "pgx_alleles.json"
    build_pgx_alleles.write_output(
        target,
        build_pgx_alleles.build_records(
            {"TPMT": _fake_cpic({"*1": {"rs1": "A"}})},
            when="2026-08-04T00:00:00Z"),
        build_pgx_alleles.reconcile(unverified_entries(), {}),
        when="2026-08-04T00:00:00Z",
        genes=["TPMT"])
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload["_meta"]["spdx"] == "CC0-1.0"
    assert payload["_meta"]["excluded_columns"] == list(
        build_pgx_alleles.FORBIDDEN_COLUMNS)
    assert len(payload["reconciliation"]["entries"]) == 17


def test_target_genes_defaults_to_the_genes_diplotype_calls():
    from backend.diplotype import GENES
    assert set(build_pgx_alleles.target_genes()) == set(GENES)


def test_target_genes_accepts_an_explicit_list():
    assert build_pgx_alleles.target_genes("cyp2c19, tpmt") == ["CYP2C19", "TPMT"]


# ===========================================================================
# 11. The manifest, which is what a licence auditor reads first
# ===========================================================================

@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def test_manifest_parses_as_json(manifest):
    assert isinstance(manifest, dict)


def test_manifest_names_its_source_of_truth(manifest):
    assert manifest["generated_from"] == "backend/external.py"


def test_manifest_carries_a_policy_string(manifest):
    policy = manifest["policy"]
    assert isinstance(policy, str) and len(policy) > 200
    assert "MIT" in policy
    assert "CC0" in policy


def test_manifest_tool_ids_match_external_exactly(manifest):
    assert {t["id"] for t in manifest["tools"]} == set(external.TOOLS)


def test_manifest_blocked_ids_match_external_exactly(manifest):
    assert {b["id"] for b in manifest["blocked"]} == set(external.BLOCKED)


def test_manifest_panel_ids_match_external_exactly(manifest):
    assert {p["id"] for p in manifest["panels"]} == set(external.PANELS)


def test_manifest_has_no_extra_or_missing_entries(manifest):
    assert manifest["counts"] == {
        "tools": len(external.TOOLS),
        "blocked": len(external.BLOCKED),
        "panels": len(external.PANELS),
    }


def test_manifest_tool_fields_agree_with_external(manifest):
    for entry in manifest["tools"]:
        source = external.TOOLS[entry["id"]]
        assert entry["name"] == source["name"]
        assert entry["purpose"] == source["purpose"]
        assert entry["capability"] == source["capability"]
        assert entry["licence"] == source["licence"]
        assert entry["spdx"] == source.get("spdx")
        assert entry["homepage"] == source["homepage"]
        assert entry["composable_with_mit"] == bool(source.get("composable", True))
        assert entry["commercial_ok"] == bool(source.get("commercial_ok"))
        assert entry["redistributable"] == bool(source.get("redistributable"))
        assert entry["verified"] == source["verified"]


def test_manifest_blocked_entries_agree_with_external(manifest):
    for entry in manifest["blocked"]:
        source = external.BLOCKED[entry["id"]]
        assert entry["name"] == source["name"]
        assert entry["licence"] == source["licence"]
        assert entry["spdx"] == source.get("spdx")
        assert entry["reason"] == source["reason"]
        assert entry["replacement"] == source.get("replacement")
        assert entry["verified"] == source["verified"]


def test_manifest_blocked_entries_all_state_a_reason_and_a_replacement(manifest):
    for entry in manifest["blocked"]:
        assert entry["reason"].strip()
        assert entry["replacement"]


def test_manifest_blocked_entries_are_never_commercial_or_redistributable(manifest):
    for entry in manifest["blocked"]:
        assert entry["commercial_ok"] is False
        assert entry["redistributable"] is False


def test_manifest_panel_entries_agree_with_external(manifest):
    for entry in manifest["panels"]:
        source = external.PANELS[entry["id"]]
        assert entry["name"] == source["name"]
        assert entry["purpose"] == source["purpose"]
        assert entry["licence"] == source["licence"]
        assert entry["spdx"] == source.get("spdx")
        assert entry["commercial_ok"] == bool(source.get("commercial_ok"))
        assert entry["verified"] == source["verified"]
        assert entry["files"] == list(source.get("files", []))
        assert entry["excluded"] == dict(source.get("excluded", {}))
        assert entry["ethics_gate"] == bool(source.get("ethics_gate", False))


def test_manifest_records_the_hgdp_ethics_gate(manifest):
    hgdp = [p for p in manifest["panels"] if p["id"] == "hgdp_optional"][0]
    assert hgdp["ethics_gate"] is True


def test_manifest_records_the_onekg_panel_exclusions(manifest):
    panel = [p for p in manifest["panels"] if p["id"] == "onekg_sgdp"][0]
    assert set(panel["excluded"]) == {"sgdp_restricted", "hgdp", "aadr"}


def test_manifest_names_the_builder_for_every_panel(manifest):
    for entry in manifest["panels"]:
        assert entry["builder"] == "data/build_panel.py"


def test_manifest_every_entry_has_the_audit_fields(manifest):
    required = ("id", "name", "purpose", "capability", "licence", "spdx",
                "homepage", "composable_with_mit", "commercial_ok",
                "redistributable", "verified")
    for group in ("tools", "blocked", "panels"):
        for entry in manifest[group]:
            for field in required:
                assert field in entry, f"{group}/{entry['id']} missing {field}"


def test_manifest_roots_are_outside_the_repository(manifest):
    assert manifest["roots"]["tools"].startswith("~/.dnainsight")
    assert manifest["roots"]["panels"].startswith("~/.dnainsight")


def test_manifest_gpl_tools_are_marked_not_composable(manifest):
    for entry in manifest["tools"]:
        if str(entry["spdx"] or "").startswith("GPL-3.0"):
            assert entry["composable_with_mit"] is False


# ===========================================================================
# 12. House style
# ===========================================================================

# Built from its code point rather than written literally, because this file is
# one of the files being checked and a literal here would fail its own test.
EM_DASH = chr(0x2014)


@pytest.mark.parametrize("path", NEW_FILES, ids=lambda p: p.name)
def test_no_em_dashes_in_the_new_files(path):
    assert EM_DASH not in path.read_text(encoding="utf-8")


@pytest.mark.parametrize("path", NEW_FILES[:2], ids=lambda p: p.name)
def test_builders_do_not_import_pandas_or_numpy(path):
    text = path.read_text(encoding="utf-8")
    assert "import pandas" not in text
    assert "import numpy" not in text
