"""Tests for backend.ancestry, the global and local ancestry layer.

No network, no external tools, no reference panel files. Every degraded path is
constructed rather than relying on fastmixture happening to be absent, so these
behave the same on a developer machine that has the tools installed.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

from backend import ancestry, external
from backend.ancestry import (
    CAVEAT_ARRAY_DENSITY,
    CAVEAT_MODEL_DEPENDENT,
    CAVEAT_PANEL_BIAS,
    CI_METHOD_LABELS,
    MANDATORY_CAVEATS,
    MIN_RESOLVABLE_FRACTION,
    MIN_RESOLVABLE_MARKERS,
    SUPERPOP_COLOURS,
    AncestryError,
    ancestry_caveats,
    chromosome_painting,
    global_ancestry,
    haplogroup_note,
    inverse_normal_cdf,
    is_phased_vcf,
    local_ancestry,
    marker_coverage,
    normal_cdf,
    panel_manifest,
    panel_unavailable,
    parse_fam,
    parse_population_map,
    parse_q_file,
    percentile_interval,
    resolvable,
    segments_from_calls,
    wilson_interval,
    write_plink,
    z_for_level,
)


# ---------------------------------------------------------------------------
# Isolation
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """Point ~/.dnainsight at a fresh temp directory for every test."""
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    external.reset_cache()
    yield
    external.reset_cache()


def build_panel(files: dict | None = None) -> Path:
    """Create a stub onekg_sgdp panel so the panel gate opens."""
    base = Path(external.panel_root()) / "onekg_sgdp"
    base.mkdir(parents=True, exist_ok=True)
    contents = files or {}
    for name in external.PANELS["onekg_sgdp"]["files"]:
        (base / name).write_text(contents.get(name, "stub\n"), encoding="utf-8")
    return base


Q_TEXT = "0.6000 0.3000 0.1000\n0.1000 0.1000 0.8000\n"
FAM_TEXT = "DNAINSIGHT SAMPLE 0 0 0 -9\n"
POP_TEXT = (
    "sample\tpopulation\tsuperpop\n"
    "HG001\tCEU\tEUR\n"
    "HG002\tCEU\tEUR\n"
    "HG003\tYRI\tAFR\n"
    "HG004\tCHB\tEAS\n"
)


# ---------------------------------------------------------------------------
# .Q parsing
# ---------------------------------------------------------------------------

def test_parse_q_file_reads_a_matrix_of_floats():
    rows = parse_q_file(Q_TEXT)
    assert rows == [[0.6, 0.3, 0.1], [0.1, 0.1, 0.8]]


def test_parse_q_file_reads_from_a_path(tmp_path):
    target = tmp_path / "ancestry.Q"
    target.write_text(Q_TEXT, encoding="utf-8")
    assert parse_q_file(target) == [[0.6, 0.3, 0.1], [0.1, 0.1, 0.8]]


def test_parse_q_file_skips_blank_and_comment_lines():
    rows = parse_q_file("# header\n\n0.5 0.5\n\n")
    assert rows == [[0.5, 0.5]]


def test_parse_q_file_rejects_a_ragged_matrix():
    with pytest.raises(AncestryError):
        parse_q_file("0.5 0.5\n0.2 0.3 0.5\n")


def test_parse_q_file_rejects_a_non_numeric_field():
    with pytest.raises(AncestryError):
        parse_q_file("0.5 EUR\n")


def test_parse_q_file_rows_sum_to_one_within_tolerance():
    for row in parse_q_file(Q_TEXT):
        assert sum(row) == pytest.approx(1.0, abs=1e-6)


def test_parse_q_file_on_empty_input_returns_no_rows():
    assert parse_q_file("") == []


# ---------------------------------------------------------------------------
# .fam and population map parsing
# ---------------------------------------------------------------------------

def test_parse_fam_reads_the_six_plink_columns():
    rows = parse_fam(FAM_TEXT)
    assert rows[0]["fid"] == "DNAINSIGHT"
    assert rows[0]["iid"] == "SAMPLE"


def test_parse_fam_pads_a_short_row_rather_than_failing():
    rows = parse_fam("FAM SAMPLE\n")
    assert rows[0]["iid"] == "SAMPLE"
    assert rows[0]["phenotype"] == "0"


def test_parse_population_map_counts_samples_per_population():
    parsed = parse_population_map(POP_TEXT)
    counts = {p["code"]: p["samples"] for p in parsed["populations"]}
    assert counts == {"CEU": 2, "YRI": 1, "CHB": 1}


def test_parse_population_map_preserves_first_seen_order():
    # The column order of a .Q matrix is fixed by the panel, so reordering
    # here would silently relabel somebody's ancestry.
    assert parse_population_map(POP_TEXT)["order"] == ["CEU", "YRI", "CHB"]


def test_parse_population_map_labels_populations_from_the_frequency_table():
    populations = {p["code"]: p for p in parse_population_map(POP_TEXT)["populations"]}
    assert "Yoruba" in populations["YRI"]["label"]
    assert populations["YRI"]["superpop"] == "AFR"


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def test_normal_cdf_matches_known_values():
    assert normal_cdf(0.0) == pytest.approx(0.5)
    assert normal_cdf(1.96) == pytest.approx(0.975, abs=1e-3)


def test_inverse_normal_cdf_inverts_the_cdf():
    for p in (0.025, 0.25, 0.5, 0.75, 0.975):
        assert normal_cdf(inverse_normal_cdf(p)) == pytest.approx(p, abs=1e-6)


def test_inverse_normal_cdf_rejects_a_probability_outside_the_open_interval():
    with pytest.raises(AncestryError):
        inverse_normal_cdf(0.0)
    with pytest.raises(AncestryError):
        inverse_normal_cdf(1.0)


def test_z_for_level_gives_the_familiar_multiplier():
    assert z_for_level(0.95) == pytest.approx(1.96, abs=1e-3)


def test_wilson_interval_brackets_the_estimate():
    low, high = wilson_interval(0.4, 500)
    assert low < 0.4 < high


def test_wilson_interval_never_leaves_the_zero_to_one_range():
    low, high = wilson_interval(0.0, 30)
    assert low >= 0.0 and high <= 1.0
    low, high = wilson_interval(1.0, 30)
    assert low >= 0.0 and high <= 1.0


def test_wilson_interval_with_no_markers_is_completely_unconstrained():
    # No informative markers means the estimate rests on nothing, and a tight
    # interval around nothing would be the dishonest answer.
    assert wilson_interval(0.5, 0) == (0.0, 1.0)


def test_wilson_interval_narrows_as_markers_increase():
    narrow = wilson_interval(0.4, 10000)
    wide = wilson_interval(0.4, 100)
    assert (narrow[1] - narrow[0]) < (wide[1] - wide[0])


def test_percentile_interval_brackets_the_replicates():
    low, high = percentile_interval([0.1 * i for i in range(11)])
    assert low <= 0.1 and high >= 0.9


def test_percentile_interval_refuses_to_invent_an_interval_from_one_draw():
    assert percentile_interval([0.5]) is None
    assert percentile_interval([]) is None


# ---------------------------------------------------------------------------
# Marker coverage and resolvability
# ---------------------------------------------------------------------------

def markers(count: int, prefix: str) -> list[str]:
    return [f"{prefix}{i}" for i in range(count)]


def test_marker_coverage_counts_what_the_array_actually_reads():
    informative = {"CEU": markers(100, "rs")}
    typed = set(markers(80, "rs"))
    coverage = marker_coverage(typed, informative)
    assert coverage["CEU"]["markers_read"] == 80
    assert coverage["CEU"]["informative_markers"] == 100
    assert coverage["CEU"]["coverage"] == pytest.approx(0.8)


def test_marker_coverage_accepts_a_genotype_map_as_input():
    informative = {"CEU": ["rs1", "rs2", "rs3", "rs4"]}
    coverage = marker_coverage({"rs1": "AG", "rs3": "CC"}, informative)
    assert coverage["CEU"]["markers_read"] == 2


def test_marker_coverage_is_case_insensitive_on_rsids():
    coverage = marker_coverage(["RS1", "Rs2"], {"CEU": ["rs1", "rs2"]})
    assert coverage["CEU"]["markers_read"] == 2


def test_a_well_covered_population_is_resolvable():
    coverage = marker_coverage(markers(900, "rs"), {"CEU": markers(1000, "rs")})
    assert coverage["CEU"]["resolvable"] is True
    assert coverage["CEU"]["state"] == "resolvable"


def test_a_low_coverage_population_is_not_resolvable_rather_than_zero():
    # The heart of it: zero percent is a measurement, not resolvable is an
    # admission. They are different claims and must not render the same.
    coverage = marker_coverage(markers(11, "rs"), {"CEU": markers(900, "rs")})
    entry = coverage["CEU"]
    assert entry["resolvable"] is False
    assert entry["state"] == "not_resolvable"
    assert entry["coverage"] != 0.0
    assert "informative markers" in entry["reason"]


def test_an_absolute_marker_floor_applies_even_at_a_high_fraction():
    # 25 markers out of 30 is 83 percent and still cannot separate anything.
    coverage = marker_coverage(markers(25, "rs"), {"CEU": markers(30, "rs")})
    assert coverage["CEU"]["resolvable"] is False
    assert str(MIN_RESOLVABLE_MARKERS) in coverage["CEU"]["reason"]


def test_a_fractional_floor_applies_even_at_a_large_absolute_count():
    # 200 markers sounds impressive and is a rounding error against 90,000.
    coverage = marker_coverage(markers(200, "rs"), {"CEU": markers(90000, "rs")})
    assert coverage["CEU"]["resolvable"] is False


def test_a_population_with_no_informative_markers_is_not_resolvable():
    coverage = marker_coverage(["rs1"], {"XXX": []})
    assert coverage["XXX"]["resolvable"] is False
    assert "no ancestry-informative markers" in coverage["XXX"]["reason"]


def test_resolvable_requires_both_thresholds():
    assert resolvable(MIN_RESOLVABLE_MARKERS, 100) is True
    assert resolvable(MIN_RESOLVABLE_MARKERS - 1, 100) is False
    assert resolvable(1000, int(1000 / MIN_RESOLVABLE_FRACTION) + 1000) is False


# ---------------------------------------------------------------------------
# Panel manifest
# ---------------------------------------------------------------------------

def test_panel_manifest_reports_unbuilt_state_without_raising():
    manifest = panel_manifest("onekg_sgdp")
    assert manifest["available"] is False
    assert manifest["state"] == "not_built"
    assert manifest["populations"] == []
    assert manifest["content_hash"] is None


def test_panel_manifest_uses_none_not_zero_for_unbuilt_counts():
    # Zero populations and "not built" are different states.
    manifest = panel_manifest("onekg_sgdp")
    assert manifest["population_count"] is None
    assert manifest["sample_count"] is None
    assert manifest["marker_count"] is None


def test_panel_manifest_always_carries_the_licence_and_verification_date():
    manifest = panel_manifest("onekg_sgdp")
    assert "1000 Genomes" in manifest["licence"]
    assert manifest["licence_verified"]
    assert manifest["commercial_ok"] is True


def test_panel_manifest_lists_the_files_it_expects_and_their_state():
    manifest = panel_manifest("onekg_sgdp")
    names = {f["name"] for f in manifest["files"]}
    assert names == set(external.PANELS["onekg_sgdp"]["files"])
    assert all(f["present"] is False for f in manifest["files"])


def test_panel_manifest_hashes_the_files_once_they_exist():
    build_panel({"populations.tsv": POP_TEXT, "panel.map": "1\nrs2\nrs3\n"})
    manifest = panel_manifest("onekg_sgdp")
    assert manifest["available"] is True
    assert manifest["content_hash"] is not None
    assert all(f["sha256"] for f in manifest["files"])


def test_panel_manifest_reports_populations_and_sample_counts_when_built():
    build_panel({"populations.tsv": POP_TEXT})
    manifest = panel_manifest("onekg_sgdp")
    assert manifest["population_count"] == 3
    assert manifest["sample_count"] == 4


def test_panel_manifest_counts_markers_from_the_map_file():
    build_panel({"panel.map": "rs1\nrs2\nrs3\nrs4\n"})
    assert panel_manifest("onekg_sgdp")["marker_count"] == 4


def test_panel_manifest_content_hash_is_stable_across_calls():
    build_panel()
    assert panel_manifest("onekg_sgdp")["content_hash"] == \
        panel_manifest("onekg_sgdp")["content_hash"]


def test_panel_manifest_content_hash_changes_when_the_panel_changes():
    build_panel()
    before = panel_manifest("onekg_sgdp")["content_hash"]
    build_panel({"panel.map": "different\n"})
    assert panel_manifest("onekg_sgdp")["content_hash"] != before


def test_panel_manifest_refuses_to_describe_an_unknown_panel():
    manifest = panel_manifest("not_a_panel")
    assert manifest["available"] is False
    assert manifest["state"] == "unknown"
    assert "Unknown reference panel" in manifest["reason"]


def test_panel_manifest_records_the_excluded_sources_and_why():
    build_panel()
    excluded = panel_manifest("onekg_sgdp")["excluded"]
    assert "hgdp" in excluded and "sgdp_restricted" in excluded


# ---------------------------------------------------------------------------
# PLINK writing
# ---------------------------------------------------------------------------

def test_write_plink_emits_the_bed_magic_bytes(tmp_path):
    calls = [{"rsid": "rs1", "chromosome": "1", "position": 100,
              "allele1": "A", "allele2": "G"}]
    result = write_plink(calls, tmp_path / "study")
    assert Path(result["bed"]).read_bytes()[:3] == bytes([0x6C, 0x1B, 0x01])


def test_write_plink_writes_one_bim_line_per_marker(tmp_path):
    calls = [
        {"rsid": "rs1", "chromosome": "1", "position": 100,
         "allele1": "A", "allele2": "G"},
        {"rsid": "rs2", "chromosome": "1", "position": 200,
         "allele1": "C", "allele2": "C"},
    ]
    result = write_plink(calls, tmp_path / "study")
    lines = Path(result["bim"]).read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    assert result["written"] == 2


def test_write_plink_fam_round_trips_through_the_fam_parser(tmp_path):
    result = write_plink([{"rsid": "rs1", "chromosome": "1", "position": 100,
                           "allele1": "A", "allele2": "G"}],
                         tmp_path / "study", sample="TESTER")
    assert parse_fam(Path(result["fam"]))[0]["iid"] == "TESTER"


def test_write_plink_skips_markers_with_no_coordinates(tmp_path):
    result = write_plink({"rs1": "AG"}, tmp_path / "study")
    assert result["written"] == 0
    assert result["skipped"][0]["reason"] == "no chromosome or position"


# ---------------------------------------------------------------------------
# Degraded payloads
# ---------------------------------------------------------------------------

def test_global_ancestry_without_the_tool_is_tagged_tool_missing():
    result = global_ancestry({"rs1": "AG"})
    assert result["available"] is False
    assert result["problem"] == "tool_missing"
    assert result["proportions"] == []
    assert result["not_attempted"] is True


def test_global_ancestry_with_the_tool_but_no_panel_is_tagged_panel_missing(monkeypatch):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    result = global_ancestry({"rs1": "AG"})
    assert result["available"] is False
    assert result["problem"] == "panel_missing"
    assert "DATA problem" in result["reason"]


def test_the_two_global_ancestry_failures_do_not_share_a_message(monkeypatch):
    tool_missing = global_ancestry({"rs1": "AG"})
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    panel_missing = global_ancestry({"rs1": "AG"})
    assert tool_missing["problem"] != panel_missing["problem"]
    assert tool_missing["reason"] != panel_missing["reason"]


def test_global_ancestry_rejects_an_unknown_mode():
    result = global_ancestry({"rs1": "AG"}, mode="magic")
    assert result["problem"] == "bad_mode"
    assert "projection" in result["reason"]


def test_global_ancestry_degraded_payloads_still_carry_the_caveats(monkeypatch):
    tool_missing = global_ancestry({"rs1": "AG"})
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    panel_missing = global_ancestry({"rs1": "AG"})
    for payload in (tool_missing, panel_missing):
        for caveat in MANDATORY_CAVEATS:
            assert caveat in payload["caveats"]


def test_global_ancestry_reports_no_usable_calls_distinctly(monkeypatch, tmp_path):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    build_panel()
    result = global_ancestry({"rs1": "AG"}, workdir=tmp_path / "work")
    assert result["problem"] == "no_input"


def test_local_ancestry_without_the_tool_is_tagged_tool_missing(tmp_path):
    result = local_ancestry(tmp_path / "nothing.vcf")
    assert result["available"] is False
    assert result["problem"] == "tool_missing"
    assert result["segments"] == []


def test_local_ancestry_with_the_tool_but_no_panel_is_tagged_panel_missing(
        monkeypatch, tmp_path):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    result = local_ancestry(tmp_path / "nothing.vcf")
    assert result["problem"] == "panel_missing"


def test_local_ancestry_refuses_unphased_input(monkeypatch, tmp_path):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    build_panel()
    vcf = tmp_path / "unphased.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    result = local_ancestry(vcf)
    assert result["problem"] == "input_not_phased"
    assert "Phase the file first" in result["reason"]


def test_local_ancestry_reports_a_missing_input_file_distinctly(
        monkeypatch, tmp_path):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    build_panel()
    result = local_ancestry(tmp_path / "absent.vcf")
    assert result["problem"] == "no_input"


def test_is_phased_vcf_detects_the_phased_separator():
    phased = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
              "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n")
    assert is_phased_vcf(phased) is True
    assert is_phased_vcf(phased.replace("0|1", "0/1")) is False


def test_is_phased_vcf_on_an_empty_file_is_false():
    assert is_phased_vcf("") is False


def test_panel_unavailable_distinguishes_partial_from_unbuilt():
    unbuilt = panel_unavailable("onekg_sgdp", "ancestry_global")
    assert unbuilt["panel_state"] == "not_built"
    base = Path(external.panel_root()) / "onekg_sgdp"
    base.mkdir(parents=True, exist_ok=True)
    (base / "panel.map").write_text("x\n", encoding="utf-8")
    partial = panel_unavailable("onekg_sgdp", "ancestry_global")
    assert partial["panel_state"] == "partial"
    assert partial["reason"] != unbuilt["reason"]


# ---------------------------------------------------------------------------
# Proportions, intervals and the not-resolvable state
# ---------------------------------------------------------------------------

def build_records(row, labels, coverage, level=0.95):
    return ancestry._proportion_records(
        row, labels, coverage, level=level, replicates=None, tool_intervals=None)


def test_proportions_carry_an_interval_and_the_method_that_produced_it():
    coverage = marker_coverage(markers(900, "rs"), {"CEU": markers(1000, "rs")})
    reported, withheld = build_records([0.6], ["CEU"], coverage)
    assert withheld == []
    record = reported[0]
    assert record["ci_low"] is not None and record["ci_high"] is not None
    assert record["ci_method"] == "wilson_marker_count"
    assert record["ci_method_label"] == CI_METHOD_LABELS["wilson_marker_count"]


def test_the_approximate_interval_says_so_in_plain_words():
    assert "APPROXIMATE" in CI_METHOD_LABELS["wilson_marker_count"]
    assert "not a bootstrap" in CI_METHOD_LABELS["wilson_marker_count"]


def test_a_bootstrap_interval_is_labelled_as_a_bootstrap():
    coverage = marker_coverage(markers(900, "rs"), {"CEU": markers(1000, "rs")})
    reported, _ = ancestry._proportion_records(
        [0.6], ["CEU"], coverage, level=0.95,
        replicates={"CEU": [0.5, 0.55, 0.6, 0.65, 0.7]}, tool_intervals=None)
    assert reported[0]["ci_method"] == "bootstrap_marker_subsets"


def test_a_tool_supplied_interval_wins_and_is_labelled_as_such():
    coverage = marker_coverage(markers(900, "rs"), {"CEU": markers(1000, "rs")})
    reported, _ = ancestry._proportion_records(
        [0.6], ["CEU"], coverage, level=0.95, replicates=None,
        tool_intervals={"CEU": (0.55, 0.65)})
    assert reported[0]["ci_method"] == "tool_reported"
    assert reported[0]["ci_low"] == pytest.approx(0.55)


def test_the_interval_brackets_the_point_estimate():
    coverage = marker_coverage(markers(900, "rs"), {"CEU": markers(1000, "rs")})
    record = build_records([0.6], ["CEU"], coverage)[0][0]
    assert record["ci_low"] <= record["proportion"] <= record["ci_high"]


def test_an_unresolvable_population_reports_none_not_zero():
    coverage = marker_coverage(markers(11, "rs"), {"CEU": markers(900, "rs")})
    reported, withheld = build_records([0.02], ["CEU"], coverage)
    assert reported == []
    assert withheld[0]["proportion"] is None
    assert withheld[0]["percent"] is None
    assert withheld[0]["state"] == "not_resolvable"
    assert "NOT RESOLVABLE" in withheld[0]["note"]


def test_an_unresolvable_population_keeps_its_raw_estimate_visible():
    # The mass does not vanish. It is withheld and the arithmetic stays
    # auditable rather than quietly disappearing from the total.
    coverage = marker_coverage(markers(11, "rs"), {"CEU": markers(900, "rs")})
    _, withheld = build_records([0.02], ["CEU"], coverage)
    assert withheld[0]["raw_estimate"] == pytest.approx(0.02)


def test_an_unresolvable_population_gets_no_confidence_interval():
    coverage = marker_coverage(markers(11, "rs"), {"CEU": markers(900, "rs")})
    _, withheld = build_records([0.02], ["CEU"], coverage)
    assert withheld[0]["ci_method"] == "none"
    assert withheld[0]["ci_low"] is None


def test_reported_proportions_sum_to_one_within_tolerance():
    informative = {code: markers(1000, f"{code.lower()}_rs")
                   for code in ("CEU", "YRI", "CHB")}
    typed = [m for group in informative.values() for m in group[:900]]
    coverage = marker_coverage(typed, informative)
    reported, withheld = build_records([0.6, 0.3, 0.1], ["CEU", "YRI", "CHB"],
                                       coverage)
    assert withheld == []
    assert sum(r["proportion"] for r in reported) == pytest.approx(1.0, abs=1e-6)


def test_proportions_are_sorted_largest_first():
    reported, _ = build_records([0.1, 0.6, 0.3], ["CEU", "YRI", "CHB"], {})
    assert [r["population"] for r in reported] == ["YRI", "CHB", "CEU"]


def test_proportions_carry_percentages_alongside_fractions():
    record = build_records([0.1234], ["CEU"], {})[0][0]
    assert record["percent"] == pytest.approx(12.34)


# ---------------------------------------------------------------------------
# Segments and chromosome painting
# ---------------------------------------------------------------------------

def painted_calls() -> list[dict]:
    return [
        {"chromosome": "1", "position": 1000, "haplotype1": "CEU", "haplotype2": "YRI"},
        {"chromosome": "1", "position": 2000, "haplotype1": "CEU", "haplotype2": "YRI"},
        {"chromosome": "1", "position": 3000, "haplotype1": "YRI", "haplotype2": "YRI"},
        {"chromosome": "2", "position": 5000, "haplotype1": "CHB", "haplotype2": "CEU"},
    ]


def test_segments_merge_consecutive_markers_of_the_same_ancestry():
    segments = segments_from_calls(painted_calls())
    hap1_chr1 = [s for s in segments if s["haplotype"] == 1 and s["chromosome"] == "1"]
    assert len(hap1_chr1) == 2
    assert hap1_chr1[0]["ancestry"] == "CEU"
    assert hap1_chr1[0]["start"] == 1000 and hap1_chr1[0]["end"] == 2000
    assert hap1_chr1[0]["markers"] == 2


def test_segments_are_built_per_haplotype_not_averaged():
    segments = segments_from_calls(painted_calls())
    assert {s["haplotype"] for s in segments} == {1, 2}
    hap2_chr1 = [s for s in segments if s["haplotype"] == 2 and s["chromosome"] == "1"]
    assert len(hap2_chr1) == 1


def test_segments_do_not_run_across_a_chromosome_boundary():
    calls = [
        {"chromosome": "1", "position": 1000, "haplotype1": "CEU", "haplotype2": "CEU"},
        {"chromosome": "2", "position": 1000, "haplotype1": "CEU", "haplotype2": "CEU"},
    ]
    segments = [s for s in segments_from_calls(calls) if s["haplotype"] == 1]
    assert {s["chromosome"] for s in segments} == {"1", "2"}
    assert len(segments) == 2


def test_segments_carry_a_span():
    segments = segments_from_calls(painted_calls())
    assert all(s["span"] >= 0 for s in segments)


def test_chromosome_painting_produces_the_documented_segment_shape():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    assert result["available"] is True
    segment = result["segments"][0]
    for key in ("chromosome", "start", "end", "ancestry", "label", "colour",
                "haplotype", "x1", "x2", "span"):
        assert key in segment


def test_chromosome_painting_emits_a_colour_key():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    key = {entry["ancestry"]: entry["colour"] for entry in result["colour_key"]}
    assert set(key) == {"CEU", "YRI", "CHB"}
    assert key["CEU"] == SUPERPOP_COLOURS["EUR"]
    assert key["YRI"] == SUPERPOP_COLOURS["AFR"]


def test_chromosome_painting_colours_are_deterministic():
    segments = segments_from_calls(painted_calls())
    first = chromosome_painting({"available": True, "segments": segments})
    second = chromosome_painting({"available": True, "segments": segments})
    assert first["colour_key"] == second["colour_key"]


def test_chromosome_painting_scales_segments_to_the_chromosome_length():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    for segment in result["segments"]:
        assert 0.0 <= segment["x1"] <= 1.0
        assert 0.0 <= segment["x2"] <= 1.0


def test_chromosome_painting_refuses_to_scale_an_unknown_chromosome():
    calls = [{"chromosome": "ZZ", "position": 10, "haplotype1": "CEU",
              "haplotype2": "CEU"}]
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(calls)})
    assert result["segments"][0]["x1"] is None


def test_chromosome_painting_groups_by_chromosome_in_genomic_order():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    assert [c["chromosome"] for c in result["chromosomes"]] == ["1", "2"]


def test_chromosome_painting_totals_span_per_ancestry():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    codes = {t["ancestry"] for t in result["totals"]}
    assert codes == {"CEU", "YRI", "CHB"}


def test_chromosome_painting_passes_a_degraded_result_straight_through():
    degraded = local_ancestry("nowhere.vcf")
    painting = chromosome_painting(degraded)
    assert painting["available"] is False
    assert painting["segments"] == []
    assert painting["problem"] == "tool_missing"


def test_chromosome_painting_handles_junk_input_without_raising():
    painting = chromosome_painting(None)
    assert painting["available"] is False
    assert painting["caveats"]


def test_chromosome_painting_warns_that_boundaries_are_marker_positions():
    result = chromosome_painting({"available": True,
                                  "segments": segments_from_calls(painted_calls())})
    assert any("outermost marker" in c for c in result["caveats"])


# ---------------------------------------------------------------------------
# Caveats
# ---------------------------------------------------------------------------

def test_ancestry_caveats_always_include_every_mandatory_caveat():
    caveats = ancestry_caveats()
    for caveat in MANDATORY_CAVEATS:
        assert caveat in caveats


def test_ancestry_caveats_state_the_panel_under_representation():
    assert CAVEAT_PANEL_BIAS in ancestry_caveats()
    assert "under-represent" in CAVEAT_PANEL_BIAS


def test_ancestry_caveats_state_that_proportions_are_model_dependent():
    assert CAVEAT_MODEL_DEPENDENT in ancestry_caveats()
    assert "not a fact about you" in CAVEAT_MODEL_DEPENDENT


def test_ancestry_caveats_state_the_array_density_limit():
    assert CAVEAT_ARRAY_DENSITY in ancestry_caveats()
    assert "consumer array" in CAVEAT_ARRAY_DENSITY


def test_ancestry_caveats_append_the_panels_own_recorded_note():
    assert external.panel_status("onekg_sgdp")["note"] in ancestry_caveats()


# ---------------------------------------------------------------------------
# Haplogroups belong to another module
# ---------------------------------------------------------------------------

def test_haplogroup_note_degrades_when_the_module_is_absent():
    note = haplogroup_note()
    assert isinstance(note["available"], bool)
    assert note["note"]


def test_haplogroup_note_never_raises_import_error():
    # backend/haplogroups.py is owned by another agent and may land later.
    assert haplogroup_note() is not None
