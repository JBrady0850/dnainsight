"""Tests for backend.imputation, the imputation quality layer.

No network, no external tools, no reference panel. Every degraded path is
exercised by construction rather than by having Beagle absent by luck, so these
pass identically on a developer machine that happens to have Beagle installed.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest

from backend import external, imputation
from backend.imputation import (
    CAP_STEP_PREFIX,
    CAVEAT_PANEL_BIAS,
    CAVEAT_RARE_VARIANTS,
    DEFAULT_DR2_THRESHOLD,
    FILTER_TOKEN_KEYS,
    IMPUTED_MAGNITUDE_CAP,
    IMPUTED_MAGNITUDE_CEILING,
    MANDATORY_CAVEATS,
    QUALITY_BANDS,
    TYPED_MAGNITUDE_CEILING,
    ImputedQualityViolation,
    af_from_info,
    apply_imputation_cap,
    apply_imputation_cap_all,
    assert_no_imputed_pathogenic_without_quality,
    build_caveats,
    coverage_report,
    dr2_from_info,
    filter_flag_pattern,
    filter_op_pattern,
    filter_tokens,
    impute,
    max_imputed_magnitude,
    minor_allele_frequency,
    panel_unavailable,
    parse_info,
    parse_vcf,
    parse_vcf_line,
    quality_band,
    write_study_vcf,
)


# ---------------------------------------------------------------------------
# Isolation
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """Point ~/.dnainsight at a fresh temp directory for every test.

    Without this a developer with Beagle actually installed, or with a consent
    file on disk, would get different results from CI, and the degraded-path
    tests would silently stop testing anything.
    """
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    external.reset_cache()
    yield
    external.reset_cache()


def imputed_vcf_text() -> str:
    """A four-line Beagle-shaped output VCF covering every quality band."""
    return (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "1\t1000\trs1\tA\tG\t.\tPASS\tDR2=0.99;AF=0.40;IMP\tGT\t0|1\n"
        "1\t2000\trs2\tC\tT\t.\tPASS\tDR2=0.85;AF=0.20;IMP\tGT\t0|0\n"
        "1\t3000\trs3\tG\tA\t.\tPASS\tDR2=0.50;AF=0.005;IMP\tGT\t0|1\n"
        "1\t4000\trs4\tT\tC\t.\tPASS\tDR2=0.10;AF=0.002;IMP\tGT\t1|1\n"
        "1\t5000\trs5\tA\tC\t.\tPASS\tAF=0.30\tGT\t0|1\n"
    )


# ---------------------------------------------------------------------------
# INFO parsing
# ---------------------------------------------------------------------------

def test_parse_info_reads_key_value_pairs():
    assert parse_info("DR2=0.9;AF=0.25") == {"DR2": "0.9", "AF": "0.25"}


def test_parse_info_treats_a_bare_flag_as_true():
    fields = parse_info("DR2=0.9;IMP")
    assert fields["IMP"] is True


def test_parse_info_on_an_empty_column_returns_empty_dict():
    assert parse_info(".") == {}
    assert parse_info("") == {}
    assert parse_info(None) == {}


def test_parse_info_tolerates_stray_semicolons():
    assert parse_info(";DR2=0.5;;") == {"DR2": "0.5"}


def test_dr2_from_info_reads_the_dr2_key():
    assert dr2_from_info("DR2=0.87;AF=0.1") == pytest.approx(0.87)


def test_dr2_from_info_falls_back_to_r2_then_info():
    assert dr2_from_info("R2=0.62") == pytest.approx(0.62)
    assert dr2_from_info("INFO=0.44") == pytest.approx(0.44)


def test_dr2_from_info_prefers_dr2_over_the_fallbacks():
    assert dr2_from_info("R2=0.10;DR2=0.90") == pytest.approx(0.90)


def test_dr2_from_info_returns_none_when_no_quality_field_exists():
    assert dr2_from_info("AF=0.3;IMP") is None


def test_dr2_from_info_takes_the_first_value_at_a_multiallelic_site():
    assert dr2_from_info("DR2=0.80,0.20") == pytest.approx(0.80)


def test_dr2_from_info_accepts_an_already_parsed_dict():
    assert dr2_from_info({"DR2": "0.75"}) == pytest.approx(0.75)


def test_af_from_info_reads_the_alternate_allele_frequency():
    assert af_from_info("DR2=0.9;AF=0.03") == pytest.approx(0.03)


def test_minor_allele_frequency_folds_onto_the_minor_allele():
    assert minor_allele_frequency(0.9) == pytest.approx(0.1)
    assert minor_allele_frequency(0.1) == pytest.approx(0.1)
    assert minor_allele_frequency(None) is None


# ---------------------------------------------------------------------------
# VCF line parsing
# ---------------------------------------------------------------------------

def test_parse_vcf_line_extracts_dr2_and_marks_the_call_imputed():
    record = parse_vcf_line(
        "1\t1000\trs1\tA\tG\t.\tPASS\tDR2=0.95;AF=0.4;IMP\tGT\t0|1")
    assert record["rsid"] == "rs1"
    assert record["dr2"] == pytest.approx(0.95)
    assert record["imputed"] is True
    assert record["typed"] is False
    assert record["quality_band"] == "high"


def test_parse_vcf_line_without_the_imp_flag_is_a_typed_call():
    record = parse_vcf_line("1\t5000\trs5\tA\tC\t.\tPASS\tAF=0.3\tGT\t0|1")
    assert record["imputed"] is False
    assert record["typed"] is True
    assert record["quality_band"] == "typed"


def test_parse_vcf_line_returns_none_for_headers_and_blanks():
    assert parse_vcf_line("##fileformat=VCFv4.2") is None
    assert parse_vcf_line("#CHROM\tPOS") is None
    assert parse_vcf_line("") is None


def test_parse_vcf_line_returns_none_for_a_truncated_row():
    assert parse_vcf_line("1\t1000\trs1") is None


def test_parse_vcf_reads_every_data_line_and_no_headers():
    records = parse_vcf(imputed_vcf_text())
    assert len(records) == 5
    assert [r["rsid"] for r in records] == ["rs1", "rs2", "rs3", "rs4", "rs5"]


def test_parse_vcf_computes_minor_allele_frequency_per_record():
    records = {r["rsid"]: r for r in parse_vcf(imputed_vcf_text())}
    assert records["rs3"]["maf"] == pytest.approx(0.005)


# ---------------------------------------------------------------------------
# Quality bands
# ---------------------------------------------------------------------------

def test_quality_band_high_starts_at_the_documented_cut_point():
    assert quality_band(0.90) == "high"
    assert quality_band(1.00) == "high"


def test_quality_band_moderate_spans_its_documented_cut_points():
    assert quality_band(0.899999) == "moderate"
    assert quality_band(0.80) == "moderate"


def test_quality_band_low_spans_its_documented_cut_points():
    assert quality_band(0.799999) == "low"
    assert quality_band(0.30) == "low"


def test_quality_band_unusable_below_its_cut_point():
    assert quality_band(0.299999) == "unusable"
    assert quality_band(0.0) == "unusable"


def test_quality_band_missing_dr2_is_unknown_not_unusable():
    # Absent is not zero. The same distinction frequency.py draws between a
    # measured 0.0 and no data at all.
    assert quality_band(None) == "unknown"
    assert quality_band("") == "unknown"
    assert quality_band("not a number") == "unknown"


def test_quality_band_clamps_a_negative_rounding_artefact_to_unusable():
    assert quality_band(-0.01) == "unusable"


def test_quality_band_accepts_a_numeric_string():
    assert quality_band("0.95") == "high"


def test_every_quality_band_returned_is_a_declared_band():
    for value in (None, -1.0, 0.0, 0.29, 0.3, 0.79, 0.8, 0.89, 0.9, 1.0):
        assert quality_band(value) in QUALITY_BANDS


# ---------------------------------------------------------------------------
# The parity guarantee
# ---------------------------------------------------------------------------

def test_imputed_ceiling_is_strictly_below_the_typed_ceiling():
    assert IMPUTED_MAGNITUDE_CEILING < TYPED_MAGNITUDE_CEILING


def test_max_imputed_magnitude_hard_caps_below_the_threshold():
    assert max_imputed_magnitude(0.5) == IMPUTED_MAGNITUDE_CAP
    assert max_imputed_magnitude(0.79) == IMPUTED_MAGNITUDE_CAP


def test_max_imputed_magnitude_uses_the_ceiling_at_or_above_the_threshold():
    assert max_imputed_magnitude(DEFAULT_DR2_THRESHOLD) == IMPUTED_MAGNITUDE_CEILING
    assert max_imputed_magnitude(1.0) == IMPUTED_MAGNITUDE_CEILING


def test_max_imputed_magnitude_treats_a_missing_dr2_as_below_threshold():
    assert max_imputed_magnitude(None) == IMPUTED_MAGNITUDE_CAP


def test_no_dr2_whatsoever_can_reach_typed_parity():
    # The structural claim of this module: a prediction never ties a
    # measurement, at any quality, under any threshold.
    for step in range(0, 101):
        dr2 = step / 100.0
        assert max_imputed_magnitude(dr2) < TYPED_MAGNITUDE_CEILING


# ---------------------------------------------------------------------------
# The cap
# ---------------------------------------------------------------------------

def test_cap_reduces_a_low_quality_imputed_call_to_the_hard_cap():
    finding = {"rsid": "rs1", "magnitude": 7.2, "imputed": True, "dr2": 0.42}
    apply_imputation_cap(finding)
    assert finding["magnitude"] == IMPUTED_MAGNITUDE_CAP


def test_cap_records_a_named_step_in_the_audit_trail():
    finding = {"rsid": "rs1", "magnitude": 7.2, "imputed": True, "dr2": 0.42}
    apply_imputation_cap(finding)
    steps = [s for s in finding["magnitude_factors"] if s.startswith(CAP_STEP_PREFIX)]
    assert len(steps) == 1
    assert "7.20" in steps[0] and "3.00" in steps[0]


def test_cap_step_names_the_dr2_and_the_band():
    finding = {"rsid": "rs1", "magnitude": 7.2, "imputed": True, "dr2": 0.42}
    apply_imputation_cap(finding)
    step = finding["magnitude_factors"][-1]
    assert "0.420" in step
    assert "low" in step


def test_cap_appends_a_step_even_when_the_ceiling_does_not_bind():
    # A silent no-op would leave the reader unable to tell the rule ran.
    finding = {"rsid": "rs1", "magnitude": 1.5, "imputed": True, "dr2": 0.95}
    apply_imputation_cap(finding)
    assert finding["magnitude"] == 1.5
    assert any(s.startswith(CAP_STEP_PREFIX) for s in finding["magnitude_factors"])
    assert finding["imputation_capped"] is False


def test_cap_preserves_an_existing_audit_trail():
    finding = {"rsid": "rs1", "magnitude": 7.2, "imputed": True, "dr2": 0.42,
               "magnitude_factors": ["base 4.00 from clinvar"]}
    apply_imputation_cap(finding)
    assert finding["magnitude_factors"][0] == "base 4.00 from clinvar"
    assert finding["magnitude_factors"][-1].startswith(CAP_STEP_PREFIX)


def test_cap_sets_the_first_class_quality_fields():
    finding = {"rsid": "rs1", "magnitude": 5.0, "imputed": True, "dr2": 0.85}
    apply_imputation_cap(finding)
    assert finding["imputed"] is True
    assert finding["dr2"] == pytest.approx(0.85)
    assert finding["imputation_quality_band"] == "moderate"
    assert finding["magnitude_ceiling"] == IMPUTED_MAGNITUDE_CEILING


def test_cap_accepts_a_dr2_supplied_at_the_call_site():
    finding = {"rsid": "rs1", "magnitude": 9.0}
    apply_imputation_cap(finding, dr2=0.2)
    assert finding["imputed"] is True
    assert finding["magnitude"] == IMPUTED_MAGNITUDE_CAP


def test_cap_is_not_applied_to_a_typed_call():
    finding = {"rsid": "rs1", "magnitude": 9.0}
    apply_imputation_cap(finding)
    assert finding["magnitude"] == 9.0
    assert "imputed" not in finding
    assert "magnitude_factors" not in finding


def test_cap_respects_an_explicit_imputed_false():
    # A merged VCF can carry a DR2 on a row that was directly genotyped.
    finding = {"rsid": "rs1", "magnitude": 9.0, "imputed": False, "dr2": 0.1}
    apply_imputation_cap(finding)
    assert finding["magnitude"] == 9.0
    assert "magnitude_factors" not in finding


def test_cap_never_raises_a_magnitude():
    finding = {"rsid": "rs1", "magnitude": 0.5, "imputed": True, "dr2": 0.99}
    apply_imputation_cap(finding)
    assert finding["magnitude"] == 0.5


def test_cap_marks_an_untrusted_band_dubious():
    finding = {"rsid": "rs1", "magnitude": 6.0, "imputed": True, "dr2": 0.2}
    apply_imputation_cap(finding)
    assert finding["dubious"] is True


def test_cap_does_not_mark_a_high_quality_call_dubious():
    finding = {"rsid": "rs1", "magnitude": 6.0, "imputed": True, "dr2": 0.97}
    apply_imputation_cap(finding)
    assert finding.get("dubious") is not True


def test_cap_treats_a_missing_dr2_as_unknown_and_hard_caps_it():
    finding = {"rsid": "rs1", "magnitude": 8.0, "imputed": True, "dr2": None}
    apply_imputation_cap(finding)
    assert finding["imputation_quality_band"] == "unknown"
    assert finding["magnitude"] == IMPUTED_MAGNITUDE_CAP
    assert "not reported" in finding["magnitude_factors"][-1]


def test_cap_honours_a_custom_threshold():
    # DR2 0.6 is below the default threshold and above a 0.5 one, so lowering
    # the threshold lifts the ceiling from the hard cap to the parity ceiling.
    strict = {"rsid": "rs1", "magnitude": 8.0, "imputed": True, "dr2": 0.6}
    apply_imputation_cap(strict)
    assert strict["magnitude"] == IMPUTED_MAGNITUDE_CAP

    relaxed = {"rsid": "rs1", "magnitude": 8.0, "imputed": True, "dr2": 0.6}
    apply_imputation_cap(relaxed, threshold=0.5)
    assert relaxed["magnitude"] == 8.0
    assert relaxed["magnitude_ceiling"] == IMPUTED_MAGNITUDE_CEILING


def test_cap_all_processes_a_whole_list():
    findings = [
        {"rsid": "rs1", "magnitude": 9.0, "imputed": True, "dr2": 0.4},
        {"rsid": "rs2", "magnitude": 9.0},
        {"rsid": "rs3", "magnitude": 9.9, "imputed": True, "dr2": 0.99},
    ]
    apply_imputation_cap_all(findings)
    assert findings[0]["magnitude"] == IMPUTED_MAGNITUDE_CAP
    assert findings[1]["magnitude"] == 9.0
    assert findings[2]["magnitude"] == IMPUTED_MAGNITUDE_CEILING


def test_typed_versus_imputed_parity_after_capping():
    # Identical evidence, identical starting magnitude. The measurement must
    # end up strictly above the prediction.
    typed = {"rsid": "rs1", "magnitude": TYPED_MAGNITUDE_CEILING}
    imputed = {"rsid": "rs1", "magnitude": TYPED_MAGNITUDE_CEILING,
               "imputed": True, "dr2": 1.0}
    apply_imputation_cap(typed)
    apply_imputation_cap(imputed)
    assert imputed["magnitude"] < typed["magnitude"]


def test_cap_returns_non_dict_input_untouched():
    assert apply_imputation_cap("not a finding") == "not a finding"


# ---------------------------------------------------------------------------
# Degraded payloads: two failures, two messages
# ---------------------------------------------------------------------------

def test_missing_tool_payload_is_tagged_tool_missing():
    result = impute({"rs1": "AG"})
    assert result["available"] is False
    assert result["problem"] == "tool_missing"
    assert result["not_attempted"] is True
    assert result["results"] == []


def test_missing_tool_payload_explains_it_was_not_attempted():
    result = impute({"rs1": "AG"})
    assert "not attempted" in result["reason"].lower() or \
        "not installed" in result["reason"].lower()


def test_missing_panel_payload_is_tagged_panel_missing(monkeypatch):
    # Pretend Beagle is ready so the panel check is the thing that fires.
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    result = impute({"rs1": "AG"})
    assert result["available"] is False
    assert result["problem"] == "panel_missing"
    assert result["panel"] == "onekg_sgdp"


def test_missing_panel_and_missing_tool_are_different_messages(monkeypatch):
    tool_missing = impute({"rs1": "AG"})
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    panel_missing = impute({"rs1": "AG"})
    assert tool_missing["problem"] != panel_missing["problem"]
    assert tool_missing["reason"] != panel_missing["reason"]


def test_missing_panel_message_says_it_is_a_data_problem(monkeypatch):
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    result = impute({"rs1": "AG"})
    assert "DATA problem" in result["reason"]
    assert "not a tool problem" in result["reason"]


def test_panel_unavailable_names_an_unknown_panel():
    payload = panel_unavailable("no_such_panel", "imputation")
    assert payload["problem"] == "panel_missing"
    assert "Unknown reference panel" in payload["reason"]


def test_panel_unavailable_reports_a_partial_build(tmp_path, monkeypatch):
    base = Path(external.panel_root()) / "onekg_sgdp"
    base.mkdir(parents=True, exist_ok=True)
    (base / "panel.map").write_text("1 rs1 0 1000\n", encoding="utf-8")
    payload = panel_unavailable("onekg_sgdp", "imputation")
    assert payload["panel_state"] == "partial"
    assert "panel.vcf.gz" in payload["files_missing"]


def test_both_degraded_payloads_still_carry_the_mandatory_caveats(monkeypatch):
    tool_missing = impute({"rs1": "AG"})
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    panel_missing = impute({"rs1": "AG"})
    for payload in (tool_missing, panel_missing):
        for caveat in MANDATORY_CAVEATS:
            assert caveat in payload["caveats"]


# ---------------------------------------------------------------------------
# Filter tokens
# ---------------------------------------------------------------------------

def test_filter_tokens_publishes_exactly_the_three_documented_tokens():
    tokens = {t["token"] for t in filter_tokens()}
    assert tokens == {"/imputed", "/typed", "/dr2>=N"}


def test_every_filter_token_carries_every_declared_key():
    for token in filter_tokens():
        for key in FILTER_TOKEN_KEYS:
            assert key in token, f"{token['token']} is missing {key}"


def test_filter_token_kinds_match_the_two_families_filters_already_has():
    kinds = {t["name"]: t["kind"] for t in filter_tokens()}
    assert kinds == {"imputed": "flag", "typed": "flag", "dr2": "operator"}


def test_imputed_flag_predicate_matches_only_imputed_calls():
    predicate = next(t["predicate"] for t in filter_tokens() if t["name"] == "imputed")
    assert predicate({"imputed": True}) is True
    assert predicate({"imputed": False}) is False
    assert predicate({}) is False


def test_typed_flag_predicate_matches_a_finding_that_never_saw_imputation():
    # Absence of the key is the common case and must resolve to typed.
    predicate = next(t["predicate"] for t in filter_tokens() if t["name"] == "typed")
    assert predicate({}) is True
    assert predicate({"imputed": False}) is True
    assert predicate({"imputed": True}) is False


def test_dr2_operator_predicate_supports_every_declared_operator():
    token = next(t for t in filter_tokens() if t["name"] == "dr2")
    predicate = token["predicate"]
    finding = {"dr2": 0.8}
    assert predicate(finding, ">=", 0.8) is True
    assert predicate(finding, "<=", 0.8) is True
    assert predicate(finding, ">", 0.8) is False
    assert predicate(finding, "<", 0.9) is True
    assert predicate(finding, "=", 0.8) is True


def test_dr2_operator_never_matches_a_finding_with_no_dr2():
    predicate = next(t["predicate"] for t in filter_tokens() if t["name"] == "dr2")
    assert predicate({}, ">=", 0.5) is False
    assert predicate({}, "<", 0.5) is False
    assert predicate({"dr2": None}, "<", 0.5) is False


def test_filter_patterns_are_regex_alternation_fragments():
    assert set(filter_flag_pattern().split("|")) == {"imputed", "typed"}
    assert filter_op_pattern() == "dr2"


# ---------------------------------------------------------------------------
# Caveats
# ---------------------------------------------------------------------------

def test_caveats_always_include_every_mandatory_caveat():
    for caveat in MANDATORY_CAVEATS:
        assert caveat in build_caveats(None)


def test_caveats_always_include_the_reference_panel_ancestry_bias_caveat():
    assert CAVEAT_PANEL_BIAS in build_caveats(None)
    assert CAVEAT_PANEL_BIAS in build_caveats({"imputed": 10, "above_threshold": 9})


def test_caveats_always_include_the_rare_variant_caveat():
    assert CAVEAT_RARE_VARIANTS in build_caveats(None)
    assert "1 percent" in CAVEAT_RARE_VARIANTS


def test_caveats_carry_the_real_coverage_numbers():
    coverage = coverage_report(parse_vcf(imputed_vcf_text()))
    text = " ".join(build_caveats(coverage))
    assert "4 imputed variants" in text
    assert "2 (50.0 percent)" in text


def test_caveats_append_the_panels_own_recorded_note():
    caveats = build_caveats(None)
    note = external.panel_status("onekg_sgdp")["note"]
    assert note in caveats


# ---------------------------------------------------------------------------
# Coverage report
# ---------------------------------------------------------------------------

def test_coverage_report_counts_typed_and_imputed_separately():
    report = coverage_report(parse_vcf(imputed_vcf_text()))
    assert report["typed"] == 1
    assert report["imputed"] == 4
    assert report["total"] == 5


def test_coverage_report_splits_on_the_dr2_threshold():
    report = coverage_report(parse_vcf(imputed_vcf_text()))
    assert report["above_threshold"] == 2
    assert report["below_threshold"] == 2
    assert report["usable_fraction"] == pytest.approx(0.5)


def test_coverage_report_counts_every_band():
    report = coverage_report(parse_vcf(imputed_vcf_text()))
    assert report["bands"]["high"] == 1
    assert report["bands"]["moderate"] == 1
    assert report["bands"]["low"] == 1
    assert report["bands"]["unusable"] == 1
    assert sum(report["bands"].values()) == report["imputed"]


def test_coverage_report_counts_rare_variants_separately():
    report = coverage_report(parse_vcf(imputed_vcf_text()))
    assert report["rare_variants"] == 2
    assert report["rare_above_threshold"] == 0


def test_coverage_report_counts_a_missing_dr2_as_unknown():
    records = parse_vcf(
        "1\t1\trs1\tA\tG\t.\tPASS\tAF=0.4;IMP\tGT\t0|1\n")
    report = coverage_report(records)
    assert report["unknown_dr2"] == 1
    assert report["bands"]["unknown"] == 1
    assert report["mean_dr2"] is None


def test_coverage_report_on_nothing_is_all_zeroes_and_does_not_raise():
    report = coverage_report([])
    assert report["imputed"] == 0
    assert report["usable_fraction"] == 0.0
    assert report["panel_available"] is False


def test_coverage_report_computes_mean_and_median_dr2():
    report = coverage_report(parse_vcf(imputed_vcf_text()))
    assert report["mean_dr2"] == pytest.approx((0.99 + 0.85 + 0.50 + 0.10) / 4)
    assert report["median_dr2"] == pytest.approx(0.675)


def test_coverage_report_honours_a_custom_threshold():
    report = coverage_report(parse_vcf(imputed_vcf_text()), dr2_threshold=0.4)
    assert report["above_threshold"] == 3
    assert report["dr2_threshold"] == pytest.approx(0.4)


# ---------------------------------------------------------------------------
# The safety invariant
# ---------------------------------------------------------------------------

def constructed_violation() -> dict:
    """An imputed pathogenic call carrying no quality evidence at all."""
    return {
        "rsid": "rs28897696",
        "gene": "BRCA1",
        "clinical_sig": "Pathogenic",
        "clinvar_sig_code": 5,
        "imputed": True,
        "magnitude": 9.5,
    }


def test_safety_check_raises_on_a_constructed_violation():
    with pytest.raises(ImputedQualityViolation):
        assert_no_imputed_pathogenic_without_quality([constructed_violation()])


def test_safety_check_carries_the_offending_findings_on_the_exception():
    with pytest.raises(ImputedQualityViolation) as caught:
        assert_no_imputed_pathogenic_without_quality([constructed_violation()])
    assert caught.value.violations[0]["rsid"] == "rs28897696"
    assert "no DR2" in caught.value.violations[0]["reason"]


def test_safety_check_can_audit_without_raising():
    violations = assert_no_imputed_pathogenic_without_quality(
        [constructed_violation()], strict=False)
    assert len(violations) == 1


def test_safety_check_passes_a_properly_capped_high_quality_call():
    finding = dict(constructed_violation(), dr2=0.97)
    apply_imputation_cap(finding)
    assert assert_no_imputed_pathogenic_without_quality([finding]) == []


def test_safety_check_catches_a_capped_but_low_quality_pathogenic_call():
    finding = dict(constructed_violation(), dr2=0.35)
    apply_imputation_cap(finding)
    violations = assert_no_imputed_pathogenic_without_quality(
        [finding], strict=False)
    assert violations and "below the" in violations[0]["reason"]


def test_safety_check_catches_an_uncapped_call_even_with_a_good_dr2():
    # A good DR2 that never went through the cap can still outrank a typed
    # call, which is the parity failure the cap exists to prevent.
    finding = dict(constructed_violation(), dr2=0.99)
    violations = assert_no_imputed_pathogenic_without_quality(
        [finding], strict=False)
    assert violations and "cap never ran" in violations[0]["reason"]


def test_safety_check_ignores_typed_pathogenic_calls():
    finding = dict(constructed_violation())
    finding.pop("imputed")
    assert assert_no_imputed_pathogenic_without_quality([finding]) == []


def test_safety_check_ignores_imputed_benign_calls():
    finding = dict(constructed_violation(),
                   clinical_sig="Benign", clinvar_sig_code=2)
    assert assert_no_imputed_pathogenic_without_quality([finding]) == []


def test_safety_check_catches_pathogenic_declared_only_in_free_text():
    finding = {"rsid": "rs1", "imputed": True, "clinical_sig": "Likely pathogenic",
               "magnitude": 8.0}
    violations = assert_no_imputed_pathogenic_without_quality(
        [finding], strict=False)
    assert len(violations) == 1


def test_safety_check_on_an_empty_list_is_clean():
    assert assert_no_imputed_pathogenic_without_quality([]) == []
    assert assert_no_imputed_pathogenic_without_quality(None) == []


# ---------------------------------------------------------------------------
# Study VCF writing
# ---------------------------------------------------------------------------

def test_write_study_vcf_writes_one_row_per_placed_call(tmp_path):
    calls = [
        {"rsid": "rs1", "chromosome": "1", "position": 1000,
         "allele1": "A", "allele2": "G", "ref": "A"},
        {"rsid": "rs2", "chromosome": "1", "position": 2000,
         "allele1": "C", "allele2": "C", "ref": "C"},
    ]
    result = write_study_vcf(calls, tmp_path / "study.vcf")
    assert result["written"] == 2
    body = Path(result["path"]).read_text(encoding="utf-8")
    assert "rs1" in body and "0/1" in body
    assert body.startswith("##fileformat=VCFv4.2")


def test_write_study_vcf_flags_an_assumed_reference_allele(tmp_path):
    calls = [{"rsid": "rs1", "chromosome": "1", "position": 1000,
              "allele1": "A", "allele2": "G"}]
    result = write_study_vcf(calls, tmp_path / "study.vcf")
    assert result["ref_assumed"] == 1
    assert "DNAI_REF_ASSUMED" in Path(result["path"]).read_text(encoding="utf-8")
    assert any("reference allele" in w for w in result["warnings"])


def test_write_study_vcf_skips_calls_with_no_coordinates(tmp_path):
    result = write_study_vcf({"rs1": "AG"}, tmp_path / "study.vcf")
    assert result["written"] == 0
    assert result["skipped"][0]["reason"] == "no chromosome or position"


def test_write_study_vcf_skips_no_calls(tmp_path):
    calls = [{"rsid": "rs1", "chromosome": "1", "position": 1000,
              "allele1": "-", "allele2": "-"}]
    result = write_study_vcf(calls, tmp_path / "study.vcf")
    assert result["written"] == 0
    assert result["skipped"][0]["reason"] == "no call"


def test_write_study_vcf_output_round_trips_through_the_parser(tmp_path):
    calls = [{"rsid": "rs1", "chromosome": "1", "position": 1000,
              "allele1": "A", "allele2": "G", "ref": "A"}]
    result = write_study_vcf(calls, tmp_path / "study.vcf")
    records = parse_vcf(Path(result["path"]))
    assert len(records) == 1
    assert records[0]["rsid"] == "rs1"
    assert records[0]["imputed"] is False


def test_impute_reports_no_usable_calls_distinctly(monkeypatch, tmp_path):
    # Tool ready, panel ready, but the genotype map carries no coordinates.
    monkeypatch.setattr(external, "is_available", lambda tool_id: True)
    base = Path(external.panel_root()) / "onekg_sgdp"
    base.mkdir(parents=True, exist_ok=True)
    for name in external.PANELS["onekg_sgdp"]["files"]:
        (base / name).write_text("stub\n", encoding="utf-8")
    result = impute({"rs1": "AG"}, workdir=tmp_path / "work")
    assert result["problem"] == "no_input"
    assert result["not_attempted"] is True
