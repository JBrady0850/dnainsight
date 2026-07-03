"""Unit tests for backend/scanner.py — offline scanner logic only."""

import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.scanner import (
    _conditions_risk,
    classify_silo,
    lookup_bundled,
    zygosity_of,
    carries_allele,
    _extract_alt,
)


# ---------------------------------------------------------------------------
# _conditions_risk
# ---------------------------------------------------------------------------

class TestConditionsRisk:
    def test_risk_detected_poor_metabolizer(self):
        assert _conditions_risk("poor metabolizer for clopidogrel", "AG") == "risk"

    def test_risk_detected_contraindicated(self):
        assert _conditions_risk("drug is contraindicated in this genotype", "GG") == "risk"

    def test_safe_detected_normal_metabolizer(self):
        assert _conditions_risk("normal metabolizer, standard response expected", "CC") == "safe"

    def test_safe_detected_protective(self):
        assert _conditions_risk("protective allele, lower risk of disease", "TT") == "safe"

    def test_unknown_returned_for_neutral(self):
        assert _conditions_risk("variant of uncertain significance", "AG") == "unknown"

    def test_risk_beats_safe_when_both_present(self):
        # "no function" (risk) + "no effect" (safe) — risk count should win or tie
        result = _conditions_risk("no function allele, no significant adverse event reported", "AG")
        # Depends on score balance; acceptable values are 'risk' or 'unknown'
        assert result in ("risk", "unknown")


# ---------------------------------------------------------------------------
# classify_silo
# ---------------------------------------------------------------------------

class TestClassifySilo:
    def test_pre_prescription_warfarin(self):
        assert classify_silo(
            "VKORC1 variant: higher warfarin sensitivity. Reduce starting dose per CPIC."
        ) == "pre_prescription"

    def test_pre_prescription_poor_metabolizer(self):
        assert classify_silo(
            "CYP2C19 poor metabolizer phenotype. Avoid clopidogrel; consider alternative antiplatelet."
        ) == "pre_prescription"

    def test_pre_prescription_statin(self):
        assert classify_silo(
            "SLCO1B1: reduced statin transport. Use lowest effective statin dose."
        ) == "pre_prescription"

    def test_actionable_risk(self):
        silo = classify_silo("FTO obesity risk allele. Associated with increased risk of obesity.")
        assert silo in ("pre_prescription", "actionable")

    def test_informational_oxtr(self):
        assert classify_silo(
            "OXTR rs53576: GG homozygotes demonstrate higher empathy and prosocial behavior."
        ) == "informational"

    def test_informational_clock(self):
        assert classify_silo(
            "CLOCK gene variant: associated with eveningness chronotype."
        ) == "informational"


# ---------------------------------------------------------------------------
# lookup_bundled — integration test against the built snp_reference.json
# ---------------------------------------------------------------------------

class TestLookupBundled:
    """These tests require data/snp_reference.json to exist (run build_reference.py first)."""

    def test_mthfr_c677t_in_reference(self):
        result = lookup_bundled("rs1801133")
        if result is None:
            pytest.skip("snp_reference.json not built — run: python data/build_reference.py")
        assert result["gene"] == "MTHFR"
        assert result["category"] == "NEURO"

    def test_slco1b1_in_reference(self):
        result = lookup_bundled("rs4149056")
        if result is None:
            pytest.skip("snp_reference.json not built")
        assert result["gene"] == "SLCO1B1"
        assert result["category"] == "PHARM"

    def test_comt_in_reference(self):
        result = lookup_bundled("rs4680")
        if result is None:
            pytest.skip("snp_reference.json not built")
        assert result["gene"] == "COMT"

    def test_dpyd_cpic_in_reference(self):
        """New v1.1 CPIC Level A addition."""
        result = lookup_bundled("rs3918290")
        if result is None:
            pytest.skip("snp_reference.json not built")
        assert result["gene"] == "DPYD"

    def test_nudt15_in_reference(self):
        """New v1.1 addition."""
        result = lookup_bundled("rs116855232")
        if result is None:
            pytest.skip("snp_reference.json not built")
        assert result["gene"] == "NUDT15"

    def test_nonexistent_rsid_returns_none(self):
        assert lookup_bundled("rs9999999999") is None

    def test_reference_exceeds_120_snps(self):
        """v1.1 roadmap requires 120+ bundled SNPs."""
        from backend.scanner import _load_bundled
        ref = _load_bundled()
        if not ref:
            pytest.skip("snp_reference.json not built")
        assert len(ref) >= 120, f"Expected 120+ SNPs, found {len(ref)}"


# ---------------------------------------------------------------------------
# zygosity_of  (v1.2 -- factual genotype classification)
# ---------------------------------------------------------------------------

class TestZygosity:
    def test_homozygous(self):
        assert zygosity_of("A", "A") == "homozygous"

    def test_heterozygous(self):
        assert zygosity_of("A", "G") == "heterozygous"

    def test_hemizygous_single_call(self):
        assert zygosity_of("A", "N") == "hemizygous"

    def test_no_call(self):
        assert zygosity_of("N", "N") == "no_call"
        assert zygosity_of("", "") == "no_call"

    def test_case_insensitive(self):
        assert zygosity_of("c", "C") == "homozygous"


# ---------------------------------------------------------------------------
# carries_allele  (v1.2 -- authoritative carrier status)
# ---------------------------------------------------------------------------

class TestCarriesAllele:
    def test_zero_copies_reference_genotype(self):
        # rs671 GG when the variant/risk allele is A -> not a carrier
        assert carries_allele("G", "G", "A") == 0

    def test_one_copy_heterozygous(self):
        assert carries_allele("A", "G", "A") == 1

    def test_two_copies_homozygous(self):
        assert carries_allele("A", "A", "A") == 2

    def test_empty_variant_allele_returns_zero(self):
        assert carries_allele("A", "A", "") == 0

    def test_case_insensitive(self):
        assert carries_allele("a", "g", "A") == 1


# ---------------------------------------------------------------------------
# _extract_alt  (v1.2 -- pull authoritative alt allele from a MyVariant hit)
# ---------------------------------------------------------------------------

class TestExtractAlt:
    def test_scalar_alt(self):
        assert _extract_alt({"dbsnp": {"alt": "A"}}) == "A"

    def test_list_alt_takes_first(self):
        assert _extract_alt({"dbsnp": {"alt": ["T", "C"]}}) == "T"

    def test_multiallelic_string(self):
        assert _extract_alt({"dbsnp": {"alt": "G,A"}}) == "G"

    def test_fallback_parses_id(self):
        assert _extract_alt({"_id": "chr1:g.11856378G>A"}) == "A"

    def test_missing_returns_empty(self):
        assert _extract_alt({}) == ""
