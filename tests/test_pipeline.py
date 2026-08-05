"""Tests for backend.pipeline, the v2 scan orchestrator.

Covers docs/API_V2.md sections 1, 2, 3.3 and 3.4. No network: every test runs
against the bundled offline reference, and anything that needs a generated data
file skips rather than fails when that file is absent.

The order of the stages is the point of this module. Strand has to be resolved
before frequency and before scoring, because a flipped genotype looks up the
wrong frequency and then scores off it.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import pipeline as pipeline_module
from backend.pipeline import (
    DEFAULT_OPTIONS,
    PHASES,
    apply_orientation,
    attach_provenance,
    available_subsystems,
    build_genotype_map,
    collect_genosets,
    collect_prs,
    collect_traits,
    compute_qc,
    compute_ranges,
    compute_summary,
    enrich_findings,
    run_full_scan,
    snps_from_merged,
)
from backend.pipeline import _ensure_contract_keys

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"
REFERENCE_FILE = DATA / "snp_reference.json"

# Phases the contract publishes for the scan status endpoint, section 3.3.
DOCUMENTED_PHASES = ("bundled", "orientation", "frequency", "genosets",
                     "traits", "prs", "api", "writing", "complete")

# Options the scan body publishes, section 3.3.
DOCUMENTED_OPTIONS = ("use_api", "population", "include_genosets",
                      "include_traits", "include_prs", "use_snpedia")

# ---------------------------------------------------------------------------
# A synthetic three-file source set: two "self" exports that disagree at one
# position, plus a mother for the comparison rows.
# ---------------------------------------------------------------------------

SELF_A = [
    {"rsid": "rs1801133", "chromosome": "1", "position": 11856378,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs1065852", "chromosome": "22", "position": 42526694,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs429358", "chromosome": "19", "position": 45411941,
     "allele1": "T", "allele2": "T"},
    {"rsid": "rs7412", "chromosome": "19", "position": 45412079,
     "allele1": "T", "allele2": "T"},
    {"rsid": "rs4988235", "chromosome": "2", "position": 136608646,
     "allele1": "T", "allele2": "T"},
    {"rsid": "rs4244285", "chromosome": "10", "position": 96541616,
     "allele1": "G", "allele2": "A"},
    {"rsid": "rs1799853", "chromosome": "10", "position": 96702047,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs7903146", "chromosome": "10", "position": 114758349,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs12255372", "chromosome": "10", "position": 114808902,
     "allele1": "G", "allele2": "T"},
    {"rsid": "rs6025", "chromosome": "1", "position": 169519049,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs3892097", "chromosome": "22", "position": 42524947,
     "allele1": "G", "allele2": "A"},
    {"rsid": "rs4680", "chromosome": "22", "position": 19951271,
     "allele1": "-", "allele2": "-"},
]

# The deliberate conflict: rs1801133 is CT in the first file and CC here.
SELF_B = [
    {"rsid": "rs1801133", "chromosome": "1", "position": 11856378,
     "allele1": "C", "allele2": "C"},
    {"rsid": "rs1065852", "chromosome": "22", "position": 42526694,
     "allele1": "C", "allele2": "T"},
    {"rsid": "rs1801131", "chromosome": "1", "position": 11854476,
     "allele1": "A", "allele2": "C"},
]

MOTHER = [
    {"rsid": "rs1801133", "chromosome": "1", "position": 11856378,
     "allele1": "T", "allele2": "C"},
    {"rsid": "rs429358", "chromosome": "19", "position": 45411941,
     "allele1": "T", "allele2": "C"},
    {"rsid": "rs4988235", "chromosome": "2", "position": 136608646,
     "allele1": "G", "allele2": "G"},
]

SOURCES = [
    {"label": "self_a.txt", "role": "self", "provider": "23andme",
     "snps": SELF_A},
    {"label": "self_b.txt", "role": "self", "provider": "ancestrydna",
     "snps": SELF_B},
    {"label": "mum.txt", "role": "mother", "provider": "23andme",
     "snps": MOTHER},
]

CONFLICT_RSID = "rs1801133"

# Genotypes that hit a genoset rule, a trait rule and a polygenic model.
GENOTYPES = {
    "rs429358": ("T", "T"), "rs7412": ("T", "T"), "rs4988235": ("T", "T"),
    "rs1801133": ("C", "T"), "rs7903146": ("C", "T"), "rs12255372": ("G", "T"),
    "rs4244285": ("G", "A"), "rs6025": ("C", "T"),
}

# Every key docs/API_V2.md section 2 promises on a finding.
SECTION_2_1 = ("rsid", "entity_type", "gene", "chromosome", "position",
               "allele1", "allele2", "genotype", "token", "zygosity")
SECTION_2_2 = ("magnitude", "magnitude_source", "repute", "summary",
               "interpretation", "confidence")
SECTION_2_3 = ("clinical_sig", "clinvar_sig_code", "review_status",
               "review_stars", "cpic_level", "pgx_level", "evidence",
               "publications", "conditions", "conditions_list", "sources")
SECTION_2_4 = ("freq", "freq_population", "freq_band", "freq_color",
               "freq_derived", "freq_method", "freq_flipped",
               "freq_ambiguous", "gmaf", "minor_allele", "population_series")
SECTION_2_5 = ("orientation", "stabilized_orientation", "flipped",
               "ambiguous", "dubious", "variant_allele", "variant_copies",
               "carrier")
SECTION_2_6 = ("count", "labels", "conflict", "calls", "comparison",
               "probability", "mendelian_ok", "topics", "medicines")
SECTION_2_7 = ("criteria", "matched_rsids", "coverage", "percentile", "band",
               "reliable", "caveats")
SECTION_2 = (SECTION_2_1 + SECTION_2_2 + SECTION_2_3 + SECTION_2_4
             + SECTION_2_5 + SECTION_2_6 + SECTION_2_7)

SUBSYSTEMS = available_subsystems()

_SCAN = {}


def require_reference():
    """Skip when the bundled reference has not been generated yet."""
    if not REFERENCE_FILE.exists():
        pytest.skip("data/snp_reference.json has not been generated")


def require_subsystem(name):
    """Skip when an optional corpus has not been generated yet."""
    if not SUBSYSTEMS.get(name):
        pytest.skip("the %s corpus has not been generated" % name)


def scan():
    """Run the synthetic end to end scan once and memoise the payload."""
    require_reference()
    if "payload" not in _SCAN:
        phases = []
        _SCAN["payload"] = run_full_scan(
            SOURCES,
            use_api=False,
            population="CEU",
            include_genosets=True,
            include_traits=True,
            include_prs=True,
            use_snpedia=False,
            progress_cb=lambda phase, done=0, total=0: phases.append(phase),
        )
        _SCAN["phases"] = phases
    return _SCAN["payload"]


def scan_phases():
    """Phase names the memoised scan reported, in order."""
    scan()
    return _SCAN["phases"]


def merged_fixture():
    """A hand-built merged genotype set, so provenance tests need no data."""
    return {
        "genotypes": {
            "rs1801133": {
                "rsid": "rs1801133", "chromosome": "1", "position": 11856378,
                "allele1": "C", "allele2": "T", "genotype": "CT", "count": 2,
                "labels": ["self_a.txt", "self_b.txt"], "conflict": True,
                "calls": [
                    {"label": "self_a.txt", "allele1": "C", "allele2": "T",
                     "genotype": "CT"},
                    {"label": "self_b.txt", "allele1": "C", "allele2": "C",
                     "genotype": "CC"},
                ],
            },
        },
        "comparison": {
            "mother": {"rs1801133": {"label": "mum.txt", "role": "mother",
                                     "allele1": "C", "allele2": "T",
                                     "genotype": "TC"}},
            "father": {"rs1801133": {"label": "dad.txt", "role": "father",
                                     "allele1": "G", "allele2": "G",
                                     "genotype": "GG"}},
        },
        "counts": {"conflicts": 1, "pooled_sources": 2, "total_positions": 14},
    }


class TestPhases:
    def test_every_documented_phase_is_present(self):
        for phase in DOCUMENTED_PHASES:
            assert phase in PHASES

    def test_the_pooling_and_scoring_stages_are_present(self):
        assert "merge" in PHASES
        assert "scoring" in PHASES

    def test_the_snpedia_overlay_stage_is_present(self):
        assert "snpedia" in PHASES

    def test_the_phase_names_are_unique(self):
        assert len(set(PHASES)) == len(PHASES)

    def test_the_phases_are_strings(self):
        assert all(isinstance(phase, str) for phase in PHASES)

    def test_complete_is_the_final_phase(self):
        assert PHASES[-1] == "complete"

    def test_merge_is_the_first_phase(self):
        assert PHASES[0] == "merge"

    def test_orientation_precedes_frequency_and_scoring(self):
        assert PHASES.index("orientation") < PHASES.index("frequency")
        assert PHASES.index("frequency") < PHASES.index("scoring")


class TestDefaultOptions:
    def test_the_documented_keys_are_present(self):
        for key in DOCUMENTED_OPTIONS:
            assert key in DEFAULT_OPTIONS

    def test_there_are_no_undocumented_options(self):
        assert set(DEFAULT_OPTIONS) == set(DOCUMENTED_OPTIONS)

    def test_the_network_path_is_off_by_default(self):
        assert DEFAULT_OPTIONS["use_api"] is False

    def test_the_default_population_is_ceu(self):
        assert DEFAULT_OPTIONS["population"] == "CEU"

    def test_the_extra_entity_types_are_on_by_default(self):
        assert DEFAULT_OPTIONS["include_genosets"] is True
        assert DEFAULT_OPTIONS["include_traits"] is True
        assert DEFAULT_OPTIONS["include_prs"] is True

    def test_the_snpedia_overlay_is_off_by_default(self):
        assert DEFAULT_OPTIONS["use_snpedia"] is False


class TestBuildGenotypeMap:
    def test_rsids_are_lowercased(self):
        result = build_genotype_map({"RS1801133": {"allele1": "C",
                                                   "allele2": "T"}})
        assert "rs1801133" in result

    def test_the_original_case_is_not_kept(self):
        result = build_genotype_map({"RS1801133": {"allele1": "C",
                                                   "allele2": "T"}})
        assert "RS1801133" not in result

    def test_an_allele_pair_is_returned(self):
        result = build_genotype_map({"rs1": {"allele1": "C", "allele2": "T"}})
        assert result["rs1"] == ("C", "T")

    def test_alleles_are_uppercased(self):
        result = build_genotype_map({"rs1": {"allele1": "c", "allele2": "t"}})
        assert result["rs1"] == ("C", "T")

    def test_allele_whitespace_is_stripped(self):
        result = build_genotype_map({"rs1": {"allele1": " C ", "allele2": "T"}})
        assert result["rs1"] == ("C", "T")

    def test_rsid_whitespace_is_stripped(self):
        result = build_genotype_map({" rs1 ": {"allele1": "C", "allele2": "T"}})
        assert "rs1" in result

    def test_a_no_call_stays_a_no_call(self):
        result = build_genotype_map({"rs1": {"allele1": "N", "allele2": "N"}})
        assert result["rs1"] == ("N", "N")

    def test_a_missing_allele_is_not_imputed(self):
        result = build_genotype_map({"rs1": {}})
        assert result["rs1"] == ("", "")

    def test_a_blank_allele_is_not_imputed(self):
        result = build_genotype_map({"rs1": {"allele1": "", "allele2": ""}})
        assert result["rs1"] == ("", "")

    def test_every_position_is_carried_over(self):
        entries = {"rs1": {"allele1": "A", "allele2": "A"},
                   "rs2": {"allele1": "C", "allele2": "G"}}
        assert len(build_genotype_map(entries)) == 2

    def test_an_empty_map_is_safe(self):
        assert build_genotype_map({}) == {}

    def test_a_none_map_is_safe(self):
        assert build_genotype_map(None) == {}


class TestSnpsFromMerged:
    def test_the_scanner_keys_are_produced(self):
        entries = {"rs1": {"chromosome": "1", "position": 10,
                           "allele1": "C", "allele2": "T"}}
        assert set(snps_from_merged(entries)[0]) == {
            "rsid", "chromosome", "position", "allele1", "allele2"}

    def test_the_values_round_trip(self):
        entries = {"rs1801133": {"chromosome": "1", "position": 11856378,
                                 "allele1": "C", "allele2": "T"}}
        snp = snps_from_merged(entries)[0]
        assert snp["rsid"] == "rs1801133"
        assert snp["chromosome"] == "1"
        assert snp["position"] == 11856378
        assert snp["allele1"] == "C"
        assert snp["allele2"] == "T"

    def test_one_row_per_position(self):
        entries = {"rs1": {"allele1": "A", "allele2": "A"},
                   "rs2": {"allele1": "C", "allele2": "C"}}
        assert len(snps_from_merged(entries)) == 2

    def test_missing_fields_become_documented_defaults(self):
        snp = snps_from_merged({"rs1": {}})[0]
        assert snp["chromosome"] == ""
        assert snp["position"] == 0
        assert snp["allele1"] == ""

    def test_a_no_call_round_trips_as_a_no_call(self):
        entries = {"rs1": {"allele1": "N", "allele2": "N"}}
        snp = snps_from_merged(entries)[0]
        assert (snp["allele1"], snp["allele2"]) == ("N", "N")

    def test_an_empty_map_gives_an_empty_list(self):
        assert snps_from_merged({}) == []

    def test_a_none_map_gives_an_empty_list(self):
        assert snps_from_merged(None) == []

    def test_the_genotype_map_and_the_snp_list_agree_on_size(self):
        entries = {"rs1": {"allele1": "A", "allele2": "G"},
                   "rs2": {"allele1": "C", "allele2": "C"}}
        assert len(snps_from_merged(entries)) == len(build_genotype_map(entries))


class TestAttachProvenance:
    def test_the_pooled_call_count_is_copied(self):
        finding = attach_provenance({"rsid": "rs1801133"}, merged_fixture())
        assert finding["count"] == 2

    def test_the_contributing_labels_are_copied(self):
        finding = attach_provenance({"rsid": "rs1801133"}, merged_fixture())
        assert finding["labels"] == ["self_a.txt", "self_b.txt"]

    def test_the_conflict_flag_is_copied(self):
        finding = attach_provenance({"rsid": "rs1801133"}, merged_fixture())
        assert finding["conflict"] is True

    def test_both_conflicting_calls_are_retained(self):
        finding = attach_provenance({"rsid": "rs1801133"}, merged_fixture())
        assert [call["genotype"] for call in finding["calls"]] == ["CT", "CC"]

    def test_nothing_is_reconciled_into_a_winner(self):
        finding = attach_provenance({"rsid": "rs1801133"}, merged_fixture())
        assert len(finding["calls"]) == 2

    def test_rsid_matching_is_case_insensitive(self):
        finding = attach_provenance({"rsid": "RS1801133"}, merged_fixture())
        assert finding["count"] == 2

    def test_an_unknown_position_gets_empty_provenance(self):
        finding = attach_provenance({"rsid": "rs404"}, merged_fixture())
        assert finding["count"] == 0
        assert finding["labels"] == []
        assert finding["conflict"] is False
        assert finding["calls"] == []

    def test_a_comparison_row_is_built_per_relative(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        assert len(finding["comparison"]) == 2

    def test_the_comparison_row_shape(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        assert set(finding["comparison"][0]) == {"label", "role", "genotype",
                                                 "shared"}

    def test_the_comparison_row_carries_the_relative_label_and_role(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        rows = {row["role"]: row["label"] for row in finding["comparison"]}
        assert rows["mother"] == "mum.txt"
        assert rows["father"] == "dad.txt"

    def test_shared_is_true_when_the_relative_matches(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        rows = {row["role"]: row["shared"] for row in finding["comparison"]}
        assert rows["mother"] is True

    def test_shared_is_allele_order_insensitive(self):
        # The mother's call is written TC and the finding's is CT.
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        mother = next(r for r in finding["comparison"] if r["role"] == "mother")
        assert mother["genotype"] == "TC"
        assert mother["shared"] is True

    def test_shared_is_false_when_the_relative_differs(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged_fixture())
        rows = {row["role"]: row["shared"] for row in finding["comparison"]}
        assert rows["father"] is False

    def test_shared_is_false_when_the_finding_has_no_genotype(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": ""},
                                    merged_fixture())
        assert all(row["shared"] is False for row in finding["comparison"])

    def test_shared_is_false_when_the_relative_has_no_genotype(self):
        merged = merged_fixture()
        merged["comparison"]["mother"]["rs1801133"]["genotype"] = ""
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged)
        mother = next(r for r in finding["comparison"] if r["role"] == "mother")
        assert mother["shared"] is False

    def test_a_relative_without_the_position_gets_no_row(self):
        merged = merged_fixture()
        del merged["comparison"]["father"]["rs1801133"]
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "CT"},
                                    merged)
        assert [row["role"] for row in finding["comparison"]] == ["mother"]

    def test_no_relatives_gives_no_comparison_rows(self):
        merged = merged_fixture()
        merged["comparison"] = {}
        finding = attach_provenance({"rsid": "rs1801133"}, merged)
        assert finding["comparison"] == []

    def test_the_genotype_case_is_normalised_before_comparing(self):
        finding = attach_provenance({"rsid": "rs1801133", "genotype": "ct"},
                                    merged_fixture())
        mother = next(r for r in finding["comparison"] if r["role"] == "mother")
        assert mother["shared"] is True


class TestApplyOrientation:
    def test_the_five_documented_keys_are_set(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T"})
        for key in ("orientation", "stabilized_orientation", "flipped",
                    "ambiguous", "token"):
            assert key in finding

    def test_orientation_defaults_to_unknown(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T"})
        assert finding["orientation"] == ""

    def test_an_existing_orientation_is_kept(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "orientation": "plus"})
        assert finding["orientation"] == "plus"

    def test_the_stabilized_orientation_is_normalised(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "stabilized_orientation": "  MINUS  "})
        assert finding["stabilized_orientation"] == "minus"

    def test_a_missing_stabilized_orientation_becomes_empty(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T"})
        assert finding["stabilized_orientation"] == ""

    def test_an_at_heterozygote_is_ambiguous(self):
        finding = apply_orientation({"allele1": "A", "allele2": "T"})
        assert finding["ambiguous"] is True

    def test_a_cg_heterozygote_is_ambiguous(self):
        finding = apply_orientation({"allele1": "G", "allele2": "C"})
        assert finding["ambiguous"] is True

    def test_an_ag_heterozygote_is_not_ambiguous(self):
        finding = apply_orientation({"allele1": "A", "allele2": "G"})
        assert finding["ambiguous"] is False

    def test_a_homozygote_is_not_ambiguous(self):
        finding = apply_orientation({"allele1": "A", "allele2": "A"})
        assert finding["ambiguous"] is False

    def test_a_no_call_is_not_ambiguous(self):
        finding = apply_orientation({"allele1": "-", "allele2": "-"})
        assert finding["ambiguous"] is False

    def test_lowercase_alleles_are_still_classified(self):
        finding = apply_orientation({"allele1": "a", "allele2": "t"})
        assert finding["ambiguous"] is True

    def test_a_minus_stabilized_orientation_can_set_flipped(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "stabilized_orientation": "minus"})
        assert finding["flipped"] is True

    def test_a_minus_strand_records_the_snpedia_token(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "stabilized_orientation": "minus"})
        assert finding["snpedia_token"] == "(A;G)"

    def test_a_plus_stabilized_orientation_does_not_flip(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "stabilized_orientation": "plus"})
        assert finding["flipped"] is False

    def test_an_absent_stabilized_orientation_does_not_flip(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T"})
        assert finding["flipped"] is False

    def test_the_token_is_built_from_the_reported_alleles(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T"})
        assert finding["token"] == "(C;T)"

    def test_an_existing_token_is_not_rebuilt(self):
        finding = apply_orientation({"allele1": "C", "allele2": "T",
                                     "token": "(A;G)"})
        assert finding["token"] == "(A;G)"

    def test_a_finding_with_no_alleles_gets_an_empty_token(self):
        finding = apply_orientation({"rsid": "dgs001"})
        assert finding["token"] == ""

    def test_a_half_called_position_gets_an_empty_token(self):
        finding = apply_orientation({"allele1": "C", "allele2": ""})
        assert finding["token"] == ""

    def test_the_finding_is_updated_in_place(self):
        finding = {"allele1": "C", "allele2": "T"}
        assert apply_orientation(finding) is finding


class TestEnsureContractKeys:
    def test_the_identity_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_1:
            if key != "rsid":
                assert key in finding

    def test_the_interest_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_2:
            assert key in finding

    def test_the_evidence_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_3:
            assert key in finding

    def test_the_strand_and_quality_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_5:
            assert key in finding

    def test_the_multi_file_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_6:
            if key != "conflict":
                assert key in finding

    def test_the_entity_extra_keys_exist(self):
        finding = _ensure_contract_keys({})
        for key in SECTION_2_7:
            assert key in finding

    def test_absent_data_is_never_a_missing_key(self):
        finding = _ensure_contract_keys({})
        assert finding["magnitude"] is None
        assert finding["repute"] == ""
        assert finding["topics"] == []
        assert finding["flipped"] is False

    def test_the_entity_type_defaults_to_snp(self):
        assert _ensure_contract_keys({})["entity_type"] == "snp"

    def test_the_confidence_defaults_to_none(self):
        assert _ensure_contract_keys({})["confidence"] == "none"

    def test_existing_values_are_not_overwritten(self):
        finding = _ensure_contract_keys({"gene": "MTHFR", "magnitude": 4.5,
                                         "entity_type": "genoset"})
        assert finding["gene"] == "MTHFR"
        assert finding["magnitude"] == 4.5
        assert finding["entity_type"] == "genoset"

    def test_it_returns_the_same_dict(self):
        finding = {"rsid": "rs1"}
        assert _ensure_contract_keys(finding) is finding

    def test_the_conditions_list_is_split_from_the_conditions_string(self):
        finding = _ensure_contract_keys({"conditions": "Alpha;Beta;Gamma"})
        assert finding["conditions_list"] == ["Alpha", "Beta", "Gamma"]

    def test_the_split_strips_whitespace(self):
        finding = _ensure_contract_keys({"conditions": "Alpha ; Beta"})
        assert finding["conditions_list"] == ["Alpha", "Beta"]

    def test_the_split_drops_empty_parts(self):
        finding = _ensure_contract_keys({"conditions": "Alpha;;Beta;"})
        assert finding["conditions_list"] == ["Alpha", "Beta"]

    def test_a_single_condition_becomes_a_one_item_list(self):
        finding = _ensure_contract_keys({"conditions": "Thrombophilia"})
        assert finding["conditions_list"] == ["Thrombophilia"]

    def test_an_existing_conditions_list_is_kept(self):
        finding = _ensure_contract_keys({"conditions": "Alpha;Beta",
                                        "conditions_list": ["Gamma"]})
        assert finding["conditions_list"] == ["Gamma"]

    def test_no_conditions_gives_an_empty_list(self):
        assert _ensure_contract_keys({})["conditions_list"] == []

    def test_carrier_is_true_for_one_copy(self):
        assert _ensure_contract_keys({"variant_copies": 1})["carrier"] is True

    def test_carrier_is_true_for_two_copies(self):
        assert _ensure_contract_keys({"variant_copies": 2})["carrier"] is True

    def test_carrier_is_false_for_zero_copies(self):
        assert _ensure_contract_keys({"variant_copies": 0})["carrier"] is False

    def test_carrier_stays_unknown_when_the_copy_number_is_unknown(self):
        assert _ensure_contract_keys({})["carrier"] is None

    def test_carrier_stays_unknown_for_a_null_copy_number(self):
        finding = _ensure_contract_keys({"variant_copies": None})
        assert finding["carrier"] is None

    def test_a_derived_carrier_flag_overrides_a_stale_one(self):
        finding = _ensure_contract_keys({"variant_copies": 0, "carrier": True})
        assert finding["carrier"] is False

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   _ensure_contract_keys does not default the section 2.4 frequency keys, nor
    #   conflict, nor rsid, so a genoset, trait or prs finding reaches /findings
    #   without them and section 2 promises every key is always present.
    def test_every_section_two_key_exists(self):
        finding = _ensure_contract_keys({})
        missing = [key for key in SECTION_2 if key not in finding]
        assert missing == []

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   The eleven section 2.4 frequency keys are only set by frequency.annotate,
    #   which enrich_findings calls for SNPs alone, so no genoset, trait or prs
    #   finding ever carries them.
    def test_the_frequency_contract_keys_exist(self):
        finding = _ensure_contract_keys({})
        assert [key for key in SECTION_2_4 if key not in finding] == []

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   conflict is documented in section 2.6 but is only set by
    #   attach_provenance, which run_full_scan applies to SNP findings alone.
    def test_the_conflict_key_exists(self):
        assert "conflict" in _ensure_contract_keys({})

    # REGRESSION GUARD. This documented a real defect, now fixed:
    #   _ensure_contract_keys tests for a None entity_type, gene or genotype and
    #   then calls setdefault, which cannot replace an existing None, so the
    #   coercion the branch was written for never happens and a None reaches the
    #   template layer.
    def test_a_null_entity_type_is_coerced_to_the_default(self):
        finding = _ensure_contract_keys({"entity_type": None, "gene": None,
                                         "genotype": None})
        assert finding["entity_type"] == "snp"
        assert finding["gene"] == ""
        assert finding["genotype"] == ""


class TestCollectGenosets:
    def test_the_documented_shape(self):
        require_subsystem("genosets")
        assert set(collect_genosets(GENOTYPES)) == {"matched", "unmatched",
                                                    "incomplete"}

    def test_matched_and_incomplete_are_disjoint(self):
        require_subsystem("genosets")
        result = collect_genosets(GENOTYPES)
        matched = {g["rsid"] for g in result["matched"]}
        incomplete = {g["rsid"] for g in result["incomplete"]}
        assert not (matched & incomplete)

    def test_matched_and_unmatched_are_disjoint(self):
        require_subsystem("genosets")
        result = collect_genosets(GENOTYPES)
        matched = {g["rsid"] for g in result["matched"]}
        unmatched = {g["rsid"] for g in result["unmatched"]}
        assert not (matched & unmatched)

    def test_a_match_on_the_available_positions_is_reported_as_matched(self):
        require_subsystem("genosets")
        result = collect_genosets(GENOTYPES)
        assert result["matched"]

    def test_every_matched_genoset_has_an_id(self):
        require_subsystem("genosets")
        for genoset in collect_genosets(GENOTYPES)["matched"]:
            assert genoset.get("rsid")

    def test_every_matched_genoset_carries_a_partial_coverage_flag(self):
        require_subsystem("genosets")
        for genoset in collect_genosets(GENOTYPES)["matched"]:
            assert isinstance(genoset["partial_coverage"], bool)

    def test_a_fully_covered_match_is_not_partial(self):
        require_subsystem("genosets")
        for genoset in collect_genosets(GENOTYPES)["matched"]:
            if genoset.get("coverage") == 1.0:
                assert genoset["partial_coverage"] is False

    def test_a_partially_covered_match_is_flagged(self):
        require_subsystem("genosets")
        for genoset in collect_genosets(GENOTYPES)["matched"]:
            coverage = genoset.get("coverage")
            if isinstance(coverage, (int, float)) and coverage < 1.0:
                assert genoset["partial_coverage"] is True

    def test_a_partial_match_is_never_also_incomplete(self):
        require_subsystem("genosets")
        result = collect_genosets(GENOTYPES)
        partial = {g["rsid"] for g in result["matched"]
                   if g.get("partial_coverage")}
        incomplete = {g["rsid"] for g in result["incomplete"]}
        assert not (partial & incomplete)

    def test_an_incomplete_genoset_is_not_reported_as_absent(self):
        require_subsystem("genosets")
        result = collect_genosets(GENOTYPES)
        unmatched = {g["rsid"] for g in result["unmatched"]}
        incomplete = {g["rsid"] for g in result["incomplete"]}
        assert not (unmatched & incomplete)

    def test_the_empty_shape_when_the_corpus_is_missing(self):
        saved = pipeline_module._genosets
        try:
            pipeline_module._genosets = None
            assert collect_genosets(GENOTYPES) == {"matched": [],
                                                   "unmatched": [],
                                                   "incomplete": []}
        finally:
            pipeline_module._genosets = saved

    def test_an_empty_genotype_map_does_not_raise(self):
        assert set(collect_genosets({})) == {"matched", "unmatched",
                                             "incomplete"}

    def test_a_none_genotype_map_does_not_raise(self):
        assert set(collect_genosets(None)) == {"matched", "unmatched",
                                               "incomplete"}


class TestCollectTraits:
    def test_the_documented_shape(self):
        require_subsystem("traits")
        assert set(collect_traits(GENOTYPES)) == {"traits", "blood_type",
                                                  "findings"}

    def test_traits_are_called(self):
        require_subsystem("traits")
        assert collect_traits(GENOTYPES)["findings"]

    def test_every_finding_has_a_blank_repute(self):
        require_subsystem("traits")
        for finding in collect_traits(GENOTYPES)["findings"]:
            assert finding["repute"] == ""

    def test_every_finding_is_a_trait_entity(self):
        require_subsystem("traits")
        for finding in collect_traits(GENOTYPES)["findings"]:
            assert finding["entity_type"] == "trait"

    def test_the_blood_type_is_a_mapping(self):
        require_subsystem("traits")
        assert isinstance(collect_traits(GENOTYPES)["blood_type"], dict)

    def test_the_empty_shape_when_the_module_is_missing(self):
        saved = pipeline_module._traits
        try:
            pipeline_module._traits = None
            assert collect_traits(GENOTYPES) == {"traits": [],
                                                 "blood_type": {},
                                                 "findings": []}
        finally:
            pipeline_module._traits = saved

    def test_an_empty_genotype_map_does_not_raise(self):
        assert set(collect_traits({})) == {"traits", "blood_type", "findings"}

    def test_a_trait_repute_is_forced_even_if_the_engine_sets_one(self):
        require_subsystem("traits")
        findings = collect_traits(GENOTYPES)["findings"]
        assert not any(finding.get("repute") for finding in findings)


class TestCollectPrs:
    def test_the_documented_shape(self):
        require_subsystem("prs")
        assert set(collect_prs(GENOTYPES)) == {"results", "findings"}

    def test_models_are_computed(self):
        require_subsystem("prs")
        assert collect_prs(GENOTYPES)["results"]

    def test_every_finding_has_a_blank_repute(self):
        require_subsystem("prs")
        for finding in collect_prs(GENOTYPES)["findings"]:
            assert finding["repute"] == ""

    def test_every_finding_is_a_prs_entity(self):
        require_subsystem("prs")
        for finding in collect_prs(GENOTYPES)["findings"]:
            assert finding["entity_type"] == "prs"

    def test_the_results_are_a_list(self):
        require_subsystem("prs")
        assert isinstance(collect_prs(GENOTYPES)["results"], list)

    def test_the_empty_shape_when_the_module_is_missing(self):
        saved = pipeline_module._prs
        try:
            pipeline_module._prs = None
            assert collect_prs(GENOTYPES) == {"results": [], "findings": []}
        finally:
            pipeline_module._prs = saved

    def test_an_empty_genotype_map_does_not_raise(self):
        assert set(collect_prs({})) == {"results", "findings"}

    def test_a_polygenic_score_never_carries_a_repute(self):
        require_subsystem("prs")
        findings = collect_prs(GENOTYPES)["findings"]
        assert not any(finding.get("repute") for finding in findings)


QC_FINDINGS = [
    {"rsid": "rs1", "entity_type": "snp", "flipped": True,
     "stabilized_orientation": "minus"},
    {"rsid": "rs2", "entity_type": "snp", "ambiguous": True,
     "stabilized_orientation": "plus"},
    {"rsid": "rs3", "entity_type": "snp", "freq_ambiguous": True,
     "stabilized_orientation": "plus"},
    {"rsid": "rs4", "entity_type": "snp", "zygosity": "no_call",
     "stabilized_orientation": ""},
    {"rsid": "rs5", "entity_type": "snp", "carrier": False,
     "stabilized_orientation": "plus"},
    {"rsid": "rs6", "entity_type": "snp", "carrier": None, "dubious": True,
     "stabilized_orientation": "plus"},
    {"rsid": "dgs001", "entity_type": "genoset", "flipped": True,
     "ambiguous": True},
    {"rsid": "lactase", "entity_type": "trait", "zygosity": "no_call"},
    {"rsid": "t2d", "entity_type": "prs", "carrier": False},
]

QC_MERGED = {"counts": {"conflicts": 3, "pooled_sources": 2,
                        "total_positions": 99}}

RANGE_FINDINGS = [
    {"rsid": "rs1", "magnitude": 7.777, "publications": 412, "freq": 99.44},
    {"rsid": "rs2", "magnitude": 1.0, "publications": 0, "freq": 0.0},
    {"rsid": "rs3", "magnitude": None, "publications": None, "freq": None},
]


class TestComputeQc:
    def test_only_snp_entities_are_counted_in_the_total(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["total"] == 6

    def test_flipped_calls_are_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["flipped"] == 1

    def test_ambiguous_calls_are_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["ambiguous"] == 2

    def test_no_calls_are_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["no_call"] == 1

    def test_non_carriers_are_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["non_carrier"] == 1

    def test_an_unknown_carrier_state_is_not_a_non_carrier(self):
        rows = [{"entity_type": "snp", "carrier": None}]
        assert compute_qc(rows, {})["non_carrier"] == 0

    def test_dubious_calls_are_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["dubious"] == 1

    def test_unknown_orientation_is_counted(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["unknown_orientation"] == 1

    def test_conflicts_come_from_the_merge_counts(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["conflicts"] == 3

    def test_pooled_sources_come_from_the_merge_counts(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["pooled_sources"] == 2

    def test_positions_come_from_the_merge_counts(self):
        assert compute_qc(QC_FINDINGS, QC_MERGED)["positions"] == 99

    def test_a_genoset_is_not_counted_as_a_flipped_snp(self):
        rows = [{"entity_type": "genoset", "flipped": True}]
        assert compute_qc(rows, {})["flipped"] == 0

    def test_a_trait_no_call_is_not_counted(self):
        rows = [{"entity_type": "trait", "zygosity": "no_call"}]
        assert compute_qc(rows, {})["no_call"] == 0

    def test_the_note_explains_why_palindromes_are_capped(self):
        note = compute_qc(QC_FINDINGS, QC_MERGED)["note"]
        assert "strand" in note
        assert "A/T or C/G" in note

    def test_an_empty_finding_set_gives_zero_counts(self):
        qc = compute_qc([], {})
        assert qc["total"] == 0
        assert qc["flipped"] == 0
        assert qc["conflicts"] == 0

    def test_a_merged_set_with_no_counts_is_safe(self):
        assert compute_qc(QC_FINDINGS, {})["conflicts"] == 0

    def test_every_count_is_an_integer(self):
        qc = compute_qc(QC_FINDINGS, QC_MERGED)
        for key, value in qc.items():
            if key != "note":
                assert isinstance(value, int)


class TestComputeRanges:
    def test_the_three_documented_ranges(self):
        assert set(compute_ranges(RANGE_FINDINGS)) == {"magnitude",
                                                      "publications",
                                                      "frequency"}

    def test_the_magnitude_bound_is_data_derived(self):
        assert compute_ranges(RANGE_FINDINGS)["magnitude"] == [0.0, 7.78]

    def test_the_publication_bound_is_data_derived(self):
        assert compute_ranges(RANGE_FINDINGS)["publications"] == [0, 412]

    def test_the_frequency_bound_is_data_derived(self):
        assert compute_ranges(RANGE_FINDINGS)["frequency"] == [0.0, 99.44]

    def test_a_null_magnitude_does_not_become_a_bound(self):
        rows = [{"magnitude": None}, {"magnitude": 3.0}]
        assert compute_ranges(rows)["magnitude"] == [0.0, 3.0]

    def test_a_null_publication_count_does_not_become_a_bound(self):
        rows = [{"publications": None}, {"publications": 7}]
        assert compute_ranges(rows)["publications"] == [0, 7]

    def test_a_null_frequency_does_not_become_a_bound(self):
        rows = [{"freq": None}, {"freq": 12.5}]
        assert compute_ranges(rows)["frequency"] == [0.0, 12.5]

    def test_the_lower_bounds_are_always_zero(self):
        ranges = compute_ranges(RANGE_FINDINGS)
        assert ranges["magnitude"][0] == 0.0
        assert ranges["publications"][0] == 0
        assert ranges["frequency"][0] == 0.0

    def test_the_bounds_are_rounded_to_two_places(self):
        rows = [{"magnitude": 6.66666, "freq": 12.34567}]
        ranges = compute_ranges(rows)
        assert ranges["magnitude"][1] == 6.67
        assert ranges["frequency"][1] == 12.35

    def test_sane_defaults_on_an_empty_set(self):
        ranges = compute_ranges([])
        assert ranges["magnitude"] == [0.0, 10.0]
        assert ranges["publications"] == [0, 0]
        assert ranges["frequency"] == [0.0, 100.0]

    def test_a_set_with_no_usable_values_falls_back_to_the_defaults(self):
        rows = [{"rsid": "rs1"}, {"rsid": "rs2"}]
        assert compute_ranges(rows)["magnitude"] == [0.0, 10.0]

    def test_an_observed_zero_frequency_is_a_real_value(self):
        rows = [{"freq": 0.0}]
        assert compute_ranges(rows)["frequency"] == [0.0, 0.0]


class TestComputeSummary:
    def test_the_documented_shape(self):
        assert set(compute_summary(QC_FINDINGS)) == {"silos", "entity_types",
                                                     "reputes"}

    def test_entity_types_are_counted(self):
        entities = compute_summary(QC_FINDINGS)["entity_types"]
        assert entities["snp"] == 6
        assert entities["genoset"] == 1
        assert entities["trait"] == 1
        assert entities["prs"] == 1

    def test_a_missing_silo_counts_as_informational(self):
        summary = compute_summary([{"rsid": "rs1"}])
        assert summary["silos"] == {"informational": 1}

    def test_silos_are_counted(self):
        rows = [{"silo": "actionable"}, {"silo": "actionable"},
                {"silo": "pre_prescription"}]
        assert compute_summary(rows)["silos"] == {"actionable": 2,
                                                 "pre_prescription": 1}

    def test_the_three_repute_keys_are_always_present(self):
        assert set(compute_summary([])["reputes"]) == {"Good", "Bad", "unset"}

    def test_reputes_are_counted(self):
        rows = [{"repute": "Good"}, {"repute": "Bad"}, {"repute": ""}]
        assert compute_summary(rows)["reputes"] == {"Good": 1, "Bad": 1,
                                                    "unset": 1}

    def test_an_unexpected_repute_counts_as_unset(self):
        rows = [{"repute": "Sideways"}]
        assert compute_summary(rows)["reputes"]["unset"] == 1

    def test_an_empty_set_gives_empty_counts(self):
        summary = compute_summary([])
        assert summary["silos"] == {}
        assert summary["entity_types"] == {}


DOCUMENTED_PAYLOAD_KEYS = ("findings", "summary", "ranges", "qc", "genosets",
                           "traits", "blood_type", "prs", "conflicts",
                           "sources", "counts", "trio", "population",
                           "subsystems")

# The v1.2 keys an existing consumer still expects on every finding.
CARRIED_OVER_KEYS = ("rsid", "entity_type", "gene", "chromosome", "position",
                     "allele1", "allele2", "genotype", "token", "zygosity",
                     "magnitude", "magnitude_source", "repute", "summary",
                     "interpretation", "confidence", "clinical_sig",
                     "clinvar_sig_code", "review_status", "review_stars",
                     "cpic_level", "pgx_level", "evidence", "publications",
                     "conditions", "conditions_list", "sources",
                     "orientation", "stabilized_orientation", "flipped",
                     "ambiguous", "dubious", "variant_allele",
                     "variant_copies", "carrier", "count", "labels", "calls",
                     "comparison", "topics", "medicines", "criteria",
                     "matched_rsids", "coverage")


class TestRunFullScan:
    def test_the_documented_top_level_keys(self):
        payload = scan()
        assert set(payload) == set(DOCUMENTED_PAYLOAD_KEYS)

    def test_findings_are_produced(self):
        assert len(scan()["findings"]) > 0

    def test_every_finding_carries_the_v1_keys(self):
        for finding in scan()["findings"]:
            for key in CARRIED_OVER_KEYS:
                assert key in finding

    def test_every_finding_is_scored(self):
        for finding in scan()["findings"]:
            assert "magnitude" in finding
            assert "confidence" in finding
            assert "repute" in finding

    def test_every_finding_records_where_its_magnitude_came_from(self):
        for finding in scan()["findings"]:
            assert finding["magnitude_source"] in ("computed", "snpedia")

    def test_every_finding_has_an_audit_trail(self):
        for finding in scan()["findings"]:
            assert finding["magnitude_factors"]

    def test_no_trait_carries_a_repute(self):
        for finding in scan()["findings"]:
            if finding["entity_type"] == "trait":
                assert finding["repute"] == ""

    def test_no_polygenic_score_carries_a_repute(self):
        for finding in scan()["findings"]:
            if finding["entity_type"] == "prs":
                assert finding["repute"] == ""

    def test_no_no_call_scores_above_zero(self):
        for finding in scan()["findings"]:
            if finding.get("zygosity") == "no_call":
                assert not (finding["magnitude"] or 0) > 0

    def test_every_magnitude_is_inside_the_documented_range(self):
        for finding in scan()["findings"]:
            magnitude = finding["magnitude"]
            if isinstance(magnitude, (int, float)):
                assert 0.0 <= magnitude <= 10.0

    def test_every_repute_is_one_of_the_three_documented_values(self):
        for finding in scan()["findings"]:
            assert finding["repute"] in ("Good", "Bad", "")

    def test_every_confidence_is_one_of_the_four_documented_words(self):
        for finding in scan()["findings"]:
            assert finding["confidence"] in ("high", "moderate", "low", "none")

    def test_the_deliberate_conflict_survives_into_conflicts(self):
        rsids = {row["rsid"] for row in scan()["conflicts"]}
        assert CONFLICT_RSID in rsids

    def test_only_the_deliberate_conflict_is_reported(self):
        assert len(scan()["conflicts"]) == 1

    def test_the_conflict_retains_both_calls(self):
        row = next(r for r in scan()["conflicts"]
                   if r["rsid"] == CONFLICT_RSID)
        genotypes = {call["genotype"] for call in row["calls"]}
        assert genotypes == {"CT", "CC"}

    def test_the_conflicting_finding_is_flagged(self):
        finding = next(f for f in scan()["findings"]
                       if f["rsid"] == CONFLICT_RSID)
        assert finding["conflict"] is True

    def test_the_conflicting_finding_names_both_source_files(self):
        finding = next(f for f in scan()["findings"]
                       if f["rsid"] == CONFLICT_RSID)
        assert set(finding["labels"]) == {"self_a.txt", "self_b.txt"}

    def test_the_progress_callback_reports_phase_names(self):
        assert scan_phases()

    def test_every_reported_phase_is_a_documented_phase(self):
        for phase in scan_phases():
            assert phase in PHASES

    def test_the_pooling_phase_is_reported(self):
        assert "merge" in scan_phases()

    def test_the_offline_annotation_phase_is_reported(self):
        assert "bundled" in scan_phases()

    def test_the_orientation_phase_is_reported(self):
        assert "orientation" in scan_phases()

    def test_the_scoring_phase_is_reported(self):
        assert "scoring" in scan_phases()

    def test_the_final_phase_is_complete(self):
        assert scan_phases()[-1] == "complete"

    def test_the_network_phase_is_not_reported_for_an_offline_scan(self):
        assert "api" not in scan_phases()

    def test_the_qc_block_counts_only_snp_entities(self):
        payload = scan()
        snps = sum(1 for f in payload["findings"]
                   if f["entity_type"] == "snp")
        assert payload["qc"]["total"] == snps

    def test_the_qc_conflict_count_matches_the_conflict_list(self):
        payload = scan()
        assert payload["qc"]["conflicts"] == len(payload["conflicts"])

    def test_the_qc_block_counts_the_failed_probe(self):
        assert scan()["qc"]["no_call"] >= 1

    def test_the_summary_counts_every_finding(self):
        payload = scan()
        counted = sum(payload["summary"]["entity_types"].values())
        assert counted == len(payload["findings"])

    def test_the_ranges_are_derived_from_the_findings(self):
        payload = scan()
        magnitudes = [f["magnitude"] for f in payload["findings"]
                      if isinstance(f["magnitude"], (int, float))]
        assert payload["ranges"]["magnitude"][1] == round(max(magnitudes), 2)

    def test_the_population_is_echoed(self):
        assert scan()["population"] == "CEU"

    def test_both_self_files_are_pooled(self):
        assert scan()["counts"]["pooled_sources"] == 2

    def test_every_source_is_described(self):
        assert len(scan()["sources"]) == 3

    def test_the_source_roles_are_reported(self):
        roles = [source["role"] for source in scan()["sources"]]
        assert roles.count("self") == 2
        assert "mother" in roles

    def test_the_mother_appears_as_a_comparison_row(self):
        finding = next(f for f in scan()["findings"]
                       if f["rsid"] == CONFLICT_RSID)
        roles = {row["role"] for row in finding["comparison"]}
        assert "mother" in roles

    def test_a_relative_never_becomes_a_pooled_call(self):
        finding = next(f for f in scan()["findings"]
                       if f["rsid"] == CONFLICT_RSID)
        assert "mum.txt" not in finding["labels"]

    def test_the_trio_is_unavailable_without_a_father(self):
        assert scan()["trio"]["trio_available"] is False

    def test_the_subsystem_report_is_included(self):
        assert set(scan()["subsystems"]) == set(available_subsystems())

    def test_the_failed_probe_stays_a_no_call(self):
        no_calls = [f["rsid"] for f in scan()["findings"]
                    if f.get("zygosity") == "no_call"]
        assert "rs4680" in no_calls

    def test_the_extra_entity_types_are_appended_to_the_same_list(self):
        payload = scan()
        present = {f["entity_type"] for f in payload["findings"]}
        assert "snp" in present
        for name, entity in (("genosets", "genoset"), ("traits", "trait"),
                             ("prs", "prs")):
            if SUBSYSTEMS.get(name):
                assert entity in present

    def test_the_genoset_block_has_the_documented_shape(self):
        assert set(scan()["genosets"]) == {"matched", "unmatched", "incomplete"}

    def test_the_blood_type_block_is_a_mapping(self):
        assert isinstance(scan()["blood_type"], dict)

    def test_switching_off_an_entity_type_removes_it(self):
        require_reference()
        payload = run_full_scan(SOURCES, use_api=False, include_genosets=False,
                                include_traits=False, include_prs=False)
        present = {f["entity_type"] for f in payload["findings"]}
        assert present == {"snp"}
        assert payload["genosets"]["matched"] == []
        assert payload["traits"] == []
        assert payload["prs"] == []


class TestEnrichFindings:
    def test_it_returns_the_same_list(self):
        findings = [{"rsid": "rs1801133", "allele1": "C", "allele2": "T",
                     "genotype": "CT", "entity_type": "snp"}]
        assert enrich_findings(findings, merged_fixture()) is findings

    def test_it_fills_the_contract_keys(self):
        findings = [{"rsid": "rs1801133", "allele1": "C", "allele2": "T",
                     "genotype": "CT", "entity_type": "snp"}]
        enrich_findings(findings, merged_fixture())
        for key in SECTION_2_1 + SECTION_2_2:
            assert key in findings[0]

    def test_it_attaches_provenance(self):
        findings = [{"rsid": "rs1801133", "allele1": "C", "allele2": "T",
                     "genotype": "CT", "entity_type": "snp"}]
        enrich_findings(findings, merged_fixture())
        assert findings[0]["count"] == 2
        assert findings[0]["conflict"] is True

    def test_it_resolves_strand_before_scoring(self):
        findings = [{"rsid": "rs1", "allele1": "A", "allele2": "T",
                     "genotype": "AT", "entity_type": "snp",
                     "cpic_level": "A"}]
        enrich_findings(findings, {"genotypes": {}})
        # The palindrome was detected first, so the cap reached the score.
        assert findings[0]["ambiguous"] is True
        assert findings[0]["magnitude"] <= 2.0
        assert findings[0]["dubious"] is True

    def test_it_scores_every_finding(self):
        findings = [{"rsid": "rs1", "entity_type": "snp"},
                    {"rsid": "dgs001", "entity_type": "genoset"}]
        enrich_findings(findings, {"genotypes": {}})
        assert all(f["magnitude_source"] == "computed" for f in findings)

    def test_a_no_call_is_scored_at_zero(self):
        findings = [{"rsid": "rs1", "entity_type": "snp", "cpic_level": "A",
                     "zygosity": "no_call", "allele1": "N", "allele2": "N"}]
        enrich_findings(findings, {"genotypes": {}})
        assert findings[0]["magnitude"] == 0.0

    def test_a_progress_callback_is_invoked(self):
        seen = []
        findings = [{"rsid": "rs1", "entity_type": "snp"}]
        enrich_findings(findings, {"genotypes": {}},
                        progress_cb=lambda phase, done=0, total=0:
                        seen.append(phase))
        assert seen == ["scoring"]

    def test_an_empty_finding_list_is_safe(self):
        assert enrich_findings([], {"genotypes": {}}) == []


class TestAvailableSubsystems:
    def test_the_documented_subsystems_are_reported(self):
        for name in ("frequency", "genosets", "traits", "prs", "snpedia"):
            assert name in available_subsystems()

    def test_every_value_is_a_bool(self):
        for value in available_subsystems().values():
            assert isinstance(value, bool)

    def test_it_reports_false_when_every_module_is_absent(self):
        saved = (pipeline_module._genosets, pipeline_module._traits,
                 pipeline_module._prs, pipeline_module._snpedia)
        try:
            pipeline_module._genosets = None
            pipeline_module._traits = None
            pipeline_module._prs = None
            pipeline_module._snpedia = None
            result = available_subsystems()
            assert result["genosets"] is False
            assert result["traits"] is False
            assert result["prs"] is False
            assert result["snpedia"] is False
        finally:
            (pipeline_module._genosets, pipeline_module._traits,
             pipeline_module._prs, pipeline_module._snpedia) = saved

    def test_it_does_not_raise_when_a_loader_throws(self):
        class Exploding:
            def load_genosets(self):
                raise RuntimeError("the corpus is corrupt")

            def load_models(self):
                raise RuntimeError("the corpus is corrupt")

            def cache_status(self):
                raise RuntimeError("the cache is locked")

        saved = (pipeline_module._genosets, pipeline_module._prs,
                 pipeline_module._snpedia)
        try:
            stub = Exploding()
            pipeline_module._genosets = stub
            pipeline_module._prs = stub
            pipeline_module._snpedia = stub
            result = available_subsystems()
            assert result["genosets"] is False
            assert result["prs"] is False
            assert result["snpedia"] is False
        finally:
            (pipeline_module._genosets, pipeline_module._prs,
             pipeline_module._snpedia) = saved

    def test_a_traits_module_with_no_corpus_reports_false(self):
        class Bare:
            TRAITS = []

        saved = pipeline_module._traits
        try:
            pipeline_module._traits = Bare()
            assert available_subsystems()["traits"] is False
        finally:
            pipeline_module._traits = saved

    def test_the_report_still_carries_the_five_v2_pipeline_subsystems(self):
        # v3.0 widened this map. The original assertion was an exact set
        # equality, which was the right test while the five were the whole
        # story and the wrong test the moment a sixth subsystem existed.
        #
        # What actually needs protecting is that the five v2 keys never
        # disappear or change name, because the frontend and every v2 client
        # read them by name. So: subset, not equality, plus an explicit check
        # that each one is still a boolean rather than a richer object.
        report = available_subsystems()
        v2_keys = {"frequency", "genosets", "traits", "prs", "snpedia"}
        assert v2_keys <= set(report)
        for key in v2_keys:
            assert isinstance(report[key], bool)

    def test_v3_subsystem_flags_are_reported_as_booleans(self):
        report = available_subsystems()
        for key in ("ledger", "provenance", "sequencing", "haplogroups",
                    "relatedness", "imputation", "ancestry", "diplotype",
                    "carrier", "assistant"):
            assert key in report
            assert isinstance(report[key], bool)

    def test_external_capabilities_are_separate_keys_from_subsystems(self):
        # A subsystem flag means DNAInsight's own code is present. A capability
        # flag means a third-party tool the user installed is present and its
        # licence is accepted. Collapsing the two would make an uninstalled
        # Beagle look like a broken DNAInsight module.
        report = available_subsystems()
        assert report.get("imputation") is True
        assert report.get("ancestry") is True
        # Tool capabilities carry a tool_ prefix. Beagle's capability is named
        # "imputation" and Ollama's is "assistant", which collide head-on with
        # the module names above. Without the prefix an uninstalled Beagle
        # overwrites the module flag and the UI reports DNAInsight's own code
        # as missing.
        assert report.get("tool_imputation") is False
        assert report.get("tool_assistant") is False
        assert "tool_ancestry_global" in report
        assert "tool_haplogroup_y" in report

    def test_every_tool_capability_appears_only_under_the_tool_prefix(self):
        # The invariant, stated precisely: a capability name may coincide with a
        # module name, and two of them do, so the guarantee cannot be "these
        # names never overlap". It has to be "the tool value is only ever
        # written under tool_<name>". Then a shared name is harmless because the
        # two values live in two keys.
        from backend import external as external_mod
        report = available_subsystems()
        capabilities = external_mod.capability_report()
        for name, ready in capabilities.items():
            assert f"tool_{name}" in report
            assert report[f"tool_{name}"] == ready

    def test_a_shared_name_keeps_the_module_value_unprefixed(self):
        # Beagle's capability and DNAInsight's module are both called
        # "imputation". The bare key must report the module, which is present,
        # and never the tool, which is not installed on a normal machine.
        from backend import external as external_mod
        report = available_subsystems()
        shared = set(external_mod.capability_report()) & {
            "ledger", "provenance", "sequencing", "haplogroups", "relatedness",
            "imputation", "ancestry", "diplotype", "carrier", "assistant",
        }
        assert shared, "Expected at least one name to be shared, or this guard is dead code."
        for name in shared:
            assert report[name] is True
            assert report[f"tool_{name}"] is False
