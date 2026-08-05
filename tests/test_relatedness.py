"""Tests for backend.relatedness: IBD segments, cM honesty, phasing and household checks."""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import external
from backend import merge
from backend import relatedness as rel


AUTOSOMES = [str(c) for c in range(1, 23)]


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """No installed tools, no installed genetic map, no cached state."""
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    external.reset_cache()
    rel.reset_map_cache()
    yield
    external.reset_cache()
    rel.reset_map_cache()


# ---------------------------------------------------------------------------
# Synthetic data builders
# ---------------------------------------------------------------------------

def make_markers(n, chrom="1", start=1_000_000, step=100_000):
    return [(chrom, start + i * step, f"{chrom}_rs{i + 1}") for i in range(n)]


def positions_of(markers):
    return {rsid: (chrom, pos) for chrom, pos, rsid in markers}


def genotypes_of(markers, pair):
    return {rsid: pair for _chrom, _pos, rsid in markers}


def snp_rows(markers, pair_for):
    rows = []
    for chrom, pos, rsid in markers:
        a1, a2 = pair_for(chrom, pos, rsid)
        rows.append({"rsid": rsid, "chromosome": chrom, "position": pos,
                     "allele1": a1, "allele2": a2})
    return rows


def family_pair(fraction, *, shared_snps=301, tail_snps=40):
    """Two kits sharing the first ``fraction`` of every autosome.

    The expected total is computed here from the module's own rate table rather
    than written down, so the fixture cannot silently drift away from the code
    it is checking.
    """
    left, right, expected_cm = [], [], 0.0
    for chrom in AUTOSOMES:
        length = rel.CHROM_LENGTHS_GRCH37[chrom]
        end = int(length * fraction)
        step = max(1, (end - 1_000_000) // (shared_snps - 1))
        for i in range(shared_snps):
            pos = 1_000_000 + i * step
            rsid = f"{chrom}_s{i}"
            row = {"rsid": rsid, "chromosome": chrom, "position": pos,
                   "allele1": "A", "allele2": "G"}
            left.append(dict(row))
            right.append(dict(row))
        span_bp = (shared_snps - 1) * step
        expected_cm += (span_bp / 1_000_000.0) * rel.AVERAGE_CM_PER_MB[chrom]
        for i in range(tail_snps):
            pos = end + 2_000_000 + i * 100_000
            rsid = f"{chrom}_t{i}"
            left.append({"rsid": rsid, "chromosome": chrom, "position": pos,
                         "allele1": "A", "allele2": "A"})
            right.append({"rsid": rsid, "chromosome": chrom, "position": pos,
                          "allele1": "G", "allele2": "G"})
    return left, right, expected_cm


# ---------------------------------------------------------------------------
# IBS
# ---------------------------------------------------------------------------

def test_identical_genotypes_are_ibs_two():
    assert rel.ibs_state(("A", "G"), ("G", "A")) == 2


def test_one_shared_allele_is_ibs_one():
    assert rel.ibs_state(("A", "A"), ("A", "G")) == 1


def test_opposite_homozygotes_are_ibs_zero():
    assert rel.ibs_state(("A", "A"), ("G", "G")) == 0


def test_a_no_call_on_either_side_gives_no_ibs_state():
    assert rel.ibs_state(("N", "N"), ("A", "G")) is None
    assert rel.ibs_state(("A", "G"), ("A", "N")) is None
    assert rel.ibs_state(None, ("A", "G")) is None


def test_ibs_accepts_the_merged_genotype_dict_shape():
    left = {"allele1": "A", "allele2": "G"}
    right = {"allele1": "G", "allele2": "A"}
    assert rel.ibs_state(left, right) == 2


def test_ibs_accepts_a_two_character_string():
    assert rel.ibs_state("AG", "AG") == 2


# ---------------------------------------------------------------------------
# Segment detection
# ---------------------------------------------------------------------------

def test_identical_samples_produce_one_segment_with_a_known_snp_count():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert result["segment_count"] == 1
    segment = result["segments"][0]
    assert segment["snps"] == 500
    assert segment["start_rsid"] == "1_rs1"
    assert segment["end_rsid"] == "1_rs500"


def test_segment_coordinates_match_the_first_and_last_marker():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    segment = result["segments"][0]
    assert segment["start_bp"] == markers[0][1]
    assert segment["end_bp"] == markers[-1][1]


def test_the_builtin_walk_reports_itself_as_approximate():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert result["method"] == "builtin_ibs"
    assert result["approximate"] is True
    assert any("identity by state" in c for c in result["caveats"])


def test_two_opposite_homozygotes_break_the_run_at_a_known_place():
    markers = make_markers(500)
    left = genotypes_of(markers, ("A", "G"))
    right = genotypes_of(markers, ("A", "G"))
    for rsid in ("1_rs200", "1_rs201", "1_rs202"):
        left[rsid] = ("A", "A")
        right[rsid] = ("G", "G")
    result = rel.shared_segments(left, right, positions=positions_of(markers),
                                 min_snps=50, min_cm=1.0)
    spans = [(s["start_rsid"], s["end_rsid"], s["snps"]) for s in result["segments"]]
    assert spans == [("1_rs1", "1_rs199", 199), ("1_rs203", "1_rs500", 298)]


def test_a_single_opposite_homozygote_is_tolerated_by_default():
    markers = make_markers(500)
    left = genotypes_of(markers, ("A", "G"))
    right = genotypes_of(markers, ("A", "G"))
    left["1_rs250"] = ("A", "A")
    right["1_rs250"] = ("G", "G")
    result = rel.shared_segments(left, right, positions=positions_of(markers))
    assert result["segment_count"] == 1
    assert result["segments"][0]["opposite_homozygotes"] == 1


def test_zero_tolerance_breaks_on_the_first_opposite_homozygote():
    markers = make_markers(500)
    left = genotypes_of(markers, ("A", "G"))
    right = genotypes_of(markers, ("A", "G"))
    left["1_rs250"] = ("A", "A")
    right["1_rs250"] = ("G", "G")
    result = rel.shared_segments(left, right, positions=positions_of(markers),
                                 min_snps=50, min_cm=1.0,
                                 max_opposite_homozygotes=0)
    assert result["segment_count"] == 2
    assert result["segments"][0]["end_rsid"] == "1_rs249"


def test_the_snp_floor_filters_a_short_run():
    markers = make_markers(100)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers), min_cm=0.0)
    assert result["segments"] == []


def test_the_cm_floor_filters_a_physically_short_run():
    markers = make_markers(400, step=1000)      # 400 kb, well under 7 cM
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers), min_snps=1)
    assert result["segments"] == []


def test_segments_are_found_on_each_chromosome_separately():
    markers = make_markers(400, chrom="1") + make_markers(400, chrom="2")
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert [s["chromosome"] for s in result["segments"]] == ["1", "2"]


def test_segments_are_returned_in_chromosome_and_coordinate_order():
    markers = make_markers(400, chrom="10") + make_markers(400, chrom="2")
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert [s["chromosome"] for s in result["segments"]] == ["2", "10"]


def test_markers_without_coordinates_are_skipped_rather_than_guessed():
    markers = make_markers(400)
    genotypes = genotypes_of(markers, ("A", "G"))
    genotypes["no_coords"] = ("A", "G")
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert result["compared_snps"] == 400


def test_no_call_positions_are_excluded_from_the_comparison():
    markers = make_markers(400)
    left = genotypes_of(markers, ("A", "G"))
    right = genotypes_of(markers, ("A", "G"))
    right["1_rs5"] = ("N", "N")
    result = rel.shared_segments(left, right, positions=positions_of(markers))
    assert result["compared_snps"] == 399


def test_a_chromosome_filter_limits_the_walk():
    markers = make_markers(400, chrom="1") + make_markers(400, chrom="2")
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers),
                                 chromosomes=["2"])
    assert result["chromosomes"] == ["2"]


def test_the_payload_says_comparison_is_limited_to_loaded_kits():
    result = rel.shared_segments({}, {}, positions={})
    assert "nothing was uploaded" in result["note"].lower()


# ---------------------------------------------------------------------------
# Centimorgans
# ---------------------------------------------------------------------------

def test_cm_without_a_map_is_flagged_as_estimated():
    result = rel.centimorgans("1", 1_000_000, 51_000_000)
    assert result["cm_estimated"] is True
    assert result["source"] == "average_rate"
    assert "not a measurement" in result["note"]


def test_cm_without_a_map_uses_the_per_chromosome_rate():
    result = rel.centimorgans("19", 0, 10_000_000)
    assert result["rate_cm_per_mb"] == rel.AVERAGE_CM_PER_MB["19"]
    assert result["cm"] == pytest.approx(10 * rel.AVERAGE_CM_PER_MB["19"])


def test_an_unknown_chromosome_falls_back_to_the_genome_average():
    result = rel.centimorgans("ZZ", 0, 10_000_000)
    assert result["rate_cm_per_mb"] == rel.GENOME_AVERAGE_CM_PER_MB


def test_cm_from_a_supplied_map_is_not_flagged_as_estimated():
    genetic_map = {"1": [(0, 0.0), (10_000_000, 10.0), (20_000_000, 30.0)]}
    result = rel.centimorgans("1", 0, 20_000_000, genetic_map)
    assert result["cm_estimated"] is False
    assert result["source"] == "genetic_map"
    assert result["cm"] == pytest.approx(30.0)


def test_cm_from_a_map_interpolates_between_bracketing_points():
    genetic_map = {"1": [(0, 0.0), (10_000_000, 10.0), (20_000_000, 30.0)]}
    result = rel.centimorgans("1", 0, 15_000_000, genetic_map)
    assert result["cm"] == pytest.approx(20.0)


def test_cm_clamps_beyond_the_ends_of_the_map():
    genetic_map = {"1": [(1_000_000, 5.0), (2_000_000, 6.0)]}
    result = rel.centimorgans("1", 0, 9_000_000, genetic_map)
    assert result["cm"] == pytest.approx(1.0)


def test_cm_handles_reversed_coordinates():
    forward = rel.centimorgans("1", 1_000_000, 5_000_000)
    reverse = rel.centimorgans("1", 5_000_000, 1_000_000)
    assert forward["cm"] == reverse["cm"]


def test_load_genetic_map_returns_none_when_no_map_is_installed():
    assert rel.load_genetic_map("1") is None


def test_load_genetic_map_reads_a_plink_map_file(tmp_path):
    (tmp_path / "chr1.map").write_text(
        "1\trs1\t0.5\t100000\n1\trs2\t1.5\t200000\nbroken line\n", encoding="utf-8")
    points = rel.load_genetic_map("1", root=tmp_path)
    assert points == [(100000, 0.5), (200000, 1.5)]


def test_estimated_cm_is_propagated_onto_every_segment():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.shared_segments(genotypes, dict(genotypes),
                                 positions=positions_of(markers))
    assert result["cm_estimated"] is True
    assert all(s["cm_estimated"] for s in result["segments"])


# ---------------------------------------------------------------------------
# Totals and relationship prediction
# ---------------------------------------------------------------------------

def test_total_and_longest_accept_a_list_or_a_payload():
    segments = [{"cm": 10.0}, {"cm": 25.5}]
    assert rel.total_shared_cm(segments) == 35.5
    assert rel.longest_segment_cm(segments) == 25.5
    assert rel.total_shared_cm({"segments": segments}) == 35.5
    assert rel.longest_segment_cm({"segments": segments}) == 25.5


def test_longest_segment_of_nothing_is_zero():
    assert rel.longest_segment_cm([]) == 0.0


def test_a_parent_child_total_returns_parent_child_among_the_candidates():
    prediction = rel.predict_relationship(3500)
    assert "parent or child" in prediction["relationships"]


def test_a_quarter_sharing_total_cannot_be_narrowed_to_one_relationship():
    prediction = rel.predict_relationship(1700)
    names = prediction["relationships"]
    assert "half sibling" in names
    assert "grandparent or grandchild" in names
    assert "aunt, uncle, niece or nephew" in names
    assert len(names) >= 3


def test_prediction_never_offers_a_single_confident_answer():
    prediction = rel.predict_relationship(1700)
    assert prediction["single_answer"] is None
    assert "cannot tell a half sibling" in prediction["caveat"]


def test_a_tiny_total_is_reported_as_unrelated_or_distant():
    prediction = rel.predict_relationship(5)
    assert prediction["state"] == "unrelated_or_distant"
    assert prediction["relationships"] == []


def test_an_impossibly_high_total_suggests_the_same_kit_loaded_twice():
    prediction = rel.predict_relationship(5000)
    assert prediction["state"] == "above_range"
    assert "same person" in prediction["summary"]


def test_prediction_reports_a_degree_range_rather_than_one_degree():
    prediction = rel.predict_relationship(1700)
    low, high = prediction["degree_range"]
    assert low < high


def test_prediction_carries_through_the_longest_segment_and_count():
    prediction = rel.predict_relationship(1700, longest_cm=120.5, segment_count=42)
    assert prediction["longest_cm"] == 120.5
    assert prediction["segment_count"] == 42


def test_prediction_tolerates_a_non_numeric_total():
    prediction = rel.predict_relationship("not a number")
    assert prediction["total_cm"] == 0.0


def test_relationship_bands_overlap_which_is_the_whole_point():
    half = next(b for b in rel.RELATIONSHIP_BANDS if b["relationship"] == "half sibling")
    gran = next(b for b in rel.RELATIONSHIP_BANDS
                if b["relationship"] == "grandparent or grandchild")
    assert half["low"] < gran["high"] and gran["low"] < half["high"]


# ---------------------------------------------------------------------------
# Parental phasing
# ---------------------------------------------------------------------------

def trio_merged(child, mother, father):
    markers = make_markers(len(child))
    sources = [
        {"label": "me", "role": "self", "snps": snp_rows(
            markers, lambda c, p, r: child[int(r.split("rs")[1]) - 1])},
    ]
    if mother is not None:
        sources.append({"label": "mum", "role": "mother", "snps": snp_rows(
            markers, lambda c, p, r: mother[int(r.split("rs")[1]) - 1])})
    if father is not None:
        sources.append({"label": "dad", "role": "father", "snps": snp_rows(
            markers, lambda c, p, r: father[int(r.split("rs")[1]) - 1])})
    return merge.merge_sources(sources)


def test_a_homozygous_mother_determines_the_transmitted_allele():
    merged = trio_merged([("A", "G")], [("A", "A")], [("G", "G")])
    phased = rel.phase_by_parents(merged)["phased"]
    assert phased["1_rs1"] == {"maternal": "A", "paternal": "G",
                               "basis": "mother_homozygous"}


def test_a_homozygous_father_determines_the_transmitted_allele():
    merged = trio_merged([("A", "G")], [("A", "G")], [("G", "G")])
    phased = rel.phase_by_parents(merged)["phased"]
    assert phased["1_rs1"]["paternal"] == "G"
    assert phased["1_rs1"]["maternal"] == "A"


def test_two_heterozygous_parents_leave_the_position_unresolved():
    merged = trio_merged([("A", "G")], [("A", "G")], [("A", "G")])
    result = rel.phase_by_parents(merged)
    assert result["phased"]["1_rs1"] is None
    assert result["counts"]["ambiguous"] == 1


def test_a_homozygous_child_is_resolved_trivially_and_counted_separately():
    merged = trio_merged([("A", "A")], [("A", "G")], [("A", "G")])
    result = rel.phase_by_parents(merged)
    assert result["phased"]["1_rs1"]["basis"] == "child_homozygous"
    assert result["counts"]["trivial_homozygous"] == 1
    assert result["counts"]["informative_heterozygous"] == 0


def test_phasing_counts_resolvable_against_ambiguous_positions():
    child = [("A", "G")] * 10
    mother = [("A", "A")] * 5 + [("A", "G")] * 5
    father = [("G", "G")] * 5 + [("A", "G")] * 5
    result = rel.phase_by_parents(trio_merged(child, mother, father))
    assert result["counts"]["resolvable"] == 5
    assert result["counts"]["ambiguous"] == 5
    assert result["rate"] == 0.5


def test_a_mendelian_inconsistency_is_recorded_and_not_phased():
    merged = trio_merged([("A", "A")], [("G", "G")], [("A", "A")])
    result = rel.phase_by_parents(merged)
    assert result["phased"]["1_rs1"] is None
    assert result["inconsistent_rsids"] == ["1_rs1"]
    assert "strand" in result["trio_note"].lower()


def test_phasing_with_only_one_parent_still_resolves_what_it_can():
    merged = trio_merged([("A", "G")], [("A", "A")], None)
    result = rel.phase_by_parents(merged)
    assert result["parents_present"] == {"mother": True, "father": False}
    assert result["phased"]["1_rs1"]["maternal"] == "A"


def test_phasing_reports_unavailable_when_no_parent_is_loaded():
    merged = trio_merged([("A", "G")], None, None)
    result = rel.phase_by_parents(merged)
    assert result["available"] is False
    assert result["counts"]["no_parent_data"] == 1


# ---------------------------------------------------------------------------
# Declared roles
# ---------------------------------------------------------------------------

def test_self_and_mother_are_expected_to_share_a_parent_child_amount():
    expectation = rel.expected_for_roles("self", "mother")
    assert expectation["relationship"] == "parent or child"
    assert expectation["low"] >= 3000


def test_role_order_does_not_matter():
    assert rel.expected_for_roles("mother", "self") == rel.expected_for_roles("self", "mother")


def test_a_mate_and_a_child_get_no_expectation_so_no_false_disagreement():
    assert rel.expected_for_roles("mate", "child") is None


def test_two_declared_siblings_are_expected_to_share_a_full_sibling_amount():
    expectation = rel.expected_for_roles("sibling", "sibling")
    assert expectation["relationship"] == "full sibling"


# ---------------------------------------------------------------------------
# Household analysis
# ---------------------------------------------------------------------------

def test_household_pairs_every_loaded_kit_including_comparison_roles():
    markers = make_markers(10)
    rows = snp_rows(markers, lambda c, p, r: ("A", "G"))
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": rows},
        {"label": "mum", "role": "mother", "snps": rows},
        {"label": "dad", "role": "father", "snps": rows},
    ])
    household = rel.analyse_household(merged, min_snps=1, min_cm=0.0)
    assert household["pair_count"] == 3
    assert {k["role"] for k in household["kits"]} == {"self", "mother", "father"}


def test_a_declared_sibling_sharing_a_full_sibling_amount_agrees():
    left, right, expected = family_pair(0.75)
    assert 2200 <= expected <= 3400
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": left},
        {"label": "sis", "role": "sibling", "snps": right},
    ])
    household = rel.analyse_household(merged)
    pair = household["pairs"][0]
    assert pair["total_cm"] == pytest.approx(expected, rel=0.01)
    assert pair["agreement"]["disagrees"] is False
    assert household["disagreement_count"] == 0


def test_a_declared_sibling_sharing_a_half_sibling_amount_is_flagged():
    left, right, expected = family_pair(0.5)
    assert 1160 <= expected < 2200
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": left},
        {"label": "sis", "role": "sibling", "snps": right},
    ])
    household = rel.analyse_household(merged)
    pair = household["pairs"][0]
    assert pair["agreement"]["disagrees"] is True
    assert household["disagreement_count"] == 1
    message = pair["agreement"]["message"]
    assert "does not match the DNA" in message
    assert "half sibling" in message
    assert "chip generation" in message


def test_household_disagreements_are_listed_separately_for_the_ui():
    left, right, _ = family_pair(0.5)
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": left},
        {"label": "sis", "role": "sibling", "snps": right},
    ])
    household = rel.analyse_household(merged)
    assert household["disagreements"][0]["declared"] == "full sibling"
    assert household["disagreements"][0]["a"]["role"] == "self"


def test_a_pair_with_no_expectation_is_not_checkable_and_never_disagrees():
    markers = make_markers(10)
    rows = snp_rows(markers, lambda c, p, r: ("A", "G"))
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": rows},
        {"label": "kid", "role": "child", "snps": rows},
        {"label": "partner", "role": "mate", "snps": rows},
    ])
    household = rel.analyse_household(merged, min_snps=1, min_cm=0.0)
    pair = next(p for p in household["pairs"]
                if {p["a"]["role"], p["b"]["role"]} == {"child", "mate"})
    assert pair["agreement"]["checkable"] is False
    assert pair["agreement"]["disagrees"] is False


def test_household_states_that_nothing_left_the_machine():
    merged = merge.merge_sources([
        {"label": "me", "role": "self",
         "snps": snp_rows(make_markers(5), lambda c, p, r: ("A", "G"))},
    ])
    household = rel.analyse_household(merged)
    assert "nothing was uploaded" in household["scope"].lower()
    assert "no database of other people" in household["scope"].lower()


def test_household_carries_both_standing_caveats():
    merged = merge.merge_sources([
        {"label": "me", "role": "self",
         "snps": snp_rows(make_markers(5), lambda c, p, r: ("A", "G"))},
    ])
    household = rel.analyse_household(merged)
    assert rel.IBD_CAVEAT in household["caveats"]
    assert rel.SHARED_CM_CAVEAT in household["caveats"]


# ---------------------------------------------------------------------------
# External adapter
# ---------------------------------------------------------------------------

def test_detect_ibd_falls_back_to_the_builtin_walk_when_ibis_is_absent():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.detect_ibd(genotypes, dict(genotypes),
                            positions=positions_of(markers))
    assert result["method"] == "builtin_ibs"
    assert result["external_available"] is False
    assert result["external"]["tool_id"] == "ibis"
    assert result["external"]["available"] is False
    assert result["external"]["not_attempted"] is True
    assert "IBIS" in result["note"]
    assert result["segment_count"] == 1


def test_the_degraded_ibis_payload_explains_how_to_enable_it():
    result = rel.detect_ibd({}, {}, positions={})
    assert result["external"]["how_to_enable"]["tool"] == "IBIS"


def test_detect_ibd_uses_ibis_output_when_the_tool_is_ready(monkeypatch, tmp_path):
    seg_line = "kitA kitB 1 1000000 51000000 IBD1 0.0 57.0 57.0 500 0 0.0"
    monkeypatch.setattr(rel.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(
        rel.external, "run",
        lambda tool_id, args, **kw: subprocess.CompletedProcess(
            args=["fake"], returncode=0, stdout=seg_line + "\n", stderr=""),
    )
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.detect_ibd(genotypes, dict(genotypes),
                            positions=positions_of(markers), workdir=tmp_path)
    assert result["method"] == "ibis"
    assert result["total_cm"] == 57.0
    assert result["cm_estimated"] is False
    assert result["builtin"]["method"] == "builtin_ibs"


def test_detect_ibd_degrades_when_the_installed_tool_fails(monkeypatch, tmp_path):
    def boom(tool_id, args, **kw):
        raise rel.external.ExternalError("IBIS exited with code 2.")

    monkeypatch.setattr(rel.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(rel.external, "run", boom)
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    result = rel.detect_ibd(genotypes, dict(genotypes),
                            positions=positions_of(markers), workdir=tmp_path)
    assert result["method"] == "builtin_ibs"
    assert result["external"]["state"] == "failed"


def test_parse_ibis_segments_reads_a_seg_row():
    segments = rel.parse_ibis_segments(
        "kitA kitB 7 1000 2000 IBD1 0.0 12.5 12.5 400 0 0.0\n")
    assert segments == [{
        "chromosome": "7", "start_bp": 1000, "end_bp": 2000,
        "start_rsid": None, "end_rsid": None, "snps": 400, "cm": 12.5,
        "cm_estimated": False, "cm_source": "ibis", "ibd_type": "IBD1",
        "opposite_homozygotes": None,
    }]


def test_parse_ibis_segments_skips_rows_it_cannot_read():
    segments = rel.parse_ibis_segments("# header\nnot enough columns\n")
    assert segments == []


# ---------------------------------------------------------------------------
# PLINK writers
# ---------------------------------------------------------------------------

def test_write_plink_map_writes_chromosome_marker_zero_and_position(tmp_path):
    markers = make_markers(3)
    path = rel.write_plink_map(markers, tmp_path / "x.map")
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].split("\t") == ["1", "1_rs1", "0", "1000000"]
    assert len(lines) == 3


def test_write_plink_ped_writes_one_row_per_sample_with_allele_pairs(tmp_path):
    markers = make_markers(2)
    samples = [{"id": "kitA", "genotypes": genotypes_of(markers, ("A", "G"))}]
    path = rel.write_plink_ped(samples, markers, tmp_path / "x.ped")
    cells = Path(path).read_text(encoding="utf-8").strip().split()
    assert cells[:6] == ["kitA", "kitA", "0", "0", "0", "-9"]
    assert cells[6:] == ["A", "G", "A", "G"]


def test_write_plink_ped_writes_a_no_call_as_the_plink_missing_code(tmp_path):
    markers = make_markers(2)
    genotypes = genotypes_of(markers, ("A", "G"))
    genotypes["1_rs2"] = ("N", "N")
    path = rel.write_plink_ped([{"id": "kitA", "genotypes": genotypes}],
                               markers, tmp_path / "x.ped")
    cells = Path(path).read_text(encoding="utf-8").strip().split()
    assert cells[8:] == ["0", "0"]


# ---------------------------------------------------------------------------
# Chromosome browser
# ---------------------------------------------------------------------------

def test_chromosome_browser_returns_every_autosome_plus_x():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    payload = rel.shared_segments(genotypes, dict(genotypes),
                                  positions=positions_of(markers))
    browser = rel.chromosome_browser_data(payload)
    assert [row["chromosome"] for row in browser["chromosomes"]] == AUTOSOMES + ["X"]


def test_chromosome_browser_positions_are_fractions_ready_for_svg():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    payload = rel.shared_segments(genotypes, dict(genotypes),
                                  positions=positions_of(markers))
    browser = rel.chromosome_browser_data(payload)
    segment = browser["chromosomes"][0]["segments"][0]
    assert 0.0 <= segment["start_fraction"] < segment["end_fraction"] <= 1.0
    assert segment["snps"] == 500


def test_chromosome_browser_keeps_empty_chromosomes_so_the_bar_still_renders():
    markers = make_markers(500)
    genotypes = genotypes_of(markers, ("A", "G"))
    payload = rel.shared_segments(genotypes, dict(genotypes),
                                  positions=positions_of(markers))
    browser = rel.chromosome_browser_data(payload)
    empty = [row for row in browser["chromosomes"] if row["segment_count"] == 0]
    assert len(empty) == len(AUTOSOMES)          # every chromosome but chr1, plus X
    assert all(row["length_bp"] > 0 for row in browser["chromosomes"])


def test_chromosome_browser_accepts_a_household_pair_and_keeps_the_labels():
    left, right, _ = family_pair(0.5, shared_snps=301, tail_snps=5)
    merged = merge.merge_sources([
        {"label": "me", "role": "self", "snps": left},
        {"label": "sis", "role": "sibling", "snps": right},
    ])
    household = rel.analyse_household(merged)
    browser = rel.chromosome_browser_data(household["pairs"][0])
    assert browser["labels"] == {"a": "me", "b": "sis"}
    assert browser["segment_count"] == len(AUTOSOMES)
    assert browser["cm_estimated"] is True


def test_chromosome_browser_reports_the_build_it_scaled_against():
    browser = rel.chromosome_browser_data({"segments": []})
    assert browser["build"] == "GRCh37"
    assert browser["total_cm"] == 0.0


# ---------------------------------------------------------------------------
# House style
# ---------------------------------------------------------------------------

def test_module_source_contains_no_em_dashes():
    source = Path(rel.__file__).read_text(encoding="utf-8")
    assert "\u2014" not in source
    assert "\u2013" not in source
