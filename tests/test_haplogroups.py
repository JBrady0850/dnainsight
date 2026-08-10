"""Tests for backend.haplogroups: backbone walking, the tri-state, and honesty fields."""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import external
from backend import haplogroups as hg


# ---------------------------------------------------------------------------
# Isolation
#
# Every test runs against an empty DNAINSIGHT_HOME so no external tool can be
# resolved and no licence can be accepted. That is the state of a fresh
# install, which is the state these adapters must degrade gracefully in.
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    monkeypatch.setenv("DNAINSIGHT_HOME", str(tmp_path / "home"))
    external.reset_cache()
    yield
    external.reset_cache()


# ---------------------------------------------------------------------------
# Helpers that build genotype maps FROM the bundled tree, so no test asserts
# a biological claim about an unverified marker.
# ---------------------------------------------------------------------------

def y_key(node):
    entry = hg.Y_BACKBONE[node]
    return (entry.get("rsid") or entry["marker"]).lower()


def derived_map(node, *, only_terminal=False):
    """Genotype map carrying the derived allele along the path to ``node``."""
    chain = hg.path_to(node)
    if only_terminal:
        chain = [node]
    out = {}
    for name in chain:
        entry = hg.Y_BACKBONE[name]
        if not entry.get("marker"):
            continue
        out[y_key(name)] = (entry["derived"], entry["derived"])
    return out


def ancestral_at(node):
    entry = hg.Y_BACKBONE[node]
    return {y_key(node): (entry["ancestral"], entry["ancestral"])}


def mt_map(node):
    """Position-keyed mtDNA map carrying every defining base along a path."""
    out = {}
    for name in hg.path_to(node, hg.MT_BACKBONE):
        for variant in hg.MT_BACKBONE[name].get("defining", []):
            out[f"mt{variant['position']}"] = (variant["base"], variant["base"])
    return out


# ---------------------------------------------------------------------------
# Tree structure
# ---------------------------------------------------------------------------

def test_tree_version_constant_is_present_and_non_empty():
    assert isinstance(hg.TREE_VERSION, str) and hg.TREE_VERSION.strip()
    assert isinstance(hg.TREE_NAME, str) and hg.TREE_NAME.strip()


def test_tree_stamp_returns_name_and_version():
    stamp = hg.tree_stamp()
    assert stamp["tree_name"] == hg.TREE_NAME
    assert stamp["tree_version"] == hg.TREE_VERSION


def test_y_backbone_size_is_in_the_expected_range():
    markers = [n for n in hg.Y_BACKBONE if n != "root"]
    assert 40 <= len(markers) <= 60


def test_mt_backbone_size_is_in_the_expected_range():
    nodes = [n for n in hg.MT_BACKBONE if n != "root"]
    assert 25 <= len(nodes) <= 40


def test_every_y_node_has_an_existing_parent_and_a_marker():
    for name, entry in hg.Y_BACKBONE.items():
        if name == "root":
            continue
        assert entry["parent"] in hg.Y_BACKBONE, name
        assert entry["marker"], name
        # Every node states its derived state SOMEWHERE, but not always as a
        # single base. M17 and M91 are indels, so they carry derived_seq and
        # their single-base fields are None rather than holding a multi-base
        # string that nothing could ever match against an array call.
        if entry.get("variant_type") == "snv":
            assert entry["derived"], name
            assert entry["derived_seq"] is None, name
        else:
            assert entry["derived_seq"], name
            assert entry["derived"] is None, name


def test_every_mt_node_has_an_existing_parent_and_defining_variants():
    for name, entry in hg.MT_BACKBONE.items():
        if name == "root":
            continue
        assert entry["parent"] in hg.MT_BACKBONE, name
        assert entry["defining"], name


def test_y_backbone_covers_the_named_major_branches():
    for node in ("A", "B", "BT", "CT", "DE", "E", "C", "F", "GHIJK", "G",
                 "HIJK", "IJK", "IJ", "I", "I1", "I2", "J", "J1", "J2", "K",
                 "LT", "NO", "N", "O", "P", "Q", "R", "R1a", "R1b", "R-M269",
                 "R-U106", "R-P312", "T", "L"):
        assert node in hg.Y_BACKBONE, node


def test_mt_backbone_covers_the_named_macro_haplogroups():
    for node in ("L0", "L1", "L2", "L3", "L4", "L5", "L6", "M", "N", "R", "U",
                 "K", "H", "HV", "V", "J", "T", "W", "X", "I", "A", "B", "C", "D"):
        assert node in hg.MT_BACKBONE, node


def test_path_to_walks_root_first():
    assert hg.path_to("R-M269")[0] == "root"
    assert hg.path_to("R-M269")[-1] == "R-M269"
    assert "R1b" in hg.path_to("R-M269")


def test_path_to_unknown_node_returns_empty():
    assert hg.path_to("not-a-node") == []


def test_marker_keys_prefers_the_rsid_then_the_marker_name():
    keys = hg.marker_keys(hg.Y_BACKBONE["R-M269"])
    assert keys[0] == "rs9786153"
    assert "m269" in keys


def test_marker_keys_falls_back_to_the_marker_name_when_no_rsid_is_recorded():
    entry = hg.Y_BACKBONE["I2"]
    assert entry["rsid"] is None
    assert hg.marker_keys(entry) == ["m438"]


# ---------------------------------------------------------------------------
# Verification honesty
# ---------------------------------------------------------------------------

def test_no_y_marker_claims_to_be_verified():
    unverified = [n for n, e in hg.Y_BACKBONE.items()
                  if n != "root" and not e["verified"]]
    assert len(unverified) == len(hg.Y_BACKBONE) - 1


def test_unverified_markers_lists_every_y_node_with_a_note():
    audit = hg.unverified_markers()
    assert audit["y_total"] == len(hg.Y_BACKBONE) - 1
    assert len(audit["y"]) == audit["y_total"]
    assert all(row["note"] for row in audit["y"])


def test_unverified_markers_reports_the_unverified_mt_nodes_only():
    audit = hg.unverified_markers()
    reported = {row["node"] for row in audit["mt"]}
    assert "L3" in reported
    assert "U" not in reported


def test_some_mt_positions_are_marked_verified():
    verified = [v for e in hg.MT_BACKBONE.values()
                for v in e.get("defining", []) if v["verified"]]
    assert len(verified) >= 10


# ---------------------------------------------------------------------------
# Marker state: derived, ancestral, not testable
# ---------------------------------------------------------------------------

def test_marker_state_reports_derived():
    entry = hg.Y_BACKBONE["R-M269"]
    state, alleles = hg.marker_state(entry, {y_key("R-M269"): (entry["derived"],) * 2})
    assert state == hg.DERIVED
    assert alleles == {entry["derived"]}


def test_marker_state_reports_ancestral():
    entry = hg.Y_BACKBONE["R-M269"]
    state, _ = hg.marker_state(entry, {y_key("R-M269"): (entry["ancestral"],) * 2})
    assert state == hg.ANCESTRAL


def test_marker_absent_from_the_array_is_not_testable_and_never_ancestral():
    entry = hg.Y_BACKBONE["R-M269"]
    state, alleles = hg.marker_state(entry, {})
    assert state == hg.NOT_ON_ARRAY
    assert state in hg.NOT_TESTABLE_STATES
    assert state != hg.ANCESTRAL
    assert alleles == set()


def test_marker_present_but_no_call_is_not_testable_and_never_ancestral():
    entry = hg.Y_BACKBONE["R-M269"]
    state, _ = hg.marker_state(entry, {y_key("R-M269"): ("N", "N")})
    assert state == hg.NO_CALL
    assert state in hg.NOT_TESTABLE_STATES
    assert state != hg.ANCESTRAL


def test_marker_with_an_unexpected_allele_is_discordant_not_ancestral():
    entry = hg.Y_BACKBONE["R-M269"]
    other = next(b for b in "ACGT" if b not in (entry["derived"], entry["ancestral"]))
    state, _ = hg.marker_state(entry, {y_key("R-M269"): (other, other)})
    assert state == hg.DISCORDANT


def test_verified_only_refuses_an_unverified_marker():
    entry = hg.Y_BACKBONE["R-M269"]
    state, _ = hg.marker_state(entry, {y_key("R-M269"): (entry["derived"],) * 2},
                               verified_only=True)
    assert state == hg.UNUSABLE
    assert state in hg.NOT_TESTABLE_STATES


# ---------------------------------------------------------------------------
# Y backbone walking
# ---------------------------------------------------------------------------

def test_call_y_backbone_walks_a_fully_derived_path():
    call = hg.call_y_backbone(derived_map("R-M269"))
    assert call["haplogroup"] == "R-M269"
    assert call["state"] == "called"
    assert call["path"][0] == "root"
    # BT is ASSUMED, not confirmed, and this is the correct answer rather than
    # a gap to close. BT is defined by M91, which the v3.3.0 dbSNP audit
    # established is a 9T to 8T length polymorphism and not a base
    # substitution. A consumer array reports two base calls per position, so no
    # array genotype can ever satisfy M91. BT sits on the path to every
    # non-A haplogroup, so this shows up on essentially every call, and saying
    # "assumed" is the honest form of it.
    assert call["assumed"] == ["BT"]


def test_call_y_backbone_stops_where_the_next_marker_is_ancestral():
    genotypes = derived_map("R1b")
    genotypes.update(ancestral_at("R-M269"))
    call = hg.call_y_backbone(genotypes)
    assert call["haplogroup"] == "R1b"
    assert "R-M269" in call["stopped_because_excluded"]
    assert "R-M269" not in call["stopped_because_not_testable"]


def test_call_y_backbone_reports_an_untested_branch_separately_from_an_excluded_one():
    call = hg.call_y_backbone(derived_map("R-M269"))
    assert set(call["stopped_because_not_testable"]) == {"R-P312", "R-U106"}
    assert call["stopped_because_excluded"] == []


def test_not_on_array_markers_never_land_in_the_ancestral_bucket():
    call = hg.call_y_backbone(derived_map("R-M269"))
    ancestral_nodes = {row["node"] for row in call["buckets"][hg.ANCESTRAL]}
    missing_nodes = {row["node"] for row in call["buckets"][hg.NOT_ON_ARRAY]}
    assert ancestral_nodes.isdisjoint(missing_nodes)
    assert call["counts"]["ancestral"] == 0
    assert call["counts"]["not_on_array"] > 0


def test_bucket_counts_add_up_to_the_whole_tree():
    call = hg.call_y_backbone(derived_map("R-M269"))
    counts = call["counts"]
    total = counts["tested"] + counts["not_testable"]
    assert total == len(hg.Y_BACKBONE) - 1


def test_call_y_backbone_skips_untested_ancestors_and_records_them_as_assumed():
    call = hg.call_y_backbone(derived_map("R-M269", only_terminal=True))
    assert call["haplogroup"] == "R-M269"
    assert "R1b" in call["assumed"]
    assert call["tested_path"] == ["R-M269"]
    assert call["confidence"] == "provisional"


def test_strict_walk_refuses_to_skip_an_untested_ancestor():
    call = hg.call_y_backbone(derived_map("R-M269", only_terminal=True),
                              skip_untestable=False)
    assert call["state"] == "unresolved"
    assert call["haplogroup"] is None


def test_two_derived_branches_are_surfaced_as_a_conflict_without_a_winner():
    genotypes = derived_map("R")
    genotypes.update(derived_map("R1", only_terminal=True))
    genotypes.update(derived_map("Q", only_terminal=True))
    call = hg.call_y_backbone(genotypes)
    assert call["conflicts"]
    conflict = call["conflicts"][0]
    assert set(conflict["derived_children"]) == {"Q", "R"}
    assert "no winner" in conflict["note"].lower()


def test_verified_only_y_call_is_unresolved_while_no_marker_is_confirmed():
    call = hg.call_y_backbone(derived_map("R-M269"), verified_only=True)
    assert call["state"] == "unresolved"
    assert call["counts"]["unusable"] == len(hg.Y_BACKBONE) - 1
    assert any("verified_only" in c for c in call["caveats"])


def test_every_y_payload_carries_the_tree_name_and_version():
    for call in (hg.call_y_backbone(derived_map("R-M269")),
                 hg.call_y_backbone({}, "female"),
                 hg.call_y(derived_map("R-M269"))):
        assert call["tree_name"] == hg.TREE_NAME
        assert call["tree_version"] == hg.TREE_VERSION


# ---------------------------------------------------------------------------
# A sample with no Y chromosome
# ---------------------------------------------------------------------------

def test_female_sample_gets_a_clear_no_y_state_rather_than_an_error():
    call = hg.call_y_backbone({"rs1234": ("A", "G")}, "female")
    assert call["state"] == "no_y_data"
    assert call["haplogroup"] is None
    assert call["available"] is False


def test_no_y_message_says_plainly_that_there_is_no_y_chromosome_data():
    call = hg.call_y_backbone({}, "female")
    assert "no y chromosome data" in call["message"].lower()
    assert "not a failure" in call["message"].lower()


def test_empty_genotypes_give_no_y_data_not_an_exception():
    call = hg.call_y_backbone({})
    assert call["state"] == "no_y_data"


def test_a_male_hint_with_no_readable_markers_is_unresolved_not_no_y_data():
    call = hg.call_y_backbone({}, "male")
    assert call["state"] == "unresolved"


def test_call_y_on_a_female_sample_does_not_try_to_run_an_external_tool():
    call = hg.call_y({"rs1234": ("A", "G")}, sex_hint="female")
    assert call["state"] == "no_y_data"
    assert "no y chromosome data" in call["note"].lower()


# ---------------------------------------------------------------------------
# mtDNA
# ---------------------------------------------------------------------------

def test_call_mt_backbone_resolves_u_from_its_three_defining_positions():
    call = hg.call_mt_backbone(mt_map("U"))
    assert call["haplogroup"] == "U"
    assert call["state"] == "called"


def test_call_mt_backbone_resolves_h_from_rcrs_state_positions():
    call = hg.call_mt_backbone(mt_map("H"))
    assert call["haplogroup"] == "H"


def test_call_mt_backbone_excludes_a_contradicted_branch():
    genotypes = mt_map("U")
    genotypes["mt9055"] = ("G", "G")          # K expects A here
    call = hg.call_mt_backbone(genotypes)
    assert call["haplogroup"] == "U"
    assert "K" in call["stopped_because_excluded"]


def test_call_mt_backbone_reports_an_unreadable_branch_as_not_testable():
    call = hg.call_mt_backbone(mt_map("U"))
    assert "K" in call["stopped_because_not_testable"]
    assert call["counts"]["not_on_array"] > 0


def test_mt_verified_only_still_resolves_a_textbook_haplogroup():
    call = hg.call_mt_backbone(mt_map("U"), verified_only=True)
    assert call["haplogroup"] == "U"


def test_mt_payload_carries_the_tree_version():
    call = hg.call_mt_backbone(mt_map("U"))
    assert call["tree_version"] == hg.TREE_VERSION
    assert call["tree_name"] == hg.TREE_NAME


def test_mt_positions_from_merged_builds_position_keys():
    merged = {"genotypes": {
        "rs1": {"chromosome": "MT", "position": 10400, "allele1": "T", "allele2": "T"},
        "rs2": {"chromosome": "1", "position": 500, "allele1": "A", "allele2": "A"},
    }}
    view = hg.mt_positions_from_merged(merged)
    assert view["mt10400"] == ("T", "T")
    assert "rs2" not in view


def test_mt_positions_from_merged_accepts_alternate_chromosome_labels():
    merged = {"genotypes": {
        "rs1": {"chromosome": "chrM", "position": 73, "allele1": "A", "allele2": "A"},
        "rs2": {"chromosome": "26", "position": 263, "allele1": "G", "allele2": "G"},
    }}
    view = hg.mt_positions_from_merged(merged)
    assert "mt73" in view and "mt263" in view


# ---------------------------------------------------------------------------
# Resolution ceiling
# ---------------------------------------------------------------------------

def test_resolution_ceiling_counts_markers_in_the_tree():
    ceiling = hg.resolution_ceiling(derived_map("R-M269"), "Y")
    assert ceiling["markers_in_tree"] == len(hg.Y_BACKBONE) - 1


def test_resolution_ceiling_counts_only_the_markers_this_array_reads():
    genotypes = derived_map("R-M269")
    ceiling = hg.resolution_ceiling(genotypes, "Y")
    # One fewer than the genotypes supplied: M91 is an indel, so a base call at
    # its position is not a usable reading of it and must not be counted as
    # coverage the array does not really have.
    untypeable_on_path = 1
    assert ceiling["markers_available"] == len(genotypes) - untypeable_on_path
    assert ceiling["markers_missing"] == (
        ceiling["markers_in_tree"] - len(genotypes) + untypeable_on_path)


def test_resolution_ceiling_coverage_is_the_computed_ratio():
    genotypes = derived_map("R-M269")
    ceiling = hg.resolution_ceiling(genotypes, "Y")
    expected = round(ceiling["markers_available"] / ceiling["markers_in_tree"], 4)
    assert ceiling["coverage"] == expected


def test_resolution_ceiling_sentence_carries_the_computed_numbers():
    genotypes = derived_map("R-M269")
    ceiling = hg.resolution_ceiling(genotypes, "Y")
    sentence = ceiling["sentence"]
    assert f"{ceiling['markers_available']:,}" in sentence
    assert f"{ceiling['markers_in_tree']:,}" in sentence
    assert "R-M269" in sentence
    assert "Big Y-700" in sentence


def test_resolution_ceiling_reports_the_deepest_node_this_array_could_reach():
    genotypes = derived_map("R-M269")
    genotypes[y_key("R-U106")] = ("N", "N")   # on the array but not called
    genotypes[y_key("R-P312")] = (hg.Y_BACKBONE["R-P312"]["ancestral"],) * 2
    ceiling = hg.resolution_ceiling(genotypes, "Y")
    assert ceiling["ceiling_node"] == "R-P312"
    assert ceiling["ceiling_depth"] > ceiling["resolved_depth"]


def test_resolution_ceiling_for_a_female_sample_states_there_is_no_y_data():
    ceiling = hg.resolution_ceiling({}, "Y", sex_hint="female")
    assert ceiling["deepest_resolvable_node"] is None
    assert "no y chromosome data" in ceiling["sentence"].lower()


def test_mt_resolution_ceiling_computes_the_genome_fraction():
    ceiling = hg.resolution_ceiling(mt_map("U"), "MT", array_positions=2500)
    assert ceiling["mt_genome_bp"] == hg.MT_GENOME_BP
    assert ceiling["mt_genome_fraction"] == round(2500 / hg.MT_GENOME_BP, 4)
    assert "percent of the 16,569 base mitochondrial genome" in ceiling["sentence"]


def test_mt_resolution_ceiling_reports_no_fraction_when_it_was_not_told_one():
    ceiling = hg.resolution_ceiling(mt_map("U"), "MT")
    assert ceiling["mt_genome_fraction"] is None
    assert "percent of the" not in ceiling["sentence"]


def test_resolution_ceiling_carries_the_tree_version():
    ceiling = hg.resolution_ceiling(derived_map("R-M269"), "Y")
    assert ceiling["tree_version"] == hg.TREE_VERSION


def test_typical_array_figures_are_labelled_as_comparison_only():
    ceiling = hg.resolution_ceiling(mt_map("U"), "MT")
    assert "comparison" in ceiling["typical_array"]["note"]
    assert ceiling["typical_array"]["count"] == hg.TYPICAL_ARRAY["mt_positions"]


# ---------------------------------------------------------------------------
# Input writers
# ---------------------------------------------------------------------------

def test_write_hsd_emits_a_header_and_one_tab_separated_sample_row(tmp_path):
    path = hg.write_hsd({"mt10400": ("T", "T")}, "kit1", workdir=tmp_path)
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    assert lines[0].split("\t") == ["SampleId", "Range", "Haplogroup", "Polymorphisms"]
    cells = lines[1].split("\t")
    assert cells[0] == "kit1"
    assert cells[2] == "?"
    assert cells[3] == "10400T"


def test_write_hsd_omits_positions_that_match_rcrs(tmp_path):
    # 2706A is the rCRS base, so it is not a polymorphism and must not be listed.
    path = hg.write_hsd({"mt2706": ("A", "A"), "mt10400": ("T", "T")}, "kit2",
                        workdir=tmp_path)
    row = Path(path).read_text(encoding="utf-8").strip().splitlines()[1]
    assert "10400T" in row
    assert "2706A" not in row


def test_write_hsd_range_covers_the_positions_actually_read(tmp_path):
    path = hg.write_hsd({"mt10400": ("T", "T"), "mt13263": ("G", "G")}, "kit3",
                        workdir=tmp_path)
    row = Path(path).read_text(encoding="utf-8").strip().splitlines()[1]
    assert row.split("\t")[1] == "10400-13263"


def test_write_hsd_with_no_readable_positions_falls_back_to_the_whole_genome(tmp_path):
    path = hg.write_hsd({}, "kit4", workdir=tmp_path)
    row = Path(path).read_text(encoding="utf-8").strip().splitlines()[1]
    assert row.split("\t")[1] == f"1-{hg.MT_GENOME_BP}"


def test_write_y_array_input_lists_only_markers_that_were_called(tmp_path):
    genotypes = derived_map("R-M269")
    genotypes[y_key("R-U106")] = ("N", "N")
    path = hg.write_y_array_input(genotypes, workdir=tmp_path)
    text = Path(path).read_text(encoding="utf-8")
    body = [ln for ln in text.splitlines() if not ln.startswith("#")]
    # Minus the root, which carries no marker, and minus M91, which is an indel
    # an array cannot call. Writing it into the tool input would hand a
    # downstream caller a base call for a length polymorphism.
    assert len(body) == len(hg.path_to("R-M269")) - 2
    assert "rs16981293" not in text
    assert hg.TREE_VERSION in text


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def test_parse_yleaf_output_reads_the_haplogroup_column():
    parsed = hg.parse_yleaf_output(
        "Sample_name\tHg\tHg_marker\nkit1\tR-M269\tM269\n")
    assert parsed["haplogroup"] == "R-M269"
    assert parsed["sample"] == "kit1"
    assert parsed["parsed"] is True


def test_parse_haplogrep_output_reads_haplogroup_and_quality():
    parsed = hg.parse_haplogrep_output(
        "SampleID\tHaplogroup\tQuality\nkit1\tH1a\t0.94\n")
    assert parsed["haplogroup"] == "H1a"
    assert parsed["quality"] == "0.94"


def test_parse_cladefinder_output_reads_the_clade_column():
    parsed = hg.parse_cladefinder_output("Sample\tClade\nkit1\tR-U106\n")
    assert parsed["haplogroup"] == "R-U106"


def test_parsers_report_failure_rather_than_raising_on_junk():
    for parser in (hg.parse_yleaf_output, hg.parse_haplogrep_output,
                   hg.parse_cladefinder_output):
        parsed = parser("not a table at all")
        assert parsed["parsed"] is False
        assert parsed["haplogroup"] is None


# ---------------------------------------------------------------------------
# Agreement between callers
# ---------------------------------------------------------------------------

def test_identical_calls_agree():
    result = hg.compare_y_calls({"haplogroup": "R-M269"}, {"haplogroup": "R-M269"})
    assert result["agree"] is True
    assert result["conflict"] is False


def test_a_shallower_call_on_the_same_branch_is_not_a_disagreement():
    result = hg.compare_y_calls({"haplogroup": "R1b"}, {"haplogroup": "R-M269"})
    assert result["relation"] == "consistent_depth"
    assert result["conflict"] is False
    assert result["deeper"] == "R-M269"


def test_divergent_calls_are_surfaced_as_a_conflict_with_no_winner():
    result = hg.compare_y_calls(
        {"haplogroup": "R-M269", "source": "yleaf"},
        {"haplogroup": "I1", "source": "cladefinder"},
    )
    assert result["conflict"] is True
    assert result["agree"] is False
    assert result["resolution"] is None
    assert {c["haplogroup"] for c in result["calls"]} == {"R-M269", "I1"}


def test_a_missing_call_is_reported_incomparable_rather_than_agreeing():
    result = hg.compare_y_calls({"haplogroup": None}, {"haplogroup": "I1"})
    assert result["comparable"] is False
    assert result["agree"] is None
    assert "not agreement" in result["note"]


# ---------------------------------------------------------------------------
# Graceful degradation with no external tools installed
# ---------------------------------------------------------------------------

def test_call_y_degrades_to_the_bundled_backbone_when_yleaf_is_absent():
    call = hg.call_y(derived_map("R-M269"))
    assert call["source"] == "bundled_backbone"
    assert call["external"]["available"] is False
    assert call["external"]["tool_id"] == "yleaf"
    assert call["haplogroup"] == "R-M269"
    assert "Yleaf" in call["note"]
    assert call["deeper_call_requires"] == "Yleaf"


def test_call_mt_degrades_to_the_bundled_backbone_when_haplogrep_is_absent():
    call = hg.call_mt(mt_map("U"))
    assert call["source"] == "bundled_backbone"
    assert call["external"]["tool_id"] == "haplogrep"
    assert call["haplogroup"] == "U"
    assert "HaploGrep" in call["note"]


def test_degraded_payloads_say_the_analysis_was_not_attempted():
    call = hg.call_y(derived_map("R-M269"))
    assert call["external"]["not_attempted"] is True
    assert call["external"]["how_to_enable"] is not None


def test_second_opinion_is_not_attempted_when_clade_finder_is_absent():
    result = hg.second_opinion_y(derived_map("R-M269"))
    assert result["state"] == "not_attempted"
    assert result["conflict"] is False
    assert result["second_opinion"] is None
    assert "independently checked" in result["note"]


def test_analyse_returns_both_systems_both_ceilings_and_the_audit_list():
    result = hg.analyse(derived_map("R-M269"), sex_hint="male",
                        mt_genotypes=mt_map("U"), array_mt_positions=2500)
    assert result["y"]["haplogroup"] == "R-M269"
    assert result["mt"]["haplogroup"] == "U"
    assert result["y_ceiling"]["system"] == "Y"
    assert result["mt_ceiling"]["system"] == "MT"
    assert result["unverified"]["y_total"] == len(hg.Y_BACKBONE) - 1
    assert result["tree_version"] == hg.TREE_VERSION


# ---------------------------------------------------------------------------
# The ready path, exercised with a stand-in for the external tool
# ---------------------------------------------------------------------------

def _fake_completed(stdout):
    return subprocess.CompletedProcess(args=["fake"], returncode=0, stdout=stdout,
                                       stderr="")


def test_call_y_uses_yleaf_output_when_the_tool_is_ready(monkeypatch, tmp_path):
    monkeypatch.setattr(hg.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(
        hg.external, "run",
        lambda tool_id, args, **kw: _fake_completed("Sample_name\tHg\nkit1\tR-U106\n"),
    )
    call = hg.call_y(derived_map("R-M269"), workdir=tmp_path)
    assert call["source"] == "yleaf"
    assert call["haplogroup"] == "R-U106"
    assert call["backbone"]["haplogroup"] == "R-M269"
    assert call["bundled_tree"]["tree_version"] == hg.TREE_VERSION


def test_call_y_degrades_when_the_installed_tool_fails(monkeypatch, tmp_path):
    def boom(tool_id, args, **kw):
        raise hg.external.ExternalError("Yleaf exited with code 1.")

    monkeypatch.setattr(hg.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(hg.external, "run", boom)
    call = hg.call_y(derived_map("R-M269"), workdir=tmp_path)
    assert call["source"] == "bundled_backbone"
    assert call["external"]["state"] == "failed"
    assert call["haplogroup"] == "R-M269"


def test_second_opinion_surfaces_a_disagreement_between_two_callers(monkeypatch, tmp_path):
    monkeypatch.setattr(hg.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(
        hg.external, "run",
        lambda tool_id, args, **kw: _fake_completed("Sample\tClade\nkit1\tI1\n"),
    )
    primary = {"haplogroup": "R-M269", "source": "yleaf"}
    result = hg.second_opinion_y(derived_map("R-M269"), primary=primary,
                                 workdir=tmp_path)
    assert result["conflict"] is True
    assert result["comparison"]["resolution"] is None
    assert result["second_opinion"]["haplogroup"] == "I1"


def test_call_mt_uses_haplogrep_output_when_the_tool_is_ready(monkeypatch, tmp_path):
    monkeypatch.setattr(hg.external, "guard", lambda tool_id, capability: None)
    monkeypatch.setattr(
        hg.external, "run",
        lambda tool_id, args, **kw: _fake_completed(
            "SampleID\tHaplogroup\tQuality\nkit1\tU5a1\t0.91\n"),
    )
    call = hg.call_mt(mt_map("U"), workdir=tmp_path)
    assert call["source"] == "haplogrep"
    assert call["haplogroup"] == "U5a1"
    assert call["backbone"]["haplogroup"] == "U"


# ---------------------------------------------------------------------------
# House style
# ---------------------------------------------------------------------------

def test_module_source_contains_no_em_dashes():
    source = Path(hg.__file__).read_text(encoding="utf-8")
    assert "\u2014" not in source
    assert "\u2013" not in source
