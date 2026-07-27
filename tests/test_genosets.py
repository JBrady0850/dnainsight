"""Tests for backend.genosets criteria parsing, evaluation and the authored corpus."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.genosets import (
    GENOSET_FILE,
    NOCALL_ALLELES,
    VALID_SILOS,
    CriteriaError,
    evaluate,
    evaluate_all,
    evaluate_all_verbose,
    get_metadata,
    load_genosets,
    parse_criteria,
    referenced_genosets,
    required_rsids,
    strip_comments,
    topological_order,
)
from data.build_reference import REFERENCE

REFERENCE_RSIDS = {str(row[0]).strip().lower() for row in REFERENCE}

REQUIRED_ENTRY_KEYS = (
    "criteria",
    "magnitude",
    "repute",
    "summary",
    "interpretation",
    "category",
    "silo",
    "evidence",
)

# rs429358 C marks e4 and T marks e2 or e3. rs7412 T marks e2 and C marks
# e3 or e4. The six common diplotypes therefore look like this.
APOE_GENOSET = {
    "e2/e2": "dgs001",
    "e2/e3": "dgs002",
    "e2/e4": "dgs003",
    "e3/e3": "dgs004",
    "e3/e4": "dgs005",
    "e4/e4": "dgs006",
}

APOE_GENOTYPES = {
    "e2/e2": {"rs429358": "TT", "rs7412": "TT"},
    "e2/e3": {"rs429358": "TT", "rs7412": "TC"},
    "e2/e4": {"rs429358": "TC", "rs7412": "TC"},
    "e3/e3": {"rs429358": "TT", "rs7412": "CC"},
    "e3/e4": {"rs429358": "TC", "rs7412": "CC"},
    "e4/e4": {"rs429358": "CC", "rs7412": "CC"},
}


def _corpus():
    """The authored corpus, or a skip when the generated file is missing."""
    corpus = load_genosets()
    if not corpus:
        pytest.skip("data/genosets.json has not been generated")
    return corpus


def _apoe_hits(diplotype):
    """Which of the six APOE genosets fire for one synthetic diplotype."""
    matched = {f["rsid"] for f in evaluate_all(APOE_GENOTYPES[diplotype], _corpus())}
    return matched & set(APOE_GENOSET.values())


class TestStripComments:
    def test_leading_space_hash_comment_is_removed(self):
        assert strip_comments("   # a note\nrs1234(A;T)") == "rs1234(A;T)"

    def test_triple_hash_comment_is_removed(self):
        assert strip_comments("### heading\nrs1234(A;T)") == "rs1234(A;T)"

    def test_html_comment_is_removed(self):
        assert strip_comments("rs1234(A;T) <!-- a note -->") == "rs1234(A;T)"

    def test_multiline_html_comment_is_removed(self):
        text = "and(\n<!-- first line\nsecond line -->\nrs1234(A;T))"
        assert strip_comments(text) == "and(\nrs1234(A;T))"

    def test_unterminated_html_comment_swallows_the_remainder(self):
        assert strip_comments("rs1234(A;T) <!-- oops\nrs5678(G;G)") == "rs1234(A;T)"

    def test_indentation_is_removed(self):
        assert strip_comments("    rs1234(A;T)   ") == "rs1234(A;T)"

    def test_blank_lines_are_dropped(self):
        assert strip_comments("\n\nrs1234(A;T)\n\n") == "rs1234(A;T)"

    def test_empty_and_none_become_empty_string(self):
        assert strip_comments("") == ""
        assert strip_comments(None) == ""

    def test_comment_only_body_becomes_empty(self):
        assert strip_comments("# nothing but a comment") == ""

    def test_stripped_body_still_parses(self):
        text = "  ### rule\n  and(rs1(A), <!-- inline --> rs2(G))\n"
        assert parse_criteria(text)["op"] == "and"


class TestParseGenotypeForms:
    def test_exact_pair(self):
        node = parse_criteria("rs1234(A;T)")
        assert node == {"op": "geno", "rsid": "rs1234",
                        "alleles": ("A", "T"), "mode": "exact"}

    def test_homozygous_pair(self):
        node = parse_criteria("rs1234(T;T)")
        assert node["alleles"] == ("T", "T")
        assert node["mode"] == "exact"

    def test_single_allele_is_at_least_one_mode(self):
        node = parse_criteria("rs1234(T)")
        assert node["alleles"] == ("T",)
        assert node["mode"] == "atleast_one"

    def test_rsid_is_normalised_to_lowercase(self):
        assert parse_criteria("RS1234(A;T)")["rsid"] == "rs1234"

    def test_alleles_are_uppercased(self):
        assert parse_criteria("rs1234(a;t)")["alleles"] == ("A", "T")

    def test_i_number_identifier_is_accepted(self):
        node = parse_criteria("i12345(A;A)")
        assert node["rsid"] == "i12345"
        assert node["alleles"] == ("A", "A")

    def test_uppercase_i_number_is_normalised(self):
        assert parse_criteria("I12345(A)")["rsid"] == "i12345"

    def test_internal_whitespace_is_tolerated(self):
        node = parse_criteria("rs1234 (  A ;  T  )")
        assert node["alleles"] == ("A", "T")
        assert node["rsid"] == "rs1234"

    def test_multi_line_expression_parses(self):
        text = "and(\n    rs1801133(A;A),\n    rs1801131(G;G)\n)"
        node = parse_criteria(text)
        assert node["op"] == "and"
        assert len(node["args"]) == 2


class TestParseFunctions:
    def test_and_is_n_ary(self):
        node = parse_criteria("and(rs1(A), rs2(G), rs3(C))")
        assert node["op"] == "and"
        assert len(node["args"]) == 3

    def test_or_is_n_ary(self):
        node = parse_criteria("or(rs1(A), rs2(G), rs3(C), rs4(T))")
        assert node["op"] == "or"
        assert len(node["args"]) == 4

    def test_not_is_n_ary(self):
        node = parse_criteria("not(rs1(A), rs2(G))")
        assert node["op"] == "not"
        assert len(node["args"]) == 2

    def test_not_accepts_a_single_argument(self):
        node = parse_criteria("not(rs1(A))")
        assert node["op"] == "not"
        assert len(node["args"]) == 1

    def test_atleast_records_its_count(self):
        node = parse_criteria("atleast(2, rs1(A), rs2(G), rs3(C))")
        assert node["op"] == "atleast"
        assert node["n"] == 2
        assert len(node["args"]) == 3

    def test_function_names_are_case_insensitive(self):
        assert parse_criteria("AND(rs1(A), rs2(G))")["op"] == "and"
        assert parse_criteria("AtLeast(2, rs1(A), rs2(G))")["op"] == "atleast"

    def test_dgs_reference(self):
        assert parse_criteria("dgs026") == {"op": "gsref", "name": "dgs026"}

    def test_gs_reference_without_the_d_prefix(self):
        assert parse_criteria("gs123") == {"op": "gsref", "name": "gs123"}

    def test_genoset_reference_is_lowercased(self):
        assert parse_criteria("DGS026")["name"] == "dgs026"

    def test_arbitrary_nesting(self):
        text = "and(or(rs1(A), not(rs2(G))), atleast(2, rs3(C), rs4(T), dgs001))"
        node = parse_criteria(text)
        assert node["op"] == "and"
        assert node["args"][0]["op"] == "or"
        assert node["args"][0]["args"][1]["op"] == "not"
        assert node["args"][1]["op"] == "atleast"
        assert node["args"][1]["args"][2] == {"op": "gsref", "name": "dgs001"}

    def test_nesting_with_comments_and_newlines(self):
        text = (
            "### top level rule\n"
            "or(\n"
            "    and(rs1801133(A;A), rs1801131(G;G)),\n"
            "    # only whole-line comments are stripped\n"
            "    <!-- alternative branch -->\n"
            "    dgs007\n"
            ")\n"
        )
        node = parse_criteria(text)
        assert node["op"] == "or"
        assert node["args"][1] == {"op": "gsref", "name": "dgs007"}

    def test_a_trailing_hash_comment_is_not_a_whole_line_comment(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234(A;T)  # trailing comments are not supported")


class TestParseErrors:
    def test_unbalanced_open_parenthesis(self):
        with pytest.raises(CriteriaError):
            parse_criteria("and(rs1(A), rs2(G)")

    def test_unbalanced_close_parenthesis(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234(A;T))")

    def test_missing_genotype_parenthesis(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234")

    def test_unknown_function_name(self):
        with pytest.raises(CriteriaError):
            parse_criteria("maybe(rs1(A), rs2(G))")

    def test_unknown_bare_identifier(self):
        with pytest.raises(CriteriaError):
            parse_criteria("banana")

    def test_non_integer_atleast_count(self):
        with pytest.raises(CriteriaError):
            parse_criteria("atleast(two, rs1(A), rs2(G))")

    def test_negative_atleast_count(self):
        with pytest.raises(CriteriaError):
            parse_criteria("atleast(-1, rs1(A), rs2(G))")

    def test_zero_atleast_count(self):
        with pytest.raises(CriteriaError):
            parse_criteria("atleast(0, rs1(A))")

    def test_atleast_without_arguments(self):
        with pytest.raises(CriteriaError):
            parse_criteria("atleast(2)")

    def test_empty_argument_list(self):
        with pytest.raises(CriteriaError):
            parse_criteria("and()")

    def test_empty_argument_after_comma(self):
        with pytest.raises(CriteriaError):
            parse_criteria("and(rs1(A),)")

    def test_empty_genotype(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234()")

    def test_more_than_two_alleles(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234(A;T;G)")

    def test_trailing_junk(self):
        with pytest.raises(CriteriaError):
            parse_criteria("rs1234(A;T) rs5678(G;G)")

    def test_empty_criteria_text(self):
        with pytest.raises(CriteriaError):
            parse_criteria("")

    def test_comment_only_criteria_text(self):
        with pytest.raises(CriteriaError):
            parse_criteria("# only a comment\n")

    def test_error_message_mentions_parentheses_when_unbalanced(self):
        with pytest.raises(CriteriaError) as info:
            parse_criteria("and(rs1(A), rs2(G)")
        assert "parenthes" in str(info.value)


class TestIntrospection:
    def test_required_rsids_on_a_nested_expression(self):
        node = parse_criteria(
            "and(or(rs1801133(A), rs1801131(G)), not(rs4988235(G;G)), dgs007)"
        )
        assert required_rsids(node) == {"rs1801133", "rs1801131", "rs4988235"}

    def test_required_rsids_ignores_genoset_references(self):
        assert required_rsids(parse_criteria("or(dgs001, dgs002)")) == set()

    def test_required_rsids_on_a_single_genotype(self):
        assert required_rsids(parse_criteria("rs6025(T;T)")) == {"rs6025"}

    def test_required_rsids_deduplicates(self):
        node = parse_criteria("or(rs6025(T;T), rs6025(C;T))")
        assert required_rsids(node) == {"rs6025"}

    def test_referenced_genosets_on_a_nested_expression(self):
        node = parse_criteria(
            "and(not(dgs010), or(dgs014, and(dgs015, rs9923231(T))))"
        )
        assert referenced_genosets(node) == {"dgs010", "dgs014", "dgs015"}

    def test_referenced_genosets_is_empty_without_references(self):
        assert referenced_genosets(parse_criteria("and(rs1(A), rs2(G))")) == set()

    def test_required_rsids_rejects_a_non_node(self):
        with pytest.raises(CriteriaError):
            required_rsids("rs1234")

    def test_referenced_genosets_rejects_a_non_node(self):
        with pytest.raises(CriteriaError):
            referenced_genosets(None)


class TestEvaluationMissingData:
    def test_missing_rsid_is_false_not_none(self):
        result = evaluate(parse_criteria("rs1234(A;T)"), {})
        assert result is False

    def test_missing_rsid_is_never_imputed_to_a_major_allele(self):
        node = parse_criteria("rs1234(C;C)")
        assert evaluate(node, {}) is False
        assert evaluate(node, {"rs9999": "CC"}) is False

    def test_missing_rsid_in_at_least_one_mode_is_false(self):
        assert evaluate(parse_criteria("rs1234(T)"), {}) is False

    def test_double_dash_no_call_is_false(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": "--"}) is False

    def test_nn_no_call_is_false(self):
        assert evaluate(parse_criteria("rs1234(N;N)"), {"rs1234": "NN"}) is False

    def test_zero_no_call_is_false(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": "00"}) is False

    def test_empty_string_genotype_is_false(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": ""}) is False

    def test_none_genotype_is_false(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": None}) is False

    def test_partial_no_call_is_false(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": ("A", "-")}) is False

    def test_nocall_alleles_constant(self):
        for token in ("", "N", "-", "0", "--"):
            assert token in NOCALL_ALLELES
        assert "A" not in NOCALL_ALLELES


class TestEvaluationGenotypeMatching:
    def test_exact_match_is_order_insensitive_for_strings(self):
        node = parse_criteria("rs1234(A;T)")
        assert evaluate(node, {"rs1234": "TA"}) is True
        assert evaluate(node, {"rs1234": "AT"}) is True

    def test_exact_match_is_order_insensitive_for_tuples(self):
        node = parse_criteria("rs1234(A;T)")
        assert evaluate(node, {"rs1234": ("T", "A")}) is True

    def test_separators_inside_a_genotype_string_are_tolerated(self):
        node = parse_criteria("rs1234(A;T)")
        assert evaluate(node, {"rs1234": "T/A"}) is True
        assert evaluate(node, {"rs1234": "T|A"}) is True

    def test_lowercase_stored_alleles_match(self):
        assert evaluate(parse_criteria("rs1234(A;T)"), {"rs1234": "ta"}) is True

    def test_exact_homozygote_does_not_match_a_heterozygote(self):
        assert evaluate(parse_criteria("rs1234(T;T)"), {"rs1234": "CT"}) is False

    def test_exact_heterozygote_does_not_match_a_homozygote(self):
        assert evaluate(parse_criteria("rs1234(C;T)"), {"rs1234": "TT"}) is False

    def test_at_least_one_matches_a_heterozygote(self):
        assert evaluate(parse_criteria("rs1234(T)"), {"rs1234": "CT"}) is True

    def test_at_least_one_matches_a_homozygote(self):
        assert evaluate(parse_criteria("rs1234(T)"), {"rs1234": "TT"}) is True

    def test_at_least_one_is_false_when_the_allele_is_absent(self):
        assert evaluate(parse_criteria("rs1234(T)"), {"rs1234": "CC"}) is False

    def test_genotype_lookup_uses_the_lowercase_rsid(self):
        assert evaluate(parse_criteria("RS1234(A;A)"), {"rs1234": "AA"}) is True


class TestEvaluationBooleanOperators:
    def test_and_is_true_only_when_every_argument_is_true(self):
        node = parse_criteria("and(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "AA", "rs2": "GG"}) is True

    def test_and_is_false_when_one_argument_is_false(self):
        node = parse_criteria("and(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "AA", "rs2": "CC"}) is False

    def test_and_is_false_when_one_argument_is_missing(self):
        node = parse_criteria("and(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "AA"}) is False

    def test_or_is_true_when_any_argument_is_true(self):
        node = parse_criteria("or(rs1(A), rs2(G), rs3(C))")
        assert evaluate(node, {"rs3": "CC"}) is True

    def test_or_is_false_when_every_argument_is_false(self):
        node = parse_criteria("or(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "GG", "rs2": "CC"}) is False

    def test_not_is_true_only_when_all_arguments_are_false(self):
        node = parse_criteria("not(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "GG", "rs2": "CC"}) is True

    def test_not_is_false_when_any_argument_is_true(self):
        node = parse_criteria("not(rs1(A), rs2(G))")
        assert evaluate(node, {"rs1": "AA", "rs2": "CC"}) is False

    def test_not_of_missing_rsids_is_true(self):
        assert evaluate(parse_criteria("not(rs1(A), rs2(G))"), {}) is True

    def test_atleast_is_true_at_exactly_n_hits(self):
        node = parse_criteria("atleast(2, rs1(A), rs2(G), rs3(C))")
        assert evaluate(node, {"rs1": "AA", "rs2": "GG", "rs3": "TT"}) is True

    def test_atleast_is_true_above_n_hits(self):
        node = parse_criteria("atleast(2, rs1(A), rs2(G), rs3(C))")
        assert evaluate(node, {"rs1": "AA", "rs2": "GG", "rs3": "CC"}) is True

    def test_atleast_is_false_below_n_hits(self):
        node = parse_criteria("atleast(2, rs1(A), rs2(G), rs3(C))")
        assert evaluate(node, {"rs1": "AA"}) is False

    def test_atleast_n_greater_than_the_argument_count_is_false(self):
        node = parse_criteria("atleast(4, rs1(A), rs2(G), rs3(C))")
        assert evaluate(node, {"rs1": "AA", "rs2": "GG", "rs3": "CC"}) is False

    def test_genoset_reference_resolves_from_prior_results(self):
        node = parse_criteria("dgs001")
        assert evaluate(node, {}, {"dgs001": True}) is True
        assert evaluate(node, {}, {"dgs001": False}) is False

    def test_unknown_genoset_reference_is_false(self):
        assert evaluate(parse_criteria("dgs999"), {}, {}) is False
        assert evaluate(parse_criteria("dgs999"), {}, None) is False

    def test_not_of_a_genoset_reference(self):
        node = parse_criteria("not(dgs001, dgs002)")
        assert evaluate(node, {}, {"dgs001": False, "dgs002": False}) is True
        assert evaluate(node, {}, {"dgs001": True, "dgs002": False}) is False

    def test_evaluate_rejects_a_non_node(self):
        with pytest.raises(CriteriaError):
            evaluate("rs1234(A;T)", {})

    def test_evaluate_rejects_an_unknown_operator(self):
        with pytest.raises(CriteriaError):
            evaluate({"op": "xor", "args": []}, {})

    def test_evaluate_rejects_an_operator_without_arguments(self):
        with pytest.raises(CriteriaError):
            evaluate({"op": "and", "args": []}, {})


class TestTopologicalOrder:
    def test_names_sort_numerically_not_lexicographically(self):
        corpus = {"dgs270": "rs1(A)", "dgs100": "rs2(G)"}
        assert topological_order(corpus) == ["dgs100", "dgs270"]

    def test_nine_sorts_before_one_hundred(self):
        corpus = {"dgs100": "rs1(A)", "dgs009": "rs2(G)"}
        assert topological_order(corpus) == ["dgs009", "dgs100"]

    def test_a_reference_is_resolved_before_its_user(self):
        corpus = {"dgs002": "dgs001", "dgs001": "rs1(A)"}
        assert topological_order(corpus) == ["dgs001", "dgs002"]

    def test_a_low_numbered_user_still_follows_its_dependency(self):
        corpus = {"dgs001": "dgs270", "dgs270": "rs1(A)"}
        order = topological_order(corpus)
        assert order.index("dgs270") < order.index("dgs001")

    def test_a_chain_of_references_is_ordered(self):
        corpus = {"dgs003": "dgs002", "dgs002": "dgs001", "dgs001": "rs1(A)"}
        assert topological_order(corpus) == ["dgs001", "dgs002", "dgs003"]

    def test_a_two_node_cycle_raises(self):
        with pytest.raises(CriteriaError):
            topological_order({"dgs001": "dgs002", "dgs002": "dgs001"})

    def test_a_self_reference_raises(self):
        with pytest.raises(CriteriaError):
            topological_order({"dgs001": "dgs001"})

    def test_cycle_message_names_the_stuck_genosets(self):
        with pytest.raises(CriteriaError) as info:
            topological_order({"dgs001": "dgs002", "dgs002": "dgs001"})
        assert "cycle" in str(info.value)
        assert "dgs001" in str(info.value)

    def test_a_reference_outside_the_corpus_is_not_a_cycle(self):
        assert topological_order({"dgs001": "dgs999"}) == ["dgs001"]

    def test_empty_corpus_gives_an_empty_order(self):
        assert topological_order({}) == []

    def test_dict_entries_are_accepted_as_well_as_strings(self):
        corpus = {"dgs001": {"criteria": "rs1(A)"}, "dgs002": {"criteria": "dgs001"}}
        assert topological_order(corpus) == ["dgs001", "dgs002"]

    def test_a_bad_entry_type_raises(self):
        with pytest.raises(CriteriaError):
            topological_order({"dgs001": 42})


class TestCorpusIntegrity:
    def test_corpus_file_exists(self):
        if not GENOSET_FILE.exists():
            pytest.skip("data/genosets.json has not been generated")
        assert GENOSET_FILE.name == "genosets.json"

    def test_corpus_holds_sixty_five_authored_genosets(self):
        assert len(_corpus()) == 65

    def test_metadata_agrees_with_the_corpus_size(self):
        _corpus()
        meta = get_metadata()
        assert meta.get("genoset_count") == 65

    def test_metadata_declares_the_mit_licence(self):
        _corpus()
        assert get_metadata().get("license") == "MIT"

    def test_load_genosets_returns_empty_for_a_missing_file(self):
        assert load_genosets(Path("no-such-genosets-file.json")) == {}

    def test_every_criteria_string_parses(self):
        for name, entry in _corpus().items():
            assert parse_criteria(entry["criteria"]) is not None, name

    def test_every_criteria_string_is_non_empty(self):
        for name, entry in _corpus().items():
            assert entry["criteria"].strip(), name

    def test_every_referenced_rsid_is_in_the_bundled_reference(self):
        corpus = _corpus()
        unknown = []
        for name, entry in corpus.items():
            for rsid in required_rsids(parse_criteria(entry["criteria"])):
                if rsid not in REFERENCE_RSIDS:
                    unknown.append((name, rsid))
        assert unknown == []

    def test_every_referenced_genoset_exists_in_the_corpus(self):
        corpus = _corpus()
        dangling = []
        for name, entry in corpus.items():
            for ref in referenced_genosets(parse_criteria(entry["criteria"])):
                if ref not in corpus:
                    dangling.append((name, ref))
        assert dangling == []

    def test_corpus_has_no_cycles(self):
        corpus = _corpus()
        order = topological_order(corpus)
        assert len(order) == len(corpus)
        assert set(order) == set(corpus)

    def test_dependencies_precede_their_users_in_the_corpus_order(self):
        corpus = _corpus()
        order = topological_order(corpus)
        position = {name: index for index, name in enumerate(order)}
        for name, entry in corpus.items():
            for ref in referenced_genosets(parse_criteria(entry["criteria"])):
                assert position[ref] < position[name], f"{name} -> {ref}"

    def test_every_entry_has_the_required_keys(self):
        for name, entry in _corpus().items():
            for key in REQUIRED_ENTRY_KEYS:
                assert key in entry, f"{name} is missing {key}"

    def test_every_entry_has_a_non_empty_summary(self):
        for name, entry in _corpus().items():
            assert entry["summary"].strip(), name

    def test_every_entry_has_a_non_empty_interpretation(self):
        for name, entry in _corpus().items():
            assert entry["interpretation"].strip(), name

    def test_every_magnitude_is_a_float_between_zero_and_ten(self):
        for name, entry in _corpus().items():
            magnitude = entry["magnitude"]
            assert isinstance(magnitude, float), name
            assert 0.0 <= magnitude <= 10.0, name

    def test_every_repute_is_good_bad_or_empty(self):
        for name, entry in _corpus().items():
            assert entry["repute"] in ("Good", "Bad", ""), name

    def test_every_silo_is_valid(self):
        for name, entry in _corpus().items():
            assert entry["silo"] in VALID_SILOS, name

    def test_topics_and_medicines_are_lists(self):
        for name, entry in _corpus().items():
            assert isinstance(entry["topics"], list), name
            assert isinstance(entry["medicines"], list), name


class TestApoeDiplotypes:
    def test_all_six_apoe_genosets_are_present(self):
        corpus = _corpus()
        for diplotype, name in APOE_GENOSET.items():
            assert name in corpus, diplotype
            assert corpus[name]["aka"] == f"APOE {diplotype}"

    def test_apoe_genosets_read_only_the_two_apoe_snps(self):
        corpus = _corpus()
        for name in APOE_GENOSET.values():
            rsids = required_rsids(parse_criteria(corpus[name]["criteria"]))
            assert rsids == {"rs429358", "rs7412"}, name

    def test_e3_e3_matches_only_dgs004(self):
        assert _apoe_hits("e3/e3") == {"dgs004"}

    def test_e4_e4_matches_only_dgs006(self):
        assert _apoe_hits("e4/e4") == {"dgs006"}

    def test_e2_e2_matches_only_dgs001(self):
        assert _apoe_hits("e2/e2") == {"dgs001"}

    def test_e3_e4_matches_only_dgs005(self):
        assert _apoe_hits("e3/e4") == {"dgs005"}

    def test_e2_e4_matches_only_dgs003(self):
        assert _apoe_hits("e2/e4") == {"dgs003"}

    def test_e2_e3_matches_only_dgs002(self):
        assert _apoe_hits("e2/e3") == {"dgs002"}

    def test_each_diplotype_excludes_the_other_five(self):
        for diplotype, name in APOE_GENOSET.items():
            hits = _apoe_hits(diplotype)
            assert hits == {name}, diplotype
            for other, other_name in APOE_GENOSET.items():
                if other != diplotype:
                    assert other_name not in hits, f"{diplotype} also fired {other_name}"

    def test_allele_order_does_not_change_the_apoe_call(self):
        corpus = _corpus()
        swapped = {"rs429358": "CT", "rs7412": "CT"}
        matched = {f["rsid"] for f in evaluate_all(swapped, corpus)}
        assert "dgs003" in matched

    def test_e4_e4_carries_the_highest_apoe_magnitude(self):
        corpus = _corpus()
        magnitudes = {n: corpus[n]["magnitude"] for n in APOE_GENOSET.values()}
        assert magnitudes["dgs006"] == max(magnitudes.values())

    def test_e4_bearing_genosets_are_reputed_bad(self):
        corpus = _corpus()
        assert corpus["dgs005"]["repute"] == "Bad"
        assert corpus["dgs006"]["repute"] == "Bad"

    def test_a_no_call_at_either_apoe_snp_matches_nothing(self):
        corpus = _corpus()
        for genotypes in (
            {"rs429358": "--", "rs7412": "CC"},
            {"rs429358": "TT", "rs7412": "--"},
            {"rs429358": "TT"},
            {"rs7412": "CC"},
        ):
            matched = {f["rsid"] for f in evaluate_all(genotypes, corpus)}
            assert matched & set(APOE_GENOSET.values()) == set()


FINDING_KEYS = {
    "aka", "allele1", "allele2", "category", "chromosome", "clinical_sig",
    "conditions", "coverage", "criteria", "entity_type", "evidence", "gene",
    "genotype", "interpretation", "magnitude", "matched", "matched_rsids",
    "medicines", "position", "repute", "rsid", "silo", "sources", "summary",
    "topics", "zygosity",
}


def _all_called():
    """Every rsID the corpus tests, each given an arbitrary real call."""
    corpus = _corpus()
    rsids = set()
    for entry in corpus.values():
        rsids |= required_rsids(parse_criteria(entry["criteria"]))
    return {rsid: "AA" for rsid in rsids}


class TestEvaluateAll:
    def test_only_matches_are_returned(self):
        findings = evaluate_all(APOE_GENOTYPES["e3/e3"], _corpus())
        assert findings
        for finding in findings:
            assert finding["matched"] is True

    def test_findings_carry_the_documented_keys(self):
        findings = evaluate_all(APOE_GENOTYPES["e3/e3"], _corpus())
        for finding in findings:
            assert set(finding.keys()) == FINDING_KEYS

    def test_entity_type_is_genoset(self):
        for finding in evaluate_all(APOE_GENOTYPES["e4/e4"], _corpus()):
            assert finding["entity_type"] == "genoset"

    def test_genosets_carry_no_snp_level_fields(self):
        for finding in evaluate_all(APOE_GENOTYPES["e4/e4"], _corpus()):
            assert finding["gene"] == ""
            assert finding["chromosome"] == ""
            assert finding["position"] == 0
            assert finding["allele1"] == ""
            assert finding["allele2"] == ""
            assert finding["genotype"] == ""
            assert finding["zygosity"] == ""

    def test_coverage_is_populated_as_a_fraction(self):
        for finding in evaluate_all(APOE_GENOTYPES["e3/e4"], _corpus()):
            assert isinstance(finding["coverage"], float)
            assert 0.0 <= finding["coverage"] <= 1.0

    def test_coverage_is_one_when_every_required_rsid_is_called(self):
        findings = {f["rsid"]: f for f in evaluate_all(APOE_GENOTYPES["e3/e3"], _corpus())}
        assert findings["dgs004"]["coverage"] == 1.0
        assert findings["dgs004"]["matched_rsids"] == ["rs429358", "rs7412"]

    def test_coverage_is_zero_when_nothing_required_was_called(self):
        findings = {f["rsid"]: f for f in evaluate_all(APOE_GENOTYPES["e3/e3"], _corpus())}
        assert findings["dgs065"]["coverage"] == 0.0

    def test_matched_rsids_are_sorted(self):
        for finding in evaluate_all(_all_called(), _corpus()):
            assert finding["matched_rsids"] == sorted(finding["matched_rsids"])

    def test_results_are_sorted_by_descending_magnitude(self):
        findings = evaluate_all(_all_called(), _corpus())
        magnitudes = [f["magnitude"] or 0.0 for f in findings]
        assert magnitudes == sorted(magnitudes, reverse=True)

    def test_criteria_is_echoed_back(self):
        corpus = _corpus()
        for finding in evaluate_all(APOE_GENOTYPES["e2/e2"], corpus):
            assert finding["criteria"] == corpus[finding["rsid"]]["criteria"]

    def test_no_genotypes_at_all_still_evaluates(self):
        names = {f["rsid"] for f in evaluate_all({}, _corpus())}
        assert names == {"dgs065"}

    def test_repeated_runs_return_the_same_order(self):
        corpus = _corpus()
        first = [f["rsid"] for f in evaluate_all(_all_called(), corpus)]
        second = [f["rsid"] for f in evaluate_all(_all_called(), corpus)]
        assert first == second


class TestEvaluateAllVerbose:
    def test_returns_the_three_documented_buckets(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        assert set(result.keys()) == {"matched", "unmatched", "incomplete"}

    def test_matched_agrees_with_evaluate_all(self):
        corpus = _corpus()
        genotypes = APOE_GENOTYPES["e2/e4"]
        verbose = [f["rsid"] for f in evaluate_all_verbose(genotypes, corpus)["matched"]]
        plain = [f["rsid"] for f in evaluate_all(genotypes, corpus)]
        assert verbose == plain

    def test_every_genoset_lands_in_at_least_one_bucket(self):
        corpus = _corpus()
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], corpus)
        seen = set()
        for bucket in ("matched", "unmatched", "incomplete"):
            seen |= {f["rsid"] for f in result[bucket]}
        assert seen == set(corpus)

    def test_an_ungenotyped_genoset_is_incomplete_not_unmatched(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        incomplete = {f["rsid"] for f in result["incomplete"]}
        unmatched = {f["rsid"] for f in result["unmatched"]}
        assert "dgs007" in incomplete
        assert "dgs007" not in unmatched

    def test_unmatched_only_holds_fully_covered_genosets(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        assert result["unmatched"]
        for finding in result["unmatched"]:
            assert finding["coverage"] == 1.0
            assert finding["matched"] is False

    def test_unmatched_holds_the_five_other_apoe_genosets(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        unmatched = {f["rsid"] for f in result["unmatched"]}
        assert unmatched == {"dgs001", "dgs002", "dgs003", "dgs005", "dgs006"}

    def test_incomplete_entries_are_all_partially_covered(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        assert result["incomplete"]
        for finding in result["incomplete"]:
            assert finding["coverage"] < 1.0

    def test_a_partially_covered_hit_appears_in_matched_and_incomplete(self):
        result = evaluate_all_verbose(APOE_GENOTYPES["e3/e3"], _corpus())
        matched = {f["rsid"] for f in result["matched"]}
        incomplete = {f["rsid"] for f in result["incomplete"]}
        assert "dgs065" in matched
        assert "dgs065" in incomplete

    def test_full_coverage_empties_the_incomplete_bucket(self):
        corpus = _corpus()
        result = evaluate_all_verbose(_all_called(), corpus)
        assert result["incomplete"] == []
        assert len(result["matched"]) + len(result["unmatched"]) == len(corpus)

    def test_no_genotypes_puts_everything_in_incomplete(self):
        corpus = _corpus()
        result = evaluate_all_verbose({}, corpus)
        assert len(result["incomplete"]) == len(corpus)
        assert result["unmatched"] == []
