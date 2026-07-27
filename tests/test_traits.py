"""Tests for backend.traits blood group prediction and the neutral trait table."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.traits import (
    ABO_AB_TAGS,
    ABO_DECISIVE_TAG,
    ABO_LD_TAG,
    ABO_SUBGROUP_TAG,
    ABO_TAGS,
    RH_CAVEAT_BASE,
    RH_PRIMARY_TAG,
    RH_SUPPORT_TAG,
    TRAIT_KEYS,
    TRAITS,
    evaluate_trait,
    predict_abo,
    predict_blood_type,
    predict_rh,
    predict_traits,
    to_findings,
)

CONFIDENCE_ORDER = ("none", "low", "moderate", "high")

TRAIT_REQUIRED_KEYS = (
    "key",
    "name",
    "category",
    "rsids",
    "rules",
    "default_call",
    "magnitude",
    "evidence",
    "caveat",
)

# Synthetic ABO inputs. rs8176719 is an indel: D is the exon 6 deletion that
# gives group O, I is the intact allele. rs8176746 T and rs8176747 C ride with
# the B transferase, so G/G at both means no B allele.
ABO_GROUP_O = {ABO_DECISIVE_TAG: "DD"}
ABO_GROUP_A_HET = {ABO_DECISIVE_TAG: "DI", "rs8176746": "GG", "rs8176747": "GG"}
ABO_GROUP_A_HOM = {ABO_DECISIVE_TAG: "II", "rs8176746": "GG", "rs8176747": "GG"}
ABO_GROUP_B_HET = {ABO_DECISIVE_TAG: "DI", "rs8176746": "TT", "rs8176747": "CC"}
ABO_GROUP_B_HOM = {ABO_DECISIVE_TAG: "II", "rs8176746": "TT", "rs8176747": "CC"}
ABO_GROUP_AB = {ABO_DECISIVE_TAG: "II", "rs8176746": "GT", "rs8176747": "CG"}

RH_POSITIVE = {RH_PRIMARY_TAG: "TT", RH_SUPPORT_TAG: "CC"}
RH_NEGATIVE = {RH_PRIMARY_TAG: "CC", RH_SUPPORT_TAG: "TT"}


def trait(key):
    """One entry from the declarative TRAITS table."""
    for entry in TRAITS:
        if entry["key"] == key:
            return entry
    raise AssertionError(f"trait {key!r} is not in TRAITS")


def called(trait_key, genotypes):
    """Evaluate one named trait and fail loudly when it declines to answer."""
    result = evaluate_trait(trait(trait_key), genotypes)
    assert result is not None, f"{trait_key} returned None for {genotypes}"
    return result


class TestAboTags:
    def test_the_decisive_tag_is_the_exon_six_deletion(self):
        assert ABO_DECISIVE_TAG == "rs8176719"

    def test_the_a_versus_b_discriminators(self):
        assert ABO_AB_TAGS == ("rs8176746", "rs8176747")

    def test_the_linkage_and_subgroup_tags(self):
        assert ABO_LD_TAG == "rs505922"
        assert ABO_SUBGROUP_TAG == "rs1053878"

    def test_abo_tags_lists_the_tags_actually_used(self):
        assert ABO_DECISIVE_TAG in ABO_TAGS
        assert ABO_LD_TAG in ABO_TAGS
        for tag in ABO_AB_TAGS:
            assert tag in ABO_TAGS


class TestPredictAboWithoutTheDecisiveTag:
    def test_an_empty_file_is_unknown(self):
        result = predict_abo({})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_a_missing_decisive_tag_is_unknown_even_with_the_others(self):
        result = predict_abo({"rs8176746": "GG", "rs8176747": "GG", ABO_LD_TAG: "CC"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_an_nn_no_call_is_unknown(self):
        result = predict_abo({ABO_DECISIVE_TAG: "NN"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_a_zero_no_call_is_unknown(self):
        result = predict_abo({ABO_DECISIVE_TAG: "00"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_an_empty_genotype_is_unknown(self):
        result = predict_abo({ABO_DECISIVE_TAG: ""})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_an_uninterpretable_token_is_unknown(self):
        result = predict_abo({ABO_DECISIVE_TAG: "AA"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_the_caveat_names_the_missing_tag_and_why(self):
        caveat = predict_abo({})["caveat"]
        assert ABO_DECISIVE_TAG in caveat
        assert "deletion" in caveat
        assert "no-call" in caveat

    def test_the_caveat_says_no_group_can_be_assigned(self):
        assert "no blood group can be assigned" in predict_abo({})["caveat"]

    def test_the_missing_list_names_the_decisive_tag(self):
        assert ABO_DECISIVE_TAG in predict_abo({})["missing"]

    def test_the_missing_list_names_every_absent_tag(self):
        missing = predict_abo({})["missing"]
        for tag in (ABO_DECISIVE_TAG, ABO_LD_TAG, ABO_SUBGROUP_TAG) + ABO_AB_TAGS:
            assert tag in missing

    def test_nothing_is_decided_and_no_genotype_is_reported(self):
        result = predict_abo({})
        assert result["deciding"] == []
        assert result["genotype_call"] == ""

    def test_every_group_stays_possible_when_no_tag_resolved(self):
        alternatives = {a["abo"] for a in predict_abo({})["alternatives"]}
        assert alternatives == {"O", "A", "B", "AB"}

    def test_no_b_allele_narrows_the_alternatives_to_a_and_o(self):
        result = predict_abo({"rs8176746": "GG", "rs8176747": "GG"})
        assert {a["abo"] for a in result["alternatives"]} == {"A", "O"}

    def test_one_b_allele_narrows_the_alternatives_to_ab_and_b(self):
        result = predict_abo({"rs8176746": "GT", "rs8176747": "CG"})
        assert {a["abo"] for a in result["alternatives"]} == {"AB", "B"}

    def test_two_b_alleles_narrow_the_alternatives_to_b_and_o(self):
        result = predict_abo({"rs8176746": "TT", "rs8176747": "CC"})
        assert {a["abo"] for a in result["alternatives"]} == {"B", "O"}

    def test_the_result_carries_the_documented_keys(self):
        assert set(predict_abo({})) == {
            "abo", "genotype_call", "confidence", "deciding",
            "missing", "caveat", "alternatives",
        }

    def test_a_double_dash_no_call_is_not_a_confident_group_o(self):
        # Regression guard. "--" is the token 23andMe and AncestryDNA write for
        # a probe that failed. It was previously read as two "-" deletion
        # alleles, which made a failed read of the one decisive ABO tag come
        # back as a confident group O. A no-call must never support a call.
        result = predict_abo({ABO_DECISIVE_TAG: "--"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_a_single_dash_is_also_a_no_call(self):
        result = predict_abo({ABO_DECISIVE_TAG: "-"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_a_half_called_indel_is_refused(self):
        # ("D", "-") is one deletion plus one failed allele. Counting it as a
        # single deletion copy would silently promote a half read to a call.
        result = predict_abo({ABO_DECISIVE_TAG: ("D", "-")})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "none"

    def test_every_nocall_spelling_degrades_the_same_way(self):
        for token in ("--", "-", "NN", "00", "", "??", "N", "0"):
            result = predict_abo({ABO_DECISIVE_TAG: token})
            assert result["abo"] == "unknown", f"{token!r} produced {result['abo']}"
            assert result["confidence"] == "none", f"{token!r} was confident"

    def test_a_failed_probe_does_not_yield_a_confident_blood_type(self):
        # The end to end version of the same defect: a full blood type must not
        # come back as O positive at high confidence off a failed ABO probe.
        result = predict_blood_type({
            ABO_DECISIVE_TAG: "--", "rs590787": "TT", "rs586178": "CC",
        })
        assert result["abo"]["abo"] == "unknown"
        assert result["confidence"] == "none"
        assert result["blood_type"] == "unknown"


class TestPredictAboCalls:
    def test_a_deletion_homozygote_is_group_o(self):
        result = predict_abo(ABO_GROUP_O)
        assert result["abo"] == "O"
        assert result["genotype_call"] == "OO"
        assert result["confidence"] == "high"

    def test_group_o_does_not_need_the_a_versus_b_tags(self):
        caveat = predict_abo(ABO_GROUP_O)["caveat"]
        assert "transferase activity" in caveat
        assert "not needed" in caveat

    def test_a_multi_character_indel_token_is_not_interpreted(self):
        assert predict_abo({ABO_DECISIVE_TAG: "DEL"})["abo"] == "unknown"

    def test_one_deletion_and_no_b_allele_is_group_a(self):
        result = predict_abo(ABO_GROUP_A_HET)
        assert result["abo"] == "A"
        assert result["genotype_call"] == "AO"

    def test_one_deletion_and_a_b_allele_is_group_b(self):
        result = predict_abo(ABO_GROUP_B_HET)
        assert result["abo"] == "B"
        assert result["genotype_call"] == "BO"

    def test_two_intact_alleles_and_no_b_allele_is_group_a(self):
        result = predict_abo(ABO_GROUP_A_HOM)
        assert result["abo"] == "A"
        assert result["genotype_call"] == "AA"

    def test_two_intact_alleles_and_two_b_alleles_is_group_b(self):
        result = predict_abo(ABO_GROUP_B_HOM)
        assert result["abo"] == "B"
        assert result["genotype_call"] == "BB"

    def test_one_b_allele_on_an_intact_backbone_is_group_ab(self):
        result = predict_abo(ABO_GROUP_AB)
        assert result["abo"] == "AB"
        assert result["genotype_call"] == "AB"

    def test_both_discriminators_agreeing_gives_high_confidence(self):
        assert predict_abo(ABO_GROUP_A_HET)["confidence"] == "high"
        assert predict_abo(ABO_GROUP_B_HOM)["confidence"] == "high"

    def test_the_discriminators_are_listed_as_deciding(self):
        rsids = [d["rsid"] for d in predict_abo(ABO_GROUP_A_HET)["deciding"]]
        assert rsids[0] == ABO_DECISIVE_TAG
        for tag in ABO_AB_TAGS:
            assert tag in rsids

    def test_one_discriminator_alone_gives_moderate_confidence(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DI", "rs8176746": "GG"})
        assert result["abo"] == "A"
        assert result["confidence"] == "moderate"
        assert "rs8176746" in result["caveat"]

    def test_the_second_discriminator_alone_still_calls_a_group(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DI", "rs8176747": "CC"})
        assert result["abo"] == "B"
        assert result["confidence"] == "moderate"

    def test_disagreeing_discriminators_drop_to_low_confidence(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DI", "rs8176746": "TT", "rs8176747": "GG"})
        assert result["confidence"] == "low"
        assert "opposite strand" in result["caveat"]
        assert {a["abo"] for a in result["alternatives"]} == {"A", "B"}

    def test_an_intact_allele_without_discriminators_cannot_separate_a_from_b(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DI"})
        assert result["abo"] == "unknown"
        assert result["confidence"] == "low"
        assert {a["abo"] for a in result["alternatives"]} == {"A", "B"}
        assert "cannot be" in result["caveat"]

    def test_two_intact_alleles_without_discriminators_leave_three_options(self):
        result = predict_abo({ABO_DECISIVE_TAG: "II"})
        assert result["abo"] == "unknown"
        assert {a["abo"] for a in result["alternatives"]} == {"A", "AB", "B"}

    def test_a_missing_linkage_tag_is_reported_but_does_not_downgrade(self):
        result = predict_abo(ABO_GROUP_A_HET)
        assert result["confidence"] == "high"
        assert ABO_LD_TAG in result["caveat"]

    def test_an_inconsistent_linkage_tag_downgrades_the_call(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DD", ABO_LD_TAG: "TT"})
        assert result["abo"] == "O"
        assert result["confidence"] == "moderate"
        assert "does not fit" in result["caveat"]

    def test_a_consistent_linkage_tag_is_listed_as_deciding(self):
        result = predict_abo({ABO_DECISIVE_TAG: "DD", ABO_LD_TAG: "CC"})
        assert result["confidence"] == "high"
        assert ABO_LD_TAG in [d["rsid"] for d in result["deciding"]]

    def test_the_a1_versus_a2_subgroup_is_left_unresolved(self):
        for genotypes in (ABO_GROUP_A_HET, ABO_GROUP_AB):
            caveat = predict_abo(genotypes)["caveat"]
            assert ABO_SUBGROUP_TAG in caveat
            assert "unresolved" in caveat

    def test_every_call_repeats_the_serology_warning(self):
        for genotypes in (ABO_GROUP_O, ABO_GROUP_A_HET, ABO_GROUP_B_HOM, ABO_GROUP_AB):
            caveat = predict_abo(genotypes)["caveat"]
            assert "serological typing" in caveat
            assert "transfusion" in caveat

    def test_allele_order_does_not_change_the_group(self):
        first = predict_abo({ABO_DECISIVE_TAG: "ID", "rs8176746": "TG", "rs8176747": "GC"})
        second = predict_abo({ABO_DECISIVE_TAG: "DI", "rs8176746": "GT", "rs8176747": "CG"})
        assert first["abo"] == second["abo"]

    def test_uppercase_and_lowercase_rsid_keys_both_work(self):
        upper = {ABO_DECISIVE_TAG.upper(): "DD"}
        assert predict_abo(upper)["abo"] == "O"


class TestPredictRh:
    def test_a_missing_primary_tag_is_unknown(self):
        result = predict_rh({})
        assert result["rh"] == "unknown"
        assert result["confidence"] == "none"

    def test_the_support_tag_alone_is_not_enough(self):
        result = predict_rh({RH_SUPPORT_TAG: "TT"})
        assert result["rh"] == "unknown"
        assert result["confidence"] == "none"

    def test_a_no_call_primary_tag_is_unknown(self):
        for token in ("--", "NN", "00", ""):
            result = predict_rh({RH_PRIMARY_TAG: token})
            assert result["rh"] == "unknown"
            assert result["confidence"] == "none"

    def test_the_unknown_caveat_names_the_missing_proxy(self):
        caveat = predict_rh({})["caveat"]
        assert RH_PRIMARY_TAG in caveat
        assert "not genotyped" in caveat

    def test_a_negative_call_is_less_trustworthy_than_a_positive_one(self):
        caveat = predict_rh(RH_NEGATIVE)["caveat"]
        assert "more trustworthy than a negative one" in caveat
        assert RH_CAVEAT_BASE in caveat

    def test_the_caveat_admits_non_european_haplotypes_are_missed(self):
        caveat = predict_rh(RH_POSITIVE)["caveat"]
        assert "Non-European" in caveat
        assert "pseudogene" in caveat

    def test_the_caveat_forbids_transfusion_and_pregnancy_use(self):
        caveat = predict_rh(RH_POSITIVE)["caveat"]
        assert "transfusion" in caveat
        assert "pregnancy" in caveat

    def test_two_deletion_tag_alleles_call_negative(self):
        assert predict_rh({RH_PRIMARY_TAG: "CC"})["rh"] == "negative"

    def test_a_negative_call_without_support_is_low_confidence(self):
        assert predict_rh({RH_PRIMARY_TAG: "CC"})["confidence"] == "low"

    def test_a_supported_negative_call_reaches_moderate_confidence(self):
        assert predict_rh(RH_NEGATIVE)["confidence"] == "moderate"

    def test_a_negative_call_never_exceeds_moderate_confidence(self):
        for support in (None, "TT", "CC", "CT"):
            genotypes = {RH_PRIMARY_TAG: "CC"}
            if support is not None:
                genotypes[RH_SUPPORT_TAG] = support
            result = predict_rh(genotypes)
            assert result["rh"] == "negative"
            assert result["confidence"] in ("low", "moderate")

    def test_a_disagreeing_support_tag_is_reported(self):
        result = predict_rh({RH_PRIMARY_TAG: "CC", RH_SUPPORT_TAG: "CC"})
        assert result["confidence"] == "low"
        assert "does not agree" in result["caveat"]

    def test_no_deletion_tag_allele_calls_positive(self):
        assert predict_rh({RH_PRIMARY_TAG: "TT"})["rh"] == "positive"

    def test_one_deletion_tag_allele_still_calls_positive(self):
        result = predict_rh({RH_PRIMARY_TAG: "CT"})
        assert result["rh"] == "positive"
        assert "single intact copy" in result["caveat"]

    def test_a_supported_positive_call_reaches_high_confidence(self):
        assert predict_rh(RH_POSITIVE)["confidence"] == "high"

    def test_a_positive_call_without_support_is_moderate(self):
        result = predict_rh({RH_PRIMARY_TAG: "TT"})
        assert result["confidence"] == "moderate"
        assert RH_SUPPORT_TAG in result["caveat"]

    def test_the_deciding_list_starts_with_the_primary_proxy(self):
        deciding = predict_rh(RH_POSITIVE)["deciding"]
        assert deciding[0]["rsid"] == RH_PRIMARY_TAG
        assert deciding[1]["rsid"] == RH_SUPPORT_TAG

    def test_the_result_carries_the_documented_keys(self):
        assert set(predict_rh({})) == {"rh", "confidence", "deciding", "caveat"}


class TestPredictBloodType:
    def test_the_result_carries_the_documented_keys(self):
        assert set(predict_blood_type({})) == {
            "blood_type", "abo", "rh", "confidence", "summary",
        }

    def test_a_positive_group_a_file(self):
        result = predict_blood_type({**ABO_GROUP_A_HET, **RH_POSITIVE})
        assert result["blood_type"] == "A+"
        assert result["confidence"] == "high"

    def test_a_negative_group_o_file(self):
        result = predict_blood_type({**ABO_GROUP_O, **RH_NEGATIVE})
        assert result["blood_type"] == "O-"
        assert result["confidence"] == "moderate"

    def test_confidence_is_the_weaker_of_abo_and_rh(self):
        genotypes = {**ABO_GROUP_A_HET, RH_PRIMARY_TAG: "TT"}
        result = predict_blood_type(genotypes)
        assert result["abo"]["confidence"] == "high"
        assert result["rh"]["confidence"] == "moderate"
        assert result["confidence"] == "moderate"

    def test_confidence_never_exceeds_either_component(self):
        cases = [
            {},
            ABO_GROUP_O,
            RH_POSITIVE,
            {**ABO_GROUP_A_HET, **RH_POSITIVE},
            {**ABO_GROUP_O, **RH_NEGATIVE},
            {**ABO_GROUP_AB, RH_PRIMARY_TAG: "CC"},
            {**ABO_GROUP_B_HET, RH_PRIMARY_TAG: "CT"},
            {ABO_DECISIVE_TAG: "DI", **RH_POSITIVE},
        ]
        for genotypes in cases:
            result = predict_blood_type(genotypes)
            combined = CONFIDENCE_ORDER.index(result["confidence"])
            assert combined <= CONFIDENCE_ORDER.index(result["abo"]["confidence"])
            assert combined <= CONFIDENCE_ORDER.index(result["rh"]["confidence"])

    def test_an_unknown_abo_gives_an_unknown_blood_type(self):
        result = predict_blood_type(RH_POSITIVE)
        assert result["blood_type"] == "unknown"
        assert result["confidence"] == "none"

    def test_an_unknown_rh_gives_an_unknown_blood_type(self):
        result = predict_blood_type(ABO_GROUP_O)
        assert result["blood_type"] == "unknown"
        assert result["confidence"] == "none"

    def test_an_empty_file_gives_an_unknown_blood_type(self):
        result = predict_blood_type({})
        assert result["blood_type"] == "unknown"
        assert result["confidence"] == "none"

    def test_the_unknown_summary_reports_both_components(self):
        summary = predict_blood_type({})["summary"]
        assert "could not be predicted" in summary
        assert "ABO is unknown" in summary
        assert "RhD is unknown" in summary

    def test_a_successful_summary_says_it_is_a_prediction(self):
        summary = predict_blood_type({**ABO_GROUP_A_HET, **RH_POSITIVE})["summary"]
        assert "A+" in summary
        assert "prediction from proxy markers rather than a blood test" in summary

    def test_the_component_results_are_carried_through(self):
        genotypes = {**ABO_GROUP_AB, **RH_NEGATIVE}
        result = predict_blood_type(genotypes)
        assert result["abo"] == predict_abo(genotypes)
        assert result["rh"] == predict_rh(genotypes)

    def test_group_ab_negative_renders_correctly(self):
        result = predict_blood_type({**ABO_GROUP_AB, **RH_NEGATIVE})
        assert result["blood_type"] == "AB-"


class TestTraitsTable:
    def test_the_table_is_not_empty(self):
        assert len(TRAITS) >= 10

    def test_every_entry_has_the_required_keys(self):
        for entry in TRAITS:
            for key in TRAIT_REQUIRED_KEYS:
                assert key in entry, f"{entry.get('key')} is missing {key}"

    def test_every_entry_has_a_non_empty_caveat(self):
        for entry in TRAITS:
            assert entry["caveat"].strip(), entry["key"]
            assert len(entry["caveat"]) > 40, entry["key"]

    def test_every_magnitude_is_a_number(self):
        for entry in TRAITS:
            assert isinstance(entry["magnitude"], (int, float)), entry["key"]
            assert not isinstance(entry["magnitude"], bool), entry["key"]
            assert 0.0 <= float(entry["magnitude"]) <= 10.0, entry["key"]

    def test_every_entry_names_at_least_one_rsid(self):
        for entry in TRAITS:
            assert isinstance(entry["rsids"], list), entry["key"]
            assert entry["rsids"], entry["key"]
            for rsid in entry["rsids"]:
                assert rsid == rsid.lower(), entry["key"]
                assert rsid.startswith("rs"), entry["key"]

    def test_every_entry_has_at_least_one_rule(self):
        for entry in TRAITS:
            assert entry["rules"], entry["key"]
            for rule in entry["rules"]:
                assert "when" in rule, entry["key"]
                assert rule.get("call", "").strip(), entry["key"]
                assert rule.get("detail", "").strip(), entry["key"]

    def test_every_entry_has_a_default_call(self):
        for entry in TRAITS:
            assert entry["default_call"].strip(), entry["key"]

    def test_every_entry_has_evidence_and_a_category(self):
        for entry in TRAITS:
            assert entry["evidence"].strip(), entry["key"]
            assert entry["category"].strip(), entry["key"]

    def test_no_trait_declares_a_repute(self):
        for entry in TRAITS:
            assert "repute" not in entry, entry["key"]

    def test_trait_keys_are_unique(self):
        keys = [entry["key"] for entry in TRAITS]
        assert len(keys) == len(set(keys))

    def test_trait_names_are_unique(self):
        names = [entry["name"] for entry in TRAITS]
        assert len(names) == len(set(names))

    def test_trait_keys_constant_matches_the_table(self):
        assert TRAIT_KEYS == tuple(entry["key"] for entry in TRAITS)

    def test_every_rule_only_mentions_declared_rsids(self):
        for entry in TRAITS:
            declared = set(entry["rsids"])
            for rule in entry["rules"]:
                for rsid in _spec_rsids_from(rule["when"]):
                    assert rsid in declared, f"{entry['key']} rule uses {rsid}"


def _spec_rsids_from(spec):
    """Collect rsIDs mentioned by a rule spec, without importing a private helper."""
    found = []
    if not isinstance(spec, dict):
        return found
    for combinator in ("all", "any"):
        for sub in (spec.get(combinator) or []):
            found.extend(_spec_rsids_from(sub))
    if spec.get("rsid"):
        found.append(str(spec["rsid"]).strip().lower())
    return found


class TestEvaluateTrait:
    def test_no_genotypes_returns_none(self):
        assert evaluate_trait(trait("lactase_persistence"), {}) is None

    def test_an_unrelated_genotype_returns_none(self):
        assert evaluate_trait(trait("lactase_persistence"), {"rs1234": "AA"}) is None

    def test_a_no_call_returns_none(self):
        for token in ("--", "NN", "00", "", "N"):
            assert evaluate_trait(trait("lactase_persistence"),
                                  {"rs4988235": token}) is None

    def test_a_non_dict_genotype_mapping_returns_none(self):
        assert evaluate_trait(trait("lactase_persistence"), None) is None

    def test_the_result_carries_the_documented_keys(self):
        result = called("lactase_persistence", {"rs4988235": "CC"})
        assert set(result) == {
            "key", "name", "category", "call", "detail", "genotype",
            "rsid", "magnitude", "evidence", "caveat", "coverage",
        }

    def test_the_key_and_name_are_echoed(self):
        result = called("lactase_persistence", {"rs4988235": "CC"})
        assert result["key"] == "lactase_persistence"
        assert result["name"] == "Lactase persistence"

    def test_the_genotype_is_reported_in_sorted_form(self):
        assert called("lactase_persistence", {"rs4988235": "TC"})["genotype"] == "CT"

    def test_a_two_item_sequence_is_accepted(self):
        assert called("lactase_persistence", {"rs4988235": ("C", "C")})["genotype"] == "CC"

    def test_rsid_keys_are_case_insensitive(self):
        assert called("lactase_persistence", {"RS4988235": "CC"})["genotype"] == "CC"

    def test_magnitude_is_a_float(self):
        assert isinstance(called("lactase_persistence", {"rs4988235": "CC"})["magnitude"], float)

    def test_coverage_is_one_for_a_single_snp_trait(self):
        assert called("lactase_persistence", {"rs4988235": "CC"})["coverage"] == 1.0

    def test_coverage_is_partial_when_one_of_two_snps_is_called(self):
        assert called("milk_digestion_summary", {"rs4988235": "CC"})["coverage"] == 0.5

    def test_the_deciding_rsid_is_reported(self):
        result = called("milk_digestion_summary", {"rs4988235": "CC"})
        assert result["rsid"] == "rs4988235"

    def test_an_unmatched_genotype_falls_back_to_the_default_call(self):
        result = called("lactase_persistence", {"rs4988235": "AA"})
        assert result["call"] == "Undetermined"
        assert result["detail"] == ""

    def test_the_first_matching_rule_wins(self):
        result = called("milk_digestion_summary",
                        {"rs4988235": "CC", "rs182549": "CC"})
        assert result["call"] == "Both LCT tags non-persistent"

    def test_a_later_rule_is_used_when_the_first_cannot_be_decided(self):
        result = called("milk_digestion_summary", {"rs4988235": "CC"})
        assert result["call"] == "Primary LCT tag non-persistent"


class TestLactaseRule:
    def test_cc_is_non_persistence(self):
        result = called("lactase_persistence", {"rs4988235": "CC"})
        assert result["call"] == "Lactase non-persistent"
        assert "after weaning" in result["detail"]

    def test_one_t_allele_is_persistence(self):
        result = called("lactase_persistence", {"rs4988235": "CT"})
        assert "persistent" in result["call"]
        assert "non-persistent" not in result["call"]

    def test_two_t_alleles_are_persistence(self):
        result = called("lactase_persistence", {"rs4988235": "TT"})
        assert result["call"] == "Lactase persistent"

    def test_allele_order_does_not_matter(self):
        first = called("lactase_persistence", {"rs4988235": "TC"})
        second = called("lactase_persistence", {"rs4988235": "CT"})
        assert first["call"] == second["call"]

    def test_the_caveat_limits_the_tag_to_european_ancestry(self):
        result = called("lactase_persistence", {"rs4988235": "CC"})
        assert "Europeans" in result["caveat"]
        assert "African" in result["caveat"]


class TestAlcoholFlushRule:
    def test_a_heterozygote_gives_the_flush_call(self):
        result = called("alcohol_flush", {"rs671": "AG"})
        assert result["call"] == "Alcohol flush reaction likely"
        assert "acetaldehyde" in result["detail"]

    def test_allele_order_does_not_matter(self):
        assert called("alcohol_flush", {"rs671": "GA"})["call"] == \
            "Alcohol flush reaction likely"

    def test_two_variant_copies_give_the_strong_call(self):
        assert called("alcohol_flush", {"rs671": "AA"})["call"] == \
            "Strong alcohol flush reaction"

    def test_no_variant_copy_gives_no_flush(self):
        assert called("alcohol_flush", {"rs671": "GG"})["call"] == \
            "No flush reaction expected"

    def test_the_caveat_notes_the_east_asian_frequency(self):
        assert "East Asian" in called("alcohol_flush", {"rs671": "AG"})["caveat"]


class TestCombinatorRules:
    def test_min_copies_matches_a_heterozygote(self):
        assert called("bitter_taste_ptc", {"rs713598": "CG"})["call"] == "Bitter taster"

    def test_an_exact_genotype_rule_matches_the_non_taster(self):
        assert called("bitter_taste_ptc", {"rs713598": "CC"})["call"] == \
            "Bitter non-taster"

    def test_an_all_combinator_needs_every_part(self):
        result = called("bitter_taste_ptc", {"rs713598": "GG", "rs1726866": "AA"})
        assert result["call"] == "Strong bitter taster"

    def test_an_all_combinator_is_undecided_when_one_part_is_missing(self):
        result = called("bitter_taste_ptc", {"rs713598": "GG"})
        assert result["call"] == "Bitter taster"

    def test_an_any_combinator_fires_on_one_true_part(self):
        result = called("mc1r_freckling", {"rs1805007": "TT"})
        assert result["call"] == "Two MC1R red-hair variants"

    def test_an_any_combinator_falls_through_when_all_parts_are_false(self):
        result = called("mc1r_freckling", {"rs1805007": "CC", "rs1805008": "CC"})
        assert result["call"] == "No MC1R R151C or R160W variant detected"

    def test_a_compound_carrier_is_recognised(self):
        result = called("mc1r_freckling", {"rs1805007": "CT", "rs1805008": "CT"})
        assert result["call"] == "Compound MC1R variant carrier"

    def test_a_single_variant_carrier_is_recognised(self):
        result = called("mc1r_freckling", {"rs1805007": "CT", "rs1805008": "CC"})
        assert result["call"] == "MC1R variant carrier"


class TestPredictTraits:
    def test_no_genotypes_gives_no_results(self):
        assert predict_traits({}) == []

    def test_no_calls_give_no_results(self):
        assert predict_traits({"rs4988235": "--", "rs671": "NN"}) == []

    def test_only_traits_the_file_can_answer_are_returned(self):
        keys = {r["key"] for r in predict_traits({"rs671": "AG"})}
        assert keys == {"alcohol_flush"}

    def test_a_shared_rsid_answers_every_trait_that_uses_it(self):
        keys = {r["key"] for r in predict_traits({"rs4988235": "CC"})}
        assert keys == {"lactase_persistence", "milk_digestion_summary"}

    def test_one_entry_per_called_trait(self):
        results = predict_traits({"rs671": "AG", "rs4988235": "TT", "rs12913832": "GG"})
        assert len(results) == len({r["key"] for r in results})
        assert len(results) == 4

    def test_results_follow_the_table_order(self):
        results = predict_traits({"rs671": "AG", "rs4988235": "TT"})
        keys = [r["key"] for r in results]
        assert keys == [k for k in TRAIT_KEYS if k in set(keys)]

    def test_every_result_has_a_non_empty_call(self):
        results = predict_traits({"rs671": "AG", "rs4988235": "TT", "rs1815739": "CT"})
        assert results
        for result in results:
            assert result["call"].strip()

    def test_an_ungenotyped_trait_is_absent_rather_than_undetermined(self):
        keys = {r["key"] for r in predict_traits({"rs671": "AG"})}
        assert "eye_colour_herc2" not in keys


def _every_trait_genotyped():
    """A real call at every rsID any trait reads, so no trait is skipped."""
    genotypes = {}
    for entry in TRAITS:
        for rsid in entry["rsids"]:
            genotypes[rsid] = "AG"
    return genotypes


class TestEyeColour:
    def test_two_low_expression_alleles_predict_light_eyes(self):
        assert called("eye_colour_herc2", {"rs12913832": "GG"})["call"] == \
            "Blue or light eyes likely"

    def test_a_heterozygote_is_intermediate(self):
        assert called("eye_colour_herc2", {"rs12913832": "AG"})["call"] == \
            "Intermediate or brown eyes"

    def test_two_high_expression_alleles_predict_brown_eyes(self):
        assert called("eye_colour_herc2", {"rs12913832": "AA"})["call"] == \
            "Brown eyes likely"

    def test_the_caveat_limits_the_variant_to_european_ancestry(self):
        caveat = called("eye_colour_herc2", {"rs12913832": "GG"})["caveat"]
        assert "Europeans" in caveat
        assert "outside Europe" in caveat

    def test_the_caveat_frames_it_as_probability_not_determination(self):
        caveat = called("eye_colour_herc2", {"rs12913832": "AG"})["caveat"]
        assert "probability, not a determination" in caveat

    def test_the_caveat_admits_many_genes_contribute(self):
        caveat = called("eye_colour_herc2", {"rs12913832": "AA"})["caveat"]
        assert "many genes contribute" in caveat


class TestToFindings:
    def test_no_results_give_no_findings(self):
        assert to_findings([]) == []
        assert to_findings(None) == []

    def test_every_finding_is_typed_as_a_trait(self):
        findings = to_findings(predict_traits(_every_trait_genotyped()))
        assert findings
        for finding in findings:
            assert finding["entity_type"] == "trait", finding["rsid"]

    def test_every_finding_has_an_empty_repute(self):
        findings = to_findings(predict_traits(_every_trait_genotyped()))
        assert findings
        for finding in findings:
            assert finding["repute"] == "", finding["rsid"]

    def test_no_trait_is_ever_coloured_good_or_bad(self):
        genotypes = _every_trait_genotyped()
        blood = predict_blood_type({**ABO_GROUP_O, **RH_NEGATIVE})
        findings = to_findings(predict_traits(genotypes), blood)
        assert findings
        for finding in findings:
            assert finding["repute"] not in ("Good", "Bad")

    def test_every_finding_sits_in_the_informational_silo(self):
        findings = to_findings(predict_traits(_every_trait_genotyped()))
        for finding in findings:
            assert finding["silo"] == "informational"
            assert finding["clinical_sig"] == "informational"

    def test_the_finding_rsid_is_the_trait_key(self):
        findings = to_findings(predict_traits({"rs671": "AG"}))
        assert findings[0]["rsid"] == "alcohol_flush"
        assert findings[0]["source_rsid"] == "rs671"

    def test_alleles_and_zygosity_come_from_the_genotype(self):
        findings = to_findings(predict_traits({"rs671": "AG"}))
        assert (findings[0]["allele1"], findings[0]["allele2"]) == ("A", "G")
        assert findings[0]["zygosity"] == "heterozygous"

    def test_a_homozygote_is_reported_as_homozygous(self):
        findings = to_findings(predict_traits({"rs671": "AA"}))
        assert findings[0]["zygosity"] == "homozygous"

    def test_the_summary_joins_the_call_and_the_detail(self):
        results = predict_traits({"rs671": "AG"})
        findings = to_findings(results)
        assert findings[0]["summary"].startswith(results[0]["call"])
        assert results[0]["detail"] in findings[0]["summary"]

    def test_the_interpretation_carries_the_limitation(self):
        findings = to_findings(predict_traits({"rs671": "AG"}))
        assert "Limitation:" in findings[0]["interpretation"]

    def test_the_caveat_is_carried_through_verbatim(self):
        results = predict_traits({"rs671": "AG"})
        findings = to_findings(results)
        assert findings[0]["caveat"] == results[0]["caveat"]

    def test_findings_declare_their_source(self):
        for finding in to_findings(predict_traits({"rs671": "AG"})):
            assert finding["sources"] == ["dnainsight_traits"]

    def test_magnitude_and_coverage_are_carried_through(self):
        findings = to_findings(predict_traits({"rs671": "AG"}))
        assert isinstance(findings[0]["magnitude"], float)
        assert findings[0]["coverage"] == 1.0

    def test_no_blood_row_is_added_without_a_blood_result(self):
        findings = to_findings(predict_traits({"rs671": "AG"}))
        assert [f["rsid"] for f in findings] == ["alcohol_flush"]

    def test_a_blood_result_is_appended_as_a_neutral_trait(self):
        blood = predict_blood_type({**ABO_GROUP_A_HET, **RH_POSITIVE})
        findings = to_findings(predict_traits({"rs671": "AG"}), blood)
        blood_finding = findings[-1]
        assert blood_finding["rsid"] == "blood_type"
        assert blood_finding["entity_type"] == "trait"
        assert blood_finding["repute"] == ""
        assert blood_finding["category"] == "Blood"
        assert blood_finding["genotype"] == "AO"

    def test_the_blood_row_repeats_both_limitations(self):
        blood = predict_blood_type({**ABO_GROUP_A_HET, **RH_POSITIVE})
        findings = to_findings([], blood)
        assert len(findings) == 1
        interpretation = findings[0]["interpretation"]
        assert "ABO limitation:" in interpretation
        assert "RhD limitation:" in interpretation

    def test_an_unknown_blood_type_reports_zero_coverage(self):
        findings = to_findings([], predict_blood_type({}))
        assert findings[0]["coverage"] == 0.0
        assert findings[0]["genotype"] == "NN"

    def test_a_known_blood_type_reports_full_coverage(self):
        blood = predict_blood_type({**ABO_GROUP_O, **RH_NEGATIVE})
        findings = to_findings([], blood)
        assert findings[0]["coverage"] == 1.0

    def test_the_blood_row_points_at_the_decisive_tag(self):
        findings = to_findings([], predict_blood_type({}))
        assert findings[0]["source_rsid"] == ABO_DECISIVE_TAG

    def test_findings_never_claim_a_chromosome_or_position(self):
        findings = to_findings(predict_traits(_every_trait_genotyped()),
                               predict_blood_type({}))
        for finding in findings:
            assert finding["chromosome"] == ""
            assert finding["position"] == 0
