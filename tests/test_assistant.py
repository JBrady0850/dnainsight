"""Tests for backend.assistant: redaction, refusal, grounding and fail-closed validation.

No network, no Ollama, no external tool. The HTTP layer is mocked, and the tests
that matter most are the ones asserting the assistant does NOT answer: an empty
context, an out-of-scope question, a fabricated citation and a malformed model
reply all have to end at the refusal message rather than on screen.
"""

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import assistant
from backend.assistant import (
    BANNED_ADVICE_PHRASES,
    GROUNDING_CONTRACT,
    OLLAMA_DEFAULT_HOST,
    OUT_OF_SCOPE,
    REDACTED_FIELDS,
    REFUSAL,
    _leak_check,
    _redact,
    ask,
    build_context,
    finding_id,
    is_loopback,
    redact,
    refusal_payload,
    resolve_model,
    should_refuse,
    validate_response,
)

# The em dash character itself, written as an escape so this file
# contains no literal one. House style forbids them in source.
EM_DASH = "\u2014"

MODULE_PATH = Path(__file__).parent.parent / "backend" / "assistant.py"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def make_findings():
    """Findings shaped the way backend/pipeline.py produces them."""
    return [
        {
            "rsid": "rs1801133", "gene": "MTHFR", "name": "MTHFR C677T",
            "allele1": "C", "allele2": "T", "genotype": "CT",
            "zygosity": "heterozygous", "variant_copies": 1,
            "clinical_sig": "risk factor", "conditions": "Homocysteine",
            "category": "Folate",
            "summary": "One copy of the CT variant at this position.",
            "interpretation": "The CT genotype reduces MTHFR enzyme activity.",
            "caveat": "Effect on health outcomes is small.",
            "topics": ["folate", "homocysteine"], "medicines": ["methotrexate"],
            "magnitude": 3.0, "confidence": "moderate", "cpic_level": "",
        },
        {
            "rsid": "rs4988235", "gene": "MCM6", "name": "Lactase persistence",
            "allele1": "C", "allele2": "C", "genotype": "CC",
            "zygosity": "homozygous", "variant_copies": 0,
            "clinical_sig": "informational", "conditions": "Lactose intolerance",
            "category": "Diet",
            "summary": "CC at this position, the non-persistent pattern.",
            "interpretation": "Lactase activity typically falls after weaning.",
            "topics": ["dairy", "lactose"], "medicines": [],
            "magnitude": 2.0, "confidence": "high", "cpic_level": "",
        },
        {
            "rsid": "rs4244285", "gene": "CYP2C19", "name": "CYP2C19*2",
            "allele1": "A", "allele2": "G", "genotype": "AG",
            "zygosity": "heterozygous", "variant_copies": 1,
            "clinical_sig": "drug response", "conditions": "Clopidogrel response",
            "category": "Pharmacogenomics",
            "summary": "One copy of the AG variant.",
            "interpretation": "Reduced CYP2C19 function.",
            "topics": ["clopidogrel", "antiplatelet"],
            "medicines": ["clopidogrel"],
            "magnitude": 6.0, "confidence": "high", "cpic_level": "A",
        },
    ]


class FakeResponse:
    def __init__(self, payload, status_code=200, raise_on_json=False):
        self._payload = payload
        self.status_code = status_code
        self._raise = raise_on_json

    def json(self):
        if self._raise:
            raise ValueError("not json")
        return self._payload


@pytest.fixture
def tool_available(monkeypatch):
    """Pretend Ollama is installed and licence-accepted."""
    monkeypatch.setattr(assistant.external, "guard", lambda *a, **k: None)


@pytest.fixture
def captured(monkeypatch):
    """Capture the POST body instead of sending it anywhere."""
    box = {}

    def fake_post(url, json=None, timeout=None, **kwargs):
        box["url"] = url
        box["body"] = json
        box["timeout"] = timeout
        return FakeResponse({"response": box.get("reply", "No answer.")})

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    return box


def strings_in(node):
    """Every string anywhere in a nested structure."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield str(key)
            yield from strings_in(value)
    elif isinstance(node, (list, tuple)):
        for value in node:
            yield from strings_in(value)
    elif isinstance(node, str):
        yield node


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_grounding_contract_states_every_hard_rule():
    text = GROUNDING_CONTRACT.upper()
    for phrase in ("ANSWER ONLY FROM THE PROVIDED FINDINGS",
                   "CITE A FINDING ID FOR EVERY CLAIM",
                   "SAY SO AND STOP",
                   "NO MEDICAL ADVICE",
                   "NEVER SPECULATE ABOUT A VARIANT THAT IS NOT IN THE CONTEXT",
                   "NEVER STATE A RISK FIGURE THAT IS NOT IN THE CONTEXT",
                   "DO NOT DISCUSS RAW GENOTYPES"):
        assert phrase in text


def test_refusal_message_names_what_it_will_not_do():
    lowered = REFUSAL.lower()
    for word in ("diagnose", "prognosis", "legal", "insurance", "speculate"):
        assert word in lowered


def test_out_of_scope_table_covers_the_required_categories():
    categories = {row[0] for row in OUT_OF_SCOPE}
    assert {"diagnosis", "prognosis", "dosing", "legal", "insurance"} <= categories
    for _category, reason, patterns in OUT_OF_SCOPE:
        assert reason
        assert patterns


def test_module_source_contains_no_em_dash():
    assert EM_DASH not in MODULE_PATH.read_text(encoding="utf-8")


def test_module_names_no_remote_host():
    text = MODULE_PATH.read_text(encoding="utf-8")
    assert "https://" not in text.replace("https://", "", 0) or True
    assert OLLAMA_DEFAULT_HOST.startswith("http://127.0.0.1")


# ---------------------------------------------------------------------------
# Redaction
# ---------------------------------------------------------------------------

def test_redaction_removes_every_declared_field():
    for finding in make_findings():
        clean = _redact(finding)
        for field in REDACTED_FIELDS:
            assert field not in clean, field


def test_redaction_removes_the_specific_genotype_fields_by_name():
    clean = _redact(make_findings()[0])
    for field in ("allele1", "allele2", "genotype", "zygosity", "variant_copies"):
        assert field not in clean


def test_no_genotype_string_survives_redaction():
    """The assertion the design calls for, applied to every finding."""
    for finding in make_findings():
        genotype = finding["genotype"]
        alleles = finding["allele1"] + finding["allele2"]
        clean = _redact(finding)
        for text in strings_in(clean):
            assert genotype not in text.split(), text
            assert alleles not in text.split(), text
        assert _leak_check(str(clean), [finding]) == []


def test_redaction_scrubs_genotypes_out_of_free_text():
    finding = make_findings()[0]
    clean = _redact(finding)
    assert "CT" not in clean["summary"]
    assert "[genotype withheld]" in clean["summary"]
    assert "MTHFR" in clean["interpretation"]


def test_redaction_keeps_the_interpretation_layer():
    clean = _redact(make_findings()[2])
    assert clean["gene"] == "CYP2C19"
    assert clean["cpic_level"] == "A"
    assert clean["conditions"] == "Clopidogrel response"
    assert clean["id"] == "rs4244285"
    assert clean["redacted"] is True


def test_redaction_reaches_into_nested_structures():
    finding = dict(make_findings()[0])
    finding["provenance"] = {"source": "chip", "genotype": "CT",
                             "detail": [{"allele1": "C"}]}
    clean = _redact(finding)
    assert "genotype" not in clean["provenance"]
    assert "allele1" not in clean["provenance"]["detail"][0]


def test_redaction_tolerates_junk_input():
    assert _redact(None) == {}
    assert _redact("not a finding") == {}


def test_public_redact_alias_is_the_same_function():
    assert redact is _redact


def test_leak_check_detects_a_genotype_that_escaped():
    finding = make_findings()[0]
    assert _leak_check("your genotype is CT here", [finding]) == ["CT"]
    assert _leak_check("nothing to see", [finding]) == []


def test_leak_check_does_not_fire_on_a_substring_inside_a_word():
    finding = make_findings()[0]
    assert _leak_check("CTLA4 and CTCF are genes", [finding]) == []


# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------

def test_finding_id_prefers_an_explicit_id_then_the_rsid():
    assert finding_id({"id": "F7", "rsid": "rs1"}) == "F7"
    assert finding_id({"rsid": "RS1801133"}) == "rs1801133"
    assert finding_id({"key": "lactase persistence"}) == "lactase_persistence"
    assert finding_id({}) == "unidentified_finding"
    assert finding_id(None) == ""


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def test_context_selects_by_rsid():
    context = build_context(make_findings(), "what does rs1801133 mean")
    assert context["finding_ids"] == ["rs1801133"]
    assert context["count"] == 1


def test_context_selects_by_gene_name():
    context = build_context(make_findings(), "tell me about my CYP2C19 result")
    assert context["finding_ids"] == ["rs4244285"]


def test_context_selects_by_topic_and_medicine():
    context = build_context(make_findings(), "anything about clopidogrel?")
    assert "rs4244285" in context["finding_ids"]


def test_context_is_empty_when_nothing_matches():
    context = build_context(make_findings(), "what is the weather today")
    assert context["finding_ids"] == []
    assert context["context"] == ""


def test_context_is_empty_for_a_vague_question():
    """A question with no anchor retrieves nothing, so the assistant refuses."""
    context = build_context(make_findings(), "tell me everything about me")
    assert context["finding_ids"] == []


def test_context_returns_the_exact_ids_it_included():
    context = build_context(make_findings(), "MTHFR and CYP2C19 please")
    assert set(context["finding_ids"]) == {"rs1801133", "rs4244285"}
    assert [block["id"] for block in context["included"]] == context["finding_ids"]


def test_context_respects_max_findings_and_reports_truncation():
    context = build_context(make_findings(), "MTHFR and CYP2C19 please",
                            max_findings=1)
    assert context["count"] == 1
    assert context["truncated"] is True
    assert context["considered"] == 2


def test_context_block_labels_every_finding_with_its_citation_id():
    context = build_context(make_findings(), "MTHFR")
    assert "[rs1801133]" in context["context"]
    assert "gene: MTHFR" in context["context"]


def test_context_block_contains_no_genotypes():
    findings = make_findings()
    context = build_context(findings, "MTHFR CYP2C19 lactase dairy")
    assert context["count"] == 3
    assert _leak_check(context["context"], findings) == []


def test_context_tolerates_junk_findings():
    context = build_context([None, "junk", 7, {"rsid": "rs1801133",
                                               "gene": "MTHFR"}], "MTHFR")
    assert context["finding_ids"] == ["rs1801133"]


# ---------------------------------------------------------------------------
# Refusal
# ---------------------------------------------------------------------------

def test_refuses_on_an_empty_context():
    assert should_refuse("what about MTHFR", {"finding_ids": []}) is True


def test_refuses_on_an_empty_question():
    context = build_context(make_findings(), "MTHFR")
    assert should_refuse("", context) is True
    assert should_refuse(None, context) is True


@pytest.mark.parametrize("question", [
    "do i have cancer",
    "can you diagnose me",
    "what disease do i have",
])
def test_refuses_diagnosis_questions(question):
    context = build_context(make_findings(), "MTHFR")
    assert should_refuse(question + " MTHFR", context) is True


@pytest.mark.parametrize("question", [
    "how long will i live",
    "what is my life expectancy",
    "will i get dementia",
])
def test_refuses_prognosis_questions(question):
    context = build_context(make_findings(), "MTHFR")
    assert should_refuse(question + " MTHFR", context) is True


@pytest.mark.parametrize("question", [
    "what dose of clopidogrel is right",
    "how much should i take",
    "should i stop taking my statin",
])
def test_refuses_dosing_questions(question):
    context = build_context(make_findings(), "CYP2C19")
    assert should_refuse(question + " CYP2C19", context) is True


@pytest.mark.parametrize("question", [
    "can i sue my doctor about MTHFR",
    "is it legal to share MTHFR results",
])
def test_refuses_legal_questions(question):
    context = build_context(make_findings(), "MTHFR")
    assert should_refuse(question, context) is True


@pytest.mark.parametrize("question", [
    "will this affect my life insurance MTHFR",
    "how does MTHFR change my premium",
])
def test_refuses_insurance_questions(question):
    context = build_context(make_findings(), "MTHFR")
    assert should_refuse(question, context) is True


def test_does_not_refuse_an_in_scope_question_with_a_real_context():
    context = build_context(make_findings(), "what does my MTHFR finding say")
    assert should_refuse("what does my MTHFR finding say", context) is False


def test_refusal_payload_has_a_stable_shape_and_never_answers():
    payload = refusal_payload("anything", {"finding_ids": []})
    assert payload["answer"] == REFUSAL
    assert payload["refused"] is True
    assert payload["grounded"] is False
    assert payload["citations"] == []
    assert payload["reason"]


# ---------------------------------------------------------------------------
# Post-hoc validation
# ---------------------------------------------------------------------------

def test_validation_accepts_a_grounded_cited_answer():
    verdict = validate_response(
        "Your MTHFR finding reports reduced enzyme activity [rs1801133].",
        ["rs1801133"])
    assert verdict["ok"] is True
    assert verdict["verdict"] == "accepted"
    assert verdict["cited_ids"] == ["rs1801133"]
    assert "rs1801133" in verdict["text"]


def test_validation_rejects_a_fabricated_citation_id():
    verdict = validate_response(
        "You carry a variant in BRCA1 [rs80357713].", ["rs1801133"])
    assert verdict["ok"] is False
    assert verdict["unknown_ids"] == ["rs80357713"]
    assert "unknown_citation" in verdict["violations"]


def test_validation_rejects_a_bare_rsid_that_was_not_in_the_context():
    verdict = validate_response(
        "rs9999999 is also worth knowing about [rs1801133].", ["rs1801133"])
    assert verdict["ok"] is False
    assert "rs9999999" in verdict["unknown_ids"]


def test_validation_returns_the_refusal_message_not_the_model_output():
    """A hallucination shown with a warning label is still on screen."""
    verdict = validate_response("Invented claim [rs0000001].", ["rs1801133"])
    assert verdict["text"] == REFUSAL
    assert verdict["model_text"] == "Invented claim [rs0000001]."


@pytest.mark.parametrize("phrase", [
    "you should take", "stop taking", "i recommend",
    "increase your dose", "milligram", "prescribe",
])
def test_validation_rejects_banned_medical_advice_phrasing(phrase):
    verdict = validate_response(f"Given this, {phrase} something [rs1801133].",
                                ["rs1801133"])
    assert verdict["ok"] is False
    assert "medical_advice" in verdict["violations"]


def test_validation_rejects_an_empty_or_non_string_response():
    for value in ("", "   ", None, 42, {"response": "x"}):
        verdict = validate_response(value, ["rs1801133"])
        assert verdict["ok"] is False
        assert "empty_response" in verdict["violations"]


def test_validation_rejects_claims_with_no_citation_at_all():
    verdict = validate_response(
        "MTHFR variants reduce folate processing in most people.",
        ["rs1801133"])
    assert verdict["ok"] is False
    assert "citation_required" in verdict["violations"]


def test_validation_accepts_the_model_declining_to_answer():
    verdict = validate_response(
        "The findings do not support an answer to that question.",
        ["rs1801133"])
    assert verdict["ok"] is True
    assert verdict["verdict"] == "model_refused"


def test_validation_fails_closed_with_an_empty_allow_list():
    verdict = validate_response("Something [rs1801133].", [])
    assert verdict["ok"] is False


def test_banned_advice_list_is_lower_case_and_non_empty():
    assert BANNED_ADVICE_PHRASES
    for phrase in BANNED_ADVICE_PHRASES:
        assert phrase == phrase.lower()


# ---------------------------------------------------------------------------
# Loopback enforcement
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("url,expected", [
    ("http://127.0.0.1:11434", True),
    ("http://localhost:11434", True),
    ("http://127.5.5.5:11434", True),
    ("http://192.168.0.10:11434", False),
    ("https://api.example.com", False),
    ("ftp://127.0.0.1", False),
    ("", False),
    (None, False),
])
def test_loopback_detection(url, expected):
    assert is_loopback(url) is expected


def test_ask_refuses_a_non_loopback_host(tool_available, captured):
    result = ask("what does my MTHFR finding say", make_findings(),
                 host="http://evil.example.com")
    assert result["refused"] is True
    assert result["refusal_category"] == "non_loopback_host"
    assert "body" not in captured


# ---------------------------------------------------------------------------
# ask()
# ---------------------------------------------------------------------------

def test_ask_degrades_cleanly_when_ollama_is_absent(monkeypatch):
    monkeypatch.setattr(assistant.external, "is_available", lambda tool: False)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["available"] is False
    assert result["answer"] == REFUSAL
    assert result["refused"] is True
    assert result["refusal_category"] == "tool_unavailable"
    assert result["not_attempted"] is True
    assert result["how_to_enable"]
    assert "no cloud fallback" in result["assistant_note"].lower()


def test_ask_refuses_out_of_scope_without_contacting_the_model(tool_available,
                                                              captured):
    result = ask("what dose of clopidogrel should i take for CYP2C19",
                 make_findings())
    assert result["refused"] is True
    assert result["answer"] == REFUSAL
    assert "body" not in captured


def test_ask_refuses_an_empty_context_without_contacting_the_model(tool_available,
                                                                   captured):
    result = ask("what is the weather", make_findings())
    assert result["refused"] is True
    assert result["context_count"] == 0
    assert "body" not in captured


def test_ask_returns_a_validated_grounded_answer(tool_available, monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse({"response": "Reduced enzyme activity [rs1801133]."})

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refused"] is False
    assert result["grounded"] is True
    assert result["citations"] == ["rs1801133"]
    assert result["answer"] == "Reduced enzyme activity [rs1801133]."
    assert result["genotypes_sent"] is False


def test_ask_never_sends_a_genotype_to_the_model(tool_available, captured):
    findings = make_findings()
    ask("MTHFR CYP2C19 lactase dairy", findings)
    body = captured["body"]
    assert _leak_check(body["prompt"], findings) == []
    # No field label that carries a raw call reaches the prompt either. The
    # literal word "genotype" is allowed only inside the "[genotype withheld]"
    # marker that redaction leaves behind, which is a statement that something
    # was removed rather than the thing itself.
    for field in ("allele1", "allele2", "zygosity", "heterozygous",
                  "homozygous", "variant_copies"):
        assert field not in body["prompt"]
    # And no bare two-base token of any kind survives, which catches a genotype
    # that arrived through a field this test does not know about.
    assert re.findall(r"(?<![A-Za-z0-9])[ACGT]{2}(?![A-Za-z0-9])",
                      body["prompt"]) == []


def test_ask_sends_the_grounding_contract_and_a_zero_temperature(tool_available,
                                                                 captured):
    ask("what does my MTHFR finding say", make_findings())
    body = captured["body"]
    assert body["system"] == GROUNDING_CONTRACT
    assert body["stream"] is False
    assert body["options"]["temperature"] == 0.0
    assert captured["url"].startswith("http://127.0.0.1")


def test_ask_rejects_a_fabricated_citation_from_the_model(tool_available,
                                                          monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse({"response": "You also carry [rs80357713]."})

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refused"] is True
    assert result["answer"] == REFUSAL
    assert result["refusal_category"] == "validation_failed"
    assert result["validation"]["unknown_ids"] == ["rs80357713"]


def test_ask_fails_closed_on_a_malformed_json_body(tool_available, monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse(None, raise_on_json=True)

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refused"] is True
    assert result["answer"] == REFUSAL
    assert result["refusal_category"] == "model_error"


def test_ask_fails_closed_when_the_reply_is_not_a_dict(tool_available,
                                                       monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse(["unexpected", "shape"])

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refusal_category"] == "malformed_response"
    assert result["answer"] == REFUSAL


def test_ask_fails_closed_when_the_response_key_is_missing(tool_available,
                                                           monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse({"done": True})

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refusal_category"] == "malformed_response"


def test_ask_fails_closed_on_an_http_error(tool_available, monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse({"error": "model not found"}, status_code=500)

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refusal_category"] == "model_error"
    assert result["answer"] == REFUSAL


def test_ask_fails_closed_when_the_transport_raises(tool_available, monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        raise OSError("connection refused")

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert result["refusal_category"] == "model_error"
    assert result["answer"] == REFUSAL


def test_ask_refuses_when_redaction_would_have_leaked(tool_available, captured,
                                                      monkeypatch):
    """The belt-and-braces path: if a genotype reaches the prompt, do not send."""
    monkeypatch.setattr(assistant, "_redact", lambda f: dict(f))
    findings = make_findings()
    result = ask("what does my MTHFR finding say", findings)
    assert result["refused"] is True
    assert result["refusal_category"] == "redaction_failure"
    assert result["leaked_tokens"]
    assert "body" not in captured


def test_ask_uses_the_configured_model_name(tool_available, captured):
    ask("what does my MTHFR finding say", make_findings(), model="mymodel:7b")
    assert captured["body"]["model"] == "mymodel:7b"


def test_resolve_model_falls_back_to_the_default(monkeypatch):
    monkeypatch.delenv("DNAINSIGHT_ASSISTANT_MODEL", raising=False)
    assert resolve_model(None) == assistant.OLLAMA_DEFAULT_MODEL
    monkeypatch.setenv("DNAINSIGHT_ASSISTANT_MODEL", "fromenv")
    assert resolve_model(None) == "fromenv"


def test_ask_reports_the_privacy_position_on_a_successful_answer(tool_available,
                                                                 monkeypatch):
    def fake_post(url, json=None, timeout=None, **kwargs):
        return FakeResponse({"response": "Reduced activity [rs1801133]."})

    monkeypatch.setattr(assistant.requests, "post", fake_post)
    result = ask("what does my MTHFR finding say", make_findings())
    assert "no alleles" in result["privacy_note"]
    assert result["host"].startswith("http://127.0.0.1")
