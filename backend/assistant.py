"""
assistant.py -- a refusal-first, citation-required, fully local genome assistant.

WHY THIS MODULE EXISTS
----------------------
Three products already bolt a language model onto personal genomics, and each
gets one important thing wrong:

  SelfDecode DecodyGPT and Sequencing.com's Sequencing AI are cloud services.
  Using them means handing your genome to somebody else's inference stack, which
  is the exact thing this project exists to avoid.

  OSGenome2 wired up a local Ollama with no grounding constraints at all, so the
  model is free to invent a variant, invent a risk figure, and say it with the
  same confidence it says true things.

A refusal-first, citation-required, entirely local assistant is unclaimed
territory. This is that.

THE THREAT MODEL, STATED PLAINLY
--------------------------------
Hallucination in a medical context is the worst failure mode available to this
codebase. A wrong magnitude is a bad ranking. A wrong sentence about a variant
somebody does not have, delivered in fluent prose next to their real name, is
something else. So the design inverts the usual default:

  REFUSAL IS THE DEFAULT. Answering is the exception, permitted only when the
  local evidence store contains findings that actually bear on the question.

  EVERY CLAIM CARRIES A CITATION. A response citing an ID that was not in the
  context is rejected outright, because a model that invents a citation has
  already demonstrated it is inventing.

  VALIDATION IS POST-HOC AND FAILS CLOSED. The system prompt is a request, not
  a guarantee. Anything the model returns is re-checked here, and anything that
  fails the check is replaced by the refusal message rather than shown with a
  warning label. A warning label next to a hallucination is still a
  hallucination on screen.

GENOTYPES NEVER LEAVE THIS PROCESS
----------------------------------
Not to the network, and not to the local model either. :func:`_redact` strips
allele1, allele2, genotype, zygosity, copy counts and any other raw call before
anything is assembled into a prompt, and :func:`_leak_check` re-scans the
finished prompt for the genotype strings of the findings it was built from. If
any survived, the request is not sent at all. The model receives finding text,
gene, interpretation and citation IDs, which is everything it needs to summarise
and nothing it needs to re-derive a genotype.

The reasoning is not paranoia about Ollama specifically. It is that a prompt is
a value that gets logged, cached, and occasionally shipped to a bug report, and
the cheapest way to guarantee a genotype is never in one of those is to
guarantee it was never in the prompt.

LOOPBACK ONLY
-------------
The Ollama endpoint is validated to be a loopback address before any request is
made, and access goes through ``external.guard`` first, so the assistant is
unavailable until the user has installed Ollama and accepted its licence, the
same gate every other external tool sits behind.
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, Iterable
from urllib.parse import urlparse

import requests

try:  # package import, with a flat fallback matching backend/traits.py
    from backend import external as _external
except ImportError:  # pragma: no cover - only outside the package
    import external as _external  # type: ignore

external = _external

__all__ = [
    "GROUNDING_CONTRACT", "REFUSAL", "REDACTED_FIELDS",
    "BANNED_ADVICE_PHRASES", "OUT_OF_SCOPE", "STOPWORDS",
    "OLLAMA_DEFAULT_HOST", "OLLAMA_DEFAULT_MODEL", "OLLAMA_ENDPOINT",
    "DEFAULT_TIMEOUT", "MAX_FINDINGS",
    "build_context", "ask", "validate_response", "should_refuse",
    "refusal_payload", "redact", "finding_id", "is_loopback",
    "resolve_host", "resolve_model",
]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OLLAMA_DEFAULT_HOST = "http://127.0.0.1:11434"
OLLAMA_DEFAULT_MODEL = "llama3.1:8b"
OLLAMA_ENDPOINT = "/api/generate"
DEFAULT_TIMEOUT = 120
MAX_FINDINGS = 40

REFUSAL = (
    "I cannot answer that from your data.\n\n"
    "This assistant answers only from the findings DNAInsight has already "
    "produced from your file, and only when those findings actually address the "
    "question. It does not diagnose, does not give prognosis, does not discuss "
    "medication amounts, and does not answer legal or insurance questions. It "
    "also will not speculate about a variant that is not in your results.\n\n"
    "If this is a health question, a licensed clinician, pharmacist or genetic "
    "counsellor is the right place for it. The American Board of Genetic "
    "Counseling maintains a directory at findageneticcounselor.com ."
)

GROUNDING_CONTRACT = """You are a careful assistant inside DNAInsight, a local, offline personal genomics application. You are talking to the person whose DNA file was analysed.

HARD RULES. These are not preferences. A response that breaks any of them is discarded before the user sees it, so breaking one costs the user an answer.

1. ANSWER ONLY FROM THE PROVIDED FINDINGS. The FINDINGS block below is the complete set of facts available to you. If something is not in it, you do not know it. Your own background knowledge about genetics is not evidence here and must not be used to add, extend, correct or contradict a finding.

2. CITE A FINDING ID FOR EVERY CLAIM. Write the ID in square brackets immediately after the claim, like this: [rs1801133]. A sentence that states a fact without a citation is a violation. Use only the IDs listed in the FINDINGS block, exactly as written.

3. IF THE FINDINGS DO NOT SUPPORT AN ANSWER, SAY SO AND STOP. Do not fill the gap. Do not offer a general explanation instead. Say plainly that the findings do not cover the question, and end the response. This is a correct and complete answer, not a failure.

4. NO MEDICAL ADVICE AND NO MEDICATION AMOUNTS. Do not tell the user to begin, continue or change anything. Do not name amounts, schedules or units. Do not diagnose, and do not give a prognosis. Where a finding has clinical relevance, say that it is something to discuss with a prescriber, pharmacist or genetic counsellor.

5. NEVER SPECULATE ABOUT A VARIANT THAT IS NOT IN THE CONTEXT. If the user asks about a gene or an rsID that does not appear in the FINDINGS block, say that it is not in their results. Do not describe what it usually does. Do not say what it would mean if they had it.

6. NEVER STATE A RISK FIGURE THAT IS NOT IN THE CONTEXT. No percentages, no odds ratios, no relative risks, no lifetime risks, unless that exact number appears in a finding you cite. Do not convert, round, combine or estimate a number that is not there.

7. DO NOT DISCUSS RAW GENOTYPES. You have not been given any, deliberately. If the user asks which letters they carry, tell them the assistant is not given that information and point them at the findings table in the application.

STYLE. Plain English, short sentences, no reassurance you cannot support and no alarm you cannot support. Uncertainty stated as uncertainty. If a finding carries a caveat, the caveat goes in the answer, not in a footnote."""

# Fields stripped from every finding before it goes anywhere near a prompt.
# The rule is simple and deliberately over-broad: anything that is, encodes, or
# reconstructs the person's own call at a position. Interpretation stays,
# because interpretation is what the assistant is for.
REDACTED_FIELDS: frozenset[str] = frozenset({
    "allele1", "allele2", "allele_1", "allele_2", "a1", "a2",
    "alleles", "genotype", "genotypes", "genotype_call", "genotype_key",
    "raw", "raw_genotype", "raw_call", "call", "calls", "gt",
    "zygosity", "variant_copies", "copies", "dosage", "diplotype",
    "allele_status", "detected_alleles", "user_genotype", "your_genotype",
})

# Phrasing that must never reach the user, whatever the model produced.
# Deliberately specific: a broad pattern such as "you have" fires on ordinary
# sentences and would turn the validator into a random refusal generator, which
# trains users to ignore it.
BANNED_ADVICE_PHRASES: tuple[str, ...] = (
    "you should take", "you should stop", "you should start",
    "stop taking", "start taking", "you must take", "you need to take",
    "i recommend", "we recommend", "i advise", "my advice",
    "increase your dose", "decrease your dose", "reduce your dose",
    "adjust your dose", "double the dose", "half the dose",
    "milligram", "mg twice", "mg daily", "mg per day", "mg once",
    "prescribe", "you are diagnosed", "you have been diagnosed",
    "diagnose you", "this means you have the disease",
    "no need to see a doctor", "you do not need a doctor",
    "instead of your medication",
)

# Questions this assistant declines by category rather than by judgement call.
# Each entry is (category, reason, patterns).
OUT_OF_SCOPE: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("diagnosis",
     "This asks for a diagnosis. A diagnosis comes from a clinician who can "
     "examine you, not from a genotype file.",
     ("do i have", "do i suffer", "am i sick", "diagnose", "diagnosis",
      "is this cancer", "have i got", "what disease do i", "what illness")),
    ("prognosis",
     "This asks for a prognosis. Nothing in a DNA file can tell you how long "
     "you will live or how a condition will progress.",
     ("how long will i live", "life expectancy", "when will i die",
      "will i get", "am i going to develop", "prognosis", "how long do i have")),
    ("dosing",
     "This asks about medication amounts. Nothing in this application discusses "
     "amounts, schedules or units, and that is a decision for a prescriber.",
     ("what dose", "how much should i take", "how many mg", "dosage",
      "how much of", "should i take", "should i stop taking",
      "should i start taking", "can i stop my", "titrate")),
    ("legal",
     "This is a legal question. This application cannot answer one.",
     ("sue", "lawsuit", "legal advice", "in court", "custody", "my lawyer",
      "is it legal", "paternity test for court")),
    ("insurance",
     "This is an insurance question. Genetic information and insurance is a "
     "matter of law and policy that varies by country and by product, and "
     "getting it wrong here could cost you money.",
     ("insurance", "insurer", "underwriting", "underwrite", "premium",
      "life cover", "health cover", "will my policy")),
)

# Words too common to indicate what a question is about. Matching on them would
# pull the entire findings table into every context block, which is retrieval in
# name only.
STOPWORDS: frozenset[str] = frozenset({
    "about", "after", "again", "against", "and", "any", "anything", "are",
    "because", "been", "before", "being", "between", "both", "but", "can",
    "cannot", "could", "did", "does", "doing", "done", "down", "during", "each",
    "explain", "few", "for", "from", "further", "had", "has", "have", "having",
    "her", "here", "hers", "him", "his", "how", "into", "its", "itself", "just",
    "know", "like", "make", "many", "mean", "means", "might", "more", "most",
    "much", "must", "myself", "need", "not", "now", "off", "once", "only",
    "other", "our", "ours", "out", "over", "own", "please", "result", "results",
    "same", "she", "should", "show", "some", "such", "tell", "than", "that",
    "the", "their", "theirs", "them", "then", "there", "these", "they", "this",
    "those", "through", "too", "under", "until", "very", "was", "way", "were",
    "what", "when", "where", "which", "while", "who", "whom", "why", "will",
    "with", "would", "you", "your", "yours",
})

_RSID_RE = re.compile(r"\brs\d+\b", re.IGNORECASE)
_CITATION_RE = re.compile(r"\[([^\[\]]{1,120})\]")
_WORD_RE = re.compile(r"[a-z0-9][a-z0-9_.\-]*")
_GENOTYPE_TOKEN_RE = re.compile(r"^[ACGTDI]{1,4}$")


# ---------------------------------------------------------------------------
# Identity and redaction
# ---------------------------------------------------------------------------

def finding_id(finding: Any) -> str:
    """Stable citation ID for one finding.

    rsID first, because that is what a user recognises and what a model is most
    likely to type correctly. Falling back to an explicit id, then to the trait
    key, then to a slug of the name. Never empty, because an empty citation ID
    would silently pass the "was this ID in the context" check.
    """
    if not isinstance(finding, dict):
        return ""
    for key in ("id", "finding_id"):
        value = str(finding.get(key) or "").strip()
        if value:
            return value
    rsid = str(finding.get("rsid") or "").strip()
    if rsid:
        return rsid.lower()
    for key in ("key", "gene", "name"):
        value = str(finding.get(key) or "").strip()
        if value:
            return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return "unidentified_finding"


def _genotype_tokens(finding: Any) -> set[str]:
    """Every string in this finding that IS a genotype, for the leak scan."""
    tokens: set[str] = set()
    if not isinstance(finding, dict):
        return tokens
    for key in ("genotype", "genotype_call", "raw_genotype", "call"):
        value = str(finding.get(key) or "").strip().upper()
        if value and _GENOTYPE_TOKEN_RE.match(value):
            tokens.add(value)
    a1 = str(finding.get("allele1") or "").strip().upper()
    a2 = str(finding.get("allele2") or "").strip().upper()
    if a1 and a2 and _GENOTYPE_TOKEN_RE.match(a1 + a2):
        tokens.add(a1 + a2)
        tokens.add(a2 + a1)
    return {t for t in tokens if len(t) >= 2}


def _scrub_text(text: str, tokens: Iterable[str]) -> str:
    """Remove standalone genotype tokens from free text.

    Targeted rather than aggressive: only the exact genotype strings belonging
    to THIS finding are removed, and only where they stand alone. A blanket
    two-letter-uppercase filter would mangle gene names and citations, and a
    mangled citation is worse than a preserved one.
    """
    out = text
    for token in sorted(set(tokens), key=len, reverse=True):
        if not token:
            continue
        pattern = re.compile(
            r"(?<![A-Za-z0-9])" + re.escape(token) + r"(?![A-Za-z0-9])"
        )
        out = pattern.sub("[genotype withheld]", out)
    return out


def _redact(finding: Any) -> dict:
    """Strip every raw call from a finding, leaving only interpretation.

    Removes allele1, allele2, genotype, zygosity, copy counts and every other
    field in REDACTED_FIELDS, at every level of nesting, and then scrubs the
    finding's own genotype strings out of the remaining free text. What survives
    is the rsID, the gene, the summary, the interpretation, the caveat and the
    evidence fields, which is exactly what a grounded answer needs.
    """
    if not isinstance(finding, dict):
        return {}
    tokens = _genotype_tokens(finding)

    def clean(node: Any) -> Any:
        if isinstance(node, dict):
            return {k: clean(v) for k, v in node.items()
                    if str(k).strip().lower() not in REDACTED_FIELDS}
        if isinstance(node, (list, tuple)):
            return [clean(v) for v in node]
        if isinstance(node, str):
            return _scrub_text(node, tokens)
        return node

    out = clean(finding)
    out["id"] = finding_id(finding)
    out["redacted"] = True
    return out


# Public alias. The underscore name is the one the design calls for; this exists
# so callers outside the module are not importing a private name.
redact = _redact


def _leak_check(text: str, findings: Iterable[Any]) -> list[str]:
    """Return genotype strings that survived into ``text``. Empty means clean.

    This is the belt to redaction's braces. Redaction can be defeated by a
    finding that puts its genotype somewhere unexpected, for example inside a
    nested provenance dict added by a later feature. Scanning the finished
    prompt catches that, and the caller refuses to send rather than trimming.
    """
    leaked: list[str] = []
    for finding in (findings or []):
        for token in _genotype_tokens(finding):
            pattern = re.compile(
                r"(?<![A-Za-z0-9])" + re.escape(token) + r"(?![A-Za-z0-9])"
            )
            if pattern.search(text):
                leaked.append(token)
    return sorted(set(leaked))


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def _question_terms(question: Any) -> tuple[set[str], set[str]]:
    """Split a question into (rsIDs mentioned, meaningful lower-case words)."""
    text = str(question or "").lower()
    rsids = {m.group(0).lower() for m in _RSID_RE.finditer(text)}
    words = {
        w for w in (m.group(0) for m in _WORD_RE.finditer(text))
        if len(w) >= 4 and w not in STOPWORDS and not w.startswith("rs")
    }
    return rsids, words


def _finding_terms(finding: dict) -> set[str]:
    """Every term a finding can legitimately be retrieved by."""
    terms: set[str] = set()
    for key in ("gene", "name", "category", "conditions", "key", "condition"):
        value = str(finding.get(key) or "").lower()
        for match in _WORD_RE.finditer(value):
            token = match.group(0)
            if len(token) >= 3 and token not in STOPWORDS:
                terms.add(token)
    for key in ("topics", "medicines", "sources", "drugs"):
        values = finding.get(key) or []
        if isinstance(values, str):
            values = [values]
        for value in values:
            for match in _WORD_RE.finditer(str(value).lower()):
                token = match.group(0)
                if len(token) >= 3 and token not in STOPWORDS:
                    terms.add(token)
    return terms


def _score(finding: dict, rsids: set[str], words: set[str]) -> tuple[int, list[str]]:
    """Relevance of one finding to one question, with the terms that matched.

    Scoring is intentionally strict. A finding with no term in common with the
    question scores zero and is dropped, which means a vague question retrieves
    nothing and the assistant refuses. That is the correct outcome: an assistant
    handed the whole findings table will summarise the whole findings table, and
    a summary nobody asked for is where invention starts.
    """
    matched: list[str] = []
    score = 0
    fid = finding_id(finding).lower()
    rsid = str(finding.get("rsid") or "").strip().lower()
    if rsid and rsid in rsids:
        score += 100
        matched.append(rsid)
    elif fid and fid in rsids:
        score += 100
        matched.append(fid)

    gene = str(finding.get("gene") or "").strip().lower()
    if gene and gene in words:
        score += 40
        matched.append(gene)

    terms = _finding_terms(finding)
    overlap = sorted(terms & words)
    for term in overlap:
        if term != gene:
            score += 12
            matched.append(term)
    return score, sorted(set(matched))


def build_context(findings: Iterable[Any],
                  question: Any,
                  *,
                  max_findings: int = MAX_FINDINGS) -> dict:
    """Select the findings that bear on a question and render the context block.

    Retrieval is restricted to the local evidence store. There is no embedding
    model, no vector database and no external lookup, because none of those is
    needed to match an rsID, a gene name, a topic or a drug name, and each of
    them would add a dependency or a network call this project does not accept.

    Returns a dict carrying the rendered context block AND the exact list of
    finding IDs it contains. The second half matters more than the first: it is
    the allow-list :func:`validate_response` checks citations against, so a
    model that cites anything else is caught.
    """
    rsids, words = _question_terms(question)
    scored: list[tuple[int, float, dict, list[str]]] = []
    for finding in (findings or []):
        if not isinstance(finding, dict):
            continue
        score, matched = _score(finding, rsids, words)
        if score <= 0:
            continue
        try:
            magnitude = float(finding.get("magnitude") or 0.0)
        except (TypeError, ValueError):
            magnitude = 0.0
        scored.append((score, magnitude, finding, matched))

    scored.sort(key=lambda row: (-row[0], -row[1], finding_id(row[2])))
    limit = max(0, int(max_findings))
    selected = scored[:limit]

    included: list[dict] = []
    ids: list[str] = []
    matched_terms: set[str] = set()
    for _score_value, _magnitude, finding, matched in selected:
        block = _redact(finding)
        if not isinstance(block, dict):
            block = {}
        # The ID is recomputed from the ORIGINAL finding rather than trusted from
        # the redacted copy, so a redactor that drops or rewrites it cannot
        # silently produce a context whose allow-list does not match its content.
        block["id"] = finding_id(finding)
        included.append(block)
        ids.append(block["id"])
        matched_terms.update(matched)

    lines: list[str] = []
    for block in included:
        lines.append(f"[{block['id']}]")
        for label, key in (("gene", "gene"), ("name", "name"),
                           ("category", "category"),
                           ("conditions", "conditions"),
                           ("clinical significance", "clinical_sig"),
                           ("evidence", "evidence"),
                           ("CPIC level", "cpic_level"),
                           ("confidence", "confidence"),
                           ("summary", "summary"),
                           ("interpretation", "interpretation"),
                           ("caveat", "caveat")):
            value = block.get(key)
            if isinstance(value, (list, tuple)):
                value = ", ".join(str(v) for v in value)
            text = str(value or "").strip()
            if text:
                lines.append(f"  {label}: {text}")
        lines.append("")

    context = "\n".join(lines).strip()
    return {
        "context": context,
        "finding_ids": ids,
        "included": included,
        "count": len(included),
        "considered": len(scored),
        "truncated": len(scored) > len(selected),
        "max_findings": limit,
        "matched_terms": sorted(matched_terms),
        "question_rsids": sorted(rsids),
    }


# ---------------------------------------------------------------------------
# Refusal
# ---------------------------------------------------------------------------

def _out_of_scope(question: Any) -> tuple[str, str]:
    """Return (category, reason) when a question is out of scope, else ('', '')."""
    text = " " + str(question or "").lower().strip() + " "
    for category, reason, patterns in OUT_OF_SCOPE:
        for pattern in patterns:
            if pattern in text:
                return category, reason
    return "", ""


def _context_ids(context: Any) -> list[str]:
    if isinstance(context, dict):
        return [str(i) for i in (context.get("finding_ids") or [])]
    if isinstance(context, str):
        return [context] if context.strip() else []
    if isinstance(context, (list, tuple)):
        return [str(i) for i in context]
    return []


def should_refuse(question: Any, context: Any) -> bool:
    """True when this question must not be sent to a model at all.

    Refuses when the question is empty, when the context is empty, or when the
    question asks for something out of scope: diagnosis, prognosis, medication
    amounts, or a legal or insurance question. Refusing before the request is
    made, rather than filtering afterwards, means an out-of-scope question never
    reaches the model and therefore never has an answer to leak.
    """
    if not str(question or "").strip():
        return True
    if not _context_ids(context):
        return True
    category, _reason = _out_of_scope(question)
    return bool(category)


def _refusal_reason(question: Any, context: Any) -> str:
    if not str(question or "").strip():
        return "No question was asked."
    if not _context_ids(context):
        return (
            "Nothing in your results addresses this question. The assistant "
            "answers only from findings DNAInsight already produced from your "
            "file, and none of them matched."
        )
    category, reason = _out_of_scope(question)
    if category:
        return reason
    return "Refused."


def refusal_payload(question: Any = "",
                    context: Any = None,
                    *,
                    reason: str = "",
                    category: str = "refused") -> dict:
    """The standard refusal shape, so every refusal path returns the same thing."""
    return {
        "available": True,
        "answer": REFUSAL,
        "refused": True,
        "refusal_category": category,
        "reason": reason or _refusal_reason(question, context),
        "grounded": False,
        "citations": [],
        "context_finding_ids": _context_ids(context),
        "context_count": len(_context_ids(context)),
        "model": None,
        "validation": None,
        "question": str(question or ""),
    }


# ---------------------------------------------------------------------------
# Post-hoc validation
# ---------------------------------------------------------------------------

def _extract_citations(text: str) -> list[str]:
    """Every ID the response appears to cite, bracketed or bare rsID."""
    found: list[str] = []
    for match in _CITATION_RE.finditer(text):
        inner = match.group(1)
        for part in re.split(r"[,;|]", inner):
            token = part.strip().strip("[]").strip()
            if token:
                found.append(token.lower())
    for match in _RSID_RE.finditer(text):
        found.append(match.group(0).lower())
    return sorted(set(found))


_MODEL_REFUSAL_MARKERS = (
    "do not support", "does not support", "not in your results",
    "not in the findings", "no finding", "cannot answer", "can not answer",
    "not covered by", "nothing in your results", "do not cover",
    "does not cover", "insufficient", "not enough information",
)


def validate_response(text: Any, allowed_ids: Iterable[Any]) -> dict:
    """Check a model response before any of it reaches the user. Fails closed.

    Rejects, in this order:

      * an empty or non-string response;
      * a response citing an ID that was not in the context, which is the
        signature of an invented fact wearing a citation;
      * a response containing banned medical-advice or dosing phrasing;
      * a response that makes claims and cites nothing at all.

    On rejection, ``text`` in the returned verdict is the refusal message, NOT
    the model output. Returning the model output alongside a warning would put
    the hallucination on screen, which is the outcome this function exists to
    prevent.

    A response that is itself a refusal, meaning it says the findings do not
    support an answer and cites nothing, is accepted. That is the model doing
    exactly what rule 3 of the grounding contract asks for.
    """
    allowed = {str(i).strip().lower() for i in (allowed_ids or []) if str(i).strip()}
    verdict = {
        "ok": False,
        "verdict": "rejected",
        "text": REFUSAL,
        "model_text": text if isinstance(text, str) else "",
        "cited_ids": [],
        "unknown_ids": [],
        "allowed_ids": sorted(allowed),
        "violations": [],
        "reason": "",
    }

    if not isinstance(text, str) or not text.strip():
        verdict["violations"].append("empty_response")
        verdict["reason"] = (
            "The model returned nothing usable, so the refusal message is shown "
            "instead."
        )
        return verdict

    lowered = text.lower()
    cited = _extract_citations(text)
    verdict["cited_ids"] = cited

    unknown = sorted(c for c in cited if c not in allowed)
    if unknown:
        verdict["unknown_ids"] = unknown
        verdict["violations"].append("unknown_citation")
        verdict["reason"] = (
            "The response cited "
            + ", ".join(unknown)
            + ", which was not in the context it was given. A citation to "
              "something that was never supplied means the content was invented, "
              "so the whole response is discarded."
        )
        return verdict

    banned = [p for p in BANNED_ADVICE_PHRASES if p in lowered]
    if banned:
        verdict["violations"].append("medical_advice")
        verdict["reason"] = (
            "The response contained phrasing this application does not permit ("
            + ", ".join(banned)
            + "), so it was discarded rather than shown with a warning."
        )
        return verdict

    if not cited:
        if any(marker in lowered for marker in _MODEL_REFUSAL_MARKERS):
            verdict.update({
                "ok": True, "verdict": "model_refused", "text": text,
                "reason": "The model declined to answer from the supplied "
                          "findings, which is a permitted and complete response.",
            })
            return verdict
        verdict["violations"].append("citation_required")
        verdict["reason"] = (
            "The response made claims without citing a single finding ID. Every "
            "claim has to be traceable to a finding, so an uncited response is "
            "discarded."
        )
        return verdict

    verdict.update({
        "ok": True, "verdict": "accepted", "text": text,
        "reason": f"Grounded in {len(cited)} cited finding(s), all of which were "
                  f"in the supplied context.",
    })
    return verdict


# ---------------------------------------------------------------------------
# Local model access
# ---------------------------------------------------------------------------

def is_loopback(url: Any) -> bool:
    """True only for a loopback host. Anything else is not a local model.

    An assistant that can be pointed at a remote host by an environment variable
    is a cloud assistant with extra steps, so the host is validated rather than
    trusted, and an unexpected value fails closed.
    """
    try:
        parsed = urlparse(str(url or ""))
    except ValueError:
        return False
    if parsed.scheme not in ("http", "https"):
        return False
    host = (parsed.hostname or "").strip().lower()
    if host in ("localhost", "127.0.0.1", "::1", "0:0:0:0:0:0:0:1"):
        return True
    return host.startswith("127.")


def resolve_host(host: Any = None) -> str:
    """The Ollama base URL, from the argument, then the environment, then default."""
    if host:
        return str(host)
    return os.environ.get("DNAINSIGHT_OLLAMA_HOST") or OLLAMA_DEFAULT_HOST


def resolve_model(model: Any = None) -> str:
    if model:
        return str(model)
    return os.environ.get("DNAINSIGHT_ASSISTANT_MODEL") or OLLAMA_DEFAULT_MODEL


def _build_prompt(question: Any, context_block: str) -> str:
    return (
        "FINDINGS. These are the only facts you may use. Each one begins with "
        "its citation ID in square brackets.\n\n"
        f"{context_block}\n\n"
        "END OF FINDINGS.\n\n"
        f"QUESTION: {str(question or '').strip()}\n\n"
        "Answer using only the findings above, citing an ID in square brackets "
        "for every claim. If they do not answer the question, say so and stop."
    )


def ask(question: Any,
        findings: Iterable[Any],
        *,
        model: Any = None,
        host: Any = None,
        timeout: int = DEFAULT_TIMEOUT,
        max_findings: int = MAX_FINDINGS) -> dict:
    """Answer a question from local findings using a local model, or refuse.

    The order of operations is the security design and is not negotiable:

      1. ``external.guard("ollama", "assistant")``. Ollama is an external tool
         like every other one, so it is unavailable until installed and
         licence-accepted, and absence degrades to the standard payload rather
         than raising.
      2. Build the context from the local evidence store only.
      3. Refuse if the context is empty or the question is out of scope. Nothing
         is sent in that case.
      4. Redact, assemble the prompt, then RE-SCAN the finished prompt for
         genotype strings. If any survived, do not send.
      5. POST to the loopback Ollama endpoint.
      6. Validate the response post-hoc and fail closed.

    Returns a dict. ``answer`` is either a validated grounded response or the
    refusal message, and never anything in between.
    """
    blocked = external.guard("ollama", "assistant")
    if blocked is not None:
        payload = dict(blocked)
        payload.update({
            "answer": REFUSAL,
            "refused": True,
            "refusal_category": "tool_unavailable",
            "grounded": False,
            "citations": [],
            "context_finding_ids": [],
            "context_count": 0,
            "model": None,
            "validation": None,
            "question": str(question or ""),
            "assistant_note": (
                "The local assistant needs Ollama installed on this machine and "
                "its licence accepted. Nothing was sent anywhere. No cloud "
                "fallback exists, deliberately: the alternative to a local model "
                "is no model, not somebody else's."
            ),
        })
        return payload

    context = build_context(findings, question, max_findings=max_findings)

    if should_refuse(question, context):
        return refusal_payload(question, context)

    prompt = _build_prompt(question, context["context"])

    leaked = _leak_check(prompt, findings)
    if leaked:
        return {
            **refusal_payload(question, context,
                              reason=(
                                  "A genotype string survived redaction and "
                                  "reached the prompt, so the request was not "
                                  "sent. This is a bug in DNAInsight, not a "
                                  "problem with your question, and it fails "
                                  "closed on purpose."),
                              category="redaction_failure"),
            "leaked_tokens": leaked,
        }

    base = resolve_host(host)
    if not is_loopback(base):
        return refusal_payload(
            question, context,
            reason=(f"The configured model host {base!r} is not a loopback "
                    f"address. This assistant talks to a model on this machine "
                    f"and nowhere else, so the request was not sent."),
            category="non_loopback_host")

    chosen = resolve_model(model)
    url = base.rstrip("/") + OLLAMA_ENDPOINT
    body = {
        "model": chosen,
        "system": GROUNDING_CONTRACT,
        "prompt": prompt,
        "stream": False,
        # Temperature zero because creativity is the failure mode here.
        "options": {"temperature": 0.0},
    }

    try:
        response = requests.post(url, json=body, timeout=timeout)
        status = getattr(response, "status_code", 0)
        if status and int(status) >= 400:
            raise ValueError(f"model returned HTTP {status}")
        data = response.json()
    except Exception as exc:                     # noqa: BLE001 - fail closed
        return refusal_payload(
            question, context,
            reason=(f"The local model could not be reached or returned "
                    f"something unusable ({type(exc).__name__}). Nothing is "
                    f"shown rather than something unverified."),
            category="model_error")

    if not isinstance(data, dict):
        return refusal_payload(
            question, context,
            reason="The local model returned a response this application could "
                   "not parse, so nothing is shown.",
            category="malformed_response")

    text = data.get("response")
    if not isinstance(text, str):
        return refusal_payload(
            question, context,
            reason="The local model's reply contained no response text, so "
                   "nothing is shown.",
            category="malformed_response")

    verdict = validate_response(text, context["finding_ids"])
    return {
        "available": True,
        "answer": verdict["text"],
        "refused": not verdict["ok"],
        "refusal_category": "" if verdict["ok"] else "validation_failed",
        "reason": verdict["reason"],
        "grounded": bool(verdict["ok"]),
        "citations": verdict["cited_ids"] if verdict["ok"] else [],
        "context_finding_ids": context["finding_ids"],
        "context_count": context["count"],
        "context_truncated": context["truncated"],
        "model": chosen,
        "host": base,
        "validation": verdict,
        "question": str(question or ""),
        "genotypes_sent": False,
        "privacy_note": (
            "The prompt contained finding text, gene names, interpretations and "
            "citation IDs. It contained no alleles, no genotypes and no "
            "zygosity, and it was sent to a model on this machine over loopback "
            "only."
        ),
    }
