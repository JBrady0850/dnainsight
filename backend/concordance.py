"""
concordance.py -- how much your own kits agree with each other, and where they do not.

Derived, read-only and additive. Nothing here mutates the merged set, and no
result of this module ever changes which genotype the rest of the application
uses. merge.py already pools kits, already detects disagreement and already
refuses to reconcile it. What it never said was WHICH two files disagreed, and
out of how many positions they could have disagreed over.

THE TRAP THIS MODULE EXISTS TO AVOID
------------------------------------
Most apparent vendor disagreement is not error, it is strand orientation. One
company reports a SNP on the plus strand and another reports it on the minus
strand, so the same person reads AA in one file and TT in the other. Nothing is
wrong with either file and nothing is wrong with the person.

Publishing that as a vendor error is a false accusation about a named company,
expressed as a number, which is the most credible form a false accusation can
take. So every disagreement is CLASSIFIED before it is COUNTED:

  orientation_artifact  one call is the exact complement of the other, and
                        neither call is an irreducible heterozygote, so the
                        complement explains the difference completely
  indeterminate         one of the calls sits on a palindromic A/T or C/G
                        heterozygote, where a strand flip and a real difference
                        are indistinguishable, so no verdict is available
  genuine               neither of the above

The indeterminate bucket is never folded into either neighbour. Folding it into
genuine overstates vendor disagreement. Folding it into artifact hides real
disagreement. It is reported as its own number and the caller may decide what to
do with it, which is the same treatment "not testable on your array" gets in
genosets.py and "strand ambiguous" gets in scoring.py.

FIVE RULES THE IMPLEMENTATION OBEYS
-----------------------------------
1. All strand logic is delegated to orientation.py, which already owns the
   complement table and the ambiguity rule. A second copy would drift, and drift
   in strand handling produces results that are internally consistent and
   externally backwards.

2. Genotypes are validated against ACGT specifically, NOT against the keys of
   orientation.COMPLEMENT. That table deliberately tolerates no-call and indel
   symbols so that flipping a whole file does not destroy them. Inheriting that
   tolerance here would let "NN" become evidence that two companies disagree.

3. No rate is ever reported without its denominator. Every rate travels with
   ``shared``. A pair with nothing in common reports ``comparable: False`` and a
   rate of None, never 0.0 and never 100.0, because "they never agreed" and "we
   never compared them" are different claims.

4. ``findings_covered`` is None when no findings were supplied. None means
   nobody asked. Zero means someone asked and the answer was none.

5. A single kit returns ``available: False`` with ``not_attempted: True``. One
   kit is an absent comparison, not a failed one. Same degradation contract as
   every other v3 capability.

6. Same-company pairs are compared and flagged with ``same_provider``, never
   dropped. Two kits from one vendor years apart ran on different chips, and
   their agreement with each other is a real number worth seeing.

KNOWN PROPERTY, STATED RATHER THAN HIDDEN
-----------------------------------------
A palindromic heterozygote reads the same on either strand, so two kits will
always agree there whatever their orientations. Those positions therefore count
towards agreement, and the agreement rate is very slightly optimistic as a
result. They are left in because excluding positions from a denominator for
being too easy is its own distortion, and because the indeterminate count next
to the rate already tells the reader that palindromic sites are in play.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, FrozenSet, List, Optional, Tuple

from . import orientation as _orientation

# The four real bases, and nothing else. See rule 2 above: this is deliberately
# NOT orientation.COMPLEMENT.keys().
BASES: FrozenSet[str] = frozenset({"A", "C", "G", "T"})

AGREEMENT = "agreement"
GENUINE = "genuine"
ORIENTATION_ARTIFACT = "orientation_artifact"
INDETERMINATE = "indeterminate"
NOT_COMPARABLE = "not_comparable"

# The three buckets that a disagreement can land in. Agreement is not a conflict
# and not_comparable is not a disagreement, so neither belongs here, and the sum
# of these three is always the conflict total.
CONFLICT_CLASSES: Tuple[str, str, str] = (GENUINE, ORIENTATION_ARTIFACT, INDETERMINATE)

_REASONS = {
    AGREEMENT: "Both kits report the same genotype.",
    ORIENTATION_ARTIFACT: (
        "One call is the exact complement of the other. This is a strand "
        "reporting difference between the two files, not a disagreement about "
        "the person."
    ),
    INDETERMINATE: (
        "At least one call is a palindromic A/T or C/G heterozygote, which "
        "reads the same on either strand. A strand flip and a real difference "
        "cannot be told apart here, so no verdict is available."
    ),
    GENUINE: (
        "The two calls differ and no strand explanation accounts for it. This "
        "is a real disagreement between the two files."
    ),
    NOT_COMPARABLE: (
        "At least one side is not a pair of ACGT base calls, so there is "
        "nothing to compare. A no-call is not evidence of disagreement."
    ),
}

NOTE = (
    "Disagreement between two of your files is usually a strand reporting "
    "difference, not an error by either company. Those are counted separately "
    "as orientation artifacts. Positions where the two cannot be told apart are "
    "counted separately again as indeterminate, and are never added to the "
    "genuine total."
)

UNAVAILABLE_REASON = (
    "Cross-vendor concordance needs at least two DNA files loaded into this "
    "profile. With one kit there is nothing to compare it against, which is an "
    "absent comparison rather than a failed one."
)


# ---------------------------------------------------------------------------
# Genotype validation and strand helpers
# ---------------------------------------------------------------------------

def _bases(genotype: Any) -> Optional[Tuple[str, str]]:
    """Return the two ACGT alleles of a genotype, or None if it is not one.

    Rule 2 lives here. Anything that is not exactly two real base calls is
    refused: no-calls, indel tokens, single alleles, empty strings, None and
    anything that is not a string at all.
    """
    if genotype is None:
        return None
    text = str(genotype).strip().upper()
    if len(text) != 2:
        return None
    if text[0] not in BASES or text[1] not in BASES:
        return None
    return text[0], text[1]


def canonical_genotype(genotype: Any) -> Optional[str]:
    """Return the genotype as a sorted two-base string, or None if unparsable.

    "GA" and "AG" are the same unordered call and both come back as "AG", which
    is the same key merge.genotype_key produces, so the two never disagree.
    """
    pair = _bases(genotype)
    if pair is None:
        return None
    first, second = _orientation.sort_alleles(pair[0], pair[1])
    if first not in BASES or second not in BASES:
        return None
    return first + second


def complement_genotype(genotype: Any) -> Optional[str]:
    """Return the opposite-strand reading of a genotype, or None if unparsable.

    Delegates the actual complement to orientation.py (rule 1) and validates the
    input against ACGT first (rule 2), so an indel or no-call token can never
    survive the round trip and re-enter the comparison as a real call.
    """
    pair = _bases(genotype)
    if pair is None:
        return None
    flipped = _orientation.complement_genotype(pair[0], pair[1])
    first, second = _orientation.sort_alleles(flipped[0], flipped[1])
    if first not in BASES or second not in BASES:
        return None
    return first + second


def is_palindromic(genotype: Any) -> bool:
    """True when the genotype is an irreducibly palindromic heterozygote.

    That is A/T or C/G in either order, and only heterozygous. A homozygous AA
    is not palindromic in the sense that matters here, because AA and TT are
    distinguishable from each other and a flip between them is detectable. The
    rule itself belongs to orientation.is_ambiguous_pair and is not restated.
    """
    pair = _bases(genotype)
    if pair is None:
        return False
    return bool(_orientation.is_ambiguous_pair(pair[0], pair[1]))


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_conflict(genotype_a: Any, genotype_b: Any) -> Dict[str, Any]:
    """Classify one position in one pair of kits.

    Returns ``{"a", "b", "classification", "reason"}`` where classification is
    exactly one of AGREEMENT, ORIENTATION_ARTIFACT, INDETERMINATE, GENUINE or
    NOT_COMPARABLE.

    The ORDER of the tests is the whole design:

      * unparsable first, so a no-call can never reach a verdict at all
      * agreement next, because two identical calls are not a conflict whatever
        strand they were read on
      * palindromic next, so a position where the strand cannot be recovered is
        never resolved into either of the two remaining answers
      * complement next, which explains the classic whole-file strand flip
      * genuine last, and only once every other explanation has been ruled out

    Genuine is the residue, never the default. That is deliberate: the bucket
    that becomes a public claim about a named company is the one nothing else
    accounted for.
    """
    left = canonical_genotype(genotype_a)
    right = canonical_genotype(genotype_b)

    if left is None or right is None:
        verdict = NOT_COMPARABLE
    elif left == right:
        verdict = AGREEMENT
    elif is_palindromic(left) or is_palindromic(right):
        verdict = INDETERMINATE
    elif complement_genotype(left) == right:
        verdict = ORIENTATION_ARTIFACT
    else:
        verdict = GENUINE

    return {
        "a": left if left is not None else _clean(genotype_a),
        "b": right if right is not None else _clean(genotype_b),
        "classification": verdict,
        "reason": _REASONS[verdict],
    }


def _clean(value: Any) -> str:
    """Uppercase a genotype for display without pretending it parsed."""
    if value is None:
        return ""
    return str(value).strip().upper()


# ---------------------------------------------------------------------------
# Kits, providers and coverage
# ---------------------------------------------------------------------------

def provider_by_label(merged: Any) -> Dict[str, str]:
    """Map every source label in the merge to the provider it declared.

    An undeclared provider maps to the empty string rather than to a guess. The
    label is often a filename and a filename is not a company, so inferring one
    from it would put a vendor's name on a number they may have had no part in.
    """
    out: Dict[str, str] = {}
    for src in ((merged or {}).get("sources") or []):
        if not isinstance(src, dict):
            continue
        label = str(src.get("label") or "").strip()
        if not label:
            continue
        out[label] = str(src.get("provider") or "").strip()
    return out


def _primary_labels(merged: Any) -> List[str]:
    """The user's own pooled kits, in a stable order, deduplicated.

    Falls back to the source list and then to the labels seen on calls, so a
    hand-built merge result without "primary_labels" still works.
    """
    merged = merged or {}
    labels: List[str] = []

    for label in (merged.get("primary_labels") or []):
        text = str(label or "").strip()
        if text and text not in labels:
            labels.append(text)
    if labels:
        return labels

    for src in (merged.get("sources") or []):
        if isinstance(src, dict) and str(src.get("role") or "") == "self":
            text = str(src.get("label") or "").strip()
            if text and text not in labels:
                labels.append(text)
    if labels:
        return labels

    for entry in (merged.get("genotypes") or {}).values():
        for call in (entry or {}).get("calls") or []:
            text = str((call or {}).get("label") or "").strip()
            if text and text not in labels:
                labels.append(text)
    return labels


def _call_genotype(call: Any) -> Optional[str]:
    """The canonical genotype of one call row, or None when it is not usable."""
    if not isinstance(call, dict):
        return None
    genotype = call.get("genotype")
    if genotype is None:
        genotype = f"{call.get('allele1') or ''}{call.get('allele2') or ''}"
    return canonical_genotype(genotype)


def _calls_by_label(entry: Any, wanted: FrozenSet[str]) -> Dict[str, str]:
    """Parsable calls at one position, keyed by label, restricted to ``wanted``.

    The first call for a label wins, which matches merge.dedupe_snps: a repeated
    rsID inside one file is a duplicated probe, not a disagreement.
    """
    out: Dict[str, str] = {}
    for call in ((entry or {}).get("calls") or []):
        if not isinstance(call, dict):
            continue
        label = str(call.get("label") or "").strip()
        if label not in wanted or label in out:
            continue
        genotype = _call_genotype(call)
        if genotype is None:
            continue
        out[label] = genotype
    return out


def coverage_by_provider(merged: Any) -> List[Dict[str, Any]]:
    """How many positions each provider group actually called.

    Kits that declared the same provider share a group, because two 23andMe
    exports are two 23andMe exports. Kits that declared NO provider each get
    their own group: an absent provider is not a company, and pooling every
    undeclared kit into one bucket would invent a vendor that agreed with itself.
    """
    merged = merged or {}
    labels = _primary_labels(merged)
    wanted = frozenset(labels)
    providers = provider_by_label(merged)

    by_label: Dict[str, int] = {label: 0 for label in labels}
    group_positions: Dict[str, int] = {}

    def key_for(label: str) -> str:
        provider = providers.get(label, "")
        return provider if provider else f"unknown:{label}"

    for entry in (merged.get("genotypes") or {}).values():
        calls = _calls_by_label(entry, wanted)
        if not calls:
            continue
        seen_groups = set()
        for label in calls:
            by_label[label] = by_label.get(label, 0) + 1
            seen_groups.add(key_for(label))
        for group in seen_groups:
            group_positions[group] = group_positions.get(group, 0) + 1

    groups: Dict[str, Dict[str, Any]] = {}
    for label in labels:
        provider = providers.get(label, "")
        key = key_for(label)
        group = groups.setdefault(key, {
            "key": key,
            "provider": provider,
            "declared": bool(provider),
            "labels": [],
            "positions": 0,
            "by_label": {},
        })
        group["labels"].append(label)
        group["by_label"][label] = by_label.get(label, 0)

    for key, group in groups.items():
        group["labels"] = sorted(group["labels"])
        group["positions"] = group_positions.get(key, 0)

    return [groups[key] for key in sorted(groups)]


# ---------------------------------------------------------------------------
# The pair matrix
# ---------------------------------------------------------------------------

def _rate(count: int, shared: int) -> Optional[float]:
    """A rate, or None when there is no denominator to compute it over.

    Rule 3. Returning 0.0 for a pair that was never compared would read as
    "these two never disagreed", which is a claim nobody is entitled to make.
    """
    if shared <= 0:
        return None
    return round(count / shared, 6)


def _new_pair(a: str, b: str, providers: Dict[str, str]) -> Dict[str, Any]:
    provider_a = providers.get(a, "")
    provider_b = providers.get(b, "")
    return {
        "a": a,
        "b": b,
        "provider_a": provider_a,
        "provider_b": provider_b,
        "same_provider": bool(
            provider_a and provider_b and provider_a.casefold() == provider_b.casefold()
        ),
        "shared": 0,
        "agree": 0,
        "conflicts": 0,
        GENUINE: 0,
        ORIENTATION_ARTIFACT: 0,
        INDETERMINATE: 0,
        "not_comparable": 0,
        "comparable": False,
        "agreement_rate": None,
        "conflict_rate": None,
        "genuine_conflict_rate": None,
    }


def pair_matrix(merged: Any) -> List[Dict[str, Any]]:
    """Every pair of the user's own kits, compared position by position.

    Only pooled ("self") kits are compared. A parent's kit disagreeing with the
    user's is not a vendor question, it is inheritance, and relatedness.py owns
    that. Comparison kits are therefore excluded here rather than silently
    inflating the denominators.

    Same-provider pairs are included and flagged (rule 6). Every pair is
    returned even when it shares nothing, carrying ``comparable: False`` and
    None rates, because a pair that could not be compared is a fact about the
    data the user should see.
    """
    merged = merged or {}
    labels = _primary_labels(merged)
    if len(labels) < 2:
        return []

    wanted = frozenset(labels)
    providers = provider_by_label(merged)
    ordered = sorted(labels)
    pairs: Dict[Tuple[str, str], Dict[str, Any]] = {
        (a, b): _new_pair(a, b, providers) for a, b in combinations(ordered, 2)
    }

    for entry in (merged.get("genotypes") or {}).values():
        calls = _calls_by_label(entry, wanted)
        raw = _raw_labels(entry, wanted)
        for (a, b), bucket in pairs.items():
            left, right = calls.get(a), calls.get(b)
            if left is None or right is None:
                # Both files touched the position but at least one call is not
                # a usable pair of bases. That is a coverage fact, never a
                # disagreement (rule 2).
                if a in raw and b in raw:
                    bucket["not_comparable"] += 1
                continue
            verdict = classify_conflict(left, right)["classification"]
            bucket["shared"] += 1
            if verdict == AGREEMENT:
                bucket["agree"] += 1
            else:
                bucket[verdict] += 1

    out: List[Dict[str, Any]] = []
    for (a, b) in combinations(ordered, 2):
        bucket = pairs[(a, b)]
        shared = bucket["shared"]
        bucket["conflicts"] = (
            bucket[GENUINE] + bucket[ORIENTATION_ARTIFACT] + bucket[INDETERMINATE]
        )
        bucket["comparable"] = shared > 0
        bucket["agreement_rate"] = _rate(bucket["agree"], shared)
        bucket["conflict_rate"] = _rate(bucket["conflicts"], shared)
        bucket["genuine_conflict_rate"] = _rate(bucket[GENUINE], shared)
        out.append(bucket)
    return out


def _raw_labels(entry: Any, wanted: FrozenSet[str]) -> FrozenSet[str]:
    """Labels that produced a call row at this position, parsable or not."""
    found = set()
    for call in ((entry or {}).get("calls") or []):
        if not isinstance(call, dict):
            continue
        label = str(call.get("label") or "").strip()
        if label in wanted:
            found.add(label)
    return frozenset(found)


# ---------------------------------------------------------------------------
# Findings coverage
# ---------------------------------------------------------------------------

def _finding_rsids(findings: Any) -> List[str]:
    out: List[str] = []
    for finding in (findings or []):
        if isinstance(finding, dict):
            rsid = str(finding.get("rsid") or "").strip().lower()
        else:
            rsid = str(finding or "").strip().lower()
        if rsid:
            out.append(rsid)
    return out


def _findings_summary(merged: Any, findings: Any, wanted: FrozenSet[str]) -> Dict[str, Optional[int]]:
    """Counts of findings that two or more of the user's kits actually read.

    None throughout when ``findings`` is None (rule 4). The distinction is the
    same one the degradation contract draws everywhere else in v3: not asked and
    asked-with-a-zero-answer are different states and must not share a value.
    """
    if findings is None:
        return {"findings_total": None, "findings_covered": None, "findings_conflicted": None}

    genotypes = (merged or {}).get("genotypes") or {}
    rsids = _finding_rsids(findings)
    covered = 0
    conflicted = 0
    for rsid in rsids:
        calls = _calls_by_label(genotypes.get(rsid), wanted)
        if len(calls) < 2:
            continue
        covered += 1
        values = list(calls.values())
        if any(
            classify_conflict(x, y)["classification"] == GENUINE
            for x, y in combinations(values, 2)
        ):
            conflicted += 1
    return {
        "findings_total": len(rsids),
        "findings_covered": covered,
        "findings_conflicted": conflicted,
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def analyse(merged: Any, findings: Any = None) -> Dict[str, Any]:
    """Cross-vendor concordance across the user's own kits.

    ``merged`` is a merge.merge_sources result and is never modified. ``findings``
    is optional; pass None to say nobody asked about findings, and a list to ask.

    With fewer than two pooled kits the result is ``available: False`` with
    ``not_attempted: True`` and ``totals: None`` (rule 5). Totals of zero would
    claim a comparison was run and found nothing, which did not happen.
    """
    merged = merged or {}
    labels = _primary_labels(merged)
    wanted = frozenset(labels)
    coverage = coverage_by_provider(merged)
    findings_summary = _findings_summary(merged, findings, wanted)

    if len(labels) < 2:
        result: Dict[str, Any] = {
            "available": False,
            "not_attempted": True,
            "reason": UNAVAILABLE_REASON,
            "kits": len(labels),
            "labels": list(labels),
            "coverage": coverage,
            "pairs": [],
            "totals": None,
            "note": NOTE,
        }
        result.update(findings_summary)
        return result

    pairs = pair_matrix(merged)
    keys = ("shared", "agree", "conflicts", GENUINE, ORIENTATION_ARTIFACT,
            INDETERMINATE, "not_comparable")
    totals = {key: sum(pair[key] for pair in pairs) for key in keys}
    totals["comparable_pairs"] = sum(1 for pair in pairs if pair["comparable"])
    totals["pairs"] = len(pairs)
    totals["agreement_rate"] = _rate(totals["agree"], totals["shared"])
    totals["conflict_rate"] = _rate(totals["conflicts"], totals["shared"])
    totals["genuine_conflict_rate"] = _rate(totals[GENUINE], totals["shared"])

    result = {
        "available": True,
        "kits": len(labels),
        "labels": list(labels),
        "coverage": coverage,
        "pairs": pairs,
        "totals": totals,
        "note": NOTE,
    }
    result.update(findings_summary)
    return result
