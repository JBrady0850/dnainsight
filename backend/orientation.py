"""Strand orientation helpers for matching consumer genotypes to SNPedia.

The strand problem, in plain terms:

DNA is double stranded. Every base on one strand is paired with its
complement on the other strand (A with T, C with G). A single SNP can
therefore be described in two equally correct ways depending on which
strand you chose to read. A variant written as (C;T) on one strand is
written as (G;A) on the other. Nothing about the person changed, only the
reporting convention.

Consumer raw data files from 23andMe, AncestryDNA, MyHeritage,
FamilyTreeDNA and LivingDNA all report genotypes on the PLUS strand of
the GRCh37 human reference assembly. That gives us one fixed, known
convention for the incoming side.

SNPedia is the other side, and it is not always plus. Each SNPedia SNP
page carries a field called StabilizedOrientation whose value is the
lowercase string "plus" or "minus". That field, and only that field,
governs how the genotype page titles are written (for example
Rs1801133(C;T)). SNPedia also carries a separate field called
Orientation, which tracks the orientation on the current build (GRCh38)
and can legitimately disagree with StabilizedOrientation. When matching a
consumer genotype against a genotype page title we must match
StabilizedOrientation, falling back to Orientation only when
StabilizedOrientation is missing.

So the rule is short: when StabilizedOrientation is "minus", complement
both alleles before looking up the genotype. When it is "plus", or when
we simply do not know, pass the alleles through untouched.

The unavoidable catch is the palindromic SNPs. For a heterozygous A/T or
C/G call, complementing gives back the very same unordered pair, so there
is no way to tell an unflipped read from a flipped one. Those cases are
flagged as ambiguous rather than quietly flipped and trusted, so that
downstream reporting can warn instead of asserting.
"""

from typing import Any, Dict, FrozenSet, List, Optional, Tuple

COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
    "-": "-",
    "D": "D",
    "I": "I",
    "0": "N",
    "": "N",
}

NOCALL_ALLELES: FrozenSet[str] = frozenset({"", "N", "-", "0", "D", "I"})

AMBIGUOUS_PAIRS: FrozenSet[FrozenSet[str]] = frozenset(
    {frozenset({"A", "T"}), frozenset({"C", "G"})}
)

_PLUS_TOKENS: FrozenSet[str] = frozenset({"plus", "+", "1", "+1", "fwd", "forward"})
_MINUS_TOKENS: FrozenSet[str] = frozenset({"minus", "-", "-1", "rev", "reverse"})

_NOCALL_TOKEN: str = "(-;-)"


def complement_allele(allele: str) -> str:
    """Return the complement of a single allele token.

    Case insensitive, always returns uppercase. Indel and no-call tokens
    map to themselves (a no-call has no strand). Anything unrecognised
    becomes "N" rather than raising, because consumer files do contain
    odd tokens and a report should degrade instead of crashing.
    """
    if allele is None:
        return "N"
    key = str(allele).strip().upper()
    return COMPLEMENT.get(key, "N")


def complement_genotype(a1: str, a2: str) -> Tuple[str, str]:
    """Complement both alleles of a genotype, preserving their order."""
    return complement_allele(a1), complement_allele(a2)


def is_no_call(allele: str) -> bool:
    """True when the allele token carries no usable base call."""
    if allele is None:
        return True
    return str(allele).strip().upper() in NOCALL_ALLELES


def is_ambiguous_pair(a1: str, a2: str) -> bool:
    """True when the pair is an irreducibly palindromic heterozygote.

    Only heterozygous A/T and C/G calls are ambiguous: complementing them
    yields the same unordered pair, so the strand cannot be recovered. A
    homozygous A;A is not ambiguous in this sense, because A;A and T;T are
    distinguishable. No-calls are never ambiguous.
    """
    if is_no_call(a1) or is_no_call(a2):
        return False
    first = str(a1).strip().upper()
    second = str(a2).strip().upper()
    if first == second:
        return False
    return frozenset({first, second}) in AMBIGUOUS_PAIRS


def normalize_orientation(value: Any) -> str:
    """Normalise an orientation value to "plus", "minus" or "" (unknown).

    Accepts the spellings that turn up in SNPedia dumps and in hand
    maintained overrides: "plus", "PLUS", "+", "1", "minus", "-", "-1",
    None and the empty string. Anything unrecognised is reported as
    unknown rather than guessed at.
    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return ""
    if isinstance(value, int):
        token = str(value)
    else:
        token = str(value).strip().lower()
    if not token:
        return ""
    if token in _PLUS_TOKENS:
        return "plus"
    if token in _MINUS_TOKENS:
        return "minus"
    return ""


def _clean_allele(allele: str) -> str:
    """Internal: uppercase an allele, collapsing every no-call token to "N"."""
    if is_no_call(allele):
        return "N"
    key = str(allele).strip().upper()
    if key not in COMPLEMENT:
        return "N"
    return key


def sort_alleles(a1: str, a2: str) -> Tuple[str, str]:
    """Return the pair in SNPedia title order (alphabetical, no-calls last).

    SNPedia writes genotype titles with the alleles sorted alphabetically,
    so a raw (T;C) call is titled (C;T). No-call tokens are normalised to
    "N" and sorted to the end so a partial call reads as (A;N).
    """
    first = _clean_allele(a1)
    second = _clean_allele(a2)
    pair = [first, second]
    pair.sort(key=lambda base: (base == "N", base))
    return pair[0], pair[1]


def genotype_token(a1: str, a2: str) -> str:
    """Return the SNPedia style genotype token, for example "(C;T)".

    Alleles are put in SNPedia title order first. A pair with no usable
    call at all is rendered as the conventional no-call token "(-;-)".
    """
    first, second = sort_alleles(a1, a2)
    if first == "N" and second == "N":
        return _NOCALL_TOKEN
    return "({0};{1})".format(first, second)


def orient_to_snpedia(
    a1: str,
    a2: str,
    stabilized_orientation: Optional[str] = None,
    orientation: Optional[str] = None,
) -> Dict[str, Any]:
    """Put a plus strand consumer genotype into SNPedia title orientation.

    The incoming pair is assumed to be GRCh37 plus strand, which is what
    23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA and LivingDNA export.
    If the SNP's StabilizedOrientation is "minus", both alleles are
    complemented before sorting, because that is the strand SNPedia used
    when it named the genotype pages. A "plus" or unknown orientation is
    passed through unflipped. The `orientation` argument (the GRCh38
    field) is consulted only when StabilizedOrientation is absent.

    Returns a dict with the keys: allele1, allele2, genotype, token,
    flipped, ambiguous, orientation, source_genotype.
    """
    raw1 = _clean_allele(a1)
    raw2 = _clean_allele(a2)

    used = normalize_orientation(stabilized_orientation)
    if stabilized_orientation is None or not str(stabilized_orientation).strip():
        used = normalize_orientation(orientation)

    ambiguous = is_ambiguous_pair(a1, a2)
    has_real_base = raw1 != "N" or raw2 != "N"

    if used == "minus" and has_real_base:
        out1 = "N" if raw1 == "N" else complement_allele(raw1)
        out2 = "N" if raw2 == "N" else complement_allele(raw2)
        flipped = True
    else:
        out1, out2 = raw1, raw2
        flipped = False

    allele1, allele2 = sort_alleles(out1, out2)
    return {
        "allele1": allele1,
        "allele2": allele2,
        "genotype": allele1 + allele2,
        "token": genotype_token(out1, out2),
        "flipped": flipped,
        "ambiguous": ambiguous,
        "orientation": used,
        "source_genotype": raw1 + raw2,
    }


def candidate_tokens(a1: str, a2: str) -> List[str]:
    """Return the tokens worth trying when orientation metadata is missing.

    The unflipped token comes first (plus strand is the more common case),
    then the complemented token. Duplicates are dropped, which is what
    happens for palindromic pairs such as (A;T). Tokens containing a
    no-call are never returned, since there is nothing to look up.
    """
    tokens: List[str] = []
    plus_pair = sort_alleles(a1, a2)
    minus_pair = sort_alleles(*complement_genotype(a1, a2))
    for first, second in (plus_pair, minus_pair):
        if first == "N" or second == "N":
            continue
        token = "({0};{1})".format(first, second)
        if token not in tokens:
            tokens.append(token)
    return tokens


def flip_report(findings: List[Dict[str, Any]]) -> Dict[str, int]:
    """Aggregate a QC summary of strand handling across many findings.

    Counts how many findings were complemented, how many sit on an
    irreducibly palindromic pair, and how many had no usable orientation
    metadata at all (a falsy "orientation" value, including a missing key).
    """
    total = 0
    flipped = 0
    ambiguous = 0
    unknown = 0
    for finding in findings or []:
        total += 1
        if not isinstance(finding, dict):
            unknown += 1
            continue
        if finding.get("flipped"):
            flipped += 1
        if finding.get("ambiguous"):
            ambiguous += 1
        if not finding.get("orientation"):
            unknown += 1
    return {
        "total": total,
        "flipped": flipped,
        "ambiguous": ambiguous,
        "unknown_orientation": unknown,
    }
