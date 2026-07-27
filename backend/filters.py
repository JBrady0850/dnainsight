"""
filters.py -- server-side filtering, sorting and faceting for the findings view.

Pure functions over a list of finding dicts. No Flask, no database, no network,
so every rule here is directly testable and the HTTP layer stays thin.

Implements docs/API_V2.md sections 3.4, 3.5 and 3.6.

TWO RULES THAT ARE EASY TO GET WRONG, AND WHY THEY MATTER
---------------------------------------------------------
1. Frequency and publication filters MUST NOT apply to genosets, traits or
   polygenic scores. Those entities have no single position, so they have no
   frequency and no citation count. A naive filter drops them the moment the
   user nudges a slider off zero, and the user concludes the feature is broken.
   The reference product carried the same exemption for the same reason.

2. A null magnitude sorts and filters as 1, not as 0. An unscored variant is
   "nobody has assessed this", which belongs just above "assessed and boring"
   and well below "interesting". Treating it as 0 buries genuinely unknown
   variants; treating it as 10 floods the top of the list.
"""

from __future__ import annotations

import re
from typing import Any, Callable, Iterable

__all__ = [
    "SORT_KEYS", "EXEMPT_ENTITIES", "UNSCORED_MAGNITUDE",
    "parse_csv", "parse_bool", "parse_float", "parse_int",
    "parse_query", "matches_query", "apply_filters", "sort_findings",
    "build_facets", "summarise", "paginate", "filter_and_sort",
]

# Entity types exempt from frequency and publication filtering. See rule 1.
EXEMPT_ENTITIES: frozenset = frozenset({"genoset", "trait", "prs"})

UNSCORED_MAGNITUDE = 1.0

SORT_KEYS: dict[str, Callable[[dict], Any]] = {
    "magnitude":    lambda f: _mag(f),
    "frequency":    lambda f: _numeric(f.get("freq"), -1.0),
    "publications": lambda f: _numeric(f.get("publications"), -1),
    "location":     lambda f: (_chrom_order(f.get("chromosome")),
                               _numeric(f.get("position"), 0)),
    "modified":     lambda f: str(f.get("discovered_at") or ""),
    "gmaf":         lambda f: _numeric(f.get("gmaf"), -1.0),
    "stars":        lambda f: _numeric(f.get("review_stars"), -1),
    "gene":         lambda f: str(f.get("gene") or "~"),
    "rsid":         lambda f: _rsid_order(f.get("rsid")),
    "coverage":     lambda f: _numeric(f.get("coverage"), -1.0),
}


# ---------------------------------------------------------------------------
# Coercion helpers
# ---------------------------------------------------------------------------

def _numeric(value: Any, default: Any) -> Any:
    if value is None or isinstance(value, bool):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _mag(finding: dict) -> float:
    """Magnitude for sorting and filtering, with null treated as 1."""
    value = finding.get("magnitude")
    if value is None or isinstance(value, bool):
        return UNSCORED_MAGNITUDE
    try:
        return float(value)
    except (TypeError, ValueError):
        return UNSCORED_MAGNITUDE


def _chrom_order(chrom: Any) -> tuple:
    """Sort chromosomes 1 to 22 numerically, then X, Y, MT, then anything else."""
    text = str(chrom or "").strip().upper().replace("CHR", "")
    if text.isdigit():
        return (0, int(text), "")
    order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if text in order:
        return (1, order[text], "")
    return (2, 0, text)


def _rsid_order(rsid: Any) -> tuple:
    """Sort rsIDs numerically, so rs99 precedes rs100."""
    text = str(rsid or "").strip().lower()
    match = re.match(r"^(rs|i)(\d+)$", text)
    if match:
        return (0, match.group(1), int(match.group(2)))
    return (1, text, 0)


def parse_csv(value: Any) -> list[str]:
    """Split a comma-separated parameter into a clean list."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        items = list(value)
    else:
        items = str(value).split(",")
    return [str(i).strip() for i in items if str(i).strip()]


def parse_bool(value: Any, default: bool = False) -> bool:
    """Accept the usual truthy spellings a query string can carry."""
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes", "on", "y", "t")


def parse_float(value: Any, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: Any, default: int | None = None) -> int | None:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Free-text query grammar, docs/API_V2.md section 3.5
# ---------------------------------------------------------------------------

_REGION_RE = re.compile(
    r"^chr(?:o|om|omosome)?\s*([0-9]{1,2}|x|y|m|mt)"
    r"(?::(\d+)(?:\s*([-+])\s*(\d+))?)?$",
    re.IGNORECASE,
)
_OP_RE = re.compile(r"/(clnsig|stars|mag|count|freq)\s*(>=|<=|>|<|=)\s*([\d.,]+)",
                    re.IGNORECASE)
_FLAG_RE = re.compile(r"/(dubious|flipped|ambiguous|conflict|carrier|nocall)",
                      re.IGNORECASE)


def parse_query(raw: Any) -> dict:
    """Split a free-text query into structured operators plus leftover text.

    Returns::

        {"text": str, "regex": compiled or None, "region": dict or None,
         "ops": {name: (operator, [values])}, "flags": set of str}

    An unparseable regular expression degrades to a literal substring match
    rather than raising, because a user typing an unbalanced bracket should get
    results, not a 500.
    """
    text = str(raw or "").strip()
    out: dict = {"text": "", "regex": None, "region": None, "ops": {}, "flags": set()}
    if not text:
        return out

    for match in _FLAG_RE.finditer(text):
        out["flags"].add(match.group(1).lower())
    text = _FLAG_RE.sub(" ", text)

    for match in _OP_RE.finditer(text):
        name, operator, values = match.group(1).lower(), match.group(2), match.group(3)
        out["ops"][name] = (operator, [v for v in values.split(",") if v])
    text = _OP_RE.sub(" ", text)

    text = text.strip()
    region_match = _REGION_RE.match(text)
    if region_match:
        chrom, start, sign, span = region_match.groups()
        region = {"chromosome": chrom.upper(), "start": None, "end": None}
        if start is not None:
            begin = int(start)
            if sign == "-" and span:
                region["start"], region["end"] = begin, int(span)
            elif sign == "+" and span:
                region["start"], region["end"] = begin, begin + int(span)
            else:
                region["start"] = region["end"] = begin
        out["region"] = region
        return out

    if text:
        out["text"] = text
        try:
            out["regex"] = re.compile(text, re.IGNORECASE)
        except re.error:
            out["regex"] = re.compile(re.escape(text), re.IGNORECASE)
    return out


_SEARCH_FIELDS = (
    "rsid", "gene", "summary", "interpretation", "conditions", "token",
    "genotype", "evidence", "clinical_sig", "criteria", "name", "aka",
)


def matches_query(finding: dict, parsed: dict) -> bool:
    """Test one finding against a parsed query."""
    flags = parsed.get("flags") or set()
    if "dubious" in flags and not finding.get("dubious"):
        return False
    if "flipped" in flags and not finding.get("flipped"):
        return False
    if "ambiguous" in flags and not (finding.get("ambiguous")
                                     or finding.get("freq_ambiguous")):
        return False
    if "conflict" in flags and not finding.get("conflict"):
        return False
    if "carrier" in flags and finding.get("carrier") is not True:
        return False
    if "nocall" in flags and finding.get("zygosity") != "no_call":
        return False

    for name, (operator, values) in (parsed.get("ops") or {}).items():
        if name == "clnsig":
            code = finding.get("clinvar_sig_code")
            wanted = {parse_int(v) for v in values}
            if code not in wanted:
                return False
            continue
        field = {"stars": "review_stars", "mag": "magnitude",
                 "count": "count", "freq": "freq"}[name]
        actual = _mag(finding) if field == "magnitude" else _numeric(finding.get(field), None)
        target = parse_float(values[0])
        if actual is None or target is None:
            return False
        if operator == ">=" and not actual >= target:
            return False
        if operator == "<=" and not actual <= target:
            return False
        if operator == ">" and not actual > target:
            return False
        if operator == "<" and not actual < target:
            return False
        if operator == "=" and not abs(actual - target) < 1e-9:
            return False

    region = parsed.get("region")
    if region:
        chrom = str(finding.get("chromosome") or "").strip().upper().replace("CHR", "")
        if chrom != region["chromosome"]:
            return False
        if region["start"] is not None:
            pos = _numeric(finding.get("position"), None)
            if pos is None or not (region["start"] <= pos <= region["end"]):
                return False
        return True

    regex = parsed.get("regex")
    if regex is not None:
        haystack = []
        for field in _SEARCH_FIELDS:
            value = finding.get(field)
            if value:
                haystack.append(str(value))
        for field in ("topics", "medicines", "conditions_list", "matched_rsids", "labels"):
            values = finding.get(field) or []
            if isinstance(values, (list, tuple)):
                haystack.extend(str(v) for v in values)
        if not regex.search(" ".join(haystack)):
            return False
    return True


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def apply_filters(findings: Iterable[dict], params: dict) -> list[dict]:
    """Apply every documented filter to a finding list.

    ``params`` is a plain dict of the query parameters in docs/API_V2.md 3.4.
    Missing keys mean "do not filter on this".
    """
    silo = str(params.get("silo") or "").strip()
    entity_types = set(parse_csv(params.get("entity_type")))
    min_mag = parse_float(params.get("min_magnitude"))
    max_mag = parse_float(params.get("max_magnitude"))
    reputes = {r.lower() for r in parse_csv(params.get("repute"))}
    min_pubs = parse_int(params.get("min_publications"))
    max_pubs = parse_int(params.get("max_publications"))
    min_freq = parse_float(params.get("min_freq"))
    max_freq = parse_float(params.get("max_freq"))
    require_freq = parse_bool(params.get("require_frequency"))
    clinvar_only = parse_bool(params.get("clinvar_only"))
    clinvar_sig = {parse_int(v) for v in parse_csv(params.get("clinvar_sig"))}
    min_stars = parse_int(params.get("min_stars"))
    genes = {g.upper() for g in parse_csv(params.get("gene"))}
    topics = {t.lower() for t in parse_csv(params.get("topic"))}
    medicines = {m.lower() for m in parse_csv(params.get("medicine"))}
    conditions = {c.lower() for c in parse_csv(params.get("condition"))}
    zygosity = {z.lower() for z in parse_csv(params.get("zygosity"))}
    carrier_only = parse_bool(params.get("carrier_only"))
    conflicts_only = parse_bool(params.get("conflicts_only"))
    ambiguous_only = parse_bool(params.get("ambiguous_only"))
    parsed_query = parse_query(params.get("q"))

    # When ClinVar mode is on and no explicit list was given, default to
    # pathogenic and likely pathogenic, matching the reference product.
    if clinvar_only and not clinvar_sig:
        clinvar_sig = {5, 4}

    out: list[dict] = []
    for f in findings or []:
        entity = str(f.get("entity_type") or "snp")
        exempt = entity in EXEMPT_ENTITIES

        if silo and f.get("silo") != silo:
            continue
        if entity_types and entity not in entity_types:
            continue

        mag = _mag(f)
        if min_mag is not None and mag < min_mag:
            continue
        if max_mag is not None and mag > max_mag:
            continue

        if reputes:
            rep = (f.get("repute") or "unset").lower()
            rep = rep if rep in ("good", "bad") else "unset"
            if rep not in reputes:
                continue

        # Rule 1: publication and frequency filters skip exempt entities.
        if not exempt:
            pubs = _numeric(f.get("publications"), 0)
            if min_pubs is not None and pubs < min_pubs:
                continue
            if max_pubs is not None and pubs > max_pubs:
                continue

            freq = f.get("freq")
            if require_freq and freq is None:
                continue
            if freq is None:
                if min_freq is not None and min_freq > 0:
                    continue
            else:
                if min_freq is not None and freq < min_freq:
                    continue
                if max_freq is not None and freq > max_freq:
                    continue

        if clinvar_only and f.get("clinvar_sig_code") is None:
            continue
        if clinvar_sig and f.get("clinvar_sig_code") not in clinvar_sig:
            continue
        if min_stars is not None and _numeric(f.get("review_stars"), 0) < min_stars:
            continue

        if genes and str(f.get("gene") or "").upper() not in genes:
            continue
        if topics and not ({str(t).lower() for t in (f.get("topics") or [])} & topics):
            continue
        if medicines and not ({str(m).lower() for m in (f.get("medicines") or [])} & medicines):
            continue
        if conditions:
            have = {str(c).lower() for c in (f.get("conditions_list") or [])}
            have.add(str(f.get("conditions") or "").lower())
            if not (have & conditions):
                continue

        if zygosity and str(f.get("zygosity") or "").lower() not in zygosity:
            continue
        if carrier_only and f.get("carrier") is not True:
            continue
        if conflicts_only and not f.get("conflict"):
            continue
        if ambiguous_only and not (f.get("ambiguous") or f.get("freq_ambiguous")
                                   or f.get("flipped")):
            continue

        if not matches_query(f, parsed_query):
            continue

        out.append(f)
    return out


# ---------------------------------------------------------------------------
# Sorting and pagination
# ---------------------------------------------------------------------------

def sort_findings(findings: list[dict], sort: str = "magnitude",
                  order: str = "desc") -> list[dict]:
    """Sort by any documented key in either direction.

    An unknown sort key falls back to magnitude rather than raising, so a stale
    bookmark cannot break the view.
    """
    key = SORT_KEYS.get(str(sort or "").strip().lower(), SORT_KEYS["magnitude"])
    reverse = str(order or "desc").strip().lower() != "asc"
    try:
        return sorted(findings, key=key, reverse=reverse)
    except TypeError:
        # Mixed types in one column, fall back to a string comparison.
        return sorted(findings, key=lambda f: str(key(f)), reverse=reverse)


def paginate(findings: list[dict], limit: int | None = 200,
             offset: int = 0) -> list[dict]:
    """Slice a result set. ``limit`` of 0 or None means everything."""
    start = max(0, int(offset or 0))
    if not limit:
        return findings[start:]
    return findings[start:start + int(limit)]


# ---------------------------------------------------------------------------
# Facets and summaries
# ---------------------------------------------------------------------------

def build_facets(findings: Iterable[dict]) -> dict:
    """Count every filterable value, so each dropdown can render 'Name (n)'."""
    buckets: dict[str, dict[str, int]] = {
        "genes": {}, "topics": {}, "medicines": {}, "conditions": {},
        "silos": {}, "entity_types": {}, "zygosity": {}, "categories": {},
        "reputes": {}, "clinvar_significance": {}, "review_stars": {},
        "cpic_levels": {}, "freq_bands": {}, "confidence": {},
        # clinvar_diseases is a DISTINCT facet from conditions, and the contract
        # (section 3.6) lists both. `conditions` holds whatever common-name
        # wording the source used; `clinvar_diseases` holds only the names that
        # came with an actual ClinVar record, which are formally controlled by
        # ClinVar. Merging them would let a loose common name masquerade as a
        # curated disease term. routes_v2 returns this dict verbatim, so a
        # missing bucket is a missing documented response key.
        "clinvar_diseases": {},
    }

    def bump(bucket: str, value: Any) -> None:
        text = str(value or "").strip()
        if not text:
            return
        buckets[bucket][text] = buckets[bucket].get(text, 0) + 1

    for f in findings or []:
        bump("genes", f.get("gene"))
        bump("silos", f.get("silo"))
        bump("entity_types", f.get("entity_type"))
        bump("zygosity", f.get("zygosity"))
        bump("categories", f.get("category"))
        bump("reputes", f.get("repute") or "Not Set")
        bump("clinvar_significance", f.get("clinical_sig"))
        bump("cpic_levels", f.get("cpic_level"))
        bump("freq_bands", f.get("freq_band"))
        bump("confidence", f.get("confidence"))
        stars = f.get("review_stars")
        if isinstance(stars, int):
            bump("review_stars", str(stars))
        for topic in f.get("topics") or []:
            bump("topics", topic)
        for med in f.get("medicines") or []:
            bump("medicines", med)
        for cond in f.get("conditions_list") or []:
            bump("conditions", cond)
        # Only findings that actually carry a ClinVar record contribute to the
        # ClinVar disease facet, so the two buckets stay meaningfully different.
        if f.get("clinvar_sig_code") is not None:
            for cond in f.get("clinvar_diseases") or f.get("conditions_list") or []:
                bump("clinvar_diseases", cond)
        # Only findings that actually carry a ClinVar record contribute to the
        # ClinVar disease facet, so the two buckets stay meaningfully different.
        if f.get("clinvar_sig_code") is not None:
            for cond in f.get("clinvar_diseases") or f.get("conditions_list") or []:
                bump("clinvar_diseases", cond)

    def render(bucket: dict[str, int]) -> list[dict]:
        return [{"value": k, "count": v}
                for k, v in sorted(bucket.items(), key=lambda kv: (-kv[1], kv[0]))]

    return {name: render(bucket) for name, bucket in buckets.items()}


def summarise(findings: Iterable[dict]) -> dict:
    """Counts per silo, entity type and repute for a result set."""
    silos: dict[str, int] = {}
    entities: dict[str, int] = {}
    reputes = {"Good": 0, "Bad": 0, "unset": 0}
    for f in findings or []:
        silo = f.get("silo") or "informational"
        silos[silo] = silos.get(silo, 0) + 1
        ent = f.get("entity_type") or "snp"
        entities[ent] = entities.get(ent, 0) + 1
        rep = f.get("repute") or "unset"
        reputes[rep if rep in ("Good", "Bad") else "unset"] += 1
    return {"silos": silos, "entity_types": entities, "reputes": reputes}


def filter_and_sort(findings: list[dict], params: dict) -> dict:
    """Full pipeline: filter, sort, paginate, and report the counts.

    Returns the response body described in docs/API_V2.md section 3.4, minus the
    keys only the HTTP layer knows (population, qc).
    """
    filtered = apply_filters(findings, params)
    ordered = sort_findings(
        filtered,
        sort=params.get("sort") or "magnitude",
        order=params.get("order") or "desc",
    )
    limit = parse_int(params.get("limit"), 200)
    offset = parse_int(params.get("offset"), 0) or 0
    page = paginate(ordered, limit=limit, offset=offset)
    return {
        "findings": page,
        "total": len(ordered),
        "returned": len(page),
        "offset": offset,
        "filtered_summary": summarise(ordered),
        "facets": build_facets(ordered),
    }
