"""
merge.py -- multi-file pooling, conflict retention, family roles and trio checks.

Design intent, taken from how Promethease actually behaves:

  * Two files from the SAME person (for example a 23andMe export plus an
    AncestryDNA export) are POOLED. The merged result is the UNION of their
    genotypes.
  * Where the two files disagree at a position, BOTH calls are retained and
    surfaced to the user. There is no voting, no confidence weighting and no
    automatic winner. Disagreement is information, not noise to be hidden.
  * Uploading a DIFFERENT person's file does not create a report for them.
    Their genotypes appear as comparison rows inside the primary person's
    report.

Only the role "self" is pooled. Every other role becomes a comparison set.
All returned lists are sorted by rsID so repeated runs produce identical output.
"""

from typing import Any

try:  # orientation.py may be added to the package after this module
    from backend import orientation as _orientation
except ImportError:  # pragma: no cover - only before orientation.py lands
    try:
        import orientation as _orientation  # type: ignore
    except ImportError:
        _orientation = None  # type: ignore


ROLES = ("self", "mother", "father", "mate", "child", "sibling", "other", "ignore")
PRIMARY_ROLE = "self"
NOCALL_KEY = "NN"

# Tokens consumer arrays use for a missing call, plus the D/I indel tokens that
# parsers.py already normalises away.
_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "D", "I", "?", "."}

_ROLE_ALIASES = {
    "me": "self",
    "myself": "self",
    "primary": "self",
    "proband": "self",
    "mom": "mother",
    "mum": "mother",
    "dad": "father",
    "spouse": "mate",
    "partner": "mate",
    "husband": "mate",
    "wife": "mate",
    "son": "child",
    "daughter": "child",
    "brother": "sibling",
    "sister": "sibling",
    "skip": "ignore",
    "exclude": "ignore",
    "excluded": "ignore",
}


class MergeError(Exception):
    """Raised when a list of upload sources cannot be merged."""


def _clean_allele(value: Any) -> str:
    """Normalise a single allele token to an upper-case stripped string."""
    return str(value if value is not None else "").strip().upper()


def is_real_call(a1: Any, a2: Any) -> bool:
    """Return True when both alleles are genuine base calls.

    A no-call token on either side (blank, N, -, 0, or the D/I indel
    placeholders emitted by several vendors) makes the whole genotype a
    no-call, because half a diploid call cannot be compared to another file.
    """
    left, right = _clean_allele(a1), _clean_allele(a2)
    if left in _NOCALL_ALLELES or right in _NOCALL_ALLELES:
        return False
    return True


def _sorted_pair(a1: Any, a2: Any) -> tuple[str, str]:
    """Return the two alleles in a canonical (sorted) order.

    Delegates to orientation.sort_alleles when that module is present, but
    validates its answer and falls back to a plain sort, so merging never
    depends on another module being importable.
    """
    left, right = _clean_allele(a1), _clean_allele(a2)
    if _orientation is not None:
        try:
            pair = _orientation.sort_alleles(left, right)
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                first, second = _clean_allele(pair[0]), _clean_allele(pair[1])
                if sorted((first, second)) == sorted((left, right)):
                    return first, second
        except Exception:
            pass
    ordered = sorted((left, right))
    return ordered[0], ordered[1]


def genotype_key(a1: Any, a2: Any) -> str:
    """Return the order-insensitive genotype key, or "NN" for a no-call.

    "AG" and "GA" both key to "AG", so two vendors that report the same
    heterozygote in opposite allele order are never treated as disagreeing.
    """
    if not is_real_call(a1, a2):
        return NOCALL_KEY
    first, second = _sorted_pair(a1, a2)
    return f"{first}{second}"


def normalize_role(value: Any) -> str:
    """Map free-text role input onto one of ROLES.

    Matching is case-insensitive and tolerates spaces, hyphens and a small
    alias table ("mom", "dad", "spouse"). Anything unrecognised or empty
    becomes "other", which is a comparison role, never the primary one.
    """
    role = str(value if value is not None else "").strip().lower()
    role = role.replace(" ", "_").replace("-", "_")
    role = _ROLE_ALIASES.get(role, role)
    return role if role in ROLES else "other"


def is_pooled_role(role: Any) -> bool:
    """Return True only for the primary role, the one that gets pooled."""
    return normalize_role(role) == PRIMARY_ROLE


def dedupe_snps(snps: Any) -> list[dict]:
    """Collapse repeated rsIDs inside a SINGLE file, keeping the first real call.

    Some vendor exports list the same position twice (duplicated probes, or a
    position that appears in more than one block of the file). Within one file
    those repeats are not a conflict, so the first genuine call wins and a
    later duplicate is dropped. A no-call placeholder is replaced if a real
    call for the same rsID appears further down the file.
    """
    out: list[dict] = []
    seen: dict[str, int] = {}
    for snp in (snps or []):
        if not isinstance(snp, dict):
            continue
        rsid = str(snp.get("rsid") or "").strip().lower()
        if not rsid:
            continue
        real = is_real_call(snp.get("allele1"), snp.get("allele2"))
        index = seen.get(rsid)
        if index is None:
            seen[rsid] = len(out)
            out.append(dict(snp))
            continue
        kept = out[index]
        if real and not is_real_call(kept.get("allele1"), kept.get("allele2")):
            out[index] = dict(snp)
    return out


def _prepare_sources(sources: Any) -> list[dict]:
    """Validate, normalise and de-duplicate the caller's source list."""
    if sources is None or isinstance(sources, (str, bytes, dict)):
        raise MergeError("merge_sources expects a list of source dicts")
    prepared: list[dict] = []
    used_labels: set[str] = set()
    for position, src in enumerate(sources):
        if not isinstance(src, dict):
            raise MergeError(f"source at index {position} is not a dict")
        role = normalize_role(src.get("role"))
        if role == "ignore":
            continue
        label = str(src.get("label") or "").strip() or f"file{position + 1}"
        if label in used_labels:
            suffix = 2
            while f"{label} ({suffix})" in used_labels:
                suffix += 1
            label = f"{label} ({suffix})"
        used_labels.add(label)
        prepared.append({
            "label": label,
            "role": role,
            "provider": str(src.get("provider") or "").strip(),
            "snps": dedupe_snps(src.get("snps") or []),
        })
    if not prepared:
        raise MergeError("no usable sources: the list was empty or every source was ignored")
    return prepared


def _new_entry(rsid: str, snp: dict) -> dict:
    """Build an empty pooled-genotype entry for one rsID."""
    return {
        "rsid": rsid,
        "chromosome": str(snp.get("chromosome") or "").strip(),
        "position": snp.get("position") or 0,
        "allele1": "N",
        "allele2": "N",
        "genotype": NOCALL_KEY,
        "count": 0,
        "labels": [],
        "conflict": False,
        "calls": [],
    }


def merge_sources(sources: list[dict]) -> dict:
    """Pool every "self" source and expose all other roles as comparison sets.

    Each source is
        {"label": str, "role": str, "provider": str,
         "snps": [{rsid, chromosome, position, allele1, allele2}, ...]}

    Sources whose role is "ignore" are dropped before anything else happens.
    Sources whose role is "self" are pooled into one primary genotype set:
    the union of their positions. Where two pooled files both make a real call
    at the same rsID and the sorted allele pairs differ, the position is marked
    conflict=True and BOTH calls stay in its "calls" list. A no-call in one
    file is never a conflict; the real call simply wins and "count" only ever
    counts real calls.

    Every other role becomes an entry in "comparison", keyed by role where the
    role is unambiguous and by label when the same role is uploaded twice.

    Sources are walked in the order given and every returned mapping is
    rebuilt in rsID order, so the output is deterministic.
    """
    prepared = _prepare_sources(sources)
    pooled = [s for s in prepared if s["role"] == PRIMARY_ROLE]
    others = [s for s in prepared if s["role"] != PRIMARY_ROLE]

    stats: dict[str, dict] = {
        s["label"]: {
            "label": s["label"],
            "role": s["role"],
            "provider": s["provider"],
            "snp_count": len(s["snps"]),
            "contributed": 0,
            "overlapped": 0,
            "conflicting": 0,
        }
        for s in prepared
    }

    genotypes: dict[str, dict] = {}
    for src in pooled:
        label = src["label"]
        stat = stats[label]
        for snp in src["snps"]:
            rsid = str(snp.get("rsid") or "").strip().lower()
            if not rsid:
                continue
            entry = genotypes.get(rsid)
            if entry is None:
                entry = _new_entry(rsid, snp)
                genotypes[rsid] = entry
                stat["contributed"] += 1
            else:
                stat["overlapped"] += 1
                if not entry["chromosome"]:
                    entry["chromosome"] = str(snp.get("chromosome") or "").strip()
                if not entry["position"]:
                    entry["position"] = snp.get("position") or 0
            a1, a2 = snp.get("allele1"), snp.get("allele2")
            if not is_real_call(a1, a2):
                continue
            first, second = _sorted_pair(a1, a2)
            key = f"{first}{second}"
            if entry["count"] == 0:
                entry["allele1"], entry["allele2"], entry["genotype"] = first, second, key
            elif key != entry["genotype"]:
                entry["conflict"] = True
            entry["count"] += 1
            entry["labels"].append(label)
            entry["calls"].append({
                "label": label,
                "allele1": first,
                "allele2": second,
                "genotype": key,
            })

    genotypes = {rsid: genotypes[rsid] for rsid in sorted(genotypes)}

    conflicts: list[dict] = []
    for rsid, entry in genotypes.items():
        if not entry["conflict"]:
            continue
        conflicts.append({
            "rsid": rsid,
            "chromosome": entry["chromosome"],
            "position": entry["position"],
            "calls": [{"label": c["label"], "genotype": c["genotype"]} for c in entry["calls"]],
        })
        for label in set(entry["labels"]):
            stats[label]["conflicting"] += 1

    comparison: dict[str, dict] = {}
    for src in others:
        key = src["role"] if src["role"] not in comparison else src["label"]
        if key in comparison:
            key = f"{src['role']}:{src['label']}"
        rows: dict[str, dict] = {}
        for snp in src["snps"]:
            rsid = str(snp.get("rsid") or "").strip().lower()
            if not rsid or not is_real_call(snp.get("allele1"), snp.get("allele2")):
                continue
            first, second = _sorted_pair(snp.get("allele1"), snp.get("allele2"))
            rows[rsid] = {
                "allele1": first,
                "allele2": second,
                "genotype": f"{first}{second}",
                "label": src["label"],
                "role": src["role"],
            }
        comparison[key] = {rsid: rows[rsid] for rsid in sorted(rows)}
        stat = stats[src["label"]]
        stat["contributed"] = len(rows)
        stat["overlapped"] = sum(1 for rsid in rows if rsid in genotypes)

    unique_by_label = {src["label"]: 0 for src in pooled}
    for entry in genotypes.values():
        labels = set(entry["labels"])
        if len(labels) == 1:
            only = labels.pop()
            unique_by_label[only] = unique_by_label.get(only, 0) + 1

    return {
        "genotypes": genotypes,
        "comparison": comparison,
        "conflicts": conflicts,
        "sources": [stats[s["label"]] for s in prepared],
        "counts": {
            "total_positions": len(genotypes),
            "pooled_sources": len(pooled),
            "comparison_sources": len(others),
            "union": len(genotypes),
            "shared": sum(1 for e in genotypes.values() if e["count"] >= 2),
            "conflicts": len(conflicts),
            "unique_by_label": unique_by_label,
        },
        "primary_labels": [s["label"] for s in pooled],
    }


def _as_key(value: Any) -> str:
    """Coerce a genotype string or 2-item sequence to a canonical genotype key."""
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return genotype_key(value[0], value[1])
    text = str(value if value is not None else "").strip().upper()
    if len(text) == 2:
        return genotype_key(text[0], text[1])
    return NOCALL_KEY


def transmission_probability(child_gt: Any, parent1_gt: Any, parent2_gt: Any) -> float | None:
    """Probability of the child genotype given both parents, Mendelian model.

    Each parent contributes one of its two alleles with equal probability and
    the two parental gametes are drawn independently, giving four equally
    likely combinations. The returned value is the share of those four
    combinations that produce the observed child genotype, so it is one of
    0.0, 0.25, 0.5, 0.75 or 1.0.

    Returns None when any of the three genotypes is a no-call, and 0.0 when
    the child genotype cannot arise from these parents at all.
    """
    child = _as_key(child_gt)
    p1 = _as_key(parent1_gt)
    p2 = _as_key(parent2_gt)
    if NOCALL_KEY in (child, p1, p2):
        return None
    hits = 0
    for left in (p1[0], p1[1]):
        for right in (p2[0], p2[1]):
            if genotype_key(left, right) == child:
                hits += 1
    return hits / 4.0


def mendelian_violation(child_gt: Any, parent1_gt: Any, parent2_gt: Any) -> bool:
    """Return True when the child genotype is impossible from these parents.

    A no-call anywhere is not a violation, it is simply uninformative.
    """
    probability = transmission_probability(child_gt, parent1_gt, parent2_gt)
    return probability is not None and probability == 0.0


TRIO_NOTE = (
    "A high violation rate almost never means non-paternity. It usually means a "
    "strand or array mismatch: the two files were produced on different chips or "
    "reported on opposite strands, so the same position is written as, for "
    "example, A/G in one file and T/C in the other. Check strand orientation and "
    "chip versions before drawing any conclusion about relatedness."
)


def trio_annotate(merged: dict) -> dict:
    """Annotate every primary genotype that has both a mother and a father call.

    For each such position this adds "parent1" (mother), "parent2" (father),
    "probability" and "mendelian_ok" to the genotype entry inside ``merged``,
    mutating it in place, and returns a summary.

    Interpretation warning: a high violation rate almost never means
    non-paternity. It usually means a strand or array mismatch between the
    uploaded files (different chip versions, or genotypes reported on opposite
    strands), which makes genuinely compatible calls look impossible. The same
    warning is returned in the summary's "note" field so the UI can show it
    next to the number.
    """
    genotypes = (merged or {}).get("genotypes") or {}
    comparison = (merged or {}).get("comparison") or {}
    mother = comparison.get("mother") or {}
    father = comparison.get("father") or {}
    available = bool(mother) and bool(father)

    compared = 0
    violations = 0
    violation_rsids: list[str] = []

    if available:
        for rsid in sorted(genotypes):
            entry = genotypes[rsid]
            mrow, frow = mother.get(rsid), father.get(rsid)
            if not mrow or not frow:
                continue
            child = entry.get("genotype", NOCALL_KEY)
            p1 = mrow.get("genotype", NOCALL_KEY)
            p2 = frow.get("genotype", NOCALL_KEY)
            probability = transmission_probability(child, p1, p2)
            if probability is None:
                continue
            entry["parent1"] = p1
            entry["parent2"] = p2
            entry["probability"] = probability
            entry["mendelian_ok"] = probability > 0.0
            compared += 1
            if probability == 0.0:
                violations += 1
                violation_rsids.append(rsid)

    rate = (violations / compared) if compared else 0.0
    return {
        "trio_available": available,
        "compared": compared,
        "violations": violations,
        "violation_rsids": sorted(violation_rsids),
        "violation_rate": round(rate, 4),
        "note": TRIO_NOTE,
    }


def comparison_rows(merged: dict, rsid: str) -> list[dict]:
    """Rows the UI shows underneath one genotype card.

    One row per non-self source that has a real call at ``rsid``:
    {"label", "role", "genotype", "shared"}. "shared" is True when that
    relative's call matches the primary pooled genotype, which is what makes a
    row worth showing at all.
    """
    merged = merged or {}
    key = str(rsid or "").strip().lower()
    primary = (merged.get("genotypes") or {}).get(key) or {}
    primary_gt = primary.get("genotype", NOCALL_KEY)
    rows: list[dict] = []
    for group in sorted((merged.get("comparison") or {}).keys()):
        row = (merged["comparison"][group] or {}).get(key)
        if not row:
            continue
        genotype = row.get("genotype", NOCALL_KEY)
        rows.append({
            "label": row.get("label", group),
            "role": row.get("role", group),
            "genotype": genotype,
            "shared": bool(primary_gt != NOCALL_KEY and genotype == primary_gt),
        })
    return rows
