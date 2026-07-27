"""
genosets.py -- genoset criteria language: parser, evaluator and corpus loader.

A *genoset* is a named boolean rule over several SNP genotypes. It is either
present or absent, never scored and never quantitative. Genosets carry
Magnitude, Repute and Summary but have NO frequency, GMAF, chromosome,
position, gene or ClinVar fields, and they are exempt from publication-count
and frequency filtering in the UI.

Criteria language
-----------------
    rs1234(A;T)             exact genotype, order insensitive
    rs1234(T;T)             homozygous T
    rs1234(T)               at least one T allele observed
    and(X, Y, ...)          all arguments true
    or(X, Y, ...)           any argument true
    not(X, Y, ...)          ALL arguments false (n-ary, NOT unary)
    atleast(N, X, Y, ...)   at least N arguments true
    dgs001                  reference to another genoset (lower-numbered only)

23andMe internal identifiers of the form ``i12345`` are accepted wherever an
rsID is accepted. rsIDs are case-insensitive and normalised to lowercase.

Comments: whole-line ``#`` and ``###`` comments (leading whitespace allowed)
and ``<!-- ... -->`` HTML comments are stripped. Expressions may span multiple
lines and may contain arbitrary internal whitespace.

Missing-data semantics (authoritative; matches Promethease)
-----------------------------------------------------------
A SNP that is absent from the person's data evaluates to **False**. It is
never null, never imputed to the population-major allele and never
permutation tested. A no-call genotype is likewise False. Users rely on this
behaviour, which is why :func:`evaluate_all_verbose` reports an ``incomplete``
bucket -- so the UI can say "not testable on your array" rather than "absent".

Allele orientation
------------------
Genotypes are expressed on the dbSNP plus (forward) strand, which is the
orientation used by 23andMe and AncestryDNA raw data files and by
``backend/orientation.py``. Callers must pass oriented alleles.
"""

import json
import re
from pathlib import Path

_BASE = Path(__file__).parent.parent
GENOSET_FILE = _BASE / "data" / "genosets.json"

# Alleles that mean "this position was not called".
NOCALL_ALLELES = {"", "N", "-", "0", "--"}

_FUNCS = {"and", "or", "not", "atleast"}
_RSID_RE = re.compile(r"^(?:rs|i)\d+$", re.IGNORECASE)
_GSREF_RE = re.compile(r"^d?gs\d+$", re.IGNORECASE)
_HTML_COMMENT_RE = re.compile(r"<!--.*?-->", re.DOTALL)
_HTML_OPEN_RE = re.compile(r"<!--.*\Z", re.DOTALL)
_NUM_RE = re.compile(r"(\d+)")

# Punctuation, or any run of characters that is not whitespace/punctuation.
_TOKEN_RE = re.compile(r"[(),;]|[^\s(),;]+")

VALID_SILOS = ("pre_prescription", "actionable", "informational")


class CriteriaError(Exception):
    """Raised for any malformed genoset criteria expression or corpus."""


# ---------------------------------------------------------------------------
# Comment stripping
# ---------------------------------------------------------------------------

def strip_comments(text: str) -> str:
    """Remove HTML comments and whole-line ``#`` comments, and de-indent.

    SNPedia-style criteria bodies are space-indented and may carry ``#`` or
    ``###`` comment lines plus ``<!-- ... -->`` HTML comments that can span
    several lines. Every surviving line is stripped of leading and trailing
    whitespace so the tokeniser never has to care about indentation.
    """
    if not text:
        return ""
    text = _HTML_COMMENT_RE.sub(" ", text)
    # An unterminated "<!--" swallows the remainder of the body.
    text = _HTML_OPEN_RE.sub(" ", text)
    kept = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        kept.append(line)
    return "\n".join(kept)


def _tokenise(text: str) -> list[str]:
    return _TOKEN_RE.findall(text)


# ---------------------------------------------------------------------------
# Recursive-descent parser
#
# AST nodes are plain dicts so the tree stays JSON serialisable:
#   {"op": "geno",    "rsid": "rs1234", "alleles": ("A", "T"), "mode": "exact"}
#   {"op": "geno",    "rsid": "rs1234", "alleles": ("T",), "mode": "atleast_one"}
#   {"op": "and",     "args": [...]}
#   {"op": "or",      "args": [...]}
#   {"op": "not",     "args": [...]}
#   {"op": "atleast", "n": 3, "args": [...]}
#   {"op": "gsref",   "name": "dgs026"}
# ---------------------------------------------------------------------------

_PUNCT = {"(", ")", ",", ";"}


class _Parser:
    def __init__(self, tokens: list[str]):
        self.toks = tokens
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self) -> str:
        tok = self.peek()
        if tok is None:
            raise CriteriaError(
                "unexpected end of criteria (unbalanced parentheses?)")
        self.i += 1
        return tok

    def expect(self, want: str) -> str:
        tok = self.peek()
        if tok != want:
            got = "end of criteria" if tok is None else repr(tok)
            raise CriteriaError(
                f"expected {want!r} but found {got} "
                f"(unbalanced parentheses?)")
        self.i += 1
        return want

    # -- productions ------------------------------------------------------
    def expr(self) -> dict:
        tok = self.take()
        if tok in _PUNCT:
            raise CriteriaError(f"unexpected token {tok!r} where an "
                                f"expression was expected")
        low = tok.lower()

        if low in _FUNCS:
            self.expect("(")
            if low == "atleast":
                return self._atleast()
            return {"op": low, "args": self.arglist()}

        if _RSID_RE.match(tok):
            self.expect("(")
            alleles = self.allelelist(low)
            mode = "atleast_one" if len(alleles) == 1 else "exact"
            return {"op": "geno", "rsid": low,
                    "alleles": tuple(alleles), "mode": mode}

        if _GSREF_RE.match(tok):
            return {"op": "gsref", "name": low}

        raise CriteriaError(f"unknown function or identifier {tok!r}")

    def _atleast(self) -> dict:
        ntok = self.take()
        if ntok in _PUNCT or not ntok.isdigit():
            raise CriteriaError(
                f"atleast() requires a non-negative integer as its first "
                f"argument, got {ntok!r}")
        n = int(ntok)
        if n < 1:
            raise CriteriaError("atleast() requires N >= 1")
        self.expect(",")
        args = self.arglist()
        return {"op": "atleast", "n": n, "args": args}

    def arglist(self) -> list[dict]:
        if self.peek() == ")":
            raise CriteriaError("empty argument list")
        args: list[dict] = []
        while True:
            args.append(self.expr())
            if self.peek() == ",":
                self.i += 1
                if self.peek() in (")", None):
                    raise CriteriaError("empty argument after ','")
                continue
            self.expect(")")
            return args

    def allelelist(self, rsid: str) -> list[str]:
        if self.peek() == ")":
            raise CriteriaError(
                f"{rsid}(): empty genotype, expected allele(s)")
        alleles: list[str] = []
        while True:
            tok = self.take()
            if tok in _PUNCT:
                raise CriteriaError(
                    f"{rsid}(): expected an allele but found {tok!r}")
            alleles.append(tok.upper())
            if self.peek() == ";":
                self.i += 1
                continue
            self.expect(")")
            break
        if len(alleles) > 2:
            raise CriteriaError(
                f"{rsid}(): expected 1 or 2 alleles, got {len(alleles)}")
        return alleles


def parse_criteria(text: str) -> dict:
    """Parse a genoset criteria body into an AST dict.

    Raises :class:`CriteriaError` on unbalanced parentheses, unknown function
    names, a non-integer ``atleast`` count, empty argument lists, or trailing
    junk after the end of the expression.
    """
    cleaned = strip_comments(text or "")
    tokens = _tokenise(cleaned)
    if not tokens:
        raise CriteriaError("criteria text is empty (nothing left after "
                            "comment stripping)")
    parser = _Parser(tokens)
    node = parser.expr()
    if parser.peek() is not None:
        junk = " ".join(parser.toks[parser.i:])[:60]
        hint = " (unbalanced parentheses?)" if ")" in junk else ""
        raise CriteriaError(
            f"trailing junk after end of expression: {junk!r}{hint}")
    return node


# ---------------------------------------------------------------------------
# AST introspection
# ---------------------------------------------------------------------------

def required_rsids(node: dict) -> set[str]:
    """Return every rsID this node tests directly (gsrefs are not expanded)."""
    if not isinstance(node, dict):
        raise CriteriaError(f"not an AST node: {node!r}")
    op = node.get("op")
    if op == "geno":
        return {node["rsid"]}
    if op == "gsref":
        return set()
    out: set[str] = set()
    for arg in node.get("args", []):
        out |= required_rsids(arg)
    return out


def referenced_genosets(node: dict) -> set[str]:
    """Return every genoset name this node references directly or nested."""
    if not isinstance(node, dict):
        raise CriteriaError(f"not an AST node: {node!r}")
    op = node.get("op")
    if op == "gsref":
        return {node["name"]}
    if op == "geno":
        return set()
    out: set[str] = set()
    for arg in node.get("args", []):
        out |= referenced_genosets(arg)
    return out


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def _split_genotype(value) -> list[str] | None:
    """Normalise a genotype to a list of uppercase alleles, or None if no-call.

    Accepts a 2-tuple/list of alleles or a 2-character string. Separators
    ``/``, ``|``, ``;`` and whitespace are tolerated inside strings.
    """
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        alleles = [str(a).strip().upper() for a in value]
    else:
        raw = str(value).strip().upper()
        for sep in ("/", "|", ";", " ", "\t"):
            raw = raw.replace(sep, "")
        if not raw:
            return None
        alleles = list(raw) if len(raw) > 1 else [raw]
    if not alleles:
        return None
    if any(a in NOCALL_ALLELES for a in alleles):
        return None
    return alleles


def _is_called(value) -> bool:
    """True when the array actually reported a usable genotype."""
    return _split_genotype(value) is not None


def evaluate(node: dict, genotypes: dict,
             genoset_results: dict | None = None) -> bool:
    """Evaluate an AST against a person's genotypes.

    ``genotypes`` maps lowercase rsid -> 2-tuple or 2-character string of
    ORIENTED alleles. A missing key is False. A no-call is False.
    ``genoset_results`` maps genoset name -> bool for ``gsref`` resolution; an
    unknown or unevaluated reference is False.
    """
    if not isinstance(node, dict):
        raise CriteriaError(f"not an AST node: {node!r}")
    op = node.get("op")

    if op == "geno":
        alleles = _split_genotype(genotypes.get(node["rsid"]))
        if alleles is None:
            return False
        want = [str(a).upper() for a in node["alleles"]]
        if node.get("mode") == "atleast_one":
            return want[0] in alleles
        return sorted(alleles) == sorted(want)

    if op == "gsref":
        results = genoset_results or {}
        return bool(results.get(node["name"], False))

    args = node.get("args") or []
    if not args:
        raise CriteriaError(f"{op}() has no arguments")

    if op == "and":
        return all(evaluate(a, genotypes, genoset_results) for a in args)
    if op == "or":
        return any(evaluate(a, genotypes, genoset_results) for a in args)
    if op == "not":
        return not any(evaluate(a, genotypes, genoset_results) for a in args)
    if op == "atleast":
        n = int(node.get("n", 1))
        hits = 0
        for a in args:
            if evaluate(a, genotypes, genoset_results):
                hits += 1
                if hits >= n:
                    return True
        return False

    raise CriteriaError(f"unknown AST operator {op!r}")


# ---------------------------------------------------------------------------
# Dependency ordering
# ---------------------------------------------------------------------------

def _numeric_key(name: str) -> tuple:
    """Sort key so dgs009 < dgs100 < dgs270 (numeric, not lexicographic)."""
    match = _NUM_RE.search(name or "")
    return (int(match.group(1)) if match else 0, name or "")


def _criteria_of(name: str, entry) -> str:
    if isinstance(entry, str):
        return entry
    if isinstance(entry, dict):
        return entry.get("criteria", "") or ""
    raise CriteriaError(f"{name}: corpus entry must be a dict or a string, "
                        f"got {type(entry).__name__}")


def _parse_corpus(corpus: dict) -> dict:
    """Parse every criteria string in a corpus into {name: ast}."""
    asts = {}
    for name, entry in (corpus or {}).items():
        asts[name] = parse_criteria(_criteria_of(name, entry))
    return asts


def topological_order(corpus: dict) -> list[str]:
    """Return genoset names ordered so every gsref dependency comes first.

    Ties are broken numerically by the integer part of the name, so dgs100
    is emitted before dgs270. Raises :class:`CriteriaError` on a cycle.
    """
    asts = _parse_corpus(corpus)
    deps = {
        name: {d for d in referenced_genosets(node) if d in asts}
        for name, node in asts.items()
    }
    remaining = dict(deps)
    order: list[str] = []
    while remaining:
        ready = sorted((n for n, d in remaining.items() if not d),
                       key=_numeric_key)
        if not ready:
            stuck = ", ".join(sorted(remaining, key=_numeric_key))
            raise CriteriaError(
                f"cycle detected in genoset references among: {stuck}")
        for name in ready:
            order.append(name)
            del remaining[name]
        done = set(order)
        for pending in remaining:
            remaining[pending] = remaining[pending] - done
    return order


# ---------------------------------------------------------------------------
# Corpus loading (mirrors backend/scanner.py's snp_reference.json handling)
# ---------------------------------------------------------------------------

_corpus_cache: dict | None = None
_corpus_meta: dict = {}


def load_genosets(path=None) -> dict:
    """Load ``data/genosets.json``. Returns ``{}`` when the file is missing.

    Supports both a flat ``{name: {...}}`` dict and the versioned
    ``{"_meta": ..., "genosets": ...}`` shape. Result is cached module-level.
    """
    global _corpus_cache, _corpus_meta
    use_default = path is None
    if use_default and _corpus_cache is not None:
        return _corpus_cache
    target = Path(path) if path is not None else GENOSET_FILE
    if not target.exists():
        if use_default:
            _corpus_cache, _corpus_meta = {}, {}
        return {}
    with open(target, "r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if isinstance(raw, dict) and "genosets" in raw and "_meta" in raw:
        _corpus_meta = raw.get("_meta") or {}
        _corpus_cache = raw.get("genosets") or {}
    else:
        _corpus_meta = {}
        _corpus_cache = raw or {}
    return _corpus_cache


def get_metadata() -> dict:
    """Return the ``_meta`` block from genosets.json (version, count, ...)."""
    load_genosets()
    return _corpus_meta


def _reset_cache() -> None:
    """Drop the cached corpus. Used by tests and by the updater."""
    global _corpus_cache, _corpus_meta
    _corpus_cache, _corpus_meta = None, {}


# ---------------------------------------------------------------------------
# Whole-corpus evaluation
# ---------------------------------------------------------------------------

def _transitive_rsids(name: str, asts: dict, seen=None) -> set[str]:
    """rsIDs needed by ``name`` including those reached through gsrefs."""
    seen = set() if seen is None else seen
    if name in seen:
        return set()
    seen.add(name)
    node = asts.get(name)
    if node is None:
        return set()
    out = set(required_rsids(node))
    for dep in referenced_genosets(node):
        out |= _transitive_rsids(dep, asts, seen)
    return out


def _magnitude_of(entry: dict):
    mag = entry.get("magnitude")
    if mag in (None, ""):
        return None
    try:
        return float(mag)
    except (TypeError, ValueError):
        return None


def _finding(name: str, entry: dict, asts: dict, genotypes: dict,
             matched: bool) -> dict:
    """Build a finding dict for one genoset, matching the project's schema."""
    required = sorted(_transitive_rsids(name, asts))
    present = [r for r in required if _is_called(genotypes.get(r))]
    coverage = 1.0 if not required else round(len(present) / len(required), 4)

    repute = (entry.get("repute") or "").strip()
    if repute not in ("Good", "Bad"):
        repute = ""
    silo = entry.get("silo") or "informational"
    if silo not in VALID_SILOS:
        silo = "informational"

    return {
        "rsid":           name,
        "entity_type":    "genoset",
        "gene":           "",
        "chromosome":     "",
        "position":       0,
        "allele1":        "",
        "allele2":        "",
        "genotype":       "",
        "zygosity":       "",
        "magnitude":      _magnitude_of(entry),
        "repute":         repute,
        "summary":        entry.get("summary", "") or "",
        "clinical_sig":   entry.get("clinical_sig") or entry.get("evidence", "") or "",
        "conditions":     entry.get("conditions") or entry.get("summary", "") or "",
        "interpretation": entry.get("interpretation", "") or "",
        "category":       entry.get("category", "") or "",
        "silo":           silo,
        "criteria":       _criteria_of(name, entry),
        "matched_rsids":  required,
        "coverage":       coverage,
        "sources":        ["genoset"],
        "topics":         list(entry.get("topics") or []),
        "medicines":      list(entry.get("medicines") or []),
        "aka":            entry.get("aka", "") or "",
        "evidence":       entry.get("evidence", "") or "",
        "matched":        bool(matched),
    }


def _evaluate_corpus(genotypes: dict, corpus: dict | None):
    """Shared engine: returns (corpus, asts, results dict, order)."""
    corpus = load_genosets() if corpus is None else corpus
    corpus = corpus or {}
    asts = _parse_corpus(corpus)
    order = topological_order(corpus)
    results: dict[str, bool] = {}
    for name in order:
        results[name] = evaluate(asts[name], genotypes, results)
    return corpus, asts, results, order


def evaluate_all(genotypes: dict, corpus: dict | None = None) -> list[dict]:
    """Evaluate every genoset in dependency order; return matches only."""
    corpus, asts, results, order = _evaluate_corpus(genotypes, corpus)
    out = []
    for name in order:
        if not results[name]:
            continue
        entry = corpus[name] if isinstance(corpus[name], dict) else {}
        out.append(_finding(name, entry, asts, genotypes, True))
    out.sort(key=lambda f: (-(f["magnitude"] or 0.0), _numeric_key(f["rsid"])))
    return out


def evaluate_all_verbose(genotypes: dict, corpus: dict | None = None) -> dict:
    """Split every genoset into matched / unmatched / incomplete buckets.

    ``incomplete`` holds genosets whose required rsIDs were not all genotyped
    (coverage < 1.0) so the UI can say "not testable on your array" instead of
    wrongly reporting the genoset as absent. A genoset may appear in both
    ``matched`` and ``incomplete`` when a partially covered rule still fired.
    ``unmatched`` only contains fully covered, genuinely absent genosets.
    """
    corpus, asts, results, order = _evaluate_corpus(genotypes, corpus)
    matched, unmatched, incomplete = [], [], []
    for name in order:
        entry = corpus[name] if isinstance(corpus[name], dict) else {}
        hit = results[name]
        finding = _finding(name, entry, asts, genotypes, hit)
        if hit:
            matched.append(finding)
        if finding["coverage"] < 1.0:
            incomplete.append(finding)
        elif not hit:
            unmatched.append(finding)
    matched.sort(key=lambda f: (-(f["magnitude"] or 0.0),
                                _numeric_key(f["rsid"])))
    return {"matched": matched, "unmatched": unmatched,
            "incomplete": incomplete}
