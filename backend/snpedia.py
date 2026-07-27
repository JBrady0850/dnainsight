"""
snpedia.py -- opt-in, local-only SNPedia harvester and lookup layer.

LICENCE POSITION, read this before changing anything in this file
=================================================================
SNPedia content is licensed CC-BY-NC-SA-3.0-US. DNAInsight is MIT and must
stay freely redistributable, therefore:

  * This module ships CODE ONLY. No SNPedia-derived data is bundled, vendored
    or committed to this repository, ever.
  * The harvested cache is written OUTSIDE the repository, under the user's
    own home directory (see CACHE_DIR and CACHE_PATH).
  * Harvesting requires explicit opt-in. Every entry point that would fetch
    content raises PermissionError carrying NOTICE until accept_license=True.
  * NOTICE is written into the cache's meta table when the cache is created,
    so the licence travels with the data, and export_cache() writes it out as
    a sibling LICENSE-SNPEDIA.txt.

PUBLISHER REQUIREMENTS honoured here
====================================
  * The target host is bots.snpedia.com. The www host sits behind a bot
    challenge and returns a small stub to scripts, so it is never used.
  * SNPs are ENUMERATED from Category:Is_a_snp (or the on-chip categories)
    before anything is fetched. Probing arbitrary rs numbers is explicitly
    forbidden and gets clients banned.
  * Page history and per-revision fetches are never requested.
  * Requests are rate limited, 2 per second by default and configurable, sent
    with a polite identifying User-Agent, and backed off exponentially on 429
    and 5xx responses while honouring Retry-After.

No network I/O happens at import time, so this module is fully testable
offline and every lookup degrades quietly when no cache exists.
"""

import json
import shutil
import sqlite3
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Callable

try:  # orientation.py may be added to the package after this module
    from backend import orientation as _orientation
except ImportError:  # pragma: no cover - only before orientation.py lands
    try:
        import orientation as _orientation  # type: ignore
    except ImportError:
        _orientation = None  # type: ignore


API_BASE = "https://bots.snpedia.com/api.php"
INDEX_BASE = "https://bots.snpedia.com/index.php"
PAGE_BASE = "https://www.snpedia.com/index.php"
USER_AGENT = (
    "DNAInsight/2.0 SNPedia harvester "
    "(personal, non-commercial; +https://github.com/JBrady0850/dnainsight)"
)

LICENSE_NAME = (
    "Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License"
)
LICENSE_URL = "http://creativecommons.org/licenses/by-nc-sa/3.0/us/"

SNP_CATEGORY = "Category:Is_a_snp"
GENOSET_CATEGORY = "Category:Is_a_genoset"
CHIP_CATEGORIES = {
    "23andme_v4": "Category:On_chip_23andMe_v4",
    "23andme_v5": "Category:On_chip_23andMe_v5",
    "ancestry_v2": "Category:On_chip_Ancestry_v2",
}
SCOPES = ("restricted", "chip_23andme_v5", "chip_ancestry_v2", "all")

DEFAULT_RATE_LIMIT = 2.0
CATEGORY_PAGE_LIMIT = 500
COMMIT_EVERY = 25
MAX_BACKOFF_SECONDS = 60.0
MAX_RETRIES = 5
SCHEMA_VERSION = "1"

NOTICE = """SNPedia content notice and licence
=================================
Content harvested by this module comes from SNPedia, https://www.snpedia.com .

SNPedia content is licensed under the Creative Commons Attribution-Noncommercial-Share
Alike 3.0 United States License, http://creativecommons.org/licenses/by-nc-sa/3.0/us/ .

By opting in you confirm that:

  * This cache is built for your own personal, non-commercial use.
  * The cache is stored outside the DNAInsight repository, in your home
    directory, and is not part of DNAInsight. DNAInsight itself is MIT
    licensed and ships no SNPedia data.
  * You will not redistribute this cache or any part of it as if it were
    yours, and you will not use it commercially.
  * If you do choose to share a database derived from this content, that
    database must itself carry this same Attribution-Noncommercial-Share
    Alike 3.0 United States licence and credit SNPedia.

Harvesting is rate limited and enumerates SNPedia's own category listings.
Do not modify this module to probe arbitrary rs numbers or to fetch page
histories; both are explicitly disallowed by the publisher.
"""

# The cache lives in the user's home directory, deliberately OUTSIDE the
# repository, so that no SNPedia-derived data can ever be committed.
CACHE_DIR = Path.home() / ".dnainsight"
CACHE_PATH = CACHE_DIR / "snpedia_cache.db"


def cache_path() -> Path:
    """Return the cache file path, outside the repository, under the home directory."""
    return Path(CACHE_PATH)


def _cache_dir() -> Path:
    """Return the directory holding the cache, following any redirected CACHE_PATH."""
    return cache_path().parent


# ---------------------------------------------------------------------------
# Cache schema
# ---------------------------------------------------------------------------
# "frequencies" is not in the required table list but is needed because the
# population diversity template has nowhere else to live; the required columns
# of the required tables are unchanged.
SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);
CREATE TABLE IF NOT EXISTS snps (
    rsid                    TEXT PRIMARY KEY,
    gene                    TEXT,
    chromosome              TEXT,
    position                INTEGER,
    assembly                TEXT,
    orientation             TEXT,
    stabilized_orientation  TEXT,
    gmaf                    REAL,
    max_magnitude           REAL,
    publications            INTEGER,
    summary                 TEXT,
    clinvar_sig             TEXT,
    clinvar_disease         TEXT,
    trait                   TEXT,
    topics                  TEXT,
    medicines               TEXT,
    conditions              TEXT,
    on_microarray           TEXT,
    modified                TEXT
);
CREATE TABLE IF NOT EXISTS genotypes (
    rsid      TEXT,
    token     TEXT,
    allele1   TEXT,
    allele2   TEXT,
    magnitude REAL,
    repute    TEXT,
    summary   TEXT,
    modified  TEXT,
    PRIMARY KEY (rsid, token)
);
CREATE TABLE IF NOT EXISTS genosets (
    name      TEXT PRIMARY KEY,
    magnitude REAL,
    repute    TEXT,
    summary   TEXT,
    criteria  TEXT,
    modified  TEXT
);
CREATE TABLE IF NOT EXISTS harvest_log (
    id       INTEGER PRIMARY KEY AUTOINCREMENT,
    started  TEXT,
    finished TEXT,
    scope    TEXT,
    pages    INTEGER,
    errors   INTEGER,
    note     TEXT
);
CREATE TABLE IF NOT EXISTS frequencies (
    rsid      TEXT,
    population TEXT,
    geno1_pct REAL,
    geno2_pct REAL,
    geno3_pct REAL,
    revision  TEXT,
    PRIMARY KEY (rsid, population)
);
CREATE INDEX IF NOT EXISTS idx_genotypes_rsid ON genotypes (rsid);
CREATE INDEX IF NOT EXISTS idx_snps_gene ON snps (gene);
"""


def init_cache() -> None:
    """Create the cache file, its tables and its indices, and store NOTICE.

    Safe to call repeatedly. The licence notice is written into the meta table
    at creation time so that the licence terms travel with the data even if
    the file is copied somewhere else.
    """
    directory = _cache_dir()
    directory.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(cache_path()))
    try:
        conn.executescript(SCHEMA)
        stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        rows = (
            ("schema_version", SCHEMA_VERSION),
            ("created", stamp),
            ("source", "SNPedia, https://www.snpedia.com"),
            ("license", LICENSE_NAME),
            ("license_url", LICENSE_URL),
            ("notice", NOTICE),
            ("redistribution", "Do not redistribute. Personal, non-commercial use only."),
        )
        conn.executemany("INSERT OR IGNORE INTO meta (key, value) VALUES (?, ?)", rows)
        conn.commit()
    finally:
        conn.close()


def _connect(create: bool = False) -> sqlite3.Connection | None:
    """Open the cache, or return None when it does not exist and create is False."""
    path = cache_path()
    if not path.exists():
        if not create:
            return None
        init_cache()
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    return conn


def _count(conn: sqlite3.Connection, table: str) -> int:
    """Return the row count of one table, or 0 when the table is missing."""
    try:
        row = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()
        return int(row["n"]) if row else 0
    except sqlite3.Error:
        return 0


def is_available() -> bool:
    """True when the cache exists and holds at least one harvested SNP row."""
    conn = _connect()
    if conn is None:
        return False
    try:
        return _count(conn, "snps") > 0
    finally:
        conn.close()


def cache_status() -> dict:
    """Describe the local cache without touching the network.

    A freshly initialised but unharvested cache reports available False,
    because an empty cache is of no use to the annotator.
    """
    path = cache_path()
    status = {
        "available": False,
        "path": str(path),
        "snps": 0,
        "genotypes": 0,
        "genosets": 0,
        "last_harvest": None,
        "scope": None,
        "license": LICENSE_NAME,
        "notice": NOTICE,
    }
    conn = _connect()
    if conn is None:
        return status
    try:
        status["snps"] = _count(conn, "snps")
        status["genotypes"] = _count(conn, "genotypes")
        status["genosets"] = _count(conn, "genosets")
        status["available"] = status["snps"] > 0
        try:
            row = conn.execute(
                "SELECT started, finished, scope FROM harvest_log "
                "ORDER BY id DESC LIMIT 1"
            ).fetchone()
        except sqlite3.Error:
            row = None
        if row is not None:
            status["last_harvest"] = row["finished"] or row["started"]
            status["scope"] = row["scope"]
    finally:
        conn.close()
    return status


# ---------------------------------------------------------------------------
# Polite HTTP layer
# ---------------------------------------------------------------------------
HTTP_TIMEOUT = 30


class _RateLimiter:
    """Simple wall-clock rate limiter, conservative by default."""

    def __init__(self, rate_per_second: float = DEFAULT_RATE_LIMIT) -> None:
        self.min_interval = (1.0 / rate_per_second) if rate_per_second and rate_per_second > 0 else 0.0
        self._last = 0.0

    def wait(self) -> None:
        """Block until enough time has passed since the previous request."""
        if self.min_interval <= 0:
            return
        elapsed = time.monotonic() - self._last
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last = time.monotonic()


_DEFAULT_LIMITER = _RateLimiter(DEFAULT_RATE_LIMIT)


def make_session() -> Any:
    """Return a requests Session with the polite User-Agent, or None.

    None means requests is unavailable and the urllib fallback will be used,
    which sends the same User-Agent header.
    """
    try:
        import requests  # imported lazily so the module works without it
    except ImportError:  # pragma: no cover - requests is in requirements.txt
        return None
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    return session


def _retry_after_seconds(value: Any) -> float | None:
    """Parse a Retry-After header value expressed in seconds."""
    if value is None:
        return None
    try:
        return max(0.0, float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def _http_get(url: str, params: dict | None, session: Any) -> tuple[int, str, Any]:
    """Perform one GET and return (status_code, body_text, retry_after_header)."""
    params = params or {}
    if session is not None:
        response = session.get(url, params=params, timeout=HTTP_TIMEOUT)
        return response.status_code, response.text, response.headers.get("Retry-After")
    target = url
    if params:
        target = f"{url}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(target, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT) as response:
            body = response.read().decode("utf-8", errors="replace")
            return int(getattr(response, "status", 200) or 200), body, response.headers.get("Retry-After")
    except urllib.error.HTTPError as exc:  # type: ignore[attr-defined]
        headers = getattr(exc, "headers", None)
        return int(exc.code), "", headers.get("Retry-After") if headers else None


def _request(
    url: str,
    params: dict | None = None,
    session: Any = None,
    limiter: _RateLimiter | None = None,
) -> str:
    """Rate-limited GET with exponential backoff on 429 and 5xx responses.

    Retry-After is honoured when the server sends it. A 4xx other than 429 is
    treated as fatal, because retrying a bad request is just rude.
    """
    limiter = limiter or _DEFAULT_LIMITER
    delay = 1.0
    last_error = ""
    for _attempt in range(MAX_RETRIES):
        limiter.wait()
        try:
            status, text, retry_after = _http_get(url, params, session)
        except Exception as exc:
            status, text, retry_after = 0, "", None
            last_error = str(exc)
        if status == 200:
            return text
        if status == 429 or status >= 500 or status == 0:
            pause = _retry_after_seconds(retry_after)
            if pause is None:
                pause = delay
            time.sleep(min(pause, MAX_BACKOFF_SECONDS))
            delay = min(delay * 2, MAX_BACKOFF_SECONDS)
            last_error = last_error or f"HTTP {status}"
            continue
        raise RuntimeError(f"SNPedia request failed with HTTP {status} for {url}")
    raise RuntimeError(
        f"SNPedia request gave up after {MAX_RETRIES} attempts for {url}: {last_error}"
    )


def _api_get(params: dict, session: Any = None, limiter: _RateLimiter | None = None) -> dict:
    """Call the bots.snpedia.com API and return the decoded JSON payload."""
    merged = dict(params)
    merged["format"] = "json"
    body = _request(API_BASE, merged, session=session, limiter=limiter)
    try:
        data = json.loads(body)
    except ValueError:
        if len(body) < 1000:
            raise RuntimeError(
                "SNPedia API returned a non-JSON stub. This is what the www host does to "
                f"scripts; make sure API_BASE is still {API_BASE}."
            )
        raise RuntimeError("SNPedia API returned a body that is not JSON.")
    if not isinstance(data, dict):
        raise RuntimeError("SNPedia API returned an unexpected JSON shape.")
    return data


# urllib.request already pulls urllib.error in as a side effect, but import it
# explicitly so the HTTPError branch in _http_get can never fail.
import urllib.error  # noqa: E402


# ---------------------------------------------------------------------------
# Sanctioned enumeration
# ---------------------------------------------------------------------------

def _enumerate_category(
    category: str,
    session: Any = None,
    limit: int | None = None,
    progress_cb: Callable[[int, int | None], None] | None = None,
    limiter: _RateLimiter | None = None,
) -> list[str]:
    """Walk a SNPedia category with continuation and return its page titles.

    This is the only sanctioned way to discover what exists. Enumerating a
    category is cheap for the server; probing rs numbers one by one is not,
    and is explicitly forbidden.
    """
    titles: list[str] = []
    continuation: str | None = None
    while True:
        params = {
            "action": "query",
            "list": "categorymembers",
            "cmtitle": category,
            "cmlimit": str(CATEGORY_PAGE_LIMIT),
        }
        if continuation:
            params["cmcontinue"] = continuation
        data = _api_get(params, session=session, limiter=limiter)
        members = ((data.get("query") or {}).get("categorymembers") or [])
        for member in members:
            title = str(member.get("title") or "").strip()
            if not title:
                continue
            titles.append(title)
            if limit and len(titles) >= limit:
                if progress_cb:
                    progress_cb(len(titles), limit)
                return titles[:limit]
        if progress_cb:
            progress_cb(len(titles), limit)
        continuation = (data.get("continue") or {}).get("cmcontinue")
        if not continuation:
            legacy = (data.get("query-continue") or {}).get("categorymembers") or {}
            continuation = legacy.get("cmcontinue")
        if not continuation:
            break
    return titles


def _titles_to_rsids(titles: list[str]) -> list[str]:
    """Reduce category page titles to a de-duplicated list of lower-case rsIDs."""
    out: list[str] = []
    seen: set[str] = set()
    for title in titles:
        rsid = title.strip().lower()
        if not rsid.startswith("rs") or "(" in rsid or ":" in rsid:
            continue
        if rsid in seen:
            continue
        seen.add(rsid)
        out.append(rsid)
    return out


def enumerate_snps(
    session: Any = None,
    limit: int | None = None,
    progress_cb: Callable[[int, int | None], None] | None = None,
) -> list[str]:
    """List every rsID SNPedia says it has, from Category:Is_a_snp.

    Always call this (or enumerate_chip_snps) before harvesting: the harvester
    must know which SNPs exist rather than guessing rs numbers.
    """
    titles = _enumerate_category(SNP_CATEGORY, session=session, limit=limit, progress_cb=progress_cb)
    return _titles_to_rsids(titles)


def enumerate_chip_snps(
    chip: str = "23andMe_v5",
    session: Any = None,
    limit: int | None = None,
    progress_cb: Callable[[int, int | None], None] | None = None,
) -> list[str]:
    """List the rsIDs SNPedia marks as present on one consumer chip.

    Accepted chips: "23andMe_v4", "23andMe_v5", "Ancestry_v2".
    """
    key = str(chip or "").strip().lower().replace("-", "_")
    category = CHIP_CATEGORIES.get(key)
    if category is None:
        raise ValueError(
            f"Unknown chip {chip!r}. Known chips: {', '.join(sorted(CHIP_CATEGORIES))}"
        )
    titles = _enumerate_category(category, session=session, limit=limit, progress_cb=progress_cb)
    return _titles_to_rsids(titles)


def enumerate_genosets(
    session: Any = None,
    limit: int | None = None,
    progress_cb: Callable[[int, int | None], None] | None = None,
) -> list[str]:
    """List genoset page titles from Category:Is_a_genoset."""
    titles = _enumerate_category(
        GENOSET_CATEGORY, session=session, limit=limit, progress_cb=progress_cb
    )
    return [t.strip() for t in titles if t.strip()]


# ---------------------------------------------------------------------------
# Page and property fetching
# ---------------------------------------------------------------------------

def page_title(rsid: str) -> str:
    """Return the SNPedia page title for an rsID ("rs53576" -> "Rs53576")."""
    text = str(rsid or "").strip()
    return text[:1].upper() + text[1:] if text else text


def _clean_dataitem(value: Any) -> str:
    """Strip Semantic MediaWiki dataitem decoration from one property value."""
    text = str(value if value is not None else "").strip()
    if "#" in text:
        head, _, tail = text.partition("#")
        marker = tail.replace("#", "").strip()
        if marker == "" or marker.isdigit():
            text = head
    return text.replace("_", " ").strip()


def fetch_subject(subject: str, session: Any = None) -> dict:
    """Fetch one page's semantic properties as a flat {property: value} dict.

    Uses action=browsebysubject, which returns the page's stored properties in
    a single request. Page history is never requested.
    """
    data = _api_get(
        {"action": "browsebysubject", "subject": str(subject or "").strip()},
        session=session,
    )
    properties: dict[str, Any] = {}
    for item in ((data.get("query") or {}).get("data") or []):
        name = str(item.get("property") or "").strip()
        if not name:
            continue
        values = [
            _clean_dataitem(entry.get("item"))
            for entry in (item.get("dataitem") or [])
            if entry.get("item") is not None
        ]
        values = [v for v in values if v != ""]
        if not values:
            continue
        properties[name] = values[0] if len(values) == 1 else values
    return properties


def fetch_wikitext(title: str, session: Any = None) -> str:
    """Fetch one page's raw wikitext (action=raw), never its history."""
    return _request(
        INDEX_BASE,
        {"title": str(title or "").strip(), "action": "raw"},
        session=session,
    )


def _as_percentages(values: list[str]) -> list[float] | None:
    """Coerce exactly three tokens to floats, or return None."""
    if len(values) != 3:
        return None
    out: list[float] = []
    for value in values:
        try:
            out.append(float(str(value).strip().rstrip("%")))
        except (TypeError, ValueError):
            return None
    return out


def parse_population_diversity(wikitext: str) -> dict:
    """Parse the {{population diversity}} template out of a SNP page.

    The template is positional after its three named genotype parameters:
    geno1/geno2/geno3, then repeating groups of a population code followed by
    three percentages, closing with a named HapMapRevision. Duplicate-looking
    codes such as HCB and CHB for the same Han Chinese sample are both kept as
    written, because callers need to see exactly what the page said.

    Returns {"genotypes": ["(C;C)", ...], "populations": {"CEU": [46.9, 44.2, 8.8]},
    "revision": str}. A page without the template yields empty containers.
    """
    result: dict[str, Any] = {"genotypes": [], "populations": {}, "revision": ""}
    text = str(wikitext or "")
    marker = text.lower().find("population diversity")
    if marker == -1:
        return result
    end = text.find("}}", marker)
    body = text[marker + len("population diversity"): end if end != -1 else len(text)]

    named: dict[str, str] = {}
    positional: list[str] = []
    for token in body.split("|"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            name, _, value = token.partition("=")
            named[name.strip().lower()] = value.strip()
            continue
        positional.append(token)

    result["genotypes"] = [named[key] for key in ("geno1", "geno2", "geno3") if named.get(key)]
    result["revision"] = named.get("hapmaprevision", "")

    index = 0
    while index < len(positional):
        code = positional[index]
        percentages = _as_percentages(positional[index + 1: index + 4])
        if percentages is None:
            index += 1
            continue
        result["populations"][code] = percentages
        index += 4
    return result


# ---------------------------------------------------------------------------
# Property mapping helpers
# ---------------------------------------------------------------------------

def _norm_key(name: str) -> str:
    """Normalise a property name for tolerant matching."""
    return str(name or "").strip().lower().replace(" ", "").replace("_", "")


def _first(properties: dict, *names: str) -> Any:
    """Return the first present property under any of the given names."""
    lookup = {_norm_key(k): v for k, v in (properties or {}).items()}
    for name in names:
        value = lookup.get(_norm_key(name))
        if value not in (None, "", []):
            return value
    return None


def _as_text(value: Any) -> str:
    """Flatten a property value to a single display string."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "; ".join(str(v).strip() for v in value if str(v).strip())
    return str(value).strip()


def _as_float(value: Any) -> float | None:
    """Best-effort float conversion, returning None instead of raising."""
    text = _as_text(value).split(";")[0].strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _as_int(value: Any) -> int | None:
    """Best-effort int conversion, returning None instead of raising."""
    number = _as_float(value)
    return int(number) if number is not None else None


def normalize_token(token: Any) -> str:
    """Normalise a SNPedia genotype token to the "(A;G)" form."""
    text = str(token or "").strip().upper().replace(" ", "")
    if text.startswith("(") and not text.endswith(")"):
        text = f"{text})"
    return text


def _token_alleles(token: str) -> tuple[str, str]:
    """Split a "(A;G)" token into its two allele strings."""
    inner = normalize_token(token).strip("()")
    parts = inner.split(";")
    if len(parts) != 2:
        return "", ""
    return parts[0].strip(), parts[1].strip()


def _geno_tokens(properties: dict) -> list[str]:
    """Return the genotype tokens a SNP page declares, in page order."""
    tokens: list[str] = []
    for name in ("Geno1", "Geno2", "Geno3"):
        value = _first(properties, name)
        if value is None:
            continue
        candidates = value if isinstance(value, (list, tuple)) else [value]
        for candidate in candidates:
            text = str(candidate).strip()
            if "(" in text:
                text = text[text.index("("):]
            token = normalize_token(text)
            if token.startswith("(") and token.endswith(")") and token not in tokens:
                tokens.append(token)
    return tokens


def _snp_row(rsid: str, properties: dict) -> tuple:
    """Build the snps table row for one SNP page's properties."""
    return (
        rsid,
        _as_text(_first(properties, "Gene", "Genes")),
        _as_text(_first(properties, "Chromosome", "Chr")),
        _as_int(_first(properties, "Position", "Chromosome position")),
        _as_text(_first(properties, "Assembly", "Reference")),
        _as_text(_first(properties, "Orientation")),
        _as_text(_first(properties, "StabilizedOrientation", "Stabilized orientation")),
        _as_float(_first(properties, "Gmaf", "GMAF")),
        _as_float(_first(properties, "MaxMagnitude", "Max magnitude")),
        _as_int(_first(properties, "Publications", "Publication count")),
        _as_text(_first(properties, "Summary", "Description")),
        _as_text(_first(properties, "ClinvarSig", "Clinvar significance", "Significance")),
        _as_text(_first(properties, "ClinvarDisease", "Clinvar disease", "Disease")),
        _as_text(_first(properties, "Trait", "Traits")),
        _as_text(_first(properties, "Topics", "Is a topic")),
        _as_text(_first(properties, "Medicines", "Medicine")),
        _as_text(_first(properties, "Conditions", "Condition")),
        _as_text(_first(properties, "OnMicroarray", "On chip", "Microarray")),
        _as_text(_first(properties, "Modification date", "Modified")),
    )


_SNP_INSERT = (
    "INSERT OR REPLACE INTO snps (rsid, gene, chromosome, position, assembly, "
    "orientation, stabilized_orientation, gmaf, max_magnitude, publications, summary, "
    "clinvar_sig, clinvar_disease, trait, topics, medicines, conditions, on_microarray, "
    "modified) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
)

_GENOTYPE_INSERT = (
    "INSERT OR REPLACE INTO genotypes (rsid, token, allele1, allele2, magnitude, "
    "repute, summary, modified) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
)

_FREQUENCY_INSERT = (
    "INSERT OR REPLACE INTO frequencies (rsid, population, geno1_pct, geno2_pct, "
    "geno3_pct, revision) VALUES (?, ?, ?, ?, ?, ?)"
)


def _genotype_row(rsid: str, token: str, properties: dict) -> tuple:
    """Build the genotypes table row for one genotype page's properties."""
    allele1, allele2 = _token_alleles(token)
    return (
        rsid,
        normalize_token(token),
        allele1,
        allele2,
        _as_float(_first(properties, "Magnitude")),
        _as_text(_first(properties, "Repute")),
        _as_text(_first(properties, "Summary", "Description")),
        _as_text(_first(properties, "Modification date", "Modified")),
    )


# ---------------------------------------------------------------------------
# Harvesting
# ---------------------------------------------------------------------------

def _now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _resolve_targets(
    scope: str,
    rsids: list[str] | None,
    session: Any,
    limiter: _RateLimiter,
    progress_cb: Callable[[int, int | None], None] | None,
) -> list[str]:
    """Turn a scope into the list of rsIDs to harvest, using enumeration only."""
    if scope == "restricted":
        if not rsids:
            raise ValueError(
                'scope "restricted" needs an explicit rsids list, normally the union of '
                "the bundled reference and the user's own findings."
            )
        return _titles_to_rsids([str(r) for r in rsids])
    if scope == "chip_23andme_v5":
        return enumerate_chip_snps("23andMe_v5", session=session, progress_cb=progress_cb)
    if scope == "chip_ancestry_v2":
        return enumerate_chip_snps("Ancestry_v2", session=session, progress_cb=progress_cb)
    return enumerate_snps(session=session, progress_cb=progress_cb)


def harvest(
    scope: str = "restricted",
    rsids: list[str] | None = None,
    accept_license: bool = False,
    max_requests: int | None = None,
    rate_limit: float = DEFAULT_RATE_LIMIT,
    progress_cb: Callable[[int, int | None], None] | None = None,
    force: bool = False,
    session: Any = None,
) -> dict:
    """Harvest SNPedia content into the local, out-of-repository cache.

    Scopes:
      "restricted"       only the rsIDs the caller passes. This is what
                         DNAInsight uses: the union of the bundled reference
                         and the user's own findings.
      "chip_23andme_v5"  everything SNPedia marks as on the 23andMe v5 chip.
      "chip_ancestry_v2" everything on the AncestryDNA v2 chip.
      "all"              every page in Category:Is_a_snp.

    Raises PermissionError carrying NOTICE unless accept_license is True; the
    licence has to be accepted explicitly every time, by the caller, in code
    the user can see.

    For each rsID this fetches the SNP page properties, then each declared
    genotype page, then the raw wikitext for the population frequency table.
    Rows are written transactionally and committed every 25 rsIDs, so an
    interrupted harvest resumes rather than starting over. rsIDs already in the
    cache are skipped unless force is True. max_requests bounds the number of
    content fetches, not the category enumeration that precedes them.
    """
    if accept_license is not True:
        raise PermissionError(NOTICE)
    if scope not in SCOPES:
        raise ValueError(f"Unknown scope {scope!r}. Known scopes: {', '.join(SCOPES)}")

    started = _now()
    clock = time.monotonic()
    limiter = _RateLimiter(rate_limit)
    if session is None:
        session = make_session()
    init_cache()

    targets = _resolve_targets(scope, rsids, session, limiter, progress_cb)
    requests_made = 0

    def budget_left() -> bool:
        return max_requests is None or requests_made < max_requests

    fetched = skipped = errors = genotype_rows = 0
    notes: list[str] = []
    conn = _connect(create=True)
    if conn is None:  # pragma: no cover - _connect(create=True) always returns
        raise RuntimeError("could not open the SNPedia cache")
    try:
        for index, rsid in enumerate(targets, start=1):
            if not budget_left():
                notes.append(f"stopped after {requests_made} requests: max_requests reached")
                break
            if not force:
                existing = conn.execute(
                    "SELECT 1 FROM snps WHERE rsid = ?", (rsid,)
                ).fetchone()
                if existing is not None:
                    skipped += 1
                    if progress_cb:
                        progress_cb(index, len(targets))
                    continue
            title = page_title(rsid)
            try:
                properties = fetch_subject(title, session=session)
                requests_made += 1
                conn.execute(_SNP_INSERT, _snp_row(rsid, properties))
                for token in _geno_tokens(properties):
                    if not budget_left():
                        break
                    geno_props = fetch_subject(f"{title}{token}", session=session)
                    requests_made += 1
                    conn.execute(_GENOTYPE_INSERT, _genotype_row(rsid, token, geno_props))
                    genotype_rows += 1
                if budget_left():
                    diversity = parse_population_diversity(fetch_wikitext(title, session=session))
                    requests_made += 1
                    for population, percentages in diversity["populations"].items():
                        conn.execute(_FREQUENCY_INSERT, (
                            rsid, population, percentages[0], percentages[1],
                            percentages[2], diversity["revision"],
                        ))
                fetched += 1
            except Exception as exc:
                errors += 1
                notes.append(f"{rsid}: {exc}")
                conn.rollback()
            if index % COMMIT_EVERY == 0:
                conn.commit()
            if progress_cb:
                progress_cb(index, len(targets))
        conn.commit()
        note = "; ".join(notes[:10]) if notes else "harvest completed"
        conn.execute(
            "INSERT INTO harvest_log (started, finished, scope, pages, errors, note) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (started, _now(), scope, requests_made, errors, note),
        )
        conn.commit()
    finally:
        conn.close()

    return {
        "requested": len(targets),
        "fetched": fetched,
        "skipped": skipped,
        "errors": errors,
        "genotypes": genotype_rows,
        "genosets": 0,
        "elapsed_seconds": round(time.monotonic() - clock, 3),
        "note": note,
    }


def harvest_genosets(
    accept_license: bool = False,
    max_requests: int | None = None,
    rate_limit: float = DEFAULT_RATE_LIMIT,
    progress_cb: Callable[[int, int | None], None] | None = None,
    force: bool = False,
    session: Any = None,
) -> dict:
    """Harvest genosets from Category:Is_a_genoset plus each /criteria subpage.

    Raises PermissionError carrying NOTICE unless accept_license is True.
    Genosets are stored in the same out-of-repository cache as SNPs.
    """
    if accept_license is not True:
        raise PermissionError(NOTICE)

    started = _now()
    clock = time.monotonic()
    limiter = _RateLimiter(rate_limit)
    if session is None:
        session = make_session()
    init_cache()

    names = enumerate_genosets(session=session, progress_cb=progress_cb)
    requests_made = 0
    fetched = skipped = errors = 0
    notes: list[str] = []

    conn = _connect(create=True)
    if conn is None:  # pragma: no cover
        raise RuntimeError("could not open the SNPedia cache")
    try:
        for index, name in enumerate(names, start=1):
            if max_requests is not None and requests_made >= max_requests:
                notes.append(f"stopped after {requests_made} requests: max_requests reached")
                break
            if not force:
                existing = conn.execute(
                    "SELECT 1 FROM genosets WHERE name = ?", (name,)
                ).fetchone()
                if existing is not None:
                    skipped += 1
                    continue
            try:
                properties = fetch_subject(name, session=session)
                requests_made += 1
                criteria = ""
                if max_requests is None or requests_made < max_requests:
                    try:
                        criteria = fetch_wikitext(f"{name}/criteria", session=session)
                    except Exception:
                        criteria = ""
                    requests_made += 1
                conn.execute(
                    "INSERT OR REPLACE INTO genosets (name, magnitude, repute, summary, "
                    "criteria, modified) VALUES (?, ?, ?, ?, ?, ?)",
                    (
                        name,
                        _as_float(_first(properties, "Magnitude")),
                        _as_text(_first(properties, "Repute")),
                        _as_text(_first(properties, "Summary", "Description")),
                        criteria,
                        _as_text(_first(properties, "Modification date", "Modified")),
                    ),
                )
                fetched += 1
            except Exception as exc:
                errors += 1
                notes.append(f"{name}: {exc}")
                conn.rollback()
            if index % COMMIT_EVERY == 0:
                conn.commit()
            if progress_cb:
                progress_cb(index, len(names))
        conn.commit()
        note = "; ".join(notes[:10]) if notes else "genoset harvest completed"
        conn.execute(
            "INSERT INTO harvest_log (started, finished, scope, pages, errors, note) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (started, _now(), "genosets", requests_made, errors, note),
        )
        conn.commit()
    finally:
        conn.close()

    return {
        "requested": len(names),
        "fetched": fetched,
        "skipped": skipped,
        "errors": errors,
        "genotypes": 0,
        "genosets": fetched,
        "elapsed_seconds": round(time.monotonic() - clock, 3),
        "note": note,
    }


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _fallback_tokens(a1: Any, a2: Any) -> list[str]:
    """Both allele orders and both strands, used when orientation.py is absent."""
    left = str(a1 or "").strip().upper()
    right = str(a2 or "").strip().upper()
    if not left or not right:
        return []
    pairs = [(left, right), (right, left)]
    comp_left = _COMPLEMENT.get(left, left)
    comp_right = _COMPLEMENT.get(right, right)
    pairs += [(comp_left, comp_right), (comp_right, comp_left)]
    tokens: list[str] = []
    for first, second in pairs:
        token = normalize_token(f"({first};{second})")
        if token not in tokens:
            tokens.append(token)
    return tokens


def _fallback_ambiguous(a1: Any, a2: Any) -> bool:
    """True for the strand-ambiguous A/T and C/G heterozygotes."""
    pair = {str(a1 or "").strip().upper(), str(a2 or "").strip().upper()}
    return pair in ({"A", "T"}, {"C", "G"})


def _token_candidates(
    a1: Any, a2: Any, stabilized: str | None = None, orientation: str | None = None
) -> tuple[list[str], bool, bool]:
    """Return (candidate_tokens, flipped, ambiguous), preferring orientation.py."""
    tokens: list[str] = []
    flipped = False
    ambiguous = _fallback_ambiguous(a1, a2)
    if _orientation is not None:
        try:
            result = _orientation.orient_to_snpedia(
                a1, a2, stabilized_orientation=stabilized, orientation=orientation
            )
            if isinstance(result, dict):
                token = normalize_token(result.get("token") or "")
                if token:
                    tokens.append(token)
                flipped = bool(result.get("flipped"))
                ambiguous = bool(result.get("ambiguous"))
        except Exception:
            pass
        try:
            for candidate in (_orientation.candidate_tokens(a1, a2) or []):
                token = normalize_token(candidate)
                if token and token not in tokens:
                    tokens.append(token)
        except Exception:
            pass
    for token in _fallback_tokens(a1, a2):
        if token not in tokens:
            tokens.append(token)
    return tokens, flipped, ambiguous


def lookup(rsid: str) -> dict | None:
    """Return the cached SNP-level row for one rsID, or None."""
    key = str(rsid or "").strip().lower()
    if not key:
        return None
    conn = _connect()
    if conn is None:
        return None
    try:
        row = conn.execute("SELECT * FROM snps WHERE rsid = ?", (key,)).fetchone()
        return dict(row) if row is not None else None
    except sqlite3.Error:
        return None
    finally:
        conn.close()


def lookup_genotype(rsid: str, a1: Any, a2: Any) -> dict | None:
    """Return the cached genotype row for one call, orientation aware.

    The cached stabilized_orientation is used to build the preferred token via
    orientation.orient_to_snpedia, then orientation.candidate_tokens and a
    local strand-flip fallback are tried in turn. "flipped" reports whether the
    matched token needed the complementary strand, and "ambiguous" flags the
    A/T and C/G pairs where strand cannot be resolved from the alleles alone.
    """
    key = str(rsid or "").strip().lower()
    if not key:
        return None
    conn = _connect()
    if conn is None:
        return None
    try:
        snp = conn.execute(
            "SELECT orientation, stabilized_orientation FROM snps WHERE rsid = ?", (key,)
        ).fetchone()
        stabilized = snp["stabilized_orientation"] if snp is not None else None
        orientation = snp["orientation"] if snp is not None else None
        tokens, flipped, ambiguous = _token_candidates(a1, a2, stabilized, orientation)
        given = sorted([str(a1 or "").strip().upper(), str(a2 or "").strip().upper()])
        for token in tokens:
            row = conn.execute(
                "SELECT * FROM genotypes WHERE rsid = ? AND token = ?", (key, token)
            ).fetchone()
            if row is None:
                continue
            matched = sorted([str(row["allele1"] or "").upper(), str(row["allele2"] or "").upper()])
            return {
                "rsid": key,
                "token": row["token"],
                "allele1": row["allele1"],
                "allele2": row["allele2"],
                "magnitude": row["magnitude"],
                "repute": row["repute"] or "",
                "summary": row["summary"] or "",
                "modified": row["modified"] or "",
                "flipped": bool(flipped or matched != given),
                "ambiguous": bool(ambiguous),
            }
        return None
    except sqlite3.Error:
        return None
    finally:
        conn.close()


def lookup_genoset(name: str) -> dict | None:
    """Return the cached genoset row for one genoset name, or None."""
    key = str(name or "").strip()
    if not key:
        return None
    conn = _connect()
    if conn is None:
        return None
    try:
        row = conn.execute("SELECT * FROM genosets WHERE name = ?", (key,)).fetchone()
        if row is None:
            row = conn.execute(
                "SELECT * FROM genosets WHERE LOWER(name) = ?", (key.lower(),)
            ).fetchone()
        return dict(row) if row is not None else None
    except sqlite3.Error:
        return None
    finally:
        conn.close()


ANNOTATION_KEYS = (
    "magnitude", "repute", "summary", "max_magnitude", "publications", "gmaf",
    "snpedia_topics", "snpedia_medicines", "snpedia_conditions", "orientation",
    "stabilized_orientation", "flipped", "ambiguous", "snpedia_url",
)


def _split_list(value: Any) -> list[str]:
    """Split a stored semicolon-separated property into a clean list."""
    text = str(value or "").strip()
    if not text:
        return []
    return [part.strip() for part in text.split(";") if part.strip()]


def annotate(finding: dict) -> dict:
    """Add SNPedia fields to a finding dict, in place, and never raise.

    Every key in ANNOTATION_KEYS is guaranteed to exist afterwards. When no
    cache is present, or the rsID is not cached, the keys are left as None or
    empty rather than being dropped, so callers can render a finding without
    checking whether the user ever opted in. Existing values on the finding are
    preserved; only fields the cache can fill are overwritten.
    """
    if not isinstance(finding, dict):
        return {key: None for key in ANNOTATION_KEYS}
    defaults = {
        "magnitude": None,
        "repute": "",
        "summary": "",
        "max_magnitude": None,
        "publications": None,
        "gmaf": None,
        "snpedia_topics": [],
        "snpedia_medicines": [],
        "snpedia_conditions": [],
        "orientation": None,
        "stabilized_orientation": None,
        "flipped": False,
        "ambiguous": False,
        "snpedia_url": None,
    }
    for key, value in defaults.items():
        finding.setdefault(key, value)
    try:
        rsid = str(finding.get("rsid") or "").strip().lower()
        if rsid.startswith("rs"):
            finding["snpedia_url"] = f"{PAGE_BASE}/{page_title(rsid)}"
        if not rsid:
            return finding
        row = lookup(rsid)
        if row is not None:
            finding["max_magnitude"] = row.get("max_magnitude")
            finding["publications"] = row.get("publications")
            finding["gmaf"] = row.get("gmaf")
            finding["snpedia_topics"] = _split_list(row.get("topics"))
            finding["snpedia_medicines"] = _split_list(row.get("medicines"))
            finding["snpedia_conditions"] = _split_list(row.get("conditions"))
            finding["orientation"] = row.get("orientation") or None
            finding["stabilized_orientation"] = row.get("stabilized_orientation") or None
            if not finding.get("summary"):
                finding["summary"] = row.get("summary") or ""
        genotype = lookup_genotype(rsid, finding.get("allele1"), finding.get("allele2"))
        if genotype is not None:
            finding["magnitude"] = genotype.get("magnitude")
            finding["repute"] = genotype.get("repute") or ""
            finding["summary"] = genotype.get("summary") or finding.get("summary") or ""
            finding["flipped"] = bool(genotype.get("flipped"))
            finding["ambiguous"] = bool(genotype.get("ambiguous"))
    except Exception:
        return finding
    return finding


# ---------------------------------------------------------------------------
# Cache management
# ---------------------------------------------------------------------------

def purge_cache() -> dict:
    """Delete the local cache file and report what was removed."""
    path = cache_path()
    before = cache_status()
    removed = False
    if path.exists():
        path.unlink()
        removed = True
    return {
        "removed": removed,
        "path": str(path),
        "snps": before["snps"],
        "genotypes": before["genotypes"],
        "genosets": before["genosets"],
    }


def export_cache(dest: Any) -> dict:
    """Copy the cache to ``dest`` and write a sibling LICENSE-SNPEDIA.txt.

    The licence file is not optional: share-alike means anyone the user passes
    the cache to must receive the same terms, so the copy always travels with
    NOTICE next to it.
    """
    source = cache_path()
    if not source.exists():
        raise FileNotFoundError(f"no SNPedia cache to export at {source}")
    target = Path(dest)
    if target.is_dir():
        target = target / source.name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(target))
    license_file = target.parent / "LICENSE-SNPEDIA.txt"
    license_file.write_text(NOTICE, encoding="utf-8")
    return {
        "exported": True,
        "source": str(source),
        "destination": str(target),
        "license_file": str(license_file),
        "bytes": target.stat().st_size,
        "license": LICENSE_NAME,
        "license_url": LICENSE_URL,
    }
