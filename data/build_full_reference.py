"""
build_full_reference.py -- Tier 2 full-array variant reference builder.

WHY THIS FILE EXISTS
--------------------
Tier 1 is the curated in-repo reference: the hand-written table in
data/build_reference.py plus the evidence overlay in data/evidence_overlay.py.
It is small, reviewed, and shipped inside the repository. Nothing here touches
it.

Tier 2 is this: every variant a consumer array actually genotypes that also has
a real record in a public database. It is built locally, it is large, and it is
gitignored. GitHub rejects files over 100 MB and data/reference.db goes well
past that, so the database is never committed. It is fully reproducible by
re-running this script, which is why not committing it costs nothing.

SOURCES AND LICENCES
--------------------
Only CC0 and US public domain material is used, so the repository itself stays
redistributable under MIT. data/DATA_SOURCES.md carries the full record.

    ClinVar         NCBI, US public domain, attribution requested not required
    GWAS Catalog    EMBL-EBI terms, open for research and commercial reuse
    CPIC            CC0-1.0

Deliberately NOT used: SNPedia, which is CC-BY-NC-SA-3.0-US and would
relicense this repository, and the PharmGKB or ClinPGx bulk downloads, whose
terms add a no-commercial-sale clause on top of an otherwise open licence.

MEMORY DISCIPLINE
-----------------
variant_summary.txt.gz is hundreds of megabytes decompressed, so every download
is streamed through gzip.GzipFile and parsed line by line. Nothing is ever read
whole into memory. Accumulation happens in SQLite itself through an upsert, not
in a Python dict, so peak memory stays flat no matter how large the source is.
The one exception is the GWAS pass, which has to count distinct studies per
rsID before it can decide whether to keep the rsID at all.

USAGE
-----
    python data/build_full_reference.py --array-file uploads/mine.txt
    python data/build_full_reference.py --limit 20000 --skip-gwas
    python data/build_full_reference.py --stats

READING IT BACK
---------------
The scanner consumes this through the three helpers at the bottom of this
module, so nothing in backend/ needs to know how the database was built:

    from data.build_full_reference import open_reference_db, lookup, coverage_stats
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sqlite3
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# The star mapping is NOT redefined here. backend/scoring.py owns it, and it is
# the single authoritative copy for the whole project. Two copies of a
# vocabulary this fiddly would drift within one release.
#
# Worth restating because it bites every time: the published ClinVar
# documentation lists "no classification for the individual variant", but the
# string that actually appears in variant_summary.txt.gz is "no classification
# for the single variant". REVIEW_STATUS_STARS already carries both spellings,
# along with the somatic and submitted-record vocabularies, so this builder can
# hand it a verbatim data string and get the right star count.
from backend.scoring import (  # noqa: E402
    REVIEW_STATUS_STARS,
    clinvar_sig_code,
    review_stars,
)

try:  # pragma: no cover - exercised only when requests is absent
    import requests
except ImportError:  # pragma: no cover
    requests = None  # type: ignore[assignment]

__all__ = [
    "CLINVAR_URL", "GWAS_URL", "GWAS_FALLBACK_URLS", "CPIC_URL",
    "COL_RSID", "COL_GENE", "COL_PHENOTYPE_IDS", "COL_POSITION_VCF",
    "COL_ASSEMBLY", "COL_CHROMOSOME", "COL_CLINICAL_SIG", "COL_REVIEW_STATUS",
    "COL_PHENOTYPE_LIST", "CLINVAR_COLUMNS", "CLINVAR_MANDATORY",
    "GWAS_COLUMNS", "GWAS_MANDATORY", "GWAS_MLOG_THRESHOLD",
    "GWAS_MIN_STUDIES", "TARGET_ASSEMBLY", "SCHEMA_DDL", "VARIANT_COLUMNS",
    "SOURCE_LICENCES", "DEFAULT_DB_PATH", "DEFAULT_UPLOADS_DIR",
    "DEFAULT_CACHE_DIR",
    "build_parser", "create_schema", "load_array_rsids", "resolve_columns",
    "stream_gzip_lines", "stream_text_lines", "stream_zip_lines",
    "gwas_source_lines", "build_clinvar", "build_gwas",
    "build_cpic", "write_meta", "open_reference_db", "lookup",
    "coverage_stats", "main",
]


# ---------------------------------------------------------------------------
# Source endpoints
# ---------------------------------------------------------------------------
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
# The documented endpoint, tried first. VERIFIED 404 ON 2026-07-27: EBI retired
# the /api/search/downloads/alternative route and now publishes the same
# ontology-annotated association table as a zip on the release FTP mirror. The
# documented URL is kept as the primary candidate so that a restored service is
# picked up automatically and so the recorded endpoint stays visible, with the
# live mirror behind it. GWAS is the one source here that is a moving target, so
# it gets a candidate list rather than a single address.
GWAS_URL = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
GWAS_FALLBACK_URLS: tuple[str, ...] = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations-full.zip",
)

# CPIC exposes its pair table through a PostgREST API. The select list is part
# of the request on purpose: it means the columns this project must not copy are
# never even transferred, let alone written.
#
# cpiclevel is CPIC's own assignment and is CC0-1.0 along with the rest of the
# CPIC database. clinpgxlevel and pgxtesting are NOT. They arrive inside the
# same CC0 dump, which makes them look safe, but they are PharmGKB-sourced and
# carry PharmGKB's no-commercial-sale clause. Copying either one into
# reference.db would attach that clause to a database built by an MIT project,
# so they are excluded at the wire level rather than filtered later.
CPIC_URL = "https://api.cpicpgx.org/v1/pair_view"
CPIC_SELECT = "genesymbol,drugname,cpiclevel"

USER_AGENT = "DNAInsight/2.0 (+https://github.com/dnainsight) reference builder"
_HEADERS = {"User-Agent": USER_AGENT}

# ---------------------------------------------------------------------------
# ClinVar column names, verbatim
#
# These are literal header strings from variant_summary.txt.gz and they are
# deliberately ugly. Every one of them has been mistyped at least once during
# development, so they are constants and the tests assert the exact spelling:
#
#   "RS# (dbSNP)"   has a space and a hash. It is not "RS", "rsid" or "RS#".
#   "#Symbol"       is gene_specific_summary.txt's gene column, hash included.
#                   variant_summary.txt.gz calls the same field "GeneSymbol".
#   "PhenotypeIDS"  ends in a capital S. It is not "PhenotypeIDs".
#
# Column POSITIONS are never hard-coded. ClinVar adds columns between releases
# and has done so twice in the somatic-classification rollout, so the header
# line is parsed and resolved by name every run.
# ---------------------------------------------------------------------------
COL_RSID = "RS# (dbSNP)"
COL_GENE = "#Symbol"
COL_PHENOTYPE_IDS = "PhenotypeIDS"
COL_PHENOTYPE_LIST = "PhenotypeList"
COL_POSITION_VCF = "PositionVCF"
COL_ASSEMBLY = "Assembly"
COL_CHROMOSOME = "Chromosome"
COL_CLINICAL_SIG = "ClinicalSignificance"
COL_REVIEW_STATUS = "ReviewStatus"

# Accepted spellings per logical field, most-expected first.
#
# CORRECTION, VERIFIED AGAINST THE LIVE FILE ON 2026-07-27. The handoff note in
# the "Tier 2 reference builder" traps list in CONTRIBUTING.md recorded the
# gene column as "#Symbol" before this was checked against the live file. That is
# true of ClinVar's gene_specific_summary.txt, but variant_summary.txt.gz
# publishes the gene as "GeneSymbol" in a 43 column header whose only hashed
# name is the leading "#AlleleID". COL_GENE keeps the documented spelling so the
# recorded trap stays visible, and "GeneSymbol" is listed alongside it, so this
# builder resolves the column correctly against either file. That is also why
# resolution is by name and never by position: two ClinVar tables that describe
# the same genes disagree about what the column is called.
CLINVAR_COLUMNS: dict[str, tuple[str, ...]] = {
    "rsid":          (COL_RSID, "RS# (dbSNP)"),
    "gene":          (COL_GENE, "GeneSymbol", "Symbol"),
    "clinical_sig":  (COL_CLINICAL_SIG, "GermlineClassification"),
    "review_status": (COL_REVIEW_STATUS, "GermlineReviewStatus"),
    "assembly":      (COL_ASSEMBLY,),
    "chromosome":    (COL_CHROMOSOME,),
    "position":      (COL_POSITION_VCF,),
    "condition":     (COL_PHENOTYPE_LIST,),
    "phenotype_ids": (COL_PHENOTYPE_IDS,),
}
CLINVAR_MANDATORY: tuple[str, ...] = ("rsid", "assembly", "clinical_sig", "review_status")

# Start and Stop are right-shifted relative to the VCF convention while
# PositionVCF is left-shifted, so the two disagree by one for any indel. Only
# PositionVCF is used for coordinates. A one-base error is enough to break a
# join against an array manifest.
TARGET_ASSEMBLY = "GRCh37"

# An rsID of -1 means ClinVar has no dbSNP mapping for that allele. It is not a
# real identifier and thousands of rows carry it, so those rows are dropped.
NO_DBSNP_MAPPING = "-1"

# ---------------------------------------------------------------------------
# GWAS Catalog column names
#
# P-VALUE underflows to 0.0 for thousands of rows because the published value is
# below the smallest double, so significance is read from PVALUE_MLOG, the
# negative log10, which does not underflow. 7.3 is -log10(5e-8), the field's
# genome-wide significance line.
#
# The sign of an effect size is not in OR or BETA. It lives as free text in
# "95% CI (TEXT)", which reads like "[1.01-1.05] unit decrease". A build that
# trusts the numeric column alone reports protective alleles as risk alleles.
#
# Roughly 40 percent of MAPPED_TRAIT_URI values are not EFO_ terms; they are
# Orphanet, MONDO, HP and NCIT URIs. Filtering on an EFO_ prefix would silently
# discard about two fifths of the catalogue, so the trait text is kept and the
# URI is not used as a gate.
# ---------------------------------------------------------------------------
GWAS_COLUMNS: dict[str, tuple[str, ...]] = {
    "rsids":        ("SNPS",),
    "mlog":         ("PVALUE_MLOG",),
    "study":        ("STUDY ACCESSION", "STUDY_ACCESSION"),
    "trait":        ("MAPPED_TRAIT",),
    "trait_raw":    ("DISEASE/TRAIT",),
    "trait_uri":    ("MAPPED_TRAIT_URI",),
    "risk_allele":  ("STRONGEST SNP-RISK ALLELE",),
    "ci_text":      ("95% CI (TEXT)",),
    "pubmed":       ("PUBMEDID",),
}
GWAS_MANDATORY: tuple[str, ...] = ("rsids", "mlog", "study")
GWAS_MLOG_THRESHOLD = 7.3
GWAS_MIN_STUDIES = 2

# ---------------------------------------------------------------------------
# Output location and schema
# ---------------------------------------------------------------------------
DEFAULT_DB_PATH = ROOT / "data" / "reference.db"
DEFAULT_UPLOADS_DIR = ROOT / "uploads"

# Scratch space for archives that cannot be streamed from a socket. Already
# gitignored as data/_cache/, and nothing in here is ever required to exist.
DEFAULT_CACHE_DIR = ROOT / "data" / "_cache"

# Field names match docs/API_V2.md section 2.3 so a row can be merged into a
# finding without a translation layer: clinical_sig, clinvar_sig_code,
# review_status, review_stars, cpic_level, publications and risk_allele all
# carry the documented meanings and ranges.
SCHEMA_DDL: tuple[str, ...] = (
    """
    CREATE TABLE IF NOT EXISTS variants (
        rsid             TEXT PRIMARY KEY,
        gene             TEXT,
        chromosome       TEXT,
        position         INTEGER,
        clinical_sig     TEXT,
        clinvar_sig_code INTEGER,
        review_status    TEXT,
        review_stars     INTEGER DEFAULT 0,
        condition        TEXT,
        cpic_level       TEXT,
        gwas_traits      TEXT,
        gwas_studies     INTEGER DEFAULT 0,
        publications     INTEGER DEFAULT 0,
        risk_allele      TEXT,
        source           TEXT
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene)",
    "CREATE INDEX IF NOT EXISTS idx_variants_review_stars ON variants(review_stars)",
    """
    CREATE TABLE IF NOT EXISTS meta (
        key   TEXT PRIMARY KEY,
        value TEXT
    )
    """,
)

VARIANT_COLUMNS: tuple[str, ...] = (
    "rsid", "gene", "chromosome", "position", "clinical_sig",
    "clinvar_sig_code", "review_status", "review_stars", "condition",
    "cpic_level", "gwas_traits", "gwas_studies", "publications",
    "risk_allele", "source",
)

# Recorded verbatim into the meta table on every build, so a database found on
# disk months later still states what may be done with its contents.
SOURCE_LICENCES: dict[str, str] = {
    "clinvar": (
        "ClinVar (NCBI). US Government work, public domain. NCBI requests "
        "but does not require attribution. See data/DATA_SOURCES.md."
    ),
    "gwas": (
        "GWAS Catalog (EMBL-EBI and NHGRI). EMBL-EBI terms of use: freely "
        "available, no restriction on research or commercial reuse, citation "
        "of the catalogue requested. See data/DATA_SOURCES.md."
    ),
    "cpic": (
        "CPIC (cpicpgx.org). CC0-1.0. Only genesymbol, drugname and cpiclevel "
        "are copied. The PharmGKB-sourced clinpgxlevel and pgxtesting columns "
        "are excluded because they carry a no-commercial-sale clause."
    ),
    "array": (
        "Array coverage set derived from the user's own raw data export. Never "
        "redistributed. No vendor manifest is downloaded or bundled."
    ),
}

BATCH_SIZE = 50000
PROGRESS_EVERY = 100000
HTTP_TIMEOUT = 180


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the Tier 2 builder.

    Exposed as a function rather than built inline in main so the test suite can
    assert every documented flag parses without running a build.
    """
    parser = argparse.ArgumentParser(
        prog="build_full_reference.py",
        description="Build the Tier 2 full-array variant reference database.",
    )
    parser.add_argument("--limit", type=int, default=0,
                        help="stop after N parsed rows per source, 0 means no limit")
    parser.add_argument("--skip-clinvar", action="store_true",
                        help="do not download or parse ClinVar")
    parser.add_argument("--skip-gwas", action="store_true",
                        help="do not download or parse the GWAS Catalog")
    parser.add_argument("--skip-cpic", action="store_true",
                        help="do not download or parse the CPIC pair table")
    parser.add_argument("--array-file", default=None,
                        help="path to a raw array export whose rsID column defines "
                             "the coverage set, or a directory of such files")
    parser.add_argument("--out", default=None,
                        help=f"output database path, default {DEFAULT_DB_PATH}")
    parser.add_argument("--stats", action="store_true",
                        help="print the meta row counts when the build finishes")
    parser.add_argument("--resume", action="store_true",
                        help="keep sources already recorded in the target database")
    return parser


# ---------------------------------------------------------------------------
# Schema and column resolution
# ---------------------------------------------------------------------------

def create_schema(conn: sqlite3.Connection) -> None:
    """Create the variants and meta tables and their indices if absent.

    WAL is set here rather than at connect time because the journal mode is a
    property of the database file, not of the connection, and a build that is
    interrupted mid-transaction should leave a readable file behind.
    """
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    for statement in SCHEMA_DDL:
        conn.execute(statement)
    conn.commit()


def resolve_columns(header: Iterable[str],
                    wanted: dict[str, tuple[str, ...]],
                    mandatory: Iterable[str] = ()) -> dict[str, int]:
    """Map logical field names to positions in a parsed header line.

    Every candidate spelling for a field is tried in order, and the first one
    present in the header wins. Fields with no matching column are simply absent
    from the returned dict, unless they are listed in mandatory, in which case a
    KeyError names them. Positions are never assumed.
    """
    index_of = {name: i for i, name in enumerate(header)}
    resolved: dict[str, int] = {}
    for field, candidates in wanted.items():
        for candidate in candidates:
            if candidate in index_of:
                resolved[field] = index_of[candidate]
                break
    missing = [field for field in mandatory if field not in resolved]
    if missing:
        raise KeyError(
            "header is missing required column(s) for: " + ", ".join(sorted(missing))
        )
    return resolved


def _cell(row: list[str], columns: dict[str, int], field: str) -> str:
    """Return one trimmed cell, or an empty string when the column is absent."""
    idx = columns.get(field)
    if idx is None or idx >= len(row):
        return ""
    value = row[idx].strip()
    return "" if value in ("na", "NA", "-", "NULL") else value


def _int_or_none(value: str) -> int | None:
    """Coerce a cell to int, returning None for blanks and junk."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float_or_none(value: str) -> float | None:
    """Coerce a cell to float, returning None for blanks and junk."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Array coverage set
# ---------------------------------------------------------------------------

def _rsids_from_file(path: Path) -> set[str]:
    """Extract the rsID column from one raw array export.

    Handles the 23andMe tab-delimited layout, the AncestryDNA layout with its
    split allele columns, and comma-separated variants of either. Header and
    comment lines are skipped. Vendor-internal identifiers such as i3000001 are
    not rsIDs and are dropped.
    """
    found: set[str] = set()
    try:
        handle = path.open("r", encoding="utf-8", errors="replace")
    except OSError:
        return found
    with handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                first = line.split("\t", 1)[0]
            elif "," in line:
                first = line.split(",", 1)[0]
            else:
                first = line.split(None, 1)[0]
            first = first.strip().strip('"').lower()
            if first.startswith("rs") and first[2:].isdigit():
                found.add(first)
    return found


def load_array_rsids(array_file: str | None = None,
                     uploads_dir: Path | None = None) -> set[str]:
    """Return the rsID set a consumer array genotypes, from the user's own data.

    No vendor manifest is downloaded. Illumina and Affymetrix manifests are not
    redistributable and fetching one would put a non-free file inside a build
    that is meant to be reproducible from open sources alone.

    With --array-file the given path is used, and a directory is accepted as
    well as a single file. Without it, every .txt in uploads/ is unioned, which
    is the practical superset of the arrays this installation has actually seen.
    An empty set means build without a coverage filter.
    """
    targets: list[Path] = []
    if array_file:
        candidate = Path(array_file)
        if candidate.is_dir():
            targets = sorted(candidate.glob("*.txt"))
        elif candidate.exists():
            targets = [candidate]
        else:
            print(f"  WARNING: --array-file {candidate} does not exist")
    else:
        base = uploads_dir if uploads_dir is not None else DEFAULT_UPLOADS_DIR
        if base.is_dir():
            targets = sorted(base.glob("*.txt"))

    union: set[str] = set()
    for target in targets:
        got = _rsids_from_file(target)
        if got:
            print(f"  coverage +{len(got):>9,} rsIDs from {target.name}")
        union |= got
    return union


# ---------------------------------------------------------------------------
# Streaming downloads
# ---------------------------------------------------------------------------

def stream_gzip_lines(url: str, timeout: int = HTTP_TIMEOUT) -> Iterator[str]:
    """Yield text lines from a remote gzipped file without buffering it.

    The response body is handed straight to gzip.GzipFile through the raw
    socket, so a several hundred megabyte decompressed stream costs one small
    buffer rather than several hundred megabytes of resident memory.
    """
    if requests is None:
        raise RuntimeError("the requests package is not installed")
    with requests.get(url, stream=True, timeout=timeout, headers=_HEADERS) as resp:
        resp.raise_for_status()
        resp.raw.decode_content = True
        with gzip.GzipFile(fileobj=resp.raw) as gz:
            wrapper = io.TextIOWrapper(gz, encoding="utf-8", errors="replace")
            for line in wrapper:
                yield line.rstrip("\r\n")


def stream_text_lines(url: str, timeout: int = HTTP_TIMEOUT,
                      params: dict[str, str] | None = None) -> Iterator[str]:
    """Yield text lines from a remote uncompressed file, one at a time."""
    if requests is None:
        raise RuntimeError("the requests package is not installed")
    with requests.get(url, stream=True, timeout=timeout,
                      headers=_HEADERS, params=params) as resp:
        resp.raise_for_status()
        if not resp.encoding:
            resp.encoding = "utf-8"
        for line in resp.iter_lines(decode_unicode=True):
            if line is None:
                continue
            yield str(line).rstrip("\r\n")


def stream_zip_lines(url: str, timeout: int = HTTP_TIMEOUT,
                     cache_dir: Path | None = None) -> Iterator[str]:
    """Yield text lines from the first readable member of a remote zip archive.

    A zip cannot be decompressed from a forward-only socket the way a gzip can,
    because the archive index sits at the end of the file. The download is
    therefore written to data/_cache/ in one megabyte chunks and the member is
    then read as a stream, so neither the archive nor the member is ever held in
    memory. data/_cache/ is gitignored and may be deleted at any time.
    """
    if requests is None:
        raise RuntimeError("the requests package is not installed")
    cache = cache_dir if cache_dir is not None else DEFAULT_CACHE_DIR
    cache.mkdir(parents=True, exist_ok=True)
    local = cache / url.rsplit("/", 1)[-1]

    with requests.get(url, stream=True, timeout=timeout, headers=_HEADERS) as resp:
        resp.raise_for_status()
        written = 0
        with open(local, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                if not chunk:
                    continue
                fh.write(chunk)
                written += len(chunk)
    print(f"  cached {written / 1048576:.1f} MB to {local}")

    with zipfile.ZipFile(local) as archive:
        members = [n for n in archive.namelist()
                   if n.lower().endswith((".tsv", ".txt", ".csv"))]
        if not members:
            members = [n for n in archive.namelist() if not n.endswith("/")]
        if not members:
            raise ValueError(f"{local.name} contains no readable member")
        print(f"  reading member {members[0]}")
        with archive.open(members[0]) as member:
            wrapper = io.TextIOWrapper(member, encoding="utf-8", errors="replace")
            for line in wrapper:
                yield line.rstrip("\r\n")


def gwas_source_lines(urls: Iterable[str]) -> Iterator[str]:
    """Yield lines from the first GWAS association download that responds.

    Each candidate is probed by pulling its first line, so a dead endpoint is
    detected before any parsing state is built. Zip candidates are recognised by
    extension. If every candidate fails the collected errors are raised together,
    which the caller turns into a skipped source rather than a failed build.
    """
    errors: list[str] = []
    for candidate in urls:
        reader = (stream_zip_lines(candidate) if candidate.lower().endswith(".zip")
                  else stream_text_lines(candidate))
        try:
            first = next(reader)
        except StopIteration:
            errors.append(f"{candidate} -> empty response")
            print(f"  WARNING: GWAS candidate returned nothing: {candidate}")
            continue
        except Exception as exc:
            errors.append(f"{candidate} -> {type(exc).__name__}: {exc}")
            print(f"  WARNING: GWAS candidate unavailable: {errors[-1]}")
            continue
        print(f"[gwas] using {candidate}")
        yield first
        for line in reader:
            yield line
        return
    raise RuntimeError("no GWAS association download responded. " + " | ".join(errors))


# ---------------------------------------------------------------------------
# Upsert construction
#
# Accumulation happens in SQLite, not in a Python dict, so memory stays flat.
# The consequence is that the same rsID arrives more than once and the rows have
# to be reconciled in SQL. A variant is not always two rows: an X or Y
# pseudoautosomal region variant appears four times, once per assembly per
# chromosome copy, and multi-gene submissions add more. The reconciliation rule
# is highest review stars wins, tie-broken by clinical actionability, so the
# surviving row is the best-attested one rather than whichever arrived last.
# ---------------------------------------------------------------------------

def _sig_rank_sql(alias: str) -> str:
    """Return SQL ranking a clinvar_sig_code by clinical actionability.

    Numeric CLNSIG codes are identifiers, not an ordered scale: 255 means
    "other" and must never outrank 5, which means pathogenic. This maps them
    onto a real order so the tie-break is meaningful.
    """
    return (
        "CASE {a}.clinvar_sig_code WHEN 5 THEN 6 WHEN 4 THEN 5 WHEN 6 THEN 4 "
        "WHEN 7 THEN 3 WHEN 1 THEN 2 ELSE 1 END"
    ).format(a=alias)


def _clinvar_upsert_sql() -> str:
    """Return the ClinVar upsert statement, built from VARIANT_COLUMNS.

    Columns ClinVar does not supply are left out of the update clause so a later
    ClinVar row cannot blank out GWAS or CPIC values already written for the
    same rsID.
    """
    untouched = ("rsid", "cpic_level", "gwas_traits", "gwas_studies",
                 "publications", "risk_allele", "source")
    cols = ", ".join(VARIANT_COLUMNS)
    marks = ", ".join("?" for _ in VARIANT_COLUMNS)
    sets = ", ".join(f"{c}=excluded.{c}" for c in VARIANT_COLUMNS
                     if c not in untouched)
    sets += (
        ", source=CASE WHEN COALESCE(variants.source, '') LIKE '%clinvar%' "
        "THEN variants.source "
        "WHEN COALESCE(variants.source, '') = '' THEN 'clinvar' "
        "ELSE 'clinvar+' || variants.source END"
    )
    return (
        f"INSERT INTO variants ({cols}) VALUES ({marks}) "
        f"ON CONFLICT(rsid) DO UPDATE SET {sets} "
        f"WHERE excluded.review_stars > variants.review_stars "
        f"OR (excluded.review_stars = variants.review_stars "
        f"AND {_sig_rank_sql('excluded')} > {_sig_rank_sql('variants')})"
    )


def _flush(conn: sqlite3.Connection, sql: str, batch: list[tuple]) -> None:
    """Write one batch of rows inside a transaction, then clear the batch."""
    if not batch:
        return
    conn.executemany(sql, batch)
    conn.commit()
    batch.clear()


# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------

def build_clinvar(conn: sqlite3.Connection,
                  coverage: set[str],
                  limit: int = 0,
                  url: str = CLINVAR_URL) -> int:
    """Stream ClinVar variant_summary.txt.gz into the variants table.

    Returns the number of rows written. Filters to GRCh37, drops rows with no
    dbSNP mapping, and drops rows outside the array coverage set when one was
    supplied. Coordinates come from PositionVCF because Start and Stop are
    shifted relative to it.
    """
    print(f"[clinvar] streaming {url}")
    sql = _clinvar_upsert_sql()
    batch: list[tuple] = []
    columns: dict[str, int] = {}
    seen = 0
    written = 0

    for line in stream_gzip_lines(url):
        if not line:
            continue
        row = line.split("\t")
        if not columns:
            columns = resolve_columns(row, CLINVAR_COLUMNS, CLINVAR_MANDATORY)
            print(f"[clinvar] resolved {len(columns)} of "
                  f"{len(CLINVAR_COLUMNS)} columns by name")
            continue

        seen += 1
        if seen % PROGRESS_EVERY == 0:
            print(f"[clinvar] {seen:>10,} rows read, {written:>9,} written")
        if limit and seen > limit:
            break

        if _cell(row, columns, "assembly") != TARGET_ASSEMBLY:
            continue
        raw_rsid = _cell(row, columns, "rsid")
        if not raw_rsid or raw_rsid == NO_DBSNP_MAPPING:
            continue
        rsid = raw_rsid.lower()
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"
        if coverage and rsid not in coverage:
            continue

        status = _cell(row, columns, "review_status")
        sig = _cell(row, columns, "clinical_sig")
        batch.append((
            rsid,
            _cell(row, columns, "gene") or None,
            _cell(row, columns, "chromosome") or None,
            _int_or_none(_cell(row, columns, "position")),
            sig.lower() or None,
            clinvar_sig_code(sig),
            status or None,
            review_stars(status),
            _cell(row, columns, "condition") or None,
            None, None, 0, 0, None, "clinvar",
        ))
        written += 1
        if len(batch) >= BATCH_SIZE:
            _flush(conn, sql, batch)

    _flush(conn, sql, batch)
    print(f"[clinvar] done: {seen:,} rows read, {written:,} written")
    return written


# ---------------------------------------------------------------------------
# GWAS Catalog
# ---------------------------------------------------------------------------

def _split_gwas_rsids(cell: str) -> list[str]:
    """Split a SNPS cell into individual rsIDs.

    The column holds one rsID for most rows, a semicolon list for haplotypes and
    an " x " separated pair for interaction associations. Anything that is not a
    plain rs-number, such as a chr:pos locus or a merged identifier, is dropped
    because it cannot be joined against an array export.
    """
    text = cell.replace(" x ", ";").replace(",", ";").replace(" ", ";")
    out: list[str] = []
    for token in text.split(";"):
        token = token.strip().lower()
        if token.startswith("rs") and token[2:].isdigit():
            out.append(token)
    return out


def _beta_direction(ci_text: str) -> str:
    """Return "increase", "decrease" or "" for a 95% CI (TEXT) cell.

    The numeric OR or BETA column is unsigned, so a beta of 0.12 could be a rise
    or a fall. The direction is only ever stated as free text in this column, in
    phrases such as "[0.09-0.15] unit decrease". Reading the numeric column alone
    reports protective alleles as risk alleles.
    """
    lowered = ci_text.lower()
    if "decrease" in lowered or "lower" in lowered or "reduc" in lowered:
        return "decrease"
    if "increase" in lowered or "higher" in lowered:
        return "increase"
    return ""


def _risk_allele(cell: str) -> str:
    """Extract the risk allele letter from a STRONGEST SNP-RISK ALLELE cell.

    The cell reads "rs1234567-A". A trailing "?" means the catalogue did not
    determine the allele, and that is returned as blank rather than guessed,
    because a wrong risk allele inverts carrier status.
    """
    if "-" not in cell:
        return ""
    allele = cell.rsplit("-", 1)[1].strip().upper()
    return allele if allele and set(allele) <= set("ACGT") else ""


def build_gwas(conn: sqlite3.Connection,
               coverage: set[str],
               limit: int = 0,
               url: str = GWAS_URL) -> int:
    """Stream the GWAS Catalog associations TSV into the variants table.

    Only associations at PVALUE_MLOG >= 7.3, which is -log10(5e-8), are counted,
    and only rsIDs backed by two or more distinct study accessions are kept. A
    single hit is not a replicated association and would otherwise flood the
    database with noise that scoring would then have to discount.

    Returns the number of rsIDs written.
    """
    candidates = [url, *(u for u in GWAS_FALLBACK_URLS if u != url)]
    print(f"[gwas] {len(candidates)} candidate download(s)")
    studies: dict[str, set[str]] = {}
    traits: dict[str, set[str]] = {}
    papers: dict[str, set[str]] = {}
    alleles: dict[str, str] = {}
    columns: dict[str, int] = {}
    seen = 0

    for line in gwas_source_lines(candidates):
        if not line:
            continue
        row = line.split("\t")
        if not columns:
            columns = resolve_columns(row, GWAS_COLUMNS, GWAS_MANDATORY)
            print(f"[gwas] resolved {len(columns)} of "
                  f"{len(GWAS_COLUMNS)} columns by name")
            continue

        seen += 1
        if seen % PROGRESS_EVERY == 0:
            print(f"[gwas] {seen:>10,} rows read, {len(studies):>9,} rsIDs tracked")
        if limit and seen > limit:
            break

        mlog = _float_or_none(_cell(row, columns, "mlog"))
        if mlog is None or mlog < GWAS_MLOG_THRESHOLD:
            continue
        study = _cell(row, columns, "study")
        if not study:
            continue

        trait = _cell(row, columns, "trait") or _cell(row, columns, "trait_raw")
        direction = _beta_direction(_cell(row, columns, "ci_text"))
        if trait and direction:
            trait = f"{trait} ({direction})"
        allele = _risk_allele(_cell(row, columns, "risk_allele"))
        pubmed = _cell(row, columns, "pubmed")

        for rsid in _split_gwas_rsids(_cell(row, columns, "rsids")):
            if coverage and rsid not in coverage:
                continue
            studies.setdefault(rsid, set()).add(study)
            if trait:
                traits.setdefault(rsid, set()).add(trait)
            if pubmed:
                papers.setdefault(rsid, set()).add(pubmed)
            if allele and rsid not in alleles:
                alleles[rsid] = allele

    sql = (
        "INSERT INTO variants "
        "(rsid, review_stars, gwas_traits, gwas_studies, publications, "
        " risk_allele, source) "
        "VALUES (?, 0, ?, ?, ?, ?, 'gwas') "
        "ON CONFLICT(rsid) DO UPDATE SET "
        "gwas_traits=excluded.gwas_traits, "
        "gwas_studies=excluded.gwas_studies, "
        "publications=max(COALESCE(variants.publications, 0), "
        "                 COALESCE(excluded.publications, 0)), "
        "risk_allele=COALESCE(NULLIF(variants.risk_allele, ''), excluded.risk_allele), "
        "source=CASE WHEN COALESCE(variants.source, '') LIKE '%gwas%' "
        "THEN variants.source "
        "WHEN COALESCE(variants.source, '') = '' THEN 'gwas' "
        "ELSE variants.source || '+gwas' END"
    )
    batch: list[tuple] = []
    written = 0
    for rsid, accessions in studies.items():
        if len(accessions) < GWAS_MIN_STUDIES:
            continue
        batch.append((
            rsid,
            "; ".join(sorted(traits.get(rsid, ()))) or None,
            len(accessions),
            len(papers.get(rsid, ())),
            alleles.get(rsid) or None,
        ))
        written += 1
        if len(batch) >= BATCH_SIZE:
            _flush(conn, sql, batch)
    _flush(conn, sql, batch)

    print(f"[gwas] done: {seen:,} rows read, {len(studies):,} rsIDs at "
          f"mlog>={GWAS_MLOG_THRESHOLD}, {written:,} replicated in "
          f">={GWAS_MIN_STUDIES} studies")
    return written


# ---------------------------------------------------------------------------
# CPIC
# ---------------------------------------------------------------------------

def build_cpic(conn: sqlite3.Connection, limit: int = 0,
               url: str = CPIC_URL) -> int:
    """Fetch the CPIC gene and drug pair table and apply cpiclevel by gene.

    CPIC assignments are gene-level, not variant-level, so the level is written
    to every variant already carrying that gene symbol. Only genesymbol,
    drugname and cpiclevel are requested. See the CPIC_URL comment for why
    clinpgxlevel and pgxtesting are excluded at the wire level.

    Returns the number of gene and drug pairs read.
    """
    if requests is None:
        raise RuntimeError("the requests package is not installed")
    print(f"[cpic] fetching {url}")
    params = {"select": CPIC_SELECT}
    if limit:
        params["limit"] = str(limit)
    resp = requests.get(url, timeout=HTTP_TIMEOUT, headers=_HEADERS, params=params)
    resp.raise_for_status()
    pairs = resp.json()
    if not isinstance(pairs, list):
        raise ValueError("CPIC pair table did not return a JSON array")

    best: dict[str, str] = {}
    order = {level: i for i, level in enumerate(
        ("Retired", "D", "C/D", "C", "B/C", "B", "A/B", "A"))}
    for pair in pairs:
        gene = str(pair.get("genesymbol") or "").strip()
        level = str(pair.get("cpiclevel") or "").strip()
        if not gene or not level:
            continue
        if order.get(level, -1) > order.get(best.get(gene, ""), -1):
            best[gene] = level

    conn.executemany(
        "UPDATE variants SET cpic_level=? WHERE gene=?",
        [(level, gene) for gene, level in sorted(best.items())],
    )
    conn.commit()
    print(f"[cpic] done: {len(pairs):,} pairs read, {len(best):,} genes levelled")
    return len(pairs)


# ---------------------------------------------------------------------------
# Meta
# ---------------------------------------------------------------------------

def _read_meta(conn: sqlite3.Connection) -> dict[str, str]:
    """Return the meta table as a plain dict, empty when the table is absent."""
    try:
        rows = conn.execute("SELECT key, value FROM meta").fetchall()
    except sqlite3.Error:
        return {}
    return {str(k): str(v) for k, v in rows}


def write_meta(conn: sqlite3.Connection,
               counts: dict[str, int],
               coverage_size: int,
               skipped: Iterable[str] = ()) -> None:
    """Record the build date, source versions, row counts and licences.

    Versions are recorded as the UTC date the source was fetched. ClinVar and
    the GWAS Catalog both publish rolling files with no version string in the
    payload, so the fetch date is the only honest version marker available.
    """
    now = datetime.now(timezone.utc)
    stamp = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    total = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]

    entries: dict[str, str] = {
        "build_date": stamp,
        "builder": "data/build_full_reference.py",
        "schema_version": "2",
        "assembly": TARGET_ASSEMBLY,
        "variants_rows": str(total),
        "array_coverage_rsids": str(coverage_size),
        "array_filter_applied": "yes" if coverage_size else "no",
        "gwas_mlog_threshold": str(GWAS_MLOG_THRESHOLD),
        "gwas_min_studies": str(GWAS_MIN_STUDIES),
        "skipped_sources": ", ".join(sorted(skipped)) or "none",
        "clinvar_url": CLINVAR_URL,
        "gwas_url": GWAS_URL,
        "cpic_url": CPIC_URL,
    }
    for source, count in counts.items():
        entries[f"{source}_rows"] = str(count)
        entries[f"{source}_version"] = stamp[:10]
    for source, statement in SOURCE_LICENCES.items():
        entries[f"{source}_licence"] = statement
    entries["licence_position"] = (
        "This database contains only CC0 and US public domain data plus a "
        "coverage set derived from the user's own file. It is gitignored, not "
        "redistributed, and reproducible by re-running the builder."
    )

    conn.executemany(
        "INSERT INTO meta (key, value) VALUES (?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        sorted(entries.items()),
    )
    conn.commit()


# ---------------------------------------------------------------------------
# Read-side API consumed by the scanner
# ---------------------------------------------------------------------------

def open_reference_db(path: Path | str | None = None) -> sqlite3.Connection | None:
    """Return a read-only connection to the Tier 2 database, or None if absent.

    Read-only is enforced through a URI so a scan can never write to the
    reference, and returning None rather than raising means the whole Tier 2
    layer is optional: a fresh clone with no built database still scans, using
    Tier 1 alone.
    """
    target = Path(path) if path is not None else DEFAULT_DB_PATH
    if not target.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{target.as_posix()}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error:
        return None


def lookup(rsid: str, conn: sqlite3.Connection | None = None) -> dict | None:
    """Return the Tier 2 row for one rsID as a dict, or None.

    None covers all three ways this can come up empty: no database built, no row
    for that rsID, or a database that cannot be read. None of them is an error
    worth raising, because a missing Tier 2 hit just means the scanner falls back
    to Tier 1 and the live lookups.
    """
    own = conn is None
    if conn is None:
        conn = open_reference_db()
    if conn is None:
        return None
    try:
        row = conn.execute(
            "SELECT * FROM variants WHERE rsid=?", (str(rsid).strip().lower(),)
        ).fetchone()
    except sqlite3.Error:
        return None
    finally:
        if own:
            conn.close()
    return dict(row) if row is not None else None


def coverage_stats(conn: sqlite3.Connection | None = None) -> dict[str, Any]:
    """Return the meta row counts, or an empty dict when nothing is built.

    Values that are entirely numeric come back as ints so a caller can compare
    and total them without parsing.
    """
    own = conn is None
    if conn is None:
        conn = open_reference_db()
    if conn is None:
        return {}
    try:
        meta = _read_meta(conn)
    finally:
        if own:
            conn.close()
    out: dict[str, Any] = {}
    for key, value in meta.items():
        if key.endswith(("_rows", "_rsids", "_version", "_date")) or key == "assembly":
            out[key] = int(value) if value.isdigit() else value
    return out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _run_source(name: str, fn: Any) -> tuple[int, str | None]:
    """Run one source pass, converting any failure into a warning.

    A source that cannot be downloaded must not take the build down with it. The
    partial database is still useful, and the skipped source is recorded in meta
    so nobody later mistakes a partial build for a complete one.
    """
    try:
        return int(fn()), None
    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        print(f"  WARNING: {name} source failed and was skipped. {reason}")
        return 0, reason


def _print_stats(out_path: Path) -> None:
    """Print coverage_stats() for a built database as sorted JSON."""
    print("-" * 74)
    print("coverage_stats():")
    conn = open_reference_db(out_path)
    try:
        print(json.dumps(coverage_stats(conn), indent=2, sort_keys=True))
    finally:
        if conn is not None:
            conn.close()


def main(argv: list[str] | None = None) -> int:
    """Build the Tier 2 reference database. Returns a process exit code."""
    args = build_parser().parse_args(argv)
    out_path = Path(args.out) if args.out else DEFAULT_DB_PATH
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 74)
    print("DNAInsight Tier 2 full-array reference builder")
    print("=" * 74)
    print(f"output: {out_path}")

    print("array coverage set:")
    coverage = load_array_rsids(args.array_file)
    if coverage:
        print(f"  coverage set size: {len(coverage):,} distinct rsIDs")
    else:
        bang = "!" * 66
        print(f"  {bang}")
        print("  !! NO ARRAY COVERAGE FILE FOUND. Building with NO array filter.")
        print("  !! Every rsID with a database record will be included, so the")
        print("  !! result is far larger than one array genotypes and is NOT a")
        print("  !! statement of what your chip actually measures. Pass")
        print("  !! --array-file <your raw export> for a coverage-filtered build.")
        print(f"  {bang}")

    conn = sqlite3.connect(str(out_path))
    try:
        create_schema(conn)
        existing = _read_meta(conn) if args.resume else {}
        counts: dict[str, int] = {}
        skipped: list[str] = []

        def already(source: str) -> bool:
            """True when --resume should keep an existing source pass."""
            return args.resume and existing.get(f"{source}_rows", "0") not in ("", "0")

        if args.skip_clinvar:
            skipped.append("clinvar (--skip-clinvar)")
        elif already("clinvar"):
            counts["clinvar"] = int(existing["clinvar_rows"])
            print(f"[clinvar] resumed, keeping {counts['clinvar']:,} existing rows")
        else:
            counts["clinvar"], failure = _run_source(
                "clinvar", lambda: build_clinvar(conn, coverage, args.limit))
            if failure:
                skipped.append(f"clinvar ({failure})")

        if args.skip_gwas:
            skipped.append("gwas (--skip-gwas)")
        elif already("gwas"):
            counts["gwas"] = int(existing["gwas_rows"])
            print(f"[gwas] resumed, keeping {counts['gwas']:,} existing rows")
        else:
            counts["gwas"], failure = _run_source(
                "gwas", lambda: build_gwas(conn, coverage, args.limit))
            if failure:
                skipped.append(f"gwas ({failure})")

        if args.skip_cpic:
            skipped.append("cpic (--skip-cpic)")
        elif already("cpic"):
            counts["cpic"] = int(existing["cpic_rows"])
            print(f"[cpic] resumed, keeping {counts['cpic']:,} existing pairs")
        else:
            counts["cpic"], failure = _run_source("cpic", lambda: build_cpic(conn))
            if failure:
                skipped.append(f"cpic ({failure})")

        write_meta(conn, counts, len(coverage), skipped)
        total = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    finally:
        conn.close()

    print("-" * 74)
    print(f"variants table: {total:,} rows")
    for source in sorted(counts):
        print(f"  {source:<10} {counts[source]:>12,}")
    if skipped:
        print("skipped sources:")
        for note in skipped:
            print(f"  {note}")

    if args.stats:
        _print_stats(out_path)

    print(f"WROTE {out_path}")
    print("Reminder: this file is gitignored and must never be committed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
