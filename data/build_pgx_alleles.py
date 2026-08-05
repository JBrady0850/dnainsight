"""
build_pgx_alleles.py -- CPIC star-allele definition table builder.

WHY THIS FILE EXISTS
--------------------
``backend/diplotype.py`` hand-encodes star-allele definitions from the CPIC
allele definition tables, and says so in its own docstring. Seventeen of those
entries carry ``verified: False``, meaning the rsID and the plus-strand base
were recalled rather than corroborated inside this repository. That is an
honest position and a bad permanent one: an unverified allele still gets used
for calling, and a silently flipped base inverts a diplotype.

This builder closes the loop. It fetches the live CPIC allele definition tables,
writes them to ``data/pgx_alleles.json``, and produces a reconciliation report
saying which of those seventeen it could confirm, which it found in conflict,
and which CPIC no longer defines under that name at all. It does NOT edit
``backend/diplotype.py``. A build script that rewrites a reviewed clinical table
without a human reading the diff is precisely the failure mode this project
avoids everywhere else.

SOURCE AND LICENCE
------------------
    CPIC, the Clinical Pharmacogenetics Implementation Consortium
    https://cpicpgx.org/ and the PostgREST API at https://api.cpicpgx.org/v1/
    Licence: CC0-1.0. CPIC places its database in the public domain.

Because it is CC0, the derived table IS committed to this repository, unlike
every other builder output here. ``data/DATA_SOURCES.md`` section 1 records the
verification. CC0 imposes no conditions, so nothing about it can relicense the
MIT grant this project makes.

THE COLUMN-LEVEL LICENCE TRAP, WHICH IS THE WHOLE REASON THIS FILE IS CAREFUL
-----------------------------------------------------------------------------
``data/DATA_SOURCES.md`` sections 1 and 9 record the rule that matters here:

    LICENCE CONTAMINATION TRAVELS WITH THE COLUMN, NOT WITH THE FILE.

The CPIC dump is CC0, but it carries two columns that are not. ``clinpgxlevel``
and ``pgxtesting`` are PharmGKB-sourced and inherit PharmGKB's data use
agreement, which adds a term prohibiting sale of the data or of products
containing it. That extra restriction is not part of CC-BY-SA-4.0 and cannot be
added to it, so the effective licence on those two columns is a bespoke
non-commercial one. They arrive inside an otherwise CC0 download and look safe
by association. They are not.

``data/build_full_reference.py`` excludes them AT THE WIRE LEVEL by naming its
columns in the PostgREST ``select`` parameter, so they are never transferred,
let alone written. This builder does exactly the same, and adds a belt to that
brace: :func:`assert_no_forbidden_columns` is run over every parsed row from
every endpoint and raises :class:`ForbiddenColumn` if either name ever appears.
If CPIC renames a column, or a future edit drops a select list, the build fails
loudly instead of quietly copying a non-commercial column into a file this
project commits under MIT.

Only ``genesymbol``, ``drugname`` and ``cpiclevel`` are read from the gene and
drug pair table. That is the same three columns build_full_reference.py takes,
and no more.

The allele definition endpoints are a different resource and carry different
fields (``name``, ``activityvalue``, ``clinicalfunctionalstatus``, ``dbsnpid``,
``variantallele``). They are CC0 too and none of them is PharmGKB-sourced, but
the forbidden-column assertion is applied to them anyway, because the cost of
running it is nothing and the cost of not running it is a licence breach.

WHAT ``verified`` MEANS IN THE OUTPUT
--------------------------------------
Per allele, with the source and the date recorded alongside:

    verified true    CPIC's live allele definition table defines this allele
                     under this name AND gives a concrete, unambiguous
                     plus-strand base for every one of its defining positions.

    verified false   Something was missing: the allele is not in CPIC's current
                     table under that name, or CPIC's definition leaves at
                     least one position as an IUPAC ambiguity code. The reason
                     is recorded per allele.

STRAND CONVENTION, AND A TRAP
-----------------------------
CPIC states allele definitions on the POSITIVE CHROMOSOMAL STRAND, which is the
same convention ``backend/diplotype.py`` uses and the same one 23andMe and
AncestryDNA report. That is what makes the two directly comparable, and it is
why a mismatch found by :func:`reconcile` is a real disagreement rather than a
convention difference.

Two things to know before reading a conflict as an error:

  * CPIC positions are GRCh38 while ``backend/diplotype.py`` records GRCh37.
    Coordinates therefore differ. The STRAND does not, so the base should
    match, and this builder compares bases only and never positions.
  * CPIC uses IUPAC ambiguity codes (R, Y, S, W, K, M and so on) at positions
    where an allele is defined as "either base". Those are NOT concrete alleles
    and are never reported as a match. See :func:`is_concrete_base`.

CONSENT AND OFFLINE CONTRACT
----------------------------
Nothing runs at import time. ``python data/build_pgx_alleles.py --help`` works
with no network. Without ``--accept-terms`` the run is a dry run: it prints
every URL, the licence, the excluded columns and the estimated size, and
fetches nothing. Same gate as ``backend/snpedia.py``, which returns HTTP 403
until ``accept_license`` is true, and same gate as ``data/build_panel.py``, so
there is one rule across all three.

USAGE
-----
    python data/build_pgx_alleles.py                       # dry run, default
    python data/build_pgx_alleles.py --accept-terms
    python data/build_pgx_alleles.py --accept-terms --genes CYP2C19,TPMT
    python data/build_pgx_alleles.py --accept-terms --report data/pgx_recon.md
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Read-only import. This builder reads the hand-encoded table to reconcile
# against it and never writes back to it: diplotype.py is reviewed clinical
# content and a machine must not edit it unattended.
from backend.diplotype import (                                   # noqa: E402
    ALLELE_DEFINITIONS,
    unverified_entries,
)

try:  # pragma: no cover - exercised only when requests is absent
    import requests
except ImportError:  # pragma: no cover
    requests = None  # type: ignore[assignment]

__all__ = [
    "CPIC_BASE", "CPIC_HOMEPAGE", "CPIC_LICENCE", "CPIC_LICENCE_URL",
    "CPIC_SPDX", "FORBIDDEN_COLUMNS", "FORBIDDEN_REASON",
    "PAIR_URL", "PAIR_SELECT", "ALLELE_URL", "ALLELE_SELECT",
    "ALLELE_DEFINITION_URL", "ALLELE_DEFINITION_SELECT",
    "SEQUENCE_LOCATION_URL", "SEQUENCE_LOCATION_SELECT",
    "ALLELE_LOCATION_VALUE_URL", "ALLELE_LOCATION_VALUE_SELECT",
    "IUPAC_AMBIGUITY", "DEFAULT_OUT", "ESTIMATED_BYTES",
    "ForbiddenColumn",
    "build_parser", "assert_no_forbidden_columns", "check_rows",
    "is_concrete_base", "normalise_allele_name", "match_cpic_allele",
    "target_genes", "download_plan", "print_plan", "reconcile",
    "reconciliation_markdown", "build_records", "write_output", "main",
]


# ---------------------------------------------------------------------------
# Source endpoints
#
# CPIC exposes PostgREST. The select list is part of every request on purpose:
# it means the columns this project must not copy are never even transferred.
# See the module docstring and data/DATA_SOURCES.md sections 1 and 9.
# ---------------------------------------------------------------------------

CPIC_HOMEPAGE = "https://cpicpgx.org/"
CPIC_BASE = "https://api.cpicpgx.org/v1"
CPIC_LICENCE = (
    "CPIC (cpicpgx.org). CC0-1.0: CPIC places its database in the public "
    "domain and imposes no conditions on reuse, including commercial reuse. "
    "Citation is a courtesy, not a term. Verified against the CPIC site and "
    "the database release notes; see data/DATA_SOURCES.md section 1."
)
CPIC_LICENCE_URL = "https://creativecommons.org/publicdomain/zero/1.0/"
CPIC_SPDX = "CC0-1.0"

# THE TWO COLUMNS THAT ARE NOT CC0. Never requested, never parsed, never
# written. Named here once so the assertion, the select lists and the printed
# licence notice all refer to the same tuple and cannot drift apart.
FORBIDDEN_COLUMNS: tuple[str, ...] = ("clinpgxlevel", "pgxtesting")

FORBIDDEN_REASON = (
    "clinpgxlevel and pgxtesting arrive inside the otherwise CC0 CPIC dump but "
    "originate from PharmGKB and inherit PharmGKB's data use agreement, which "
    "adds a term prohibiting sale of the data or of products containing it. "
    "That term is not part of any standard licence identifier and cannot be "
    "added to one. Copying either column into a file this MIT project commits "
    "would attach a no-commercial-sale clause to every downstream user. They "
    "are excluded at the wire level in the select list and the exclusion is "
    "re-checked on every parsed row."
)

# VERIFIED HTTP 200 ON 2026-08-04 for every endpoint below.
PAIR_URL = f"{CPIC_BASE}/pair_view"
# Exactly the three columns build_full_reference.py takes, and no more.
PAIR_SELECT = "genesymbol,drugname,cpiclevel"

ALLELE_URL = f"{CPIC_BASE}/allele"
ALLELE_SELECT = ("genesymbol,name,activityvalue,clinicalfunctionalstatus,"
                 "definitionid,strength")

ALLELE_DEFINITION_URL = f"{CPIC_BASE}/allele_definition"
ALLELE_DEFINITION_SELECT = "id,genesymbol,name,structuralvariation"

SEQUENCE_LOCATION_URL = f"{CPIC_BASE}/sequence_location"
SEQUENCE_LOCATION_SELECT = ("id,genesymbol,dbsnpid,name,chromosomelocation,"
                            "position")

ALLELE_LOCATION_VALUE_URL = f"{CPIC_BASE}/allele_location_value"
ALLELE_LOCATION_VALUE_SELECT = "alleledefinitionid,locationid,variantallele"

# IUPAC ambiguity codes. CPIC uses them where an allele is defined as "either
# base at this position". They are not concrete alleles: writing one into a
# definition table would make a caller test a genotype against a letter that
# never appears in an array export, and it would silently never match.
IUPAC_AMBIGUITY: frozenset[str] = frozenset("RYSWKMBDHVN")

USER_AGENT = "DNAInsight/3.0 (+https://github.com/dnainsight) PGx allele builder"
HTTP_TIMEOUT = 120

# data/pgx_alleles.json is committed, unlike every other builder output in this
# directory. That is a licence decision, not a size decision: CC0 content may
# be bundled, and this table is small.
DEFAULT_OUT = ROOT / "data" / "pgx_alleles.json"

# Approximate response sizes measured on 2026-08-04. Printed in the dry run so
# the user knows the download is trivial rather than being asked to accept an
# unknown.
ESTIMATED_BYTES = 6 * 1024 * 1024

BUILDER_VERSION = "3.0.0"

# Genes this project actually calls. Restricting the fetch to them keeps
# data/pgx_alleles.json small enough to commit and keeps the diff readable
# when it is rebuilt. CPIC defines 184 CYP2D6 alleles alone.
DEFAULT_GENES: tuple[str, ...] = tuple(ALLELE_DEFINITIONS)


class ForbiddenColumn(RuntimeError):
    """A non-CC0 column appeared in a parsed row. Always fatal.

    Never caught anywhere in this module. The whole point is that the build
    stops rather than writing a column whose licence this project cannot honour.
    """

    def __init__(self, column: str, where: str = "") -> None:
        message = (
            f"FORBIDDEN COLUMN {column!r} appeared in CPIC data"
            + (f" while parsing {where}" if where else "")
            + ".\n\n" + FORBIDDEN_REASON
            + "\n\nThe build is stopped. Nothing was written. Fix the select "
              "list before re-running."
        )
        super().__init__(message)
        self.column = column
        self.where = where


# ---------------------------------------------------------------------------
# The licence assertion
# ---------------------------------------------------------------------------

def assert_no_forbidden_columns(row: Any, where: str = "") -> Any:
    """Raise :class:`ForbiddenColumn` if a row carries a non-CC0 column.

    Matching is case-insensitive and ignores surrounding whitespace, because
    PostgREST column names are lower-case by convention but a hand-built
    fixture or a future CPIC release need not be, and a check that a rename
    defeats is not a check.

    Returns the row unchanged so it can be used inline in a comprehension,
    which makes it hard to forget to call.
    """
    if isinstance(row, dict):
        for key in row:
            if str(key).strip().lower() in FORBIDDEN_COLUMNS:
                raise ForbiddenColumn(str(key).strip().lower(), where)
    return row


def check_rows(rows: Iterable[Any], where: str = "") -> list[Any]:
    """Run the forbidden-column assertion over every row and return them.

    Every endpoint response in this module goes through here. There is no path
    that parses a CPIC row without passing it.
    """
    return [assert_no_forbidden_columns(row, where) for row in rows]


# ---------------------------------------------------------------------------
# Allele name and base handling
# ---------------------------------------------------------------------------

def is_concrete_base(value: Any) -> bool:
    """True when a CPIC ``variantallele`` is a single unambiguous DNA base.

    Rejects IUPAC ambiguity codes and multi-base strings. An ambiguity code
    means CPIC defines the allele as "either base here", which cannot be
    compared against the single plus-strand base backend/diplotype.py records,
    and must therefore never be reported as agreement.
    """
    text = str(value or "").strip().upper()
    if len(text) != 1:
        return False
    if text in IUPAC_AMBIGUITY:
        return False
    return text in ("A", "C", "G", "T")


def normalise_allele_name(name: Any) -> str:
    """Canonical form of an allele name for matching purposes.

    Upper-cased and stripped of whitespace. Star numbers and HGVS names are
    left otherwise intact, because "*3A" and "*3a" are the same allele while
    "*3A" and "*3B" are not, and any cleverer normalisation risks merging two
    alleles that differ by one character.
    """
    return str(name or "").strip().upper().replace(" ", "")


def match_cpic_allele(cpic_names: Iterable[str], wanted: str) -> str | None:
    """Find CPIC's name for a locally recorded allele, or None.

    Three passes, most exact first:

      1. Exact match after normalisation.
      2. Parenthetical match. CPIC names several DPYD alleles by HGVS with the
         common name in brackets, for example
         "c.1129-5923C>G, c.1236G>A (HapB3)". backend/diplotype.py records that
         allele as "HapB3", and refusing to connect the two would report a real
         definition as missing.
      3. Nothing. Returning None is a result, not a failure: CPIC retires and
         renames alleles, and "CPIC no longer defines this under this name" is
         exactly what the reconciliation report should say.

    Deliberately NOT done: fuzzy or substring matching on star numbers. "*1"
    is a substring of "*10", "*11" and "*100", and a match that loose would
    confidently reconcile the wrong allele.
    """
    target = normalise_allele_name(wanted)
    if not target:
        return None
    names = list(cpic_names)

    for name in names:
        if normalise_allele_name(name) == target:
            return name

    for name in names:
        text = str(name or "")
        if "(" not in text or ")" not in text:
            continue
        inner = text[text.rfind("(") + 1:text.rfind(")")]
        if normalise_allele_name(inner) == target:
            return name

    return None


def target_genes(selection: str = "") -> list[str]:
    """Resolve --genes into a gene list, defaulting to the genes we call."""
    raw = str(selection or "").strip()
    if not raw or raw.lower() == "all":
        return list(DEFAULT_GENES)
    wanted = [g.strip().upper() for g in raw.replace(";", ",").split(",")
              if g.strip()]
    known = {g.upper(): g for g in DEFAULT_GENES}
    # An unrecognised gene is passed through rather than dropped. CPIC defines
    # far more genes than this project calls, and somebody inspecting a new one
    # is a legitimate use of this script.
    return [known.get(g, g) for g in wanted]


# ---------------------------------------------------------------------------
# The plan, and the consent gate
# ---------------------------------------------------------------------------

def download_plan(args: argparse.Namespace) -> dict:
    """Everything this invocation would fetch, with licence and exclusions."""
    genes = target_genes(args.genes)
    return {
        "source": "CPIC",
        "homepage": CPIC_HOMEPAGE,
        "licence": CPIC_LICENCE,
        "spdx": CPIC_SPDX,
        "licence_url": CPIC_LICENCE_URL,
        "genes": genes,
        "estimated_bytes": ESTIMATED_BYTES,
        "out": str(Path(args.out) if args.out else DEFAULT_OUT),
        "requests": [
            {"url": ALLELE_DEFINITION_URL, "select": ALLELE_DEFINITION_SELECT,
             "what": "allele definition identities per gene"},
            {"url": SEQUENCE_LOCATION_URL, "select": SEQUENCE_LOCATION_SELECT,
             "what": "position to rsID mapping per gene"},
            {"url": ALLELE_LOCATION_VALUE_URL,
             "select": ALLELE_LOCATION_VALUE_SELECT,
             "what": "plus-strand base per allele per position"},
            {"url": ALLELE_URL, "select": ALLELE_SELECT,
             "what": "function and activity value per allele"},
            {"url": PAIR_URL, "select": PAIR_SELECT,
             "what": "gene and drug pairs, three columns only"},
        ],
        "forbidden": list(FORBIDDEN_COLUMNS),
        "forbidden_reason": FORBIDDEN_REASON,
    }


def print_plan(plan: dict, *, dry_run: bool) -> None:
    """Print the plan, the licence and the excluded columns.

    Printed on every run. A user passing --accept-terms is agreeing to
    something specific and is entitled to read it in the same terminal before
    the first byte moves.
    """
    bar = "=" * 74
    print(bar)
    print("DNAInsight v3.0 CPIC star-allele definition builder")
    print(bar)
    print(f"source:   {plan['source']}  {plan['homepage']}")
    print(f"licence:  {plan['spdx']}  {plan['licence_url']}")
    print(f"output:   {plan['out']}")
    print(f"genes:    {', '.join(plan['genes'])}")
    print(f"size:     about {plan['estimated_bytes'] / 1048576:.1f} MiB of JSON")
    print()
    print("LICENCE, VERBATIM")
    print("-" * 74)
    for chunk in plan["licence"].split(". "):
        if chunk.strip():
            print(f"  {chunk.strip().rstrip('.')}.")
    print()
    print("WILL REQUEST")
    print("-" * 74)
    for item in plan["requests"]:
        print(f"  {item['url']}")
        print(f"    select: {item['select']}")
        print(f"    for:    {item['what']}")
    print()
    print("WILL NOT REQUEST, AND WILL FAIL LOUDLY IF IT ARRIVES ANYWAY")
    print("-" * 74)
    print(f"  columns: {', '.join(plan['forbidden'])}")
    for chunk in plan["forbidden_reason"].split(". "):
        if chunk.strip():
            print(f"  {chunk.strip().rstrip('.')}.")
    print()
    if dry_run:
        print(bar)
        print("DRY RUN. Nothing was downloaded and nothing was written.")
        print("Re-run with --accept-terms to proceed.")
        print(bar)


# ---------------------------------------------------------------------------
# Transport
# ---------------------------------------------------------------------------

def _get(url: str, select: str, where: str,
         params: dict[str, str] | None = None) -> list[dict]:
    """One PostgREST GET, column-restricted, with the licence assertion applied.

    ``select`` is not optional and is never empty. A request without one
    returns every column, which is exactly how a PharmGKB-sourced column ends
    up in a file that claims to be CC0.
    """
    if requests is None:  # pragma: no cover - only when the dep is missing
        raise RuntimeError("the requests package is not installed")
    if not select.strip():
        raise ValueError("refusing to issue a CPIC request with no select list")
    query = {"select": select}
    if params:
        query.update(params)
    resp = requests.get(url, params=query, timeout=HTTP_TIMEOUT,
                        headers={"User-Agent": USER_AGENT})
    resp.raise_for_status()
    payload = resp.json()
    if not isinstance(payload, list):
        raise ValueError(f"{url} did not return a JSON array")
    return check_rows(payload, where)


def fetch_gene(gene: str) -> dict:
    """Fetch one gene's complete allele definition table from CPIC.

    Returns ``{allele_name: {rsid: plus_strand_base_or_ambiguity_code}}`` plus
    the per-allele function and activity rows, joined on CPIC's own identifiers
    rather than on names, because names are what change between releases.
    """
    definitions = _get(ALLELE_DEFINITION_URL, ALLELE_DEFINITION_SELECT,
                       f"allele_definition/{gene}",
                       {"genesymbol": f"eq.{gene}"})
    if not definitions:
        return {"gene": gene, "definitions": {}, "alleles": {},
                "structural": set()}

    locations = _get(SEQUENCE_LOCATION_URL, SEQUENCE_LOCATION_SELECT,
                     f"sequence_location/{gene}", {"genesymbol": f"eq.{gene}"})
    by_location = {row.get("id"): row for row in locations}

    ids = ",".join(str(row["id"]) for row in definitions if row.get("id"))
    values = _get(ALLELE_LOCATION_VALUE_URL, ALLELE_LOCATION_VALUE_SELECT,
                  f"allele_location_value/{gene}",
                  {"alleledefinitionid": f"in.({ids})"}) if ids else []

    per_definition: dict[Any, dict[str, str]] = {}
    for row in values:
        location = by_location.get(row.get("locationid")) or {}
        rsid = str(location.get("dbsnpid") or "").strip().lower()
        if not rsid.startswith("rs"):
            # CPIC records some positions with no dbSNP identifier at all.
            # Those cannot be compared against an array export and are dropped
            # rather than given a synthetic name.
            continue
        base = str(row.get("variantallele") or "").strip().upper()
        per_definition.setdefault(row.get("alleledefinitionid"), {})[rsid] = base

    alleles = _get(ALLELE_URL, ALLELE_SELECT, f"allele/{gene}",
                   {"genesymbol": f"eq.{gene}"})
    by_definition_id = {row.get("definitionid"): row for row in alleles}

    out_definitions: dict[str, dict[str, str]] = {}
    out_alleles: dict[str, dict] = {}
    structural: set[str] = set()
    for row in definitions:
        name = str(row.get("name") or "").strip()
        if not name:
            continue
        out_definitions[name] = per_definition.get(row.get("id"), {})
        if row.get("structuralvariation"):
            structural.add(name)
        meta = by_definition_id.get(row.get("id")) or {}
        out_alleles[name] = {
            "function": str(meta.get("clinicalfunctionalstatus") or ""),
            "activity": meta.get("activityvalue"),
            "strength": str(meta.get("strength") or ""),
        }
    return {"gene": gene, "definitions": out_definitions,
            "alleles": out_alleles, "structural": structural}


# ---------------------------------------------------------------------------
# Reconciliation against backend/diplotype.py
# ---------------------------------------------------------------------------

def reconcile(local_entries: Sequence[dict],
              cpic: dict[str, dict]) -> list[dict]:
    """Compare locally recorded unverified alleles against CPIC's live table.

    ``local_entries`` is the output of ``backend.diplotype.unverified_entries``.
    ``cpic`` maps a gene symbol to the dict returned by :func:`fetch_gene`.

    One of five statuses per entry:

        confirmed          Every rsID recorded locally appears in CPIC's
                           definition of that allele with the same
                           unambiguous plus-strand base.
        conflict           At least one rsID appears with a DIFFERENT base.
                           This is a real disagreement, not a strand
                           convention difference: both sources state the
                           positive chromosomal strand.
        rsid_not_in_cpic   CPIC defines the allele but its definition does not
                           include one or more of the rsIDs recorded locally.
        ambiguous          CPIC gives an IUPAC ambiguity code at a position,
                           so it defines no single base to compare against.
        allele_not_in_cpic CPIC's current table has no allele of that name for
                           that gene. Alleles do get retired and renamed.
        gene_not_fetched   The gene was not part of this run.

    Nothing here modifies backend/diplotype.py. The report is for a human.
    """
    out: list[dict] = []
    for entry in local_entries:
        gene = entry["gene"]
        allele = entry["allele"]
        record: dict[str, Any] = {
            "gene": gene,
            "allele": allele,
            "local_variants": dict(entry.get("variants") or {}),
            "cpic_allele": None,
            "cpic_variants": {},
            "status": "gene_not_fetched",
            "detail": "",
            "mismatches": [],
        }

        table = cpic.get(gene)
        if table is None:
            record["detail"] = (
                f"{gene} was not fetched in this run, so nothing was compared."
            )
            out.append(record)
            continue

        matched = match_cpic_allele(table["definitions"].keys(), allele)
        if matched is None:
            record["status"] = "allele_not_in_cpic"
            record["detail"] = (
                f"CPIC's current allele definition table for {gene} has no "
                f"allele named {allele!r}. That is not evidence the local "
                f"definition is wrong; CPIC retires and renames alleles, and "
                f"this one may have been folded into another haplotype."
            )
            out.append(record)
            continue

        record["cpic_allele"] = matched
        definition = table["definitions"].get(matched, {})
        record["cpic_variants"] = {
            rsid: base for rsid, base in definition.items()
            if rsid in record["local_variants"]
        }

        mismatches: list[dict] = []
        missing: list[str] = []
        ambiguous: list[str] = []
        for rsid, local_base in record["local_variants"].items():
            cpic_base = definition.get(rsid)
            if cpic_base is None:
                missing.append(rsid)
                continue
            if not is_concrete_base(cpic_base):
                ambiguous.append(f"{rsid}={cpic_base}")
                continue
            if cpic_base.upper() != str(local_base).strip().upper():
                mismatches.append({"rsid": rsid,
                                   "local": str(local_base).strip().upper(),
                                   "cpic": cpic_base.upper()})

        record["mismatches"] = mismatches
        if mismatches:
            record["status"] = "conflict"
            record["detail"] = (
                "CPIC gives a different plus-strand base: "
                + "; ".join(f"{m['rsid']} local {m['local']} versus CPIC "
                            f"{m['cpic']}" for m in mismatches)
                + ". Both sources state the positive chromosomal strand, so "
                  "this is a genuine disagreement and needs a human to resolve "
                  "it in backend/diplotype.py."
            )
        elif missing:
            record["status"] = "rsid_not_in_cpic"
            record["detail"] = (
                f"CPIC defines {matched} but its definition does not include "
                + ", ".join(sorted(missing))
                + ". The local entry may be using a linked tag rather than a "
                  "defining variant."
            )
        elif ambiguous:
            record["status"] = "ambiguous"
            record["detail"] = (
                "CPIC records an IUPAC ambiguity code at "
                + ", ".join(sorted(ambiguous))
                + ", which defines no single base to compare against."
            )
        else:
            record["status"] = "confirmed"
            record["detail"] = (
                f"Every locally recorded rsID appears in CPIC's definition of "
                f"{matched} with the same plus-strand base."
            )
        out.append(record)
    return out


_STATUS_ORDER: tuple[str, ...] = (
    "confirmed", "conflict", "rsid_not_in_cpic", "ambiguous",
    "allele_not_in_cpic", "gene_not_fetched",
)


def reconciliation_markdown(report: Sequence[dict], *, when: str) -> str:
    """Render the reconciliation report as plain markdown for a human to read."""
    counts = {status: 0 for status in _STATUS_ORDER}
    for row in report:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    lines: list[str] = []
    lines.append("# CPIC reconciliation of unverified star alleles")
    lines.append("")
    lines.append(f"Generated {when} by `data/build_pgx_alleles.py` against the "
                 f"live CPIC API.")
    lines.append("")
    lines.append(f"`backend.diplotype.unverified_entries()` currently lists "
                 f"**{len(report)}** alleles that this repository could not "
                 f"corroborate in-tree. This run compared every one of them "
                 f"against CPIC's allele definition tables.")
    lines.append("")
    lines.append("| status | count | meaning |")
    lines.append("|---|---|---|")
    meanings = {
        "confirmed": "CPIC agrees on every recorded rsID and base",
        "conflict": "CPIC gives a different plus-strand base",
        "rsid_not_in_cpic": "CPIC defines the allele but not with that rsID",
        "ambiguous": "CPIC records an IUPAC ambiguity code, so no comparison",
        "allele_not_in_cpic": "CPIC no longer defines an allele of that name",
        "gene_not_fetched": "gene not included in this run",
    }
    for status in _STATUS_ORDER:
        if counts.get(status):
            lines.append(f"| {status} | {counts[status]} | "
                         f"{meanings[status]} |")
    lines.append("")
    lines.append("Nothing in `backend/diplotype.py` was modified. Resolving a "
                 "conflict is a human edit to a reviewed clinical table.")
    lines.append("")

    for status in _STATUS_ORDER:
        rows = [r for r in report if r["status"] == status]
        if not rows:
            continue
        lines.append(f"## {status} ({len(rows)})")
        lines.append("")
        for row in rows:
            local = ", ".join(f"{k}={v}" for k, v in
                              sorted(row["local_variants"].items()))
            lines.append(f"- **{row['gene']} {row['allele']}** "
                         f"(local: {local or 'no variants'})")
            if row.get("cpic_allele") and row["cpic_allele"] != row["allele"]:
                lines.append(f"  - matched CPIC allele: `{row['cpic_allele']}`")
            lines.append(f"  - {row['detail']}")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def build_records(cpic: dict[str, dict], *, when: str) -> list[dict]:
    """Flatten fetched CPIC tables into the committed allele record list.

    One record per gene and allele: gene, allele name, defining variants as
    rsID to plus-strand base, activity value, function, and a verified flag
    carrying its source and date.
    """
    records: list[dict] = []
    for gene in sorted(cpic):
        table = cpic[gene]
        for allele in sorted(table["definitions"]):
            raw = table["definitions"][allele]
            concrete = {rsid: base.upper() for rsid, base in sorted(raw.items())
                        if is_concrete_base(base)}
            ambiguous = sorted(rsid for rsid, base in raw.items()
                               if not is_concrete_base(base))
            meta = table["alleles"].get(allele, {})
            verified = bool(raw) and not ambiguous
            reason = ""
            if not raw:
                reason = ("CPIC defines this allele but gives it no dbSNP-mapped "
                          "defining positions, so nothing is testable from an "
                          "array.")
            elif ambiguous:
                reason = ("CPIC records IUPAC ambiguity codes at "
                          + ", ".join(ambiguous)
                          + ", so at least one defining position has no single "
                            "plus-strand base.")
            records.append({
                "gene": gene,
                "allele": allele,
                "variants": concrete,
                "ambiguous_positions": ambiguous,
                "activity": meta.get("activity"),
                "function": meta.get("function", ""),
                "evidence_strength": meta.get("strength", ""),
                "structural_variation": allele in table["structural"],
                "verified": verified,
                "verified_source": (
                    "CPIC allele definition tables via " + ALLELE_DEFINITION_URL
                ),
                "verified_date": when[:10],
                "verified_note": reason,
            })
    return records


def write_output(path: Path, records: Sequence[dict],
                 report: Sequence[dict], *, when: str,
                 genes: Sequence[str]) -> None:
    """Write data/pgx_alleles.json, including the reconciliation report.

    The report is embedded rather than written beside the table so the two can
    never drift: a definition file whose provenance note lives in a different
    file is a definition file whose provenance note goes stale.
    """
    payload = {
        "_meta": {
            "version": BUILDER_VERSION,
            "built_at": when,
            "builder": "data/build_pgx_alleles.py",
            "source": "CPIC",
            "source_url": CPIC_BASE,
            "homepage": CPIC_HOMEPAGE,
            "licence": CPIC_LICENCE,
            "spdx": CPIC_SPDX,
            "licence_url": CPIC_LICENCE_URL,
            "strand": ("GRCh37 and GRCh38 positive chromosomal strand. CPIC "
                       "states allele definitions on the positive strand, "
                       "which is what 23andMe and AncestryDNA report and what "
                       "backend/diplotype.py records."),
            "genes": list(genes),
            "alleles": len(records),
            "excluded_columns": list(FORBIDDEN_COLUMNS),
            "excluded_columns_reason": FORBIDDEN_REASON,
            "note": ("Committed to this repository because CPIC is CC0-1.0 and "
                     "CC0 imposes no conditions. See data/DATA_SOURCES.md "
                     "sections 1 and 9."),
        },
        "alleles": list(records),
        "reconciliation": {
            "against": "backend.diplotype.unverified_entries()",
            "generated_at": when,
            "note": ("This builder does not edit backend/diplotype.py. "
                     "Resolving a conflict is a human edit to a reviewed "
                     "clinical table."),
            "entries": list(report),
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)
        fh.write("\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the CPIC allele builder."""
    parser = argparse.ArgumentParser(
        prog="build_pgx_alleles.py",
        description=("Build data/pgx_alleles.json from the CPIC allele "
                     "definition tables (CC0-1.0) and reconcile it against the "
                     "unverified entries in backend/diplotype.py."),
        epilog=("Without --accept-terms this is a dry run: it prints the "
                "endpoints, the licence and the columns it refuses to "
                "transfer, and downloads nothing."),
    )
    parser.add_argument("--accept-terms", action="store_true",
                        help="accept the printed licence and allow requests. "
                             "Without this flag the run is a dry run.")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan and exit, even with --accept-terms")
    parser.add_argument("--genes", default="",
                        help="comma-separated gene symbols, default the nine "
                             "genes backend/diplotype.py calls. 'all' is "
                             "accepted and is much larger.")
    parser.add_argument("--out", default=None,
                        help=f"output path (default {DEFAULT_OUT})")
    parser.add_argument("--report", default=None, metavar="PATH",
                        help="also write the reconciliation report as markdown "
                             "to PATH. It is embedded in the JSON either way.")
    parser.add_argument("--fail-on-conflict", action="store_true",
                        help="exit non-zero when reconciliation finds a "
                             "conflicting plus-strand base, for CI use")
    return parser


def _print_report(report: Sequence[dict]) -> dict[str, int]:
    """Print the reconciliation table and return the per-status counts."""
    counts: dict[str, int] = {}
    for row in report:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    print("-" * 74)
    print(f"RECONCILIATION against backend.diplotype.unverified_entries() "
          f"({len(report)} alleles)")
    print("-" * 74)
    for status in _STATUS_ORDER:
        rows = [r for r in report if r["status"] == status]
        if not rows:
            continue
        print(f"  {status} ({len(rows)}):")
        for row in rows:
            suffix = ""
            if row["mismatches"]:
                suffix = "  " + "; ".join(
                    f"{m['rsid']} local {m['local']} vs CPIC {m['cpic']}"
                    for m in row["mismatches"])
            print(f"    {row['gene']:<8} {row['allele']:<12}{suffix}")
    print()
    print("  backend/diplotype.py was NOT modified. Every change to a reviewed")
    print("  clinical table is a human edit.")
    return counts


def main(argv: list[str] | None = None) -> int:
    """Build data/pgx_alleles.json. Returns a process exit code.

    Exit codes: 0 for a completed build or a dry run, 1 for a failed build,
    2 when --fail-on-conflict is set and a conflict was found, 3 when a
    forbidden column arrived, which is a licence failure and is reported
    separately from an ordinary error on purpose.
    """
    args = build_parser().parse_args(argv)
    plan = download_plan(args)

    # THE GATE. Dry run is the default; --accept-terms is the only thing that
    # turns it off. Same shape as backend/snpedia.py and data/build_panel.py.
    dry_run = bool(args.dry_run or not args.accept_terms)
    print_plan(plan, dry_run=dry_run)
    if dry_run:
        return 0

    genes = plan["genes"]
    when = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    out_path = Path(args.out) if args.out else DEFAULT_OUT

    cpic: dict[str, dict] = {}
    try:
        for gene in genes:
            print(f"[cpic] fetching {gene}")
            table = fetch_gene(gene)
            cpic[gene] = table
            print(f"[cpic] {gene}: {len(table['definitions'])} allele "
                  f"definition(s)")
        # The pair table is fetched last and read for three columns only. It is
        # not written into pgx_alleles.json; it is requested here so that the
        # column exclusion is exercised on every build rather than only when
        # somebody happens to run build_full_reference.py.
        pairs = _get(PAIR_URL, PAIR_SELECT, "pair_view")
        print(f"[cpic] pair table: {len(pairs):,} row(s), "
              f"{len(FORBIDDEN_COLUMNS)} forbidden column(s) checked and absent")
    except ForbiddenColumn as exc:
        print(str(exc))
        return 3
    except Exception as exc:
        print(f"BUILD FAILED: {type(exc).__name__}: {exc}")
        print("Nothing was written.")
        return 1

    records = build_records(cpic, when=when)
    report = reconcile(unverified_entries(), cpic)
    write_output(out_path, records, report, when=when, genes=genes)

    if args.report:
        report_path = Path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(reconciliation_markdown(report, when=when),
                               encoding="utf-8")
        print(f"WROTE {report_path}")

    counts = _print_report(report)
    print("-" * 74)
    print(f"WROTE {out_path}")
    print(f"  {len(records):,} allele definition(s) across {len(genes)} gene(s)")
    print(f"  verified: {sum(1 for r in records if r['verified']):,}")
    print("  This file IS committed. CPIC is CC0-1.0 and CC0 imposes no "
          "conditions.")

    if args.fail_on_conflict and counts.get("conflict"):
        print(f"\n--fail-on-conflict: {counts['conflict']} conflicting "
              f"plus-strand base(s) found.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
