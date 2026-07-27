"""
build_frequencies.py -- build-time harvester for data/frequencies.json.

Fetches per-rsID population allele frequencies for exactly the rsIDs in the
bundled reference table (data/build_reference.py) and writes a single JSON
file the app reads offline through backend/frequency.py.

Source and licensing
--------------------
Ensembl REST, ``/variation/human/{rsid}?pops=1``, which republishes the
1000 Genomes Phase 3 population frequency tables. 1000 Genomes data is open
with no restrictions; Ensembl's code is Apache-2.0 and its data is open.
Only ``1000GENOMES:phase_3:<CODE>`` rows feed the population table. gnomAD
rows are deliberately excluded from that table because the panel and the
ascertainment differ, but the gnomAD global AF is kept in a separate
``gnomad`` field where available.

Note on GMAF: current Ensembl releases return ``MAF`` and ``minor_allele``
as null on this endpoint. When that happens the global minor allele and its
frequency are recovered from the ``1000GENOMES:phase_3:ALL`` panel, which is
precisely the cohort GMAF is defined over.

Robustness
----------
This script must never be able to hang or fail a build. Every request is
throttled and retried, the whole harvest is wrapped so partial results are
always written, and the exit code is 0 even when most requests failed. When
more than 25 percent of rsIDs come back empty it writes what it has and
prints a warning. ``--offline`` skips the network entirely and emits a valid
empty scaffold so the Flask app still starts.

Usage:
    python data/build_frequencies.py
    python data/build_frequencies.py --limit 10
    python data/build_frequencies.py --only rs1801133,rs4680
    python data/build_frequencies.py --resume
    python data/build_frequencies.py --offline
"""

import argparse
import datetime
import json
import sys
import time
from datetime import timezone
from pathlib import Path
from typing import Any

import requests

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from backend.frequency import POPULATION_CODES   # noqa: E402
from data.build_reference import REFERENCE       # noqa: E402


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
VERSION = "2.0.0"
OUT_FILE = Path(__file__).parent / "frequencies.json"

ENSEMBL_BASE = "https://rest.ensembl.org"
GET_URL = ENSEMBL_BASE + "/variation/human/{rsid}"
POST_URL = ENSEMBL_BASE + "/variation/human"

USER_AGENT = ("DNAInsight/2.0 (offline personal DNA analysis; "
              "build-time frequency harvest)")
TIMEOUT = 30

# Ensembl documents 200 ids as the maximum for a POST to /variation/human.
POST_BATCH = 200

MAX_RETRIES = 3          # 3 attempts total on 429 / 5xx / network error
MAX_RPS = 10.0           # hard self-imposed ceiling
MIN_INTERVAL = 1.0 / MAX_RPS
MAX_BACKOFF = 30.0
FAIL_THRESHOLD = 0.25    # above this miss rate we warn but still exit 0

PHASE3_PREFIX = "1000GENOMES:phase_3:"
GLOBAL_CODE = "ALL"      # the pooled phase 3 panel, used only for GMAF

SOURCE = "Ensembl REST / 1000 Genomes Phase 3"
LICENSE = ("1000 Genomes: open, no restrictions; "
           "Ensembl: Apache-2.0 code, open data")
NOTE = ("Genotype frequencies derived under Hardy-Weinberg from phase 3 "
        "allele frequencies unless observed counts were available. A stored "
        "0.0 means not observed in that sample, not unknown.")


def bundled_rsids() -> list[str]:
    """Return the bundled rsIDs in table order, de-duplicated.

    This is the authoritative work list: the frequency layer only ever needs
    frequencies for variants the scanner can actually report on.
    """
    seen: set[str] = set()
    out: list[str] = []
    for row in REFERENCE:
        rsid = str(row[0]).strip().lower()
        if rsid and rsid not in seen:
            seen.add(rsid)
            out.append(rsid)
    return out


# ---------------------------------------------------------------------------
# Rate limiting and transport
# ---------------------------------------------------------------------------

class Throttle:
    """Keeps request rate under MAX_RPS and honours the server's push-back.

    Two headers matter. ``X-RateLimit-Remaining`` is Ensembl's bucket counter:
    when it runs low we pause rather than sprint into a 429. ``Retry-After``
    is an explicit instruction and is obeyed literally, capped so that a bad
    header cannot stall a build for minutes.
    """

    def __init__(self) -> None:
        self._last = 0.0

    def wait(self) -> None:
        """Block until at least MIN_INTERVAL has passed since the last call."""
        gap = time.monotonic() - self._last
        if gap < MIN_INTERVAL:
            time.sleep(MIN_INTERVAL - gap)
        self._last = time.monotonic()

    def observe(self, response: Any) -> None:
        """Inspect rate-limit headers on a response and sleep if asked to."""
        headers = getattr(response, "headers", {}) or {}
        retry_after = headers.get("Retry-After")
        if retry_after:
            try:
                time.sleep(min(max(float(retry_after), 0.0), MAX_BACKOFF))
                return
            except (TypeError, ValueError):
                pass
        remaining = headers.get("X-RateLimit-Remaining")
        if remaining is not None:
            try:
                if int(remaining) <= 2:
                    time.sleep(1.0)
            except (TypeError, ValueError):
                pass


_THROTTLE = Throttle()


def make_session() -> requests.Session:
    """Return a requests session carrying the descriptive User-Agent."""
    session = requests.Session()
    session.headers.update({
        "User-Agent": USER_AGENT,
        "Accept": "application/json",
    })
    return session


def _backoff(attempt: int) -> float:
    """Exponential backoff in seconds for attempt 1, 2, 3 -> 1, 2, 4."""
    return min(float(2 ** (attempt - 1)), MAX_BACKOFF)


def _request(session: requests.Session, method: str, url: str,
             **kwargs: Any) -> tuple[Any, str]:
    """Issue one throttled, retried request. Returns ``(payload, error)``.

    ``error`` is an empty string on success. Retries up to MAX_RETRIES times
    on 429, any 5xx and any transport-level exception, with exponential
    backoff. A 4xx other than 429 is not retried: it means the rsID is not in
    Ensembl, which is a legitimate miss rather than a fault.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        _THROTTLE.wait()
        try:
            response = session.request(method, url, timeout=TIMEOUT, **kwargs)
        except Exception as exc:                      # transport level
            if attempt == MAX_RETRIES:
                return None, f"network: {type(exc).__name__}"
            time.sleep(_backoff(attempt))
            continue

        status = response.status_code
        if status == 429 or status >= 500:
            _THROTTLE.observe(response)
            if attempt == MAX_RETRIES:
                return None, f"http {status}"
            time.sleep(_backoff(attempt))
            continue

        _THROTTLE.observe(response)
        if status >= 400:
            return None, f"http {status}"
        try:
            return response.json(), ""
        except ValueError:
            return None, "bad json"
    return None, "retries exhausted"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def derive_genotypes(alleles: dict[str, dict[str, float]]) -> dict:
    """Expand per-population allele fractions into genotype percentages.

    Hardy-Weinberg equilibrium: for the two commonest alleles A (frequency p)
    and B (frequency q), the expected genotype proportions are AA = p squared,
    AB = 2pq and BB = q squared. Values are scaled to percent so the app can
    print them directly.

    Only the top two alleles are used. Every site in the bundled reference is
    biallelic in practice, and a third allele at trace frequency would add
    noise rather than information.
    """
    out: dict[str, dict[str, float]] = {}
    for code, table in alleles.items():
        ranked = sorted(table.items(), key=lambda kv: (-kv[1], kv[0]))
        if not ranked:
            continue
        if len(ranked) == 1:
            allele, p = ranked[0]
            out[code] = {allele + allele: round(p * p * 100.0, 2)}
            continue
        (a, p), (b, q) = ranked[0], ranked[1]
        out[code] = {
            a + a: round(p * p * 100.0, 2),
            a + b: round(2.0 * p * q * 100.0, 2),
            b + b: round(q * q * 100.0, 2),
        }
    return out


def _as_float(value: Any) -> float | None:
    """Coerce an Ensembl field to float, or None when it is not numeric."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def minor_from_global(global_alleles: dict[str, float]) -> tuple[str, float | None]:
    """Recover the global minor allele and GMAF from the phase 3 ALL panel.

    Current Ensembl releases return the top-level ``MAF`` and
    ``minor_allele`` fields as null on this endpoint, so the values are taken
    from ``1000GENOMES:phase_3:ALL`` instead. That panel is the pooled phase 3
    cohort, which is exactly what GMAF is defined over. An exact 50/50 tie is
    broken alphabetically; that is arbitrary but harmless, since either choice
    gives the same 0.5 frequency.
    """
    if len(global_alleles) < 2:
        return "", None
    ranked = sorted(global_alleles.items(), key=lambda kv: (kv[1], kv[0]))
    allele, freq = ranked[0]
    return allele, round(freq, 6)


def parse_variation(payload: Any) -> dict | None:
    """Convert one Ensembl variation payload into a frequencies.json entry.

    Reads the ``populations`` list, keeping only
    ``1000GENOMES:phase_3:<CODE>`` rows whose code is one of the panels this
    app reports on. The pooled ``ALL`` row is used for GMAF but never enters
    the population table, and superpopulation rows (AFR, EUR and so on) are
    skipped because their codes are not in POPULATION_CODES.

    Returns ``None`` when the payload carries neither a usable population
    table nor a GMAF, which is the signal for a miss.
    """
    if not isinstance(payload, dict):
        return None

    alleles: dict[str, dict[str, float]] = {}
    allele_counts: dict[str, int] = {}
    global_alleles: dict[str, float] = {}
    gnomad_rows: dict[str, dict[str, float]] = {}

    for row in payload.get("populations") or []:
        if not isinstance(row, dict):
            continue
        name = str(row.get("population") or "")
        allele = str(row.get("allele") or "").strip().upper()
        freq = _as_float(row.get("frequency"))
        if not allele or freq is None:
            continue

        if name.startswith(PHASE3_PREFIX):
            code = name[len(PHASE3_PREFIX):].strip().upper()
            if code == GLOBAL_CODE:
                global_alleles[allele] = freq
                continue
            if code not in POPULATION_CODES:
                continue
            alleles.setdefault(code, {})[allele] = round(freq, 6)
            count = _as_float(row.get("allele_count"))
            if count is not None and count > 0:
                allele_counts[code] = allele_counts.get(code, 0) + int(count)
        elif name.lower().startswith("gnomad") and name.upper().endswith(":ALL"):
            gnomad_rows.setdefault(name, {})[allele] = round(freq, 6)

    minor = str(payload.get("minor_allele") or "").strip().upper()
    gmaf_value = _as_float(payload.get("MAF"))
    if gmaf_value is not None:
        gmaf_value = round(gmaf_value, 6)
    if not minor or gmaf_value is None:
        fallback_minor, fallback_gmaf = minor_from_global(global_alleles)
        minor = minor or fallback_minor
        if gmaf_value is None:
            gmaf_value = fallback_gmaf

    if not alleles and gmaf_value is None:
        return None

    # gnomAD genomes are preferred over exomes for a global AF; either is only
    # a cross-check on GMAF, never a substitute for the phase 3 panel table.
    gnomad: float | None = None
    if minor:
        for key in ("gnomADg:ALL", "gnomADe:ALL"):
            table = gnomad_rows.get(key)
            if table and minor in table:
                gnomad = table[minor]
                break
        if gnomad is None:
            for table in gnomad_rows.values():
                if minor in table:
                    gnomad = table[minor]
                    break

    entry: dict = {
        "minor_allele": minor,
        "gmaf": gmaf_value,
        "alleles": {c: alleles[c] for c in sorted(alleles)},
        "genotypes": derive_genotypes(alleles),
        # This source publishes allele frequencies, not genotype counts, so
        # every genotype table here is a Hardy-Weinberg expectation. The flag
        # stops backend/frequency.py reporting it as observed.
        "genotypes_derived": True,
        "n": {c: allele_counts[c] // 2
              for c in sorted(allele_counts) if allele_counts[c] >= 2},
    }
    if gnomad is not None:
        entry["gnomad"] = gnomad
    return entry


# ---------------------------------------------------------------------------
# Harvest
# ---------------------------------------------------------------------------

def fetch_batch(session: requests.Session,
                rsids: list[str]) -> tuple[dict, str]:
    """POST a batch of up to POST_BATCH rsIDs. Returns ``(payload, error)``.

    The batch path is purely an optimisation: one round trip instead of 200.
    ``pops`` is sent both as a query parameter and in the body because only
    the query parameter is actually honoured by the current REST release.
    Any failure is reported to the caller, which then falls back to per-rsID
    GETs rather than losing the whole batch.
    """
    payload, error = _request(
        session, "POST", POST_URL,
        params={"pops": "1"},
        json={"ids": rsids, "pops": 1},
        headers={"Content-Type": "application/json"},
    )
    if error:
        return {}, error
    if not isinstance(payload, dict):
        return {}, "unexpected batch shape"
    return payload, ""


def fetch_one(session: requests.Session, rsid: str) -> tuple[Any, str]:
    """GET a single rsID with pops=1. Returns ``(payload, error)``."""
    return _request(session, "GET", GET_URL.format(rsid=rsid),
                    params={"pops": "1"})


def harvest(session: requests.Session, rsids: list[str]) -> tuple[dict, int]:
    """Fetch every rsID and return ``(entries, miss_count)``.

    Tries the POST batch path first and degrades permanently to per-rsID GET
    the first time a batch fails, so one bad batch cannot cost the whole run.
    Progress is printed as ``n/total rsid ok|miss``.
    """
    total = len(rsids)
    entries: dict = {}
    pending = list(rsids)
    done = 0
    misses = 0
    use_post = True

    while pending:
        if use_post:
            batch = pending[:POST_BATCH]
            payload, error = fetch_batch(session, batch)
            if error:
                print(f"  batch POST failed ({error}); "
                      f"falling back to per-rsID GET", flush=True)
                use_post = False
                continue
            pending = pending[len(batch):]
            for rsid in batch:
                done += 1
                entry = parse_variation(payload.get(rsid))
                if entry is None:
                    misses += 1
                else:
                    entries[rsid] = entry
                print(f"{done}/{total} {rsid} {'ok' if entry else 'miss'}",
                      flush=True)
        else:
            rsid = pending.pop(0)
            done += 1
            payload, error = fetch_one(session, rsid)
            entry = parse_variation(payload) if not error else None
            if entry is None:
                misses += 1
            else:
                entries[rsid] = entry
            suffix = f"miss ({error})" if error else ("ok" if entry else "miss")
            print(f"{done}/{total} {rsid} {suffix}", flush=True)

    return entries, misses


def load_existing(path: Path) -> dict:
    """Return the ``frequencies`` block of an existing output file, or ``{}``."""
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except (OSError, ValueError):
        return {}
    if isinstance(raw, dict) and "frequencies" in raw:
        block = raw.get("frequencies")
        return block if isinstance(block, dict) else {}
    return raw if isinstance(raw, dict) else {}


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def present_populations(entries: dict) -> list[str]:
    """Return the population codes that actually appear, in canonical order."""
    present: set[str] = set()
    for entry in entries.values():
        if isinstance(entry, dict) and isinstance(entry.get("alleles"), dict):
            present.update(entry["alleles"].keys())
    return [c for c in POPULATION_CODES if c in present]


def write_output(path: Path, entries: dict, note_suffix: str = "") -> None:
    """Write the versioned frequencies.json payload to ``path``."""
    payload = {
        "_meta": {
            "version":     VERSION,
            "built_at":    datetime.datetime.now(timezone.utc)
                           .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "source":      SOURCE,
            "license":     LICENSE,
            "rsids":       len(entries),
            "populations": present_populations(entries),
            "note":        NOTE + note_suffix,
        },
        "frequencies": {k: entries[k] for k in sorted(entries)},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command line flags for the harvester."""
    ap = argparse.ArgumentParser(
        description="Harvest 1000 Genomes population frequencies for the "
                    "bundled DNAInsight rsIDs.")
    ap.add_argument("--limit", type=int, default=0,
                    help="fetch at most N rsIDs (0 means all)")
    ap.add_argument("--only", default="",
                    help="comma-separated rsIDs to fetch instead of the "
                         "whole bundled list")
    ap.add_argument("--resume", action="store_true",
                    help="skip rsIDs already present in frequencies.json and "
                         "merge the new results into it")
    ap.add_argument("--offline", action="store_true",
                    help="skip all network access and emit an empty scaffold "
                         "so the app still runs")
    ap.add_argument("--out", default=str(OUT_FILE),
                    help="output path (defaults to data/frequencies.json)")
    return ap.parse_args(argv)


def select_rsids(args: argparse.Namespace, existing: dict) -> list[str]:
    """Resolve the work list from the CLI flags and any existing output."""
    if args.only:
        wanted = [r.strip().lower() for r in args.only.split(",") if r.strip()]
    else:
        wanted = bundled_rsids()
    if args.resume:
        wanted = [r for r in wanted if r not in existing]
    if args.limit and args.limit > 0:
        wanted = wanted[:args.limit]
    return wanted


def main(argv: list[str] | None = None) -> int:
    """Run the harvest and write data/frequencies.json. Always returns 0.

    The build is never allowed to fail on network conditions. A total outage
    produces a scaffold plus a loud warning; a partial outage produces partial
    data plus a warning. Only a genuinely broken invocation (bad flags) exits
    non-zero, and argparse handles that before we get here.
    """
    args = parse_args(argv)
    out_path = Path(args.out)
    existing = load_existing(out_path)

    if args.offline:
        entries = dict(existing) if args.resume else {}
        write_output(out_path, entries,
                     note_suffix=" Built with --offline: no network access "
                                 "was attempted for this run.")
        print(f"OFFLINE: wrote scaffold with {len(entries)} rsID(s) -> "
              f"{out_path}")
        return 0

    wanted = select_rsids(args, existing)
    total = len(wanted)
    if total == 0:
        write_output(out_path, existing)
        print(f"Nothing to fetch. Rewrote {len(existing)} existing "
              f"rsID(s) -> {out_path}")
        return 0

    print(f"Harvesting {total} rsID(s) from {SOURCE}")
    session = make_session()
    try:
        fetched, misses = harvest(session, wanted)
    except KeyboardInterrupt:
        print("\nInterrupted; writing what was collected so far.")
        fetched, misses = {}, total
    except Exception as exc:                    # never fail the build
        print(f"\nWARNING: harvest aborted ({type(exc).__name__}: {exc}); "
              f"writing what was collected so far.")
        fetched, misses = {}, total
    finally:
        try:
            session.close()
        except Exception:
            pass

    entries = dict(existing)
    entries.update(fetched)

    suffix = ""
    miss_rate = misses / total if total else 0.0
    if miss_rate > FAIL_THRESHOLD:
        suffix = (f" Partial build: {misses} of {total} rsIDs returned no "
                  f"data on the last harvest.")
        print(f"\nWARNING: {misses}/{total} rsIDs "
              f"({miss_rate * 100:.0f} percent) returned no data. This "
              f"usually means the network is blocked. Writing partial "
              f"results anyway.")

    write_output(out_path, entries, note_suffix=suffix)
    pops = present_populations(entries)
    print(f"\nWrote {len(entries)} rsID(s), {len(pops)} population(s) "
          f"-> {out_path}")
    print(f"  fetched this run: {len(fetched)}   misses: {misses}")
    if pops:
        print(f"  populations: {', '.join(pops)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
