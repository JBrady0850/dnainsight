"""audit_y_dbsnp.py -- check Y_BACKBONE rows against NCBI dbSNP.

WHAT THIS CAN AND CANNOT SETTLE
-------------------------------
dbSNP settles variant CLASS (substitution against indel), chromosome, GRCh38
and GRCh37 positions, and the reference and alternate allele set.

dbSNP CANNOT settle ancestral against derived. It reports REFERENCE over
ALTERNATE; Y_BACKBONE records ANCESTRAL over DERIVED. On the Y these routinely
disagree, because the reference Y descends from a lineage carrying the derived
allele at many backbone nodes. This script therefore reports which state the
reference carries and decides nothing about the assignment. That is the whole
point of the `ref_carries` field added in v3.1.1, and the first run of this
audit measured the risk it guards against: 10 of 17 determinable nodes have a
reference carrying the DERIVED allele.

WHY parse_spdi EXISTS AS ITS OWN FUNCTION
-----------------------------------------
The `spdi` field is a COMMA-SEPARATED LIST when a site is multi-allelic. A
first draft of this audit split on ":" and expected four parts, which silently
blanked the allele pair for every multi-allelic record and reported four of
them as conflicts against named markers. The draft was discarded before it was
acted on, and `tests/test_y_dbsnp_audit.py` pins the parsing so the same wrong
answer cannot be produced again.

NOTHING IS WRITTEN BACK
-----------------------
This script reads Y_BACKBONE and never modifies it. A dbSNP value may not be
copied into an `ancestral` or `derived` field on the strength of this audit
alone, because of the reference/ancestral inversion described above.

Usage:
    python tools/audit_y_dbsnp.py            print the table
    python tools/audit_y_dbsnp.py --json     also write y_backbone_audit.json
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(ROOT))

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

CONSISTENT       = "consistent"
MULTI_ALLELIC    = "consistent (multi-allelic)"
STRAND           = "strand"
CONFLICT         = "CONFLICT"
CLASS_ERROR      = "CLASS"
NOT_FOUND        = "not in dbSNP"


def parse_spdi(field) -> tuple[str, set]:
    """Return (reference allele, set of alternate alleles) across every SPDI.

    The field is a comma-separated list when the site is multi-allelic. Parsing
    only the first entry blanks the pair and turns a normal tri-allelic record
    into a false accusation of conflict.
    """
    ref, alts = "", set()
    for spdi in str(field or "").split(","):
        parts = spdi.strip().split(":")
        if len(parts) != 4:
            continue
        if not ref:
            ref = parts[2]
        if parts[3]:
            alts.add(parts[3])
    return ref, alts


def reference_state(ref: str, ancestral: str, derived: str, snp_class: str) -> str:
    """Which state the assembly reference carries, or "undetermined".

    Reported, never acted on. A builder that maps this onto `ancestral` inverts
    the tree wherever the reference carries the derived allele, which the first
    run of this audit measured at 10 of 17 determinable nodes.
    """
    if not ref or snp_class != "snv":
        return "undetermined"
    if ref == ancestral:
        return "ancestral"
    if ref == derived:
        return "derived"
    if COMPLEMENT.get(ref) == ancestral:
        return "ancestral (opposite strand)"
    if COMPLEMENT.get(ref) == derived:
        return "derived (opposite strand)"
    return "undetermined"


def classify(entry: dict, record: dict) -> dict:
    """Compare one Y_BACKBONE row against one dbSNP esummary record.

    Verdicts, in the order they are tested:

      not in dbSNP  the rsID returned nothing at all
      CLASS         dbSNP does not call it a substitution, so the row is wrong
                    in KIND rather than in value. An indel recorded as a base
                    substitution cannot be genotyped by a rule expecting bases.
      consistent    the recorded pair is exactly the observed pair
      multi-allelic the recorded pair is a subset of a larger observed set,
                    which is normal and is not a conflict
      strand        the recorded pair matches the complement of the observed
                    set, so the two sources differ in strand, not in fact
      CONFLICT      none of the above. Residue, never the default.
    """
    ref, alts = parse_spdi(record.get("spdi"))
    ancestral = str(entry.get("ancestral") or "").upper()
    derived = str(entry.get("derived") or "").upper()
    snp_class = record.get("snp_class")

    observed = ({ref} | alts) if ref else set()
    recorded = {a for a in (ancestral, derived) if a}
    flipped = {COMPLEMENT.get(a, a) for a in observed}

    if not record:
        verdict = NOT_FOUND
    elif snp_class != "snv":
        verdict = CLASS_ERROR
    elif recorded and recorded <= observed:
        verdict = CONSISTENT if len(observed) == 2 else MULTI_ALLELIC
    elif recorded and recorded <= flipped:
        verdict = STRAND
    else:
        verdict = CONFLICT

    return {
        "marker": entry.get("marker"),
        "rsid": entry.get("rsid"),
        "chr": record.get("chr"),
        "class": snp_class,
        "spdi": record.get("spdi"),
        "grch38": record.get("chrpos"),
        "grch37": record.get("chrpos_prev_assm"),
        "dbsnp_ref": ref,
        "dbsnp_alts": sorted(alts),
        "table_ancestral": ancestral,
        "table_derived": derived,
        "verdict": verdict,
        "ref_carries": reference_state(ref, ancestral, derived, snp_class),
    }


def fetch(rsids, timeout: int = 90) -> dict:
    """One batched E-utilities request for every rsID. Needs the network."""
    ids = ",".join(str(r)[2:] for r in rsids)
    url = f"{EUTILS}?db=snp&id={ids}&retmode=json"
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return json.load(fh).get("result", {})


def audit(backbone: dict, records: dict) -> list[dict]:
    """Classify every backbone row that carries an rsID. No network, no writes."""
    out = []
    for node, entry in sorted(backbone.items()):
        if not entry.get("rsid"):
            continue
        row = classify(entry, records.get(str(entry["rsid"])[2:]) or {})
        row["node"] = node
        out.append(row)
    return out


def unresolved(backbone: dict) -> list[str]:
    """Markers with no rsID, which dbSNP cannot reach.

    dbSNP does not index marker names: esearch for M91, M175 and M267 each
    returns zero hits. These rows need a name-to-rsID source before they can be
    audited at all, and are reported as unaudited rather than assumed correct.
    """
    return [str(e.get("marker")) for n, e in sorted(backbone.items())
            if n != "root" and not e.get("rsid")]


def main() -> int:
    from backend import haplogroups as H

    backbone = H.Y_BACKBONE
    rsids = [e["rsid"] for e in backbone.values() if e.get("rsid")]
    if not rsids:
        print("no rsIDs recorded in Y_BACKBONE; nothing to audit")
        return 0

    try:
        records = fetch(rsids)
    except Exception as exc:
        print(f"dbSNP unreachable: {exc}")
        return 1

    rows = audit(backbone, records)
    missing = unresolved(backbone)

    header = (f"{'node':<8}{'marker':<7}{'rsid':<12}{'class':<8}{'GRCh38':<14}"
              f"{'GRCh37':<14}{'ref/alt':<14}{'anc>der':<10}{'verdict':<28}ref_carries")
    print(header)
    print("-" * len(header))
    for r in rows:
        pair = f"{r['dbsnp_ref']}/{'/'.join(r['dbsnp_alts'])}"
        print(f"{r['node']:<8}{str(r['marker']):<7}{r['rsid']:<12}{str(r['class']):<8}"
              f"{str(r['grch38']):<14}{str(r['grch37']):<14}{pair:<14}"
              f"{r['table_ancestral']}>{r['table_derived']:<8}"
              f"{r['verdict']:<28}{r['ref_carries']}")

    print()
    for verdict in (CONSISTENT, MULTI_ALLELIC, STRAND, CLASS_ERROR, CONFLICT, NOT_FOUND):
        hits = [r["marker"] for r in rows if r["verdict"] == verdict]
        if hits:
            print(f"{verdict:<28}{len(hits):>3}   {', '.join(hits)}")

    print()
    print("what the reference carries:")
    for state in sorted({r["ref_carries"] for r in rows}):
        hits = [r["marker"] for r in rows if r["ref_carries"] == state]
        print(f"  {state:<30}{len(hits):>3}   {', '.join(hits)}")

    print()
    print(f"audited {len(rows)} of {len(rows) + len(missing)} markers.")
    print(f"{len(missing)} carry no rsID and cannot be reached by dbSNP:")
    print("  " + ", ".join(missing))
    print()
    print("dbSNP cannot settle ancestral against derived. No value above may be")
    print("copied into an ancestral or derived field on the strength of this run.")

    if "--json" in sys.argv:
        target = ROOT / "y_backbone_audit.json"
        target.write_text(json.dumps({"audited": rows, "unresolved": missing},
                                     indent=1), encoding="utf-8")
        print(f"\nwrote {target}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
