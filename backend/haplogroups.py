"""
haplogroups.py -- uniparental lineage calling with an explicit resolution ceiling.

WHAT THIS IS FOR
----------------
A haplogroup call from a consumer array is a shallow call, and every commercial
product that sells one is quiet about how shallow. A 23andMe v5 chip carries
roughly 2,000 Y-SNPs out of the tens of thousands on a modern Y tree, and about
2,500 mitochondrial positions out of 16,569 bases, which is around 15 percent of
the mitochondrial genome. FamilyTreeDNA sells depth, so it has a commercial
reason to be clear. 23andMe reports a shallow array call against an outdated
tree and does not say so.

MyHeritage added a Y-DNA haplogroup tier in 2026, sourced from the kit the user
already sent, and caps it at what FamilyTreeDNA's comparison page calls
Intermediate: Y-Adam down to the Metal Age, two notable connections and two
ancient ones, with a paid Big Y-700 upgrade for anything deeper. That tier is
capped by product decision rather than by the assay, since MyHeritage moved to
low-pass whole genome sequencing in late 2025 and low-pass WGS reads the Y.
Which is the point worth keeping in view here: the depth a vendor reports is a
commercial choice, the depth the data supports is a measurement, and only one of
those two numbers is in this module. Ours is the measurement.

Checked 2026-08-09 against the FamilyTreeDNA MyHeritage upgrade page. An earlier
version of this docstring said MyHeritage did not offer the feature at all,
which was true when written and is not now. Competitive claims in comments go
stale silently, so this one carries its date.

DNAInsight's differentiator is telling the user exactly where their data runs
out. So the resolution ceiling is a first-class returned field on every call,
computed from the person's own data, never a footnote and never hardcoded.
See :func:`resolution_ceiling`.

THE TRI-STATE, WHICH IS PROJECT INVARIANT 3
-------------------------------------------
``backend/genosets.py`` splits results into matched, unmatched and incomplete
precisely so the UI can say "not testable on your array" rather than "absent".
The same distinction is load-bearing here and is the single most important
behaviour in this module:

  derived        the marker was read and carries the derived state
  ancestral      the marker was read and carries the ancestral state
  not testable   the marker was NOT read

A marker that is not on the array must NEVER be reported as ancestral. Doing so
turns "we could not look" into "we looked and found nothing", which is how a
person gets told they are not R-U106 when nobody ever tested U106.

VERIFICATION HONESTY
--------------------
The bundled backbone below is a working offline tree, not a curated data
product. Every entry carries a ``verified`` boolean.

  * NOT ONE Y entry is verified. Every Y rsID and every ancestral/derived
    allele pair in ``Y_BACKBONE`` is a CANDIDATE recorded from literature
    recall and must be confirmed by the data builder against ISOGG or YFull
    before it is flipped to True. They are included because an unverified
    candidate that is flagged is more useful than an empty tree, and because
    ``verified_only=True`` gives any caller a way to refuse them outright.
  * mtDNA positions are rCRS coordinates, which are unambiguous, so the
    textbook macro-haplogroup positions there ARE marked verified. The
    L-lineage nodes are not; see ``MT_BACKBONE`` notes.

:func:`unverified_markers` returns the full audit list. The confidence field on
every call is capped at "provisional" while any unverified marker is on the
resolved path.

TREE VERSIONING
---------------
A haplogroup call is meaningless without the tree it was made against. R1b1a2
under one tree revision is R1b1a1b under another and the person did not change.
Every payload from this module carries ``tree_name`` and ``tree_version``, and
external results carry the external tool's own tree string when it reports one.

EXTERNAL TOOLS
--------------
Yleaf (GPL-3.0), HaploGrep 3 (MIT) and Clade Finder (MIT) refine these calls.
None of them is imported, linked or vendored. They are invoked through
``external.run``, and the subprocess boundary IS the licence boundary: a
separate process executing a program the user installed into their own home
directory does not relicense this MIT tree. Every adapter opens with
``external.guard`` and returns the standard degraded payload rather than
raising, so a user with no tools installed still gets the bundled backbone call
and an honest note about what a deeper test would add.

When Yleaf and Clade Finder both run and DISAGREE, the disagreement is surfaced
as a conflict and no winner is picked. That is deliberately the same choice
``merge.py`` makes for pooled file conflicts: disagreement between two callers
is information about reliability, and a silent merge destroys it.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

try:  # the backend package is the normal import path
    from backend import external as _external
except ImportError:  # pragma: no cover - direct-module import fallback
    import external as _external  # type: ignore

external = _external


__all__ = [
    "TREE_NAME", "TREE_VERSION", "TREE_STAMP",
    "Y_BACKBONE", "MT_BACKBONE", "MT_GENOME_BP", "TYPICAL_ARRAY",
    "DERIVED", "ANCESTRAL", "NOT_ON_ARRAY", "NO_CALL", "UNUSABLE",
    "DISCORDANT", "NOT_TESTABLE_STATES",
    "tree_stamp", "unverified_markers", "path_to", "marker_keys",
    "marker_state", "call_y_backbone", "call_mt_backbone",
    "resolution_ceiling", "mt_positions_from_merged",
    "write_hsd", "write_y_array_input",
    "parse_yleaf_output", "parse_haplogrep_output", "parse_cladefinder_output",
    "compare_y_calls", "call_y", "call_mt", "second_opinion_y", "analyse",
]


# ---------------------------------------------------------------------------
# Tree identity
# ---------------------------------------------------------------------------

TREE_NAME = "DNAInsight backbone"
# 0.2, v3.4.0: 17 Y rsIDs resolved against Karafet et al. 2008 Supplementary
# Table 1, M20's derived allele corrected, three allele pairs transposed back to
# the source orientation, and two more markers reclassified as length
# polymorphisms. A call made against 0.1 is not comparable to one made against
# 0.2, which is the entire reason this string rides on every payload.
TREE_VERSION = "0.2"

# Bumped whenever any marker, allele or edge below changes. A stored call that
# carries an older stamp must be recomputed, never re-displayed as current.
TREE_STAMP = {"tree_name": TREE_NAME, "tree_version": TREE_VERSION}


def tree_stamp() -> dict:
    """The tree identity stamped onto every payload this module returns."""
    return dict(TREE_STAMP)


# ---------------------------------------------------------------------------
# Marker states
# ---------------------------------------------------------------------------

DERIVED = "derived"
ANCESTRAL = "ancestral"
NOT_ON_ARRAY = "not_on_array"
NO_CALL = "no_call"
UNUSABLE = "unusable"
DISCORDANT = "discordant"

# Everything that is NOT a statement about the person's lineage. Grouped here
# so no caller can accidentally treat "we never read it" as "ancestral".
NOT_TESTABLE_STATES = (NOT_ON_ARRAY, NO_CALL, UNUSABLE)

# Tokens an array writes when a probe failed. Kept identical to merge.py and
# genosets.py so the three modules agree on what "no data" looks like.
_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "D", "I", "?", "."}


# ---------------------------------------------------------------------------
# Array context, used only for the "what would a deeper test add" comparison.
#
# These are approximate published figures for a 23andMe v5 chip. They are never
# presented as the user's own numbers: the user's counts are computed from the
# user's own genotype map. They exist so the ceiling sentence can say how much
# of the mitochondrial genome an array reads without hardcoding the percentage.
# ---------------------------------------------------------------------------

MT_GENOME_BP = 16569          # rCRS length, NC_012920.1
TYPICAL_ARRAY = {
    "name": "23andMe v5 style consumer array",
    "y_snps": 2000,           # approximate
    "mt_positions": 2500,     # approximate
}


# ---------------------------------------------------------------------------
# Bundled Y backbone.
#
# HONESTY RULE FOR THIS TABLE, DO NOT WEAKEN IT
# ---------------------------------------------
# `verified` is False on every single entry. The rsIDs and the ancestral and
# derived alleles below are CANDIDATES recorded from literature recall. Not one
# of them has been checked against dbSNP, ISOGG or YFull by a machine, and a
# wrong rsID here produces a confidently wrong haplogroup, which is worse than
# no haplogroup. The data builder must confirm each row and only then set
# verified=True. Until it does:
#
#   * the payload confidence is capped at "provisional",
#   * `verified_only=True` refuses to use any of them,
#   * `unverified_markers()` lists every one for audit.
#
# Where no rsID is recorded at all the marker NAME is the lookup key, which
# means the marker reads as not testable on an rsID-keyed array. That is the
# honest outcome, not a bug.
#
# Intermediate ISOGG nodes are deliberately collapsed. This is a backbone for
# telling a user roughly where they sit and, much more importantly, where their
# array stops. It is not a substitute for a full tree.
# ---------------------------------------------------------------------------

def _y(node: str, parent: str | None, label: str, marker: str,
       rsid: str | None, ancestral: str | None, derived: str | None,
       *, verified: bool = False, note: str = "",
       assembly: str | None = None, ref_carries: str | None = None,
       dbsnp_checked: bool = False) -> dict:
    """Build one Y backbone entry. ``verified`` defaults False on purpose.

    THE THREE FIELDS THAT LOOK LIKE BOOKKEEPING AND ARE NOT
    -------------------------------------------------------
    ``assembly``, ``ref_carries`` and ``dbsnp_checked`` exist because dbSNP is
    the obvious machine-readable source for confirming this table and it cannot
    confirm the thing that matters most.

    dbSNP reports REFERENCE over ALTERNATE. This table records ANCESTRAL over
    DERIVED. Those are different questions, and on the Y chromosome they
    routinely disagree, because the GRCh38 reference Y comes from a single
    non-African lineage that carries the derived allele at many backbone nodes.
    Worked example, and the reason this exists: dbSNP gives rs2032595 (M168, the
    CT node) as chrY:12702062 T>C forward, while this table gives M168 as
    ancestral C, derived T. Both are right. The reference simply carries the
    derived state there.

    A builder that mapped ``ref`` onto ``ancestral`` would invert roughly half
    this tree and every test in the suite would still pass, because the data
    would be internally consistent and externally backwards. So an entry may not
    claim ``verified`` until it also records which assembly it was checked
    against and which state that assembly's reference carries.
    ``tests/test_haplogroup_nomenclature.py`` enforces that and will fail the
    build rather than trust anyone to remember it.

    ``ref_carries`` is "ancestral", "derived", or None for not yet established.
    None is the honest default and is never treated as "ancestral".
    """
    return {
        "system": "Y",
        "node": node,
        "parent": parent,
        "label": label,
        "marker": marker,
        "rsid": rsid,
        "ancestral": (ancestral or "").upper() or None,
        "derived": (derived or "").upper() or None,
        "verified": bool(verified),
        "assembly": assembly,
        "ref_carries": ref_carries,
        "dbsnp_checked": bool(dbsnp_checked),
        # A marker is a base substitution until something establishes it is not.
        # "snv" is the honest default because every row here was recorded as a
        # substitution; where that turned out to be wrong in KIND, the audit
        # table below overwrites this and clears the single-base fields.
        # ancestral_seq and derived_seq carry whole sequences for indels, and
        # stay None for substitutions so there is one place to read each fact.
        "variant_type": "snv",
        "ancestral_seq": None,
        "derived_seq": None,
        "note": note,
    }


_UNVERIFIED_Y = (
    "rsID and allele assignment unverified; the data builder must confirm this "
    "row against ISOGG or YFull before it is trusted."
)

Y_BACKBONE: dict[str, dict] = {
    "root": {
        "system": "Y", "node": "root", "parent": None, "label": "Y-MRCA",
        "marker": None, "rsid": None, "ancestral": None, "derived": None,
        "verified": True,
        "note": "Root carries no marker. Every male is here by definition.",
    },
    # -- deepest splits ------------------------------------------------------
    "A":       _y("A", "root", "A", "M31", None, "G", "A", note=_UNVERIFIED_Y),
    "BT":      _y("BT", "root", "BT", "M91", None, "A", "T", note=_UNVERIFIED_Y),
    "B":       _y("B", "BT", "B", "M60", None, "C", "A", note=_UNVERIFIED_Y),
    "CT":      _y("CT", "BT", "CT", "M168", "rs2032595", "C", "T",
                  note="rs2032595 was supplied as an example pairing in the "
                       "Wave 3 brief and is still unconfirmed here. " + _UNVERIFIED_Y),
    # -- DE and E ------------------------------------------------------------
    "DE":      _y("DE", "CT", "DE", "M145", None, "A", "G", note=_UNVERIFIED_Y),
    "D":       _y("D", "DE", "D", "M174", None, "T", "C", note=_UNVERIFIED_Y),
    "E":       _y("E", "DE", "E", "M96", "rs9306841", "G", "C", note=_UNVERIFIED_Y),
    "E-M2":    _y("E-M2", "E", "E1b1a", "M2", None, "A", "G", note=_UNVERIFIED_Y),
    "E-M35":   _y("E-M35", "E", "E1b1b", "M35", None, "G", "C", note=_UNVERIFIED_Y),
    # -- C and F -------------------------------------------------------------
    "CF":      _y("CF", "CT", "CF", "P143", None, "C", "T", note=_UNVERIFIED_Y),
    "C":       _y("C", "CF", "C", "M130", None, "C", "T", note=_UNVERIFIED_Y),
    "C-M217":  _y("C-M217", "C", "C2", "M217", None, "A", "C", note=_UNVERIFIED_Y),
    "F":       _y("F", "CF", "F", "M89", "rs2032652", "C", "T", note=_UNVERIFIED_Y),
    # -- GHIJK spine ---------------------------------------------------------
    "GHIJK":   _y("GHIJK", "F", "GHIJK", "F1329", None, "A", "G", note=_UNVERIFIED_Y),
    "G":       _y("G", "GHIJK", "G", "M201", "rs2032636", "G", "T", note=_UNVERIFIED_Y),
    "G-P15":   _y("G-P15", "G", "G2a", "P15", None, "C", "T", note=_UNVERIFIED_Y),
    "HIJK":    _y("HIJK", "GHIJK", "HIJK", "F929", None, "A", "G", note=_UNVERIFIED_Y),
    "H":       _y("H", "HIJK", "H", "M69", None, "T", "C", note=_UNVERIFIED_Y),
    "IJK":     _y("IJK", "HIJK", "IJK", "L15", None, "G", "A", note=_UNVERIFIED_Y),
    # -- IJ, I, J ------------------------------------------------------------
    "IJ":      _y("IJ", "IJK", "IJ", "M429", None, "A", "C", note=_UNVERIFIED_Y),
    "I":       _y("I", "IJ", "I", "M170", "rs2032597", "A", "C", note=_UNVERIFIED_Y),
    "I1":      _y("I1", "I", "I1", "M253", "rs17307677", "C", "T", note=_UNVERIFIED_Y),
    "I2":      _y("I2", "I", "I2", "M438", None, "A", "G", note=_UNVERIFIED_Y),
    "I2a":     _y("I2a", "I2", "I2a", "P37.2", None, "T", "C", note=_UNVERIFIED_Y),
    "J":       _y("J", "IJ", "J", "M304", "rs13447352", "A", "C", note=_UNVERIFIED_Y),
    "J1":      _y("J1", "J", "J1", "M267", None, "C", "A", note=_UNVERIFIED_Y),
    "J2":      _y("J2", "J", "J2", "M172", None, "T", "G", note=_UNVERIFIED_Y),
    "J2a":     _y("J2a", "J2", "J2a", "M410", None, "A", "G", note=_UNVERIFIED_Y),
    # -- K and below ---------------------------------------------------------
    "K":       _y("K", "IJK", "K", "M9", "rs3900", "C", "G", note=_UNVERIFIED_Y),
    "LT":      _y("LT", "K", "LT", "L298", None, "C", "T", note=_UNVERIFIED_Y),
    "L":       _y("L", "LT", "L", "M20", "rs3911", "A", "C", note=_UNVERIFIED_Y),
    "T":       _y("T", "LT", "T", "M184", None, "G", "A", note=_UNVERIFIED_Y),
    "NO":      _y("NO", "K", "NO", "M214", None, "T", "C", note=_UNVERIFIED_Y),
    "N":       _y("N", "NO", "N", "M231", None, "G", "A", note=_UNVERIFIED_Y),
    "N-M178":  _y("N-M178", "N", "N1a1", "M178", None, "C", "T", note=_UNVERIFIED_Y),
    "O":       _y("O", "NO", "O", "M175", None, "A", "G", note=_UNVERIFIED_Y),
    "O-M122":  _y("O-M122", "O", "O2", "M122", None, "T", "C", note=_UNVERIFIED_Y),
    "K2b":     _y("K2b", "K", "K2b", "P331", None, "A", "G", note=_UNVERIFIED_Y),
    # -- P, Q, R -------------------------------------------------------------
    "P":       _y("P", "K2b", "P", "M45", "rs2032631", "G", "A", note=_UNVERIFIED_Y),
    "Q":       _y("Q", "P", "Q", "M242", "rs8179021", "C", "T", note=_UNVERIFIED_Y),
    "Q-M3":    _y("Q-M3", "Q", "Q1b1a1a", "M3", None, "C", "T", note=_UNVERIFIED_Y),
    "R":       _y("R", "P", "R", "M207", "rs2032658", "A", "G", note=_UNVERIFIED_Y),
    "R1":      _y("R1", "R", "R1", "M173", "rs2032624", "A", "C", note=_UNVERIFIED_Y),
    "R1a":     _y("R1a", "R1", "R1a", "M420", None, "T", "A", note=_UNVERIFIED_Y),
    "R1a1a":   _y("R1a1a", "R1a", "R1a1a", "M17", "rs3908", "G", "A", note=_UNVERIFIED_Y),
    "R1b":     _y("R1b", "R1", "R1b", "M343", "rs9786184", "C", "A", note=_UNVERIFIED_Y),
    "R-M269":  _y("R-M269", "R1b", "R1b1a1b", "M269", "rs9786153", "T", "C",
                  note=_UNVERIFIED_Y),
    "R-U106":  _y("R-U106", "R-M269", "R1b1a1b1a1a1", "U106", "rs16981293", "C", "T",
                  note=_UNVERIFIED_Y),
    "R-P312":  _y("R-P312", "R-M269", "R1b1a1b1a1b", "P312", "rs34276300", "A", "G",
                  note=_UNVERIFIED_Y),
}


# ---------------------------------------------------------------------------
# dbSNP audit, run 2026-08-10, recorded here rather than edited into 18 rows.
#
# Applied as a table so the provenance of every changed field is visible in one
# place and can be re-derived by `python tools/audit_y_dbsnp.py`. See
# docs/Y_BACKBONE_AUDIT.md for the full run and its limits.
#
# WHAT THIS DOES NOT DO: it does not mark anything verified. dbSNP settles
# variant class, position and the reference/alternate pair. It cannot settle
# ancestral against derived, which is the assignment `verified` is about. The
# audit measured exactly why that matters: the GRCh38 reference Y carries the
# DERIVED allele at 10 of the 17 nodes where the state is determinable, so a
# builder mapping `ref` onto `ancestral` would invert 59 percent of this table.
# ---------------------------------------------------------------------------

_DBSNP_ASSEMBLY = "GRCh38"

# marker -> which state the GRCh38 reference carries
_DBSNP_REF_CARRIES: dict[str, str] = {
    "M168": "derived",   "M96":  "derived",   "M89":  "derived",
    "M201": "ancestral", "M170": "ancestral", "M253": "derived",
    "M304": "ancestral", "M9":   "derived",   "M20":  "ancestral",
    "M45":  "derived",   "M242": "ancestral", "M207": "derived",
    "M173": "derived",   "M343": "derived",   "M269": "derived",
    "U106": "ancestral", "P312": "ancestral",
}

# Markers dbSNP reports as multi-allelic. The recorded pair is a subset of the
# observed set in every case, which is normal and is not a conflict. Recorded
# because reading only the first SPDI of the comma-separated list makes all
# four look like conflicts, which is a false accusation about real data.
_DBSNP_MULTI_ALLELIC = ("M45", "M343", "M269", "P312")

_M20_CONFLICT = (
    "CONFLICT, unresolved. dbSNP rs3911 records A/G at chrY:19571568 (GRCh38), "
    "single-allelic. This table records ancestral A, derived C. There is no C "
    "allele in dbSNP and complementing does not reconcile the two, so one of "
    "the rsID assignment and the allele pair is wrong. Not corrected here "
    "because the audit cannot say which. See docs/Y_BACKBONE_AUDIT.md."
)

# Markers that are NOT base substitutions. An indel recorded as a substitution
# is wrong in KIND, not in value, and no allele edit fixes it.
#
# marker -> (variant_type, ancestral_seq, derived_seq, note)
_INDEL_CORRECTIONS: dict[str, tuple] = {
    "M17": (
        "del", "GGGG", "GGG",
        "Single-base deletion in a four-base G homopolymer, NOT a G>A "
        "substitution. dbSNP rs3908 gives snp_class delins, SPDI "
        "NC_000024.10:19571278:GGGG:GGG, HGVS NC_000024.10:g.19571282del, "
        "SEQ=[G/-]. This row previously recorded ancestral G, derived A; there "
        "is no A allele at that site. M17 defines R1a1a, so that node cannot "
        "be called from array base calls at all.",
    ),
    "M91": (
        "del", "TTTTTTTTT", "TTTTTTTT",
        "Length polymorphism in a T homopolymer, 9T to 8T, NOT an A>T "
        "substitution. Stated as such in the Karafet et al. 2008 text. This row "
        "previously recorded ancestral A, derived T. No rsID has been "
        "established, so the sequences are from the publication rather than "
        "from dbSNP and the row stays unverified.",
    ),
}


def _apply_audit(backbone: dict[str, dict]) -> None:
    """Fold the recorded audit results into the backbone at import time.

    Deliberately separate from the table above so that a reader can see what
    was asserted from literature recall and what was later measured, rather
    than finding the two silently merged into one row.
    """
    for entry in backbone.values():
        marker = entry.get("marker")
        if not marker:
            continue

        carries = _DBSNP_REF_CARRIES.get(marker)
        if carries:
            entry["assembly"] = _DBSNP_ASSEMBLY
            entry["ref_carries"] = carries
            entry["dbsnp_checked"] = True
        if marker in _DBSNP_MULTI_ALLELIC:
            entry["multi_allelic"] = True

        correction = _INDEL_CORRECTIONS.get(marker)
        if correction:
            variant_type, ancestral_seq, derived_seq, note = correction
            entry["variant_type"] = variant_type
            entry["ancestral_seq"] = ancestral_seq
            entry["derived_seq"] = derived_seq
            # The single-base fields are cleared rather than stuffed with a
            # multi-base string. Everything that reads them, marker_state
            # included, treats them as one base, and a widened value would be
            # compared against an array call and silently never match.
            entry["ancestral"] = None
            entry["derived"] = None
            entry["verified"] = False
            entry["note"] = note
            # M17 WAS checked against dbSNP, which is how the class error was
            # found, so dbsnp_checked records that truthfully. ref_carries stays
            # None: an indel has no single reference base to carry a state.
            # M91 has no rsID, so nothing was checked and the flag stays false.
            if entry.get("rsid"):
                entry["assembly"] = _DBSNP_ASSEMBLY
                entry["dbsnp_checked"] = True

        if marker == "M20":
            entry["note"] = _M20_CONFLICT


_apply_audit(Y_BACKBONE)


# ---------------------------------------------------------------------------
# The Karafet 2008 supplement, folded in as a THIRD and separate layer.
#
# Karafet TM, Mendez FL, Meilerman MB, Underhill PA, Zegura SL, Hammer MF.
# "New binary polymorphisms reshape and increase resolution of the human Y
# chromosomal haplogroup tree." Genome Research 18:830-838, 2008.
# Supplementary Table 1, 599 markers over 12 positional columns.
#
# WHY THIS IS A THIRD PASS AND NOT MERGED INTO THE TABLE ABOVE
# -----------------------------------------------------------
# There are now three kinds of claim in this module and they have different
# strengths. The literal table is what was asserted from literature recall.
# `_apply_audit` is what dbSNP MEASURED. This pass is what a primary
# publication STATES. Merging them into one row would make it impossible to
# tell afterwards which rows are load-bearing on which source, which is the
# same reason `_apply_audit` was kept separate in v3.3.0.
#
# WHAT A PUBLICATION BUYS AND DOES NOT BUY
# ----------------------------------------
# It buys the rsID, the allele pair and the variant class. It does NOT buy
# `verified`. That flag also requires the reference orientation, and dbSNP is
# the only thing here that supplies it. So nothing below is promoted.
#
# EXTRACTION, AND WHY THE METHOD IS RECORDED
# ------------------------------------------
# Supplementary Table 1 is laid out in POSITIONAL columns. A reading-order PDF
# extractor can lift a RefSNP ID onto the wrong marker and emit a row that
# looks entirely normal. These values came from `pdftotext -layout`, and the
# column ordering was checked against pages 2, 3, 4 and 7 rendered at 170 dpi
# before any row was trusted. 583 of 599 rows parse under a strict column
# regex; the 16 that do not are line-wrap artifacts.
# ---------------------------------------------------------------------------

_KARAFET = "Karafet et al. 2008, Genome Research 18:830, Supplementary Table 1"


def _k(row: int, mutation: str, ypos: str) -> str:
    return (f"{_KARAFET} row {row}: {mutation} at Y-position {ypos} "
            f"(2008 assembly, NOT liftover-safe).")


# marker -> (row, rsID, mutation, Y-position, extra). The stored allele pair
# already matches the source in every row here; only the rsID is new.
_KARAFET_RSIDS: dict[str, tuple] = {
    "M2":   (1,   "rs3893",     "A->G", "12606580", "Listed as M2=SY81."),
    "M3":   (2,   "rs3894",     "C->T", "17605757", ""),
    "M69":  (64,  "rs2032673",  "T->C", "20353446", ""),
    "M130": (570, "rs35284970", "C->T", "2794854",  "Listed as RPS4Y711=M130."),
    "M172": (161, "rs2032604",  "T->G", "13479028", ""),
    "M174": (163, "rs2032602",  "T->C", "13463674", ""),
    "M184": (171, "rs20320",    "G->A", "13407557", "Listed as M184=USP9Y+3178."),
    "M214": (200, "rs2032674",  "T->C", "13981319", ""),
    "M217": (203, "rs2032668",  "A->C", "13946727", ""),
    "M231": (213, "rs9341278",  "G->A", "13979118", ""),
    "M438": (481, "rs17307294", "A->G", "15148198", "Listed as P215=M438."),
    "P143": (416, "rs4141886",  "G->A", "12707867",
             "The stored C->T pair is the exact reverse complement of the "
             "source G->A: the same call read on the opposite strand. The "
             "pair is left as it stands and the strand difference is recorded "
             "here rather than silently rewritten."),
}

# Markers the survey genotyped and assigned NO RefSNP ID. Recorded so the
# absence reads as measured rather than as an oversight. rsid stays None.
_KARAFET_NO_RSID: dict[str, tuple] = {
    "M35":   (32,  "G->C", "20201091"),
    "P15":   (294, "C->T", "21653414"),
    "P37.2": (315, "T->C", "13001692"),
    "M410":  (276, "A->G", "2811678"),
    "M122":  (115, "T->C", "20224062"),
}

# marker -> (row, rsID or None, ancestral, derived, Y-position, why).
# The stored pair ran the wrong way round. Not a strand difference: reversing
# the direction is required, and for M267 complementing as well.
_KARAFET_TRANSPOSED: dict[str, tuple] = {
    "M145": (135, "rs3848982", "G", "A", "20176596",
             "Was recorded A->G. Listed as M145=P205."),
    "M178": (165, None, "T", "C", "20201143",
             "Was recorded C->T. The survey assigned no RefSNP ID."),
    "M267": (220, "rs9341313", "T", "G", "21151206",
             "Was recorded C->A, which is the source pair complemented AND "
             "reversed. The source orientation is adopted verbatim."),
}

# marker -> (row, rsID, variant_type, ancestral_seq, derived_seq, position, why).
# The same class error v3.3.0 found in M17 and M91, in two more markers that
# were not flagged then. dbSNP agrees with the class on both.
_KARAFET_INDELS: dict[str, tuple] = {
    "M60": (54, "rs2032623", "ins", "", "T", "20337461-20337460",
            "Single-base insertion, NOT a C>A substitution. dbSNP rs2032623 "
            "returns snp_class ins, which agrees. This row previously recorded "
            "ancestral C, derived A; neither allele describes the variant."),
    "M175": (164, "rs2032678", "del", "CTTCTCTTCTC", "CTTCTC", "14018100-14018104",
             "Five-base deletion, NOT an A>G substitution. dbSNP rs2032678 "
             "returns snp_class delins with SPDI alleles CTTCTCTTCTC/CTTCTC, "
             "which agrees. M175 defines O, so that node cannot be called from "
             "array base calls at all."),
}

# marker -> (row, refused value, why nothing is written).
_KARAFET_HOLDS: dict[str, tuple] = {
    "M31": (28, "G->C at Y-position 20199142, no RefSNP ID assigned",
            "The stored pair is G->A. The ancestral base agrees and the "
            "derived base does not. No strand or direction operation "
            "reconciles A with C."),
    "M429": (399, "rs17306671, T->A at Y-position 12541334, listed as P125=M429",
             "The stored pair is A->C, and the reverse complement of T->A is "
             "A->T. Because the alleles conflict, the rsID that travels with "
             "them is not written either."),
}

# marker -> why this source can never resolve it. Recorded so no future
# session spends effort re-reading a supplement that predates the marker.
_KARAFET_ABSENT: dict[str, str] = {
    "F1329": "F series, Wei et al. 2013.",
    "F929":  "F series, Wei et al. 2013.",
    "L15":   "FTDNA L discovery series.",
    "L298":  "FTDNA L discovery series.",
    "P331":  "The P series in that paper stops short of P331.",
    "M420":  "Underhill et al. 2010.",
}

_M20_RESOLVED = (
    "RESOLVED. This row recorded ancestral A, derived C, and dbSNP gives "
    "rs3911 as A/G with no C allele at the site, which is why v3.2.1 could not "
    "reconcile the two. " + _KARAFET + " row 19 gives M20 as A->G at "
    "Y-position 20192842. The supplement and dbSNP agree, so the derived "
    "allele is G and the stored C was simply wrong. The rsID was never in "
    "doubt. Superseded the conflict note carried since v3.2.1."
)

_M91_RSID = (
    " " + _KARAFET + " row 84 assigns M91 rs2032651 and states the 9T to 8T "
    "contraction independently of the article body, so the sequences above are "
    "now confirmed by the supplement as well. NCBI esummary returns an empty "
    "record for rs2032651, which is a 2001-era accession and has most likely "
    "been merged forward, so dbsnp_checked stays false. The supplement also "
    "assigns M91 to haplogroup A rather than BT; that is a 2008-tree "
    "convention difference and is NOT acted on here."
)


def _apply_karafet(backbone: dict[str, dict]) -> None:
    """Fold the 2008 supplement in. Runs AFTER the dbSNP audit, by design."""
    for entry in backbone.values():
        marker = entry.get("marker")
        if not marker:
            continue

        if marker in _KARAFET_RSIDS:
            row, rsid, mutation, ypos, extra = _KARAFET_RSIDS[marker]
            entry["rsid"] = rsid
            entry["note"] = " ".join(x for x in (_k(row, mutation, ypos), extra) if x)

        elif marker in _KARAFET_NO_RSID:
            row, mutation, ypos = _KARAFET_NO_RSID[marker]
            entry["note"] = (_k(row, mutation, ypos) + " The survey genotyped "
                             "this marker and assigned it no RefSNP ID, so "
                             "rsid stays None rather than being guessed.")

        elif marker in _KARAFET_TRANSPOSED:
            row, rsid, anc, der, ypos, why = _KARAFET_TRANSPOSED[marker]
            entry["rsid"] = rsid
            entry["ancestral"] = anc
            entry["derived"] = der
            entry["note"] = _k(row, f"{anc}->{der}", ypos) + " " + why

        elif marker in _KARAFET_INDELS:
            row, rsid, vtype, anc_seq, der_seq, ypos, why = _KARAFET_INDELS[marker]
            entry["rsid"] = rsid
            entry["variant_type"] = vtype
            # An INSERTION has an empty ancestral sequence by definition.
            # That is "" and not None: None means not established, "" means
            # established as nothing there, and collapsing the two would
            # lose the distinction this module exists to preserve.
            entry["ancestral_seq"] = anc_seq
            entry["derived_seq"] = der_seq
            # Same rule as _apply_audit: clear the single-base fields rather
            # than widening them, because every reader treats them as one base.
            entry["ancestral"] = None
            entry["derived"] = None
            entry["verified"] = False
            entry["assembly"] = _DBSNP_ASSEMBLY
            entry["dbsnp_checked"] = True
            entry["note"] = _k(row, vtype, ypos) + " " + why

        elif marker in _KARAFET_HOLDS:
            row, refused, why = _KARAFET_HOLDS[marker]
            entry["note"] = (
                f"HELD. {_KARAFET} row {row} gives {refused}. {why} Exactly one "
                f"side is wrong and nothing available says which, so NO value "
                f"is written. " + _UNVERIFIED_Y)

        elif marker in _KARAFET_ABSENT:
            entry["note"] = (
                f"Absent from Karafet et al. 2008; the marker post-dates that "
                f"survey, so the 2008 supplement can never resolve it. "
                f"{_KARAFET_ABSENT[marker]} " + _UNVERIFIED_Y)

        elif marker == "M20":
            entry["derived"] = "G"
            entry["note"] = _M20_RESOLVED

        elif marker == "M91":
            entry["rsid"] = "rs2032651"
            entry["note"] = entry["note"] + _M91_RSID


_apply_karafet(Y_BACKBONE)


# ---------------------------------------------------------------------------
# Bundled mtDNA backbone.
#
# mtDNA is expressed as rCRS positions rather than rsIDs, which removes the
# rsID-mapping hazard entirely: position 10400 means one thing and only one
# thing. So the textbook macro-haplogroup rows here ARE marked verified.
#
# Direction convention, and why it is not "derived":
# rCRS is itself an H2a2a1 sequence, so several macro-haplogroup markers appear
# as back mutations relative to it. H is defined by CARRYING the rCRS base at
# 2706 and 7028 while everything outside HV carries the other base. Calling
# that "ancestral" would be flatly wrong. So each mt node records the base the
# haplogroup CARRIES plus the rCRS base, and support is "observed base equals
# the carried base". No derived/ancestral direction is claimed at all.
#
# This only works because the walk is hierarchical. Read in isolation, HV's
# 14766C would also "support" an L or M sample, since they carry the rCRS base
# there too. HV is only ever tested once R0 and R and N and L3 are established,
# and that parent constraint is what makes the marker discriminating.
#
# The L-lineage nodes are NOT verified. Their defining positions were recorded
# from recall and the direction of several is genuinely uncertain; the data
# builder must confirm them against PhyloTree.
# ---------------------------------------------------------------------------

def _v(position: int, base: str, rcrs: str, *, verified: bool = False,
       rsid: str | None = None) -> dict:
    """One defining mtDNA variant: the base this haplogroup carries at a position."""
    return {
        "position": int(position),
        "base": base.upper(),
        "rcrs": rcrs.upper(),
        "rsid": rsid,
        "verified": bool(verified),
    }


_UNVERIFIED_MT = (
    "Defining positions recorded from recall; the data builder must confirm "
    "them against PhyloTree before they are trusted."
)

MT_BACKBONE: dict[str, dict] = {
    "root": {
        "system": "MT", "node": "root", "parent": None, "label": "mt-MRCA",
        "defining": [], "verified": True,
        "note": "Root carries no defining variant. Everyone is here.",
    },
    # -- African L lineages, none verified -----------------------------------
    "L0": {
        "system": "MT", "node": "L0", "parent": "root", "label": "L0",
        "defining": [_v(1048, "T", "C")], "note": _UNVERIFIED_MT,
    },
    "L1'2'3'4'5'6": {
        "system": "MT", "node": "L1'2'3'4'5'6", "parent": "root",
        "label": "L1'2'3'4'5'6",
        "defining": [_v(3594, "T", "C"), _v(10819, "G", "A")],
        "note": "3594 is the classic L0 versus rest split. " + _UNVERIFIED_MT,
    },
    "L1": {
        "system": "MT", "node": "L1", "parent": "L1'2'3'4'5'6", "label": "L1",
        "defining": [_v(3666, "A", "G")], "note": _UNVERIFIED_MT,
    },
    "L2'3'4'5'6": {
        "system": "MT", "node": "L2'3'4'5'6", "parent": "L1'2'3'4'5'6",
        "label": "L2'3'4'5'6",
        "defining": [_v(2758, "A", "G"), _v(2885, "C", "T"), _v(7146, "G", "A")],
        "note": _UNVERIFIED_MT,
    },
    "L2": {
        "system": "MT", "node": "L2", "parent": "L2'3'4'5'6", "label": "L2",
        "defining": [_v(2416, "C", "T"), _v(8206, "A", "G")], "note": _UNVERIFIED_MT,
    },
    "L4": {
        "system": "MT", "node": "L4", "parent": "L2'3'4'5'6", "label": "L4",
        "defining": [_v(5460, "A", "G")], "note": _UNVERIFIED_MT,
    },
    "L5": {
        "system": "MT", "node": "L5", "parent": "L2'3'4'5'6", "label": "L5",
        "defining": [_v(3423, "A", "G")], "note": _UNVERIFIED_MT,
    },
    "L6": {
        "system": "MT", "node": "L6", "parent": "L2'3'4'5'6", "label": "L6",
        "defining": [_v(5836, "G", "A")], "note": _UNVERIFIED_MT,
    },
    "L3": {
        "system": "MT", "node": "L3", "parent": "L2'3'4'5'6", "label": "L3",
        "defining": [_v(769, "A", "G"), _v(1018, "A", "G"), _v(16311, "C", "T")],
        "note": _UNVERIFIED_MT,
    },
    # -- M and its two array-visible branches --------------------------------
    "M": {
        "system": "MT", "node": "M", "parent": "L3", "label": "M",
        "defining": [_v(10400, "T", "C", verified=True), _v(14783, "C", "T")],
        "note": "10400C>T is the textbook M marker and is verified.",
    },
    "C": {
        "system": "MT", "node": "C", "parent": "M", "label": "C",
        "defining": [_v(13263, "G", "A", verified=True)],
        "note": "",
    },
    "D": {
        "system": "MT", "node": "D", "parent": "M", "label": "D",
        "defining": [_v(5178, "A", "C", verified=True)],
        "note": "",
    },
    # -- N and its branches --------------------------------------------------
    "N": {
        "system": "MT", "node": "N", "parent": "L3", "label": "N",
        "defining": [_v(8701, "A", "A"), _v(9540, "T", "T"),
                     _v(10873, "T", "T"), _v(15301, "G", "G")],
        "note": "N is defined by rCRS-state back mutations relative to L3 and "
                "the direction of each is uncertain here. " + _UNVERIFIED_MT,
    },
    "A": {
        "system": "MT", "node": "A", "parent": "N", "label": "A",
        "defining": [_v(663, "G", "A", verified=True), _v(8794, "T", "C")],
        "note": "",
    },
    "I": {
        "system": "MT", "node": "I", "parent": "N", "label": "I",
        "defining": [_v(10034, "C", "T", verified=True), _v(8616, "T", "G")],
        "note": "",
    },
    "W": {
        "system": "MT", "node": "W", "parent": "N", "label": "W",
        "defining": [_v(8994, "A", "G", verified=True), _v(3505, "G", "A")],
        "note": "",
    },
    "X": {
        "system": "MT", "node": "X", "parent": "N", "label": "X",
        "defining": [_v(13966, "G", "A", verified=True), _v(6221, "C", "T")],
        "note": "",
    },
    # -- R and everything European arrays actually resolve -------------------
    "R": {
        "system": "MT", "node": "R", "parent": "N", "label": "R",
        "defining": [_v(12705, "C", "C"), _v(16223, "C", "C")],
        "note": "R is defined by back mutations to the rCRS state at 12705 and "
                "16223. " + _UNVERIFIED_MT,
    },
    "B": {
        "system": "MT", "node": "B", "parent": "R", "label": "B",
        "defining": [_v(16189, "C", "T")],
        "note": "B's real defining feature is the 8281 to 8289 nine base pair "
                "deletion, which no consumer array calls. The 16189 proxy used "
                "here is weak and is not verified. " + _UNVERIFIED_MT,
    },
    "R0": {
        "system": "MT", "node": "R0", "parent": "R", "label": "R0",
        "defining": [_v(73, "A", "A")],
        "note": _UNVERIFIED_MT,
    },
    "U": {
        "system": "MT", "node": "U", "parent": "R", "label": "U",
        "defining": [_v(11467, "G", "A", verified=True),
                     _v(12308, "G", "A", verified=True),
                     _v(12372, "A", "G", verified=True)],
        "note": "The 11467, 12308, 12372 triple is the textbook U motif.",
    },
    "K": {
        "system": "MT", "node": "K", "parent": "U", "label": "K",
        "defining": [_v(9055, "A", "G", verified=True), _v(14798, "C", "T")],
        "note": "",
    },
    "JT": {
        "system": "MT", "node": "JT", "parent": "R", "label": "JT",
        "defining": [_v(4216, "C", "T"), _v(15452, "A", "C")],
        "note": _UNVERIFIED_MT,
    },
    "J": {
        "system": "MT", "node": "J", "parent": "JT", "label": "J",
        "defining": [_v(13708, "A", "G", verified=True),
                     _v(16069, "T", "C", verified=True)],
        "note": "",
    },
    "T": {
        "system": "MT", "node": "T", "parent": "JT", "label": "T",
        "defining": [_v(4917, "G", "A", verified=True), _v(13368, "A", "G")],
        "note": "",
    },
    "HV": {
        "system": "MT", "node": "HV", "parent": "R0", "label": "HV",
        "defining": [_v(14766, "C", "C")],
        "note": "14766 is a back mutation to the rCRS state and only "
                "discriminates below R0. " + _UNVERIFIED_MT,
    },
    "H": {
        "system": "MT", "node": "H", "parent": "HV", "label": "H",
        "defining": [_v(2706, "A", "A", verified=True),
                     _v(7028, "C", "C", verified=True)],
        "note": "rCRS is itself an H sequence, so H is defined by carrying the "
                "rCRS base at both positions.",
    },
    "V": {
        "system": "MT", "node": "V", "parent": "HV", "label": "V",
        "defining": [_v(4580, "A", "G", verified=True)],
        "note": "",
    },
}


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def _children(tree: dict) -> dict[str, list[str]]:
    """Map parent node name to its children, in a stable order."""
    out: dict[str, list[str]] = {name: [] for name in tree}
    for name, entry in tree.items():
        parent = entry.get("parent")
        if parent is None:
            continue
        out.setdefault(parent, []).append(name)
    for name in out:
        out[name].sort()
    return out


def path_to(node: str, tree: dict | None = None) -> list[str]:
    """Root-first chain of node names ending at ``node``.

    Returns [] for an unknown node. A cycle in the table would hang the walk,
    so the loop is bounded by the node count and simply stops.
    """
    tree = Y_BACKBONE if tree is None else tree
    if node not in tree:
        return []
    chain: list[str] = []
    seen: set[str] = set()
    current: str | None = node
    while current is not None and current in tree and current not in seen:
        seen.add(current)
        chain.append(current)
        current = tree[current].get("parent")
    chain.reverse()
    return chain


def _depth(node: str, tree: dict) -> int:
    return max(0, len(path_to(node, tree)) - 1)


def marker_keys(entry: dict) -> list[str]:
    """Lookup keys for a Y marker, most specific first.

    The rsID is tried first when the table records one. The marker NAME is
    always tried too, because a builder-supplied genotype map may be keyed by
    marker name for rows whose rsID is not yet confirmed. Both are lowercased
    to match the project's genotype-map convention.
    """
    keys: list[str] = []
    rsid = entry.get("rsid")
    if rsid:
        keys.append(str(rsid).strip().lower())
    marker = entry.get("marker")
    if marker:
        keys.append(str(marker).strip().lower())
    return keys


def _mt_keys(variant: dict) -> list[str]:
    """Lookup keys for one mtDNA position.

    Consumer vendors name mitochondrial positions inconsistently: 23andMe uses
    internal ``i`` identifiers for many of them, others use rsIDs, and neither
    is stable across chip versions. So a position-keyed view is the reliable
    input and ``mt_positions_from_merged`` builds one. The rsID and the bare
    position string are accepted as well so a caller is never forced to
    reshape data twice.
    """
    position = int(variant.get("position") or 0)
    keys: list[str] = []
    rsid = variant.get("rsid")
    if rsid:
        keys.append(str(rsid).strip().lower())
    keys.extend([f"mt{position}", f"m{position}", f"mt-{position}", str(position)])
    return keys


_LEADING_LETTERS = re.compile(r"^[A-Za-z]+")


def snp_name(entry: dict) -> str | None:
    """The SNP-based name for a node, for example ``J-M267`` for J1.

    Two naming systems are in circulation and they describe the same nodes.
    The older letter-number form (J1, R1b1a1b) encodes a position in the tree,
    so it gets renumbered every time the tree is revised and a saved result
    silently goes stale. The SNP-based form (J-M267, R-M269) names the mutation
    instead, which does not move. FamilyTreeDNA reports the second, this
    project's labels are the first, and a user comparing the two had no way to
    know they were looking at one node rather than two.

    That confusion is not hypothetical. A widely shared claim in August 2026
    held that consumer vendors "misclassify J1 as generic J-M267". M267 is the
    SNP that defines J1, so the two strings are the same haplogroup and there
    was no misclassification to report. Emitting both names removes the
    ambiguity at the source.

    Returns None for the root, which no marker defines.
    """
    marker = (entry or {}).get("marker")
    if not marker:
        return None
    label = entry.get("label") or entry.get("node") or ""
    match = _LEADING_LETTERS.match(str(label))
    if not match:
        return None
    return f"{match.group(0)}-{marker}"


def equivalent_names(entry: dict) -> list[str]:
    """Every name this one node answers to, label first, no duplicates.

    Order is deliberate. The label leads because it is what the rest of this
    codebase keys on, and callers that take the first element get the value
    they got before this function existed.
    """
    entry = entry or {}
    names: list[str] = []
    for candidate in (entry.get("label"), snp_name(entry), entry.get("node")):
        if candidate and candidate not in names:
            names.append(candidate)
    return names


def untypeable_markers() -> dict:
    """Markers that are not base substitutions and cannot be called from an array.

    Separate from ``unverified_markers`` because these are a different failure.
    An unverified marker might be right and nobody checked. These are known to
    be wrong in KIND: a consumer array reports two base calls per position, and
    a deletion is not a base, so no genotype can ever satisfy the rule. Reporting
    them as merely unusable would hide that the ceiling is structural rather
    than a matter of coverage, which is the same distinction between "not
    present" and "never checked" that the rest of this module holds.
    """
    rows = [
        {
            "node": entry["node"],
            "marker": entry.get("marker"),
            "label": entry.get("label"),
            "variant_type": entry.get("variant_type"),
            "ancestral_seq": entry.get("ancestral_seq"),
            "derived_seq": entry.get("derived_seq"),
            "note": entry.get("note", ""),
        }
        for entry in Y_BACKBONE.values()
        if entry.get("variant_type") not in (None, "snv")
    ]
    rows.sort(key=lambda row: str(row["marker"]))
    return {
        "system": "Y",
        "count": len(rows),
        "markers": rows,
        "reason": (
            "These markers are insertions or deletions, not base substitutions. "
            "A consumer array reports two base calls at a position, so no array "
            "genotype can satisfy them and the nodes they define cannot be "
            "reached from array data at all. This is a structural ceiling, not "
            "missing coverage."
        ),
    }


def unverified_markers() -> dict:
    """Every backbone entry whose rsID or alleles are not confirmed.

    This is the audit list the data builder works from. It is returned rather
    than logged because an unverified marker that only appears in a log is an
    unverified marker nobody fixes.
    """
    y_rows = [
        {
            "node": entry["node"], "marker": entry.get("marker"),
            "rsid": entry.get("rsid"), "note": entry.get("note", ""),
            # What the builder still has to establish for this row. Carried
            # here so the audit list states the remaining work instead of just
            # naming the row: an rsID alone cannot settle which state is
            # ancestral, because dbSNP answers reference over alternate.
            "assembly": entry.get("assembly"),
            "ref_carries": entry.get("ref_carries"),
            "dbsnp_checked": bool(entry.get("dbsnp_checked")),
        }
        for name, entry in sorted(Y_BACKBONE.items())
        if name != "root" and not entry.get("verified")
    ]
    mt_rows = []
    for name, entry in sorted(MT_BACKBONE.items()):
        if name == "root":
            continue
        bad = [v for v in entry.get("defining", []) if not v.get("verified")]
        if bad:
            mt_rows.append({
                "node": entry["node"],
                "positions": [v["position"] for v in bad],
                "note": entry.get("note", ""),
            })
    return {
        "y": y_rows,
        "mt": mt_rows,
        "y_total": len([n for n in Y_BACKBONE if n != "root"]),
        "mt_total": len([n for n in MT_BACKBONE if n != "root"]),
        "policy": (
            "No entry may be marked verified until it has been confirmed "
            "against ISOGG or YFull for Y and PhyloTree for mtDNA. An "
            "unverified marker still produces a call, flagged provisional."
        ),
    }


# ---------------------------------------------------------------------------
# Genotype access
# ---------------------------------------------------------------------------

def _alleles_of(value: Any) -> set[str]:
    """Uppercase real alleles in a genotype value, empty set for a no-call.

    Accepts the project's ``(allele1, allele2)`` tuple, a 2-character string,
    and the merged-genotype dict shape, because Wave 3 gets called with all
    three depending on which side of the pipeline is asking.
    """
    if value is None:
        return set()
    if isinstance(value, dict):
        raw: Iterable[Any] = (value.get("allele1"), value.get("allele2"))
    elif isinstance(value, (list, tuple)):
        raw = value
    else:
        text = str(value).strip().upper()
        for sep in ("/", "|", ";", " ", "\t"):
            text = text.replace(sep, "")
        raw = list(text)
    out: set[str] = set()
    for allele in raw:
        token = str(allele if allele is not None else "").strip().upper()
        if token and token not in _NOCALL_ALLELES:
            out.add(token)
    return out


def _lookup(genotypes: dict, keys: Sequence[str]) -> tuple[bool, set[str]]:
    """(key_present, alleles). Present-but-no-call is NOT the same as absent."""
    genotypes = genotypes or {}
    for key in keys:
        if key in genotypes:
            return True, _alleles_of(genotypes[key])
    return False, set()


def marker_state(entry: dict, genotypes: dict, *,
                 verified_only: bool = False) -> tuple[str, set[str]]:
    """Classify one Y marker for this person.

    Returns one of DERIVED, ANCESTRAL, DISCORDANT, NO_CALL, NOT_ON_ARRAY or
    UNUSABLE, plus the observed alleles.

    NOT_ON_ARRAY and NO_CALL are kept apart from ANCESTRAL on purpose. Folding
    them together is the single mistake this module exists to avoid.
    """
    derived = entry.get("derived")
    ancestral = entry.get("ancestral")
    if not derived:
        return UNUSABLE, set()
    if verified_only and not entry.get("verified"):
        return UNUSABLE, set()
    present, alleles = _lookup(genotypes, marker_keys(entry))
    if not present:
        return NOT_ON_ARRAY, set()
    if not alleles:
        return NO_CALL, set()
    if derived in alleles:
        return DERIVED, alleles
    if ancestral and ancestral in alleles:
        return ANCESTRAL, alleles
    return DISCORDANT, alleles


def _mt_variant_state(variant: dict, genotypes: dict, *,
                      verified_only: bool = False) -> tuple[str, set[str]]:
    """Classify one mtDNA position: supported, contradicted, or not testable.

    DERIVED here means "carries the base this haplogroup is defined by". The
    word is reused only so the tri-state buckets have one vocabulary across
    both systems; no derived/ancestral direction is being claimed for mtDNA.
    """
    if verified_only and not variant.get("verified"):
        return UNUSABLE, set()
    present, alleles = _lookup(genotypes, _mt_keys(variant))
    if not present:
        return NOT_ON_ARRAY, set()
    if not alleles:
        return NO_CALL, set()
    return (DERIVED if variant["base"] in alleles else ANCESTRAL), alleles


# ---------------------------------------------------------------------------
# Bucket bookkeeping, the genosets tri-state applied to a tree
# ---------------------------------------------------------------------------

def _empty_buckets() -> dict[str, list]:
    return {DERIVED: [], ANCESTRAL: [], DISCORDANT: [],
            NO_CALL: [], NOT_ON_ARRAY: [], UNUSABLE: []}


def _bucket_summary(buckets: dict[str, list]) -> dict:
    not_testable = (len(buckets[NO_CALL]) + len(buckets[NOT_ON_ARRAY])
                    + len(buckets[UNUSABLE]))
    return {
        "derived": len(buckets[DERIVED]),
        "ancestral": len(buckets[ANCESTRAL]),
        "discordant": len(buckets[DISCORDANT]),
        "no_call": len(buckets[NO_CALL]),
        "not_on_array": len(buckets[NOT_ON_ARRAY]),
        "unusable": len(buckets[UNUSABLE]),
        "not_testable": not_testable,
        "tested": len(buckets[DERIVED]) + len(buckets[ANCESTRAL]) + len(buckets[DISCORDANT]),
    }


_TRI_STATE_NOTE = (
    "Markers that are not on your array are reported as not testable, never as "
    "ancestral. Not tested is not the same as negative."
)

_ASSUMED_NOTE = (
    "Some nodes on this path carry no marker your array reads. The walk passed "
    "through them on the strength of a deeper marker rather than testing them, "
    "and they are listed as assumed. A deeper call resting on an assumed "
    "ancestor is only as good as that assumption."
)


def _subtree_supported(node: str, children: dict[str, list[str]],
                       supported: dict[str, bool]) -> bool:
    """True when any descendant of ``node`` reads positively.

    This is what makes it safe to pass through an untested intermediate node:
    the walk only skips a node when something BELOW it is positively called, so
    it is never inventing a branch out of pure absence.
    """
    stack = list(children.get(node, []))
    seen: set[str] = set()
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        if supported.get(current):
            return True
        stack.extend(children.get(current, []))
    return False


def _advance_through_untestable(skippable: Sequence[str], chain: list[str],
                                assumed: list[str]) -> str | None:
    """Move the walk through one untested node, or refuse to choose a branch.

    Both backbones need exactly this step. A node whose own marker or defining
    positions the array never read can still be passed through when exactly one
    such branch has support below it, which is what lets a consumer chip resolve
    past a backbone marker it does not carry. Anything other than exactly one
    candidate is not a choice the data can settle, so the caller stops the walk;
    picking a branch there would invent a lineage.

    Appends the node to ``chain`` and to ``assumed`` in place, because a node
    the walk passed through must never be counted as tested. Returns the node
    moved to, or None when the walk must stop.
    """
    if len(skippable) != 1:
        return None
    node = skippable[0]
    chain.append(node)
    assumed.append(node)
    return node


def _ambiguous_skip_conflict(node: str, skippable: Sequence[str], *,
                             children_key: str, note: str,
                             markers: Sequence[Any] | None = None) -> dict:
    """Record that more than one untested branch could have been walked into.

    The key naming the branches is passed in rather than fixed, because the Y
    payload calls them ``derived_children`` and the mtDNA payload calls them
    ``candidates``. Both names are already published in docs/API_V3.md, so
    unifying them here would be an API change wearing a refactor's clothes.
    """
    record: dict = {"at": node, children_key: sorted(skippable)}
    if markers is not None:
        record["markers"] = list(markers)
    record["note"] = note
    return record


def _backbone_payload(system: str, tree: dict, *, resolved: bool,
                      stopped_at: str, chain: list[str],
                      path_markers: list[dict], assumed: list[str],
                      buckets: dict[str, list], conflicts: list[dict],
                      blocked_by_array: list[str], excluded: list[str],
                      verified_used: int, unverified_used: int,
                      caveats: list[str], extra: dict | None = None) -> dict:
    """Assemble the response both backbone walks return.

    Y and mtDNA answer the same question from different evidence, so the field
    set is deliberately identical and a caller that can read one can read the
    other. Building it in one place is what stops the two from drifting apart,
    which here would be an honesty bug rather than a cosmetic one:
    ``stopped_because_not_testable`` and ``stopped_because_excluded`` are the
    difference between "your array cannot see this branch" and "your line is
    genuinely not on it", and never confusing those two is invariant 3 of this
    project.

    ``extra`` carries fields only one system has, currently the Y sex hint. It
    is merged before ``caveats`` so the key order of both payloads is exactly
    what it was when the dict was written out twice.
    """
    payload = {
        "system": system,
        "available": bool(resolved),
        "state": "called" if resolved else "unresolved",
        "haplogroup": stopped_at if resolved else None,
        "label": tree[stopped_at]["label"] if resolved else None,
        # Both naming systems, so a result can be compared against a vendor
        # report without the reader having to know they describe one node.
        # Empty rather than guessed when nothing resolved: an unresolved call
        # has no node to name, and inventing one here would be the same class
        # of error as reporting an untested marker as ancestral.
        "equivalent_names": (equivalent_names(tree[stopped_at])
                             if resolved else []),
        "also_written_as": (snp_name(tree[stopped_at]) if resolved else None),
        "path": chain,
        "path_markers": path_markers,
        "assumed": assumed,
        "tested_path": [name for name in chain
                        if name != "root" and name not in assumed],
        "buckets": buckets,
        "counts": _bucket_summary(buckets),
        "conflicts": conflicts,
        "stopped_at": stopped_at,
        "stopped_because_not_testable": blocked_by_array,
        "stopped_because_excluded": excluded,
        "source": "bundled_backbone",
        "confidence": ("provisional" if (unverified_used or assumed)
                       else ("moderate" if resolved else "none")),
        "verified_markers_used": verified_used,
        "unverified_markers_used": unverified_used,
    }
    payload.update(extra or {})
    payload["caveats"] = caveats
    payload.update(tree_stamp())
    return payload


# ---------------------------------------------------------------------------
# Y calling
# ---------------------------------------------------------------------------

_FEMALE_HINTS = {"f", "female", "xx", "woman", "2"}
_MALE_HINTS = {"m", "male", "xy", "man", "1"}

NO_Y_MESSAGE = (
    "No Y chromosome data. A Y haplogroup describes a direct paternal line and "
    "can only be called from a Y chromosome, so there is nothing to call here. "
    "This is not a failure and not a missing result: it is the expected state "
    "for a sample with no Y. A direct paternal line can still be traced through "
    "a father, brother or paternal uncle who tests."
)


def call_y_backbone(genotypes: dict, sex_hint: Any = None, *,
                    verified_only: bool = False,
                    skip_untestable: bool = True) -> dict:
    """Walk the bundled Y backbone from the root to the deepest derived node.

    A node is entered only when its own defining marker reads DERIVED. The walk
    stops as soon as no child is derived, which means an ancestral child and a
    not-testable child produce the same STOP but are reported very differently:
    an ancestral child genuinely excludes that branch, a not-testable child
    means the array simply cannot see it and the person may well belong there.

    Two derived children at the same level cannot both be true. That is
    surfaced as a conflict and the walk stops rather than guessing, which is
    the same choice ``merge.py`` makes for pooled file conflicts.

    ``skip_untestable=True``, the default, lets the walk pass through a node
    whose own marker is absent from the array when exactly one such branch has
    a derived marker below it. Without that, a real consumer file resolves to
    nothing at all, because chips carry terminal markers like M269 and skip
    backbone markers like F1329. Every skipped node is listed in ``assumed``
    and never counted as tested. Set it False for a strict walk.

    ``verified_only=True`` refuses every unverified marker. Because no Y row is
    verified yet, that currently yields an unresolved call, which is the honest
    answer while the rsID mappings are unconfirmed.
    """
    genotypes = genotypes or {}
    hint = str(sex_hint or "").strip().lower()

    buckets = _empty_buckets()
    states: dict[str, str] = {}
    observed: dict[str, list[str]] = {}
    for name, entry in Y_BACKBONE.items():
        if name == "root":
            continue
        state, alleles = marker_state(entry, genotypes, verified_only=verified_only)
        states[name] = state
        observed[name] = sorted(alleles)
        buckets[state].append({
            "node": name, "marker": entry.get("marker"), "rsid": entry.get("rsid"),
            "state": state, "observed": sorted(alleles),
            "verified": bool(entry.get("verified")),
        })

    # Whether the file contains Y data at all is decided WITHOUT the
    # verified_only filter. "You have no Y chromosome" and "the bundled tree
    # has no confirmed markers yet" are completely different statements and a
    # woman must never be told the second one.
    has_y_data = any(
        marker_state(entry, genotypes)[0] in (DERIVED, ANCESTRAL, DISCORDANT)
        for name, entry in Y_BACKBONE.items() if name != "root"
    )

    if hint in _FEMALE_HINTS or (not has_y_data and hint not in _MALE_HINTS):
        state = "no_y_data"
        payload = {
            "system": "Y",
            "available": False,
            "state": state,
            "haplogroup": None,
            "label": None,
            "path": [],
            "path_markers": [],
            "buckets": buckets,
            "counts": _bucket_summary(buckets),
            "conflicts": [],
            "source": "bundled_backbone",
            "confidence": "none",
            "verified_markers_used": 0,
            "unverified_markers_used": 0,
            "sex_hint": hint or None,
            "message": NO_Y_MESSAGE,
            "caveats": [NO_Y_MESSAGE, _TRI_STATE_NOTE],
        }
        payload.update(tree_stamp())
        return payload

    children = _children(Y_BACKBONE)
    supported = {name: states.get(name) == DERIVED for name in states}
    current = "root"
    chain = ["root"]
    assumed: list[str] = []
    conflicts: list[dict] = []
    guard_count = 0
    while guard_count <= len(Y_BACKBONE):
        guard_count += 1
        derived_kids = [k for k in children.get(current, []) if states.get(k) == DERIVED]
        if not derived_kids and skip_untestable:
            # No child is positively derived, but a child whose own marker is
            # not on the array may still have a derived descendant. Passing
            # through it is what lets a 23andMe file resolve to R-M269 when the
            # chip carries M269 and not F1329. It is recorded as assumed, never
            # as tested, and an ambiguity here stops the walk instead of
            # picking a branch.
            skippable = [
                k for k in children.get(current, [])
                if states.get(k) in NOT_TESTABLE_STATES
                and _subtree_supported(k, children, supported)
            ]
            moved = _advance_through_untestable(skippable, chain, assumed)
            if moved is None:
                if len(skippable) > 1:
                    conflicts.append(_ambiguous_skip_conflict(
                        current, skippable,
                        children_key="derived_children",
                        markers=[Y_BACKBONE[k].get("marker")
                                 for k in sorted(skippable)],
                        note=("More than one untested branch under this node "
                              "has a derived marker below it. The walk cannot "
                              "choose between them without testing this level, "
                              "so it stops here."),
                    ))
                break
            current = moved
            continue
        if not derived_kids:
            break
        if len(derived_kids) > 1:
            # Two mutually exclusive branches cannot both be derived. Record it
            # and stop; picking one would be inventing a result.
            conflicts.append({
                "at": current,
                "derived_children": sorted(derived_kids),
                "markers": [Y_BACKBONE[k].get("marker") for k in sorted(derived_kids)],
                "note": (
                    "More than one branch under this node reads derived. Both "
                    "cannot be true. The most likely cause is a wrong rsID in "
                    "the bundled tree or a strand mismatch in the source file, "
                    "not an unusual lineage. No winner is picked."
                ),
            })
            break
        current = derived_kids[0]
        chain.append(current)

    resolved = current != "root"
    path_markers = [
        {
            "node": name,
            "marker": Y_BACKBONE[name].get("marker"),
            "rsid": Y_BACKBONE[name].get("rsid"),
            "state": states.get(name),
            "observed": observed.get(name, []),
            "verified": bool(Y_BACKBONE[name].get("verified")),
        }
        for name in chain if name != "root"
    ]
    verified_used = sum(1 for m in path_markers if m["verified"])
    unverified_used = len(path_markers) - verified_used

    # Which branches under the stopping point were genuinely excluded, and
    # which were simply invisible. This is the ceiling made local.
    stopped_at = current
    blocked_by_array = sorted(
        k for k in children.get(stopped_at, [])
        if states.get(k) in NOT_TESTABLE_STATES
    )
    excluded = sorted(
        k for k in children.get(stopped_at, [])
        if states.get(k) in (ANCESTRAL, DISCORDANT)
    )

    caveats = [_TRI_STATE_NOTE]
    if unverified_used:
        caveats.append(
            f"{unverified_used} of the {len(path_markers)} markers on this path "
            f"use an unverified rsID or allele assignment. The call is "
            f"provisional until the data builder confirms them."
        )
    if assumed:
        caveats.append(_ASSUMED_NOTE + " Assumed: " + ", ".join(assumed) + ".")
    if verified_only and buckets[UNUSABLE]:
        caveats.append(
            f"{len(buckets[UNUSABLE])} markers were refused because their rsID "
            f"or allele assignment is not yet confirmed, and verified_only was "
            f"set. Your file may well carry them."
        )
    if blocked_by_array:
        caveats.append(
            "Branches below this point were not tested because your array does "
            "not read their markers: " + ", ".join(blocked_by_array) + "."
        )

    return _backbone_payload(
        "Y", Y_BACKBONE,
        resolved=resolved, stopped_at=stopped_at, chain=chain,
        path_markers=path_markers, assumed=assumed, buckets=buckets,
        conflicts=conflicts, blocked_by_array=blocked_by_array,
        excluded=excluded, verified_used=verified_used,
        unverified_used=unverified_used, caveats=caveats,
        extra={"sex_hint": hint or None},
    )


# ---------------------------------------------------------------------------
# mtDNA calling
# ---------------------------------------------------------------------------

def mt_positions_from_merged(merged: dict) -> dict:
    """Build a position-keyed mtDNA view from a ``merge.merge_sources`` result.

    Returns ``{"mt<position>": (allele1, allele2)}``. Vendors label the
    mitochondrion "MT", "M", "chrM" or "26" depending on the file, so all four
    are accepted. A caller with a plain rsID-keyed map can pass that straight
    to :func:`call_mt_backbone` instead, but then only rows whose rsID the
    bundled tree happens to record will be readable.
    """
    out: dict[str, tuple[str, str]] = {}
    for rsid, entry in ((merged or {}).get("genotypes") or {}).items():
        chrom = str(entry.get("chromosome") or "").strip().upper()
        if chrom.startswith("CHR"):
            chrom = chrom[3:]
        if chrom not in ("MT", "M", "26"):
            continue
        position = entry.get("position") or 0
        try:
            position = int(position)
        except (TypeError, ValueError):
            continue
        if position <= 0:
            continue
        pair = (str(entry.get("allele1") or "N"), str(entry.get("allele2") or "N"))
        out[f"mt{position}"] = pair
        if rsid:
            out.setdefault(str(rsid).strip().lower(), pair)
    return out


def call_mt_backbone(genotypes: dict, *, verified_only: bool = False,
                     skip_untestable: bool = True) -> dict:
    """Walk the bundled mtDNA backbone from the root to the deepest supported node.

    A child node is entered when at least one of its defining positions is
    readable, all readable defining positions carry the haplogroup's base, and
    none contradicts. A node with no readable position at all is not entered
    and is reported as not testable, never as excluded.

    ``skip_untestable=True``, the default, lets the walk pass through a node
    none of whose positions the array reads, when exactly one such branch has a
    supported node below it. Arrays cover the mitochondrion patchily, roughly
    2,500 of 16,569 bases on a typical chip, so a strict walk usually stops at
    the root and tells the user nothing. Skipped nodes are listed in
    ``assumed`` and never counted as tested.
    """
    genotypes = genotypes or {}
    buckets = _empty_buckets()
    node_support: dict[str, dict] = {}

    for name, entry in MT_BACKBONE.items():
        if name == "root":
            continue
        support = 0
        contradict = 0
        rows = []
        for variant in entry.get("defining", []):
            state, alleles = _mt_variant_state(variant, genotypes,
                                               verified_only=verified_only)
            rows.append({
                "node": name, "position": variant["position"],
                "expected": variant["base"], "rcrs": variant["rcrs"],
                "state": state, "observed": sorted(alleles),
                "verified": bool(variant.get("verified")),
            })
            buckets[state].append(rows[-1])
            if state == DERIVED:
                support += 1
            elif state == ANCESTRAL:
                contradict += 1
        node_support[name] = {
            "support": support, "contradict": contradict,
            "testable": support + contradict, "rows": rows,
        }

    children = _children(MT_BACKBONE)
    supported_map = {
        name: (info["support"] > 0 and info["contradict"] == 0)
        for name, info in node_support.items()
    }
    current = "root"
    chain = ["root"]
    assumed: list[str] = []
    conflicts: list[dict] = []
    guard_count = 0
    while guard_count <= len(MT_BACKBONE):
        guard_count += 1
        candidates = [
            k for k in children.get(current, [])
            if node_support.get(k, {}).get("support", 0) > 0
            and node_support.get(k, {}).get("contradict", 0) == 0
        ]
        if not candidates and skip_untestable:
            skippable = [
                k for k in children.get(current, [])
                if node_support.get(k, {}).get("testable", 0) == 0
                and _subtree_supported(k, children, supported_map)
            ]
            moved = _advance_through_untestable(skippable, chain, assumed)
            if moved is None:
                if len(skippable) > 1:
                    conflicts.append(_ambiguous_skip_conflict(
                        current, skippable,
                        children_key="candidates",
                        note=("More than one untested branch under this node "
                              "has a supported node below it. The walk stops "
                              "rather than choosing between them."),
                    ))
                break
            current = moved
            continue
        if not candidates:
            break
        if len(candidates) > 1:
            best = max(candidates, key=lambda k: node_support[k]["support"])
            top = [k for k in candidates
                   if node_support[k]["support"] == node_support[best]["support"]]
            conflicts.append({
                "at": current,
                "candidates": sorted(candidates),
                "note": (
                    "More than one branch under this node is supported. Sister "
                    "haplogroups are mutually exclusive, so this points at a "
                    "wrong defining position in the bundled tree or a "
                    "heteroplasmic or miscalled position, not an unusual "
                    "lineage."
                ),
            })
            if len(top) > 1:
                break
            current = best
            chain.append(current)
            continue
        current = candidates[0]
        chain.append(current)

    resolved = current != "root"
    path_rows: list[dict] = []
    for name in chain:
        if name == "root":
            continue
        path_rows.extend(node_support[name]["rows"])
    verified_used = sum(1 for r in path_rows if r["verified"] and r["state"] == DERIVED)
    unverified_used = sum(1 for r in path_rows
                          if not r["verified"] and r["state"] == DERIVED)

    stopped_at = current
    blocked_by_array = sorted(
        k for k in children.get(stopped_at, [])
        if node_support.get(k, {}).get("testable", 0) == 0
    )
    excluded = sorted(
        k for k in children.get(stopped_at, [])
        if node_support.get(k, {}).get("contradict", 0) > 0
    )

    caveats = [_TRI_STATE_NOTE]
    if unverified_used:
        caveats.append(
            f"{unverified_used} of the supporting positions on this path are "
            f"unverified. The call is provisional until the data builder "
            f"confirms them against PhyloTree."
        )
    if assumed:
        caveats.append(_ASSUMED_NOTE + " Assumed: " + ", ".join(assumed) + ".")
    if blocked_by_array:
        caveats.append(
            "Branches below this point were not tested because your array does "
            "not read their positions: " + ", ".join(blocked_by_array) + "."
        )

    return _backbone_payload(
        "MT", MT_BACKBONE,
        resolved=resolved, stopped_at=stopped_at, chain=chain,
        path_markers=path_rows, assumed=assumed, buckets=buckets,
        conflicts=conflicts, blocked_by_array=blocked_by_array,
        excluded=excluded, verified_used=verified_used,
        unverified_used=unverified_used, caveats=caveats,
    )


# ---------------------------------------------------------------------------
# Resolution ceiling. The differentiator, computed, never hardcoded.
# ---------------------------------------------------------------------------

def _y_marker_readable(entry: dict, genotypes: dict, verified_only: bool) -> bool:
    state, _ = marker_state(entry, genotypes, verified_only=verified_only)
    return state not in NOT_TESTABLE_STATES


def _mt_node_readable(entry: dict, genotypes: dict, verified_only: bool) -> bool:
    for variant in entry.get("defining", []):
        state, _ = _mt_variant_state(variant, genotypes, verified_only=verified_only)
        if state not in NOT_TESTABLE_STATES:
            return True
    return False


def _deepest_reachable(tree: dict, readable: dict[str, bool]) -> str:
    """Deepest node this array could ever resolve to.

    A node counts when its OWN marker is readable, not when its whole ancestor
    chain is, because the walk is allowed to pass through an untested
    intermediate node and record it as assumed. So the ceiling is the deepest
    readable marker anywhere in the tree: the answer to "how deep could this
    data ever go", which is a different and more useful question than "how deep
    did it go this time".
    """
    best = "root"
    best_depth = 0
    for node, is_readable in sorted(readable.items()):
        if not is_readable:
            continue
        depth = _depth(node, tree)
        if depth > best_depth:
            best, best_depth = node, depth
    return best


def resolution_ceiling(genotypes: dict, system: str = "Y", *,
                       sex_hint: Any = None,
                       verified_only: bool = False,
                       array_positions: int | None = None) -> dict:
    """Where this person's data runs out, in numbers and in one plain sentence.

    ``markers_available`` and ``markers_in_tree`` are counted from the supplied
    genotype map and the bundled tree. Nothing here is a stored constant, which
    matters: the whole point of the field is that it describes THIS array and
    THIS person rather than a marketing average.

    ``array_positions`` is the total number of positions the array reads on the
    relevant chromosome, when the caller knows it. Supplying it lets the
    mitochondrial sentence state what fraction of the 16,569 base genome was
    read. Without it the fraction is reported as None rather than guessed.
    """
    genotypes = genotypes or {}
    system = str(system or "Y").strip().upper()
    system = "MT" if system in ("MT", "M", "MTDNA", "MITO") else "Y"

    if system == "Y":
        tree = Y_BACKBONE
        readable = {
            name: _y_marker_readable(entry, genotypes, verified_only)
            for name, entry in tree.items() if name != "root"
        }
        call = call_y_backbone(genotypes, sex_hint, verified_only=verified_only)
        unit = "Y-SNPs"
        upgrade = (
            "A Big Y-700 test sequences the whole Y chromosome and would "
            "resolve further."
        )
    else:
        tree = MT_BACKBONE
        readable = {
            name: _mt_node_readable(entry, genotypes, verified_only)
            for name, entry in tree.items() if name != "root"
        }
        call = call_mt_backbone(genotypes, verified_only=verified_only)
        unit = "mitochondrial backbone positions"
        upgrade = (
            f"A full mitochondrial sequence reads all {MT_GENOME_BP:,} bases "
            f"and would resolve further."
        )

    markers_in_tree = len(readable)
    markers_available = sum(1 for value in readable.values() if value)
    coverage = round(markers_available / markers_in_tree, 4) if markers_in_tree else 0.0
    ceiling_node = _deepest_reachable(tree, readable)
    resolved = call.get("haplogroup")

    genome_fraction = None
    if system == "MT" and array_positions:
        genome_fraction = round(float(array_positions) / MT_GENOME_BP, 4)

    typical_key = "y_snps" if system == "Y" else "mt_positions"
    typical = TYPICAL_ARRAY[typical_key]
    typical_fraction = None
    if system == "MT":
        typical_fraction = round(float(typical) / MT_GENOME_BP, 4)

    if call.get("state") == "no_y_data":
        sentence = NO_Y_MESSAGE
    elif resolved:
        head = (
            f"Your array reads {markers_available:,} of the {markers_in_tree:,} "
            f"{unit} in the {TREE_NAME} tree, version {TREE_VERSION}, and "
            f"resolves you to {resolved}."
        )
        if genome_fraction is not None:
            head += (
                f" That is {round(genome_fraction * 100, 1)} percent of the "
                f"{MT_GENOME_BP:,} base mitochondrial genome."
            )
        deeper = ""
        if ceiling_node != resolved and ceiling_node != "root":
            deeper = (
                f" The deepest this array could reach on this tree is "
                f"{ceiling_node}."
            )
        sentence = head + deeper + " " + upgrade
    else:
        sentence = (
            f"Your array reads {markers_available:,} of the {markers_in_tree:,} "
            f"{unit} in the {TREE_NAME} tree, version {TREE_VERSION}, and does "
            f"not resolve to any node below the root. " + upgrade
        )

    out = {
        "system": system,
        "markers_available": markers_available,
        "markers_in_tree": markers_in_tree,
        "markers_missing": markers_in_tree - markers_available,
        "coverage": coverage,
        "deepest_resolvable_node": resolved,
        "deepest_resolvable_label": call.get("label"),
        "ceiling_node": None if ceiling_node == "root" else ceiling_node,
        "ceiling_depth": _depth(ceiling_node, tree),
        "resolved_depth": _depth(resolved, tree) if resolved else 0,
        "array_positions": array_positions,
        "mt_genome_bp": MT_GENOME_BP if system == "MT" else None,
        "mt_genome_fraction": genome_fraction,
        "typical_array": {
            "name": TYPICAL_ARRAY["name"],
            "count": typical,
            "mt_genome_fraction": typical_fraction,
            "note": (
                "Approximate published figure for a consumer array, used only "
                "for comparison. Your own numbers above are counted from your "
                "own file."
            ),
        },
        "sentence": sentence,
        "upgrade_path": upgrade,
    }
    out.update(tree_stamp())
    return out


# ---------------------------------------------------------------------------
# Input writers for the external tools
# ---------------------------------------------------------------------------

def _workdir(workdir: str | Path | None) -> Path:
    if workdir is not None:
        path = Path(workdir)
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix="dnainsight_hap_"))


def write_y_array_input(genotypes: dict, path: str | Path | None = None, *,
                        workdir: str | Path | None = None) -> str:
    """Write the Y markers as a 23andMe-style raw array file.

    Yleaf was chosen over every alternative because it accepts SNP-array input
    directly, which is the only format a consumer array user actually has. The
    file below is that format: a comment header then rsid, chromosome, position
    and genotype, tab separated.

    Position is written as 0 where the bundled tree does not record one, since
    the backbone is a marker table and not a coordinate table. The wiring pass
    must supply real coordinates before this file is handed to a tool that
    cares about them.
    """
    target = Path(path) if path is not None else _workdir(workdir) / "y_markers.txt"
    lines = [
        "# DNAInsight Y marker export",
        f"# tree: {TREE_NAME} {TREE_VERSION}",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for name in sorted(Y_BACKBONE):
        entry = Y_BACKBONE[name]
        if name == "root" or not entry.get("marker"):
            continue
        present, alleles = _lookup(genotypes or {}, marker_keys(entry))
        if not present or not alleles:
            continue
        ordered = sorted(alleles)
        call = ordered[0] * 2 if len(ordered) == 1 else "".join(ordered[:2])
        rsid = entry.get("rsid") or entry.get("marker")
        lines.append(f"{rsid}\tY\t0\t{call}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(target)


def write_hsd(genotypes: dict, sample_id: str = "sample",
              path: str | Path | None = None, *,
              workdir: str | Path | None = None,
              range_text: str = "") -> str:
    """Write HaploGrep hsd input for one sample.

    The hsd format is deliberately simple: a header, then one row per sample
    holding the sample id, the covered range, a haplogroup placeholder and the
    list of differences from rCRS. Only positions whose observed base differs
    from the recorded rCRS base are listed, because that is what hsd means by a
    polymorphism.

    Positions the array did not read are simply absent from the row. HaploGrep
    treats an absent position as uncovered rather than as reference, which is
    the same not-testable distinction this module enforces everywhere else, so
    the range is written from the positions actually read rather than as the
    whole genome when the caller does not override it.
    """
    genotypes = genotypes or {}
    target = Path(path) if path is not None else _workdir(workdir) / f"{sample_id}.hsd"

    seen: dict[int, str] = {}
    covered: list[int] = []
    for name in sorted(MT_BACKBONE):
        entry = MT_BACKBONE[name]
        if name == "root":
            continue
        for variant in entry.get("defining", []):
            position = variant["position"]
            if position in seen or position in covered:
                continue
            present, alleles = _lookup(genotypes, _mt_keys(variant))
            if not present or not alleles:
                continue
            covered.append(position)
            observed = sorted(alleles)[0]
            if observed != variant["rcrs"]:
                seen[position] = observed

    if range_text:
        span = range_text
    elif covered:
        span = f"{min(covered)}-{max(covered)}"
    else:
        span = f"1-{MT_GENOME_BP}"

    polymorphisms = [f"{pos}{seen[pos]}" for pos in sorted(seen)]
    header = "SampleId\tRange\tHaplogroup\tPolymorphisms"
    row = "\t".join([str(sample_id), span, "?", *polymorphisms])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(header + "\n" + row + "\n", encoding="utf-8")
    return str(target)


# ---------------------------------------------------------------------------
# Output parsers.
#
# Every one of these is tolerant on purpose. The exact column layout of Yleaf,
# HaploGrep and Clade Finder varies between releases, and an adapter that
# crashes on an unexpected column is worse than one that reports it could not
# read the answer.
# ---------------------------------------------------------------------------

def _parse_table(text: str) -> list[dict]:
    """Split a whitespace or tab delimited table into dicts keyed by header."""
    rows: list[dict] = []
    lines = [ln for ln in str(text or "").splitlines() if ln.strip()]
    lines = [ln for ln in lines if not ln.lstrip().startswith("#")]
    if not lines:
        return rows
    header = [c.strip().lower() for c in _split_cells(lines[0])]
    for line in lines[1:]:
        cells = _split_cells(line)
        if not cells:
            continue
        rows.append({header[i]: cells[i] for i in range(min(len(header), len(cells)))})
    return rows


def _split_cells(line: str) -> list[str]:
    if "\t" in line:
        return [c.strip() for c in line.split("\t")]
    return line.split()


def _pick(row: dict, *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value:
            return str(value).strip()
    return ""


def parse_yleaf_output(text: str) -> dict:
    """Pull the haplogroup out of a Yleaf prediction table."""
    rows = _parse_table(text)
    for row in rows:
        call = _pick(row, "hg", "haplogroup", "predicted_haplogroup", "hg_prediction")
        if call:
            return {
                "haplogroup": call,
                "sample": _pick(row, "sample_name", "sample", "sampleid", "id"),
                "quality": _pick(row, "hg_marker", "qc-score", "score", "quality") or None,
                "rows": rows,
                "parsed": True,
            }
    return {"haplogroup": None, "sample": "", "quality": None,
            "rows": rows, "parsed": False}


def parse_haplogrep_output(text: str) -> dict:
    """Pull the haplogroup and quality score out of a HaploGrep classify table."""
    rows = _parse_table(text)
    for row in rows:
        call = _pick(row, "haplogroup", "hg")
        if call and call != "?":
            return {
                "haplogroup": call,
                "sample": _pick(row, "sampleid", "sample_id", "sample"),
                "quality": _pick(row, "quality", "rank") or None,
                "rows": rows,
                "parsed": True,
            }
    return {"haplogroup": None, "sample": "", "quality": None,
            "rows": rows, "parsed": False}


def parse_cladefinder_output(text: str) -> dict:
    """Pull the clade out of a Clade Finder result table."""
    rows = _parse_table(text)
    for row in rows:
        call = _pick(row, "clade", "haplogroup", "hg", "result")
        if call:
            return {
                "haplogroup": call,
                "sample": _pick(row, "sample", "sampleid", "id"),
                "quality": _pick(row, "score", "confidence") or None,
                "rows": rows,
                "parsed": True,
            }
    return {"haplogroup": None, "sample": "", "quality": None,
            "rows": rows, "parsed": False}


# ---------------------------------------------------------------------------
# Agreement between two independent Y callers
# ---------------------------------------------------------------------------

def _is_ancestor(candidate: str, node: str) -> bool:
    chain = path_to(node, Y_BACKBONE)
    return bool(candidate) and candidate in chain[:-1]


def compare_y_calls(call_a: dict, call_b: dict) -> dict:
    """Compare two independent Y calls without picking a winner.

    Three outcomes:

      agree        identical labels
      consistent   one is an ancestor of the other in the bundled tree, so a
                   shallow caller and a deep caller are saying the same thing
                   at different depths, which is not a disagreement
      conflict     genuinely different branches

    On a conflict both calls are retained side by side and no resolution is
    offered. That is the same rule ``merge.py`` applies to two files that
    disagree at a position: the disagreement is the finding.
    """
    name_a = (call_a or {}).get("haplogroup")
    name_b = (call_b or {}).get("haplogroup")
    tool_a = (call_a or {}).get("source") or (call_a or {}).get("tool") or "a"
    tool_b = (call_b or {}).get("source") or (call_b or {}).get("tool") or "b"
    calls = [
        {"tool": tool_a, "haplogroup": name_a,
         "tree_name": (call_a or {}).get("tree_name"),
         "tree_version": (call_a or {}).get("tree_version")},
        {"tool": tool_b, "haplogroup": name_b,
         "tree_name": (call_b or {}).get("tree_name"),
         "tree_version": (call_b or {}).get("tree_version")},
    ]

    if not name_a or not name_b:
        return {
            "comparable": False, "agree": None, "conflict": False,
            "relation": "incomparable", "calls": calls,
            "note": "One of the two callers produced no haplogroup, so there is "
                    "nothing to compare. This is not agreement.",
        }
    if name_a == name_b:
        return {"comparable": True, "agree": True, "conflict": False,
                "relation": "identical", "calls": calls,
                "note": "Both callers returned the same haplogroup."}
    if _is_ancestor(name_a, name_b) or _is_ancestor(name_b, name_a):
        deeper = name_b if _is_ancestor(name_a, name_b) else name_a
        return {
            "comparable": True, "agree": True, "conflict": False,
            "relation": "consistent_depth", "deeper": deeper, "calls": calls,
            "note": (
                "One call sits above the other on the same branch, so they do "
                "not disagree. The deeper call is simply a finer reading of the "
                "same lineage."
            ),
        }
    return {
        "comparable": True, "agree": False, "conflict": True,
        "relation": "divergent", "calls": calls,
        "resolution": None,
        "note": (
            "The two callers place this sample on different branches. Both "
            "calls are kept and neither is chosen. A disagreement between two "
            "independent callers is information about reliability, and the "
            "usual cause is a different tree revision or a marker one caller "
            "could not read, not an exotic lineage. Compare the tree versions "
            "on each call first."
        ),
    }


# ---------------------------------------------------------------------------
# External adapters.
#
# THE SUBPROCESS BOUNDARY IS THE LICENCE BOUNDARY. Yleaf is GPL-3.0. Nothing
# here imports it, links it, or copies a line of it. `external.run` starts a
# separate process running a program the user installed into their own home
# directory on explicit consent, which is why this MIT tree stays MIT. Every
# adapter degrades through `external.guard` instead of raising, so a user with
# no tools installed still gets the bundled backbone call.
# ---------------------------------------------------------------------------

# The exact flags each tool accepts have changed between releases. They are
# collected here so the wiring pass has one place to correct them rather than
# hunting through call sites.
YLEAF_ARGS = ("-r", "{input}", "-o", "{outdir}", "-rg", "hg19")
HAPLOGREP_ARGS = ("classify", "--in", "{input}", "--out", "{output}",
                  "--format", "hsd")
CLADEFINDER_ARGS = ("--input", "{input}", "--output", "{output}")


def _format_args(template: Sequence[str], **values: str) -> list[str]:
    return [str(part).format(**values) for part in template]


def _degraded(base: dict, blocked: dict, tool_name: str, upgrade: str) -> dict:
    """Attach the standard degraded payload to a bundled backbone call."""
    base = dict(base)
    base["source"] = "bundled_backbone"
    base["external"] = blocked
    base["external_available"] = False
    note = (
        f"This call came from the bundled {TREE_NAME} tree, version "
        f"{TREE_VERSION}, which is a backbone and nothing more. {upgrade} "
        f"{blocked.get('reason', '')}".strip()
    )
    base["note"] = note
    base.setdefault("caveats", []).append(note)
    base["deeper_call_requires"] = tool_name
    return base


def call_y(genotypes: dict, *, sex_hint: Any = None,
           workdir: str | Path | None = None,
           verified_only: bool = False) -> dict:
    """Y haplogroup call, refined by Yleaf when it is installed.

    With Yleaf absent this returns the bundled backbone call with ``source``
    set to "bundled_backbone" and an explicit note that a deeper call needs
    Yleaf. It never raises and never returns an empty result, because "we could
    not run the good tool" and "you have no Y haplogroup" are completely
    different statements and the user must be able to tell them apart.
    """
    fallback = call_y_backbone(genotypes, sex_hint, verified_only=verified_only)
    upgrade = (
        "Yleaf reads the whole Y-SNP panel your file contains against a current "
        "tree and would resolve further."
    )

    if fallback.get("state") == "no_y_data":
        # Checked BEFORE the tool gate on purpose. A woman must be told there
        # is no Y chromosome to call, not that Yleaf is missing. Reporting a
        # tool problem here would send her installing software that would tell
        # her exactly the same thing.
        fallback["note"] = NO_Y_MESSAGE
        return fallback

    blocked = external.guard("yleaf", "haplogroup_y")
    if blocked is not None:
        return _degraded(fallback, blocked, "Yleaf", upgrade)

    try:
        work = _workdir(workdir)
        infile = write_y_array_input(genotypes, workdir=work)
        outdir = work / "yleaf_out"
        outdir.mkdir(parents=True, exist_ok=True)
        completed = external.run(
            "yleaf",
            _format_args(YLEAF_ARGS, input=infile, outdir=str(outdir)),
        )
        parsed = parse_yleaf_output(_collect_output(completed, outdir))
    except external.ExternalError as exc:
        # A tool that is installed but fails is still a degraded result, not a
        # crash. The user gets the backbone answer and the reason.
        failure = {
            "available": False, "capability": "haplogroup_y", "tool": "Yleaf",
            "tool_id": "yleaf", "state": "failed", "reason": str(exc),
            "not_attempted": False, "results": [],
        }
        return _degraded(fallback, failure, "Yleaf", upgrade)

    if not parsed.get("haplogroup"):
        failure = {
            "available": False, "capability": "haplogroup_y", "tool": "Yleaf",
            "tool_id": "yleaf", "state": "unparsed",
            "reason": "Yleaf ran but its output table could not be read.",
            "not_attempted": False, "results": [],
        }
        return _degraded(fallback, failure, "Yleaf", upgrade)

    payload = {
        "system": "Y",
        "available": True,
        "state": "called",
        "haplogroup": parsed["haplogroup"],
        "label": parsed["haplogroup"],
        "source": "yleaf",
        "external_available": True,
        "quality": parsed.get("quality"),
        "confidence": "external",
        "backbone": fallback,
        "caveats": [
            _TRI_STATE_NOTE,
            "This call came from Yleaf against the tree Yleaf was installed "
            "with. Record that tree revision: the same sample resolves to "
            "different labels under different revisions.",
        ],
        "tree_name": "Yleaf tree",
        "tree_version": "as installed, not recorded by DNAInsight",
        "bundled_tree": tree_stamp(),
    }
    return payload


def _collect_output(completed: Any, outdir: Path) -> str:
    """Prefer a written result table, fall back to stdout.

    Yleaf writes its prediction to a file in the output directory in most
    releases and to stdout in others, so both are accepted rather than pinning
    a filename that changes.
    """
    try:
        candidates = sorted(p for p in Path(outdir).rglob("*")
                            if p.is_file() and p.suffix in (".txt", ".tsv", ".hg"))
    except OSError:
        candidates = []
    for candidate in candidates:
        try:
            text = candidate.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if text.strip():
            return text
    return getattr(completed, "stdout", "") or ""


def call_mt(genotypes: dict, *, sample_id: str = "sample",
            workdir: str | Path | None = None,
            verified_only: bool = False) -> dict:
    """mtDNA haplogroup call, refined by HaploGrep 3 when it is installed.

    HaploGrep takes VCF or hsd. hsd is written here because it needs nothing
    but the positions the array read, and because building a VCF would mean
    inventing reference bases for positions the bundled tree does not record.
    """
    fallback = call_mt_backbone(genotypes, verified_only=verified_only)
    upgrade = (
        "HaploGrep 3 classifies against a full PhyloTree catalogue and would "
        "resolve further."
    )

    blocked = external.guard("haplogrep", "haplogroup_mt")
    if blocked is not None:
        return _degraded(fallback, blocked, "HaploGrep 3", upgrade)

    try:
        work = _workdir(workdir)
        infile = write_hsd(genotypes, sample_id, workdir=work)
        outfile = work / f"{sample_id}.haplogrep.txt"
        completed = external.run(
            "haplogrep",
            _format_args(HAPLOGREP_ARGS, input=infile, output=str(outfile)),
        )
        text = ""
        if outfile.exists():
            text = outfile.read_text(encoding="utf-8", errors="replace")
        if not text.strip():
            text = getattr(completed, "stdout", "") or ""
        parsed = parse_haplogrep_output(text)
    except external.ExternalError as exc:
        failure = {
            "available": False, "capability": "haplogroup_mt",
            "tool": "HaploGrep 3", "tool_id": "haplogrep", "state": "failed",
            "reason": str(exc), "not_attempted": False, "results": [],
        }
        return _degraded(fallback, failure, "HaploGrep 3", upgrade)

    if not parsed.get("haplogroup"):
        failure = {
            "available": False, "capability": "haplogroup_mt",
            "tool": "HaploGrep 3", "tool_id": "haplogrep", "state": "unparsed",
            "reason": "HaploGrep ran but its output table could not be read.",
            "not_attempted": False, "results": [],
        }
        return _degraded(fallback, failure, "HaploGrep 3", upgrade)

    return {
        "system": "MT",
        "available": True,
        "state": "called",
        "haplogroup": parsed["haplogroup"],
        "label": parsed["haplogroup"],
        "source": "haplogrep",
        "external_available": True,
        "quality": parsed.get("quality"),
        "confidence": "external",
        "backbone": fallback,
        "caveats": [
            _TRI_STATE_NOTE,
            "HaploGrep's tree catalogue is installed separately from its "
            "binary. Record which catalogue produced this call; a haplogroup "
            "without its tree revision is not a reproducible result.",
        ],
        "tree_name": "HaploGrep PhyloTree catalogue",
        "tree_version": "as installed, not recorded by DNAInsight",
        "bundled_tree": tree_stamp(),
    }


def second_opinion_y(genotypes: dict, *, primary: dict | None = None,
                     sex_hint: Any = None,
                     workdir: str | Path | None = None) -> dict:
    """Independent Y call from Clade Finder, compared against the primary call.

    When Yleaf and Clade Finder both run and DISAGREE the disagreement is
    surfaced as a conflict and no winner is chosen. Two independent callers
    landing on different branches tells the user something real about how solid
    the call is, and quietly preferring one would throw that away.
    """
    primary = primary if primary is not None else call_y(genotypes, sex_hint=sex_hint)

    blocked = external.guard("cladefinder", "haplogroup_y_second_opinion")
    if blocked is not None:
        return {
            "system": "Y",
            "available": False,
            "state": "not_attempted",
            "primary": primary,
            "second_opinion": None,
            "comparison": None,
            "conflict": False,
            "external": blocked,
            "note": (
                "No second opinion was taken. Clade Finder is not available, so "
                "the primary call has not been independently checked. That is "
                "different from two callers agreeing."
            ),
            **tree_stamp(),
        }

    try:
        work = _workdir(workdir)
        infile = write_y_array_input(genotypes, workdir=work)
        outfile = work / "cladefinder.txt"
        completed = external.run(
            "cladefinder",
            _format_args(CLADEFINDER_ARGS, input=infile, output=str(outfile)),
        )
        text = ""
        if outfile.exists():
            text = outfile.read_text(encoding="utf-8", errors="replace")
        if not text.strip():
            text = getattr(completed, "stdout", "") or ""
        parsed = parse_cladefinder_output(text)
    except external.ExternalError as exc:
        return {
            "system": "Y", "available": False, "state": "failed",
            "primary": primary, "second_opinion": None, "comparison": None,
            "conflict": False,
            "external": {
                "available": False, "tool": "Clade Finder",
                "tool_id": "cladefinder", "state": "failed", "reason": str(exc),
                "not_attempted": False, "results": [],
            },
            "note": "Clade Finder is installed but failed, so no second "
                    "opinion was obtained.",
            **tree_stamp(),
        }

    second = {
        "haplogroup": parsed.get("haplogroup"),
        "source": "cladefinder",
        "quality": parsed.get("quality"),
        "tree_name": "Clade Finder tree",
        "tree_version": "as installed, not recorded by DNAInsight",
    }
    comparison = compare_y_calls(primary, second)
    return {
        "system": "Y",
        "available": bool(parsed.get("haplogroup")),
        "state": "compared",
        "primary": primary,
        "second_opinion": second,
        "comparison": comparison,
        "conflict": bool(comparison.get("conflict")),
        "external": {"available": True, "tool": "Clade Finder",
                     "tool_id": "cladefinder", "state": "ready"},
        "note": comparison.get("note", ""),
        **tree_stamp(),
    }


# ---------------------------------------------------------------------------
# Convenience entry point for the pipeline
# ---------------------------------------------------------------------------

def analyse(genotypes: dict, *, sex_hint: Any = None,
            mt_genotypes: dict | None = None,
            array_mt_positions: int | None = None,
            verified_only: bool = False) -> dict:
    """Both uniparental calls plus both ceilings, in one payload.

    ``mt_genotypes`` is a position-keyed mitochondrial view; see
    :func:`mt_positions_from_merged`. When it is omitted the rsID-keyed map is
    used for both systems, which works but reads far fewer mitochondrial
    positions, and the ceiling will say so.
    """
    mt_source = mt_genotypes if mt_genotypes is not None else genotypes
    y_call = call_y(genotypes, sex_hint=sex_hint, verified_only=verified_only)
    mt_call = call_mt(mt_source, verified_only=verified_only)
    return {
        "y": y_call,
        "mt": mt_call,
        "y_ceiling": resolution_ceiling(genotypes, "Y", sex_hint=sex_hint,
                                        verified_only=verified_only),
        "mt_ceiling": resolution_ceiling(mt_source, "MT",
                                         verified_only=verified_only,
                                         array_positions=array_mt_positions),
        "second_opinion": second_opinion_y(genotypes, primary=y_call,
                                           sex_hint=sex_hint),
        "unverified": unverified_markers(),
        **tree_stamp(),
    }
