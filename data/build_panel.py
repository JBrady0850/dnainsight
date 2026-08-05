"""
build_panel.py -- reference panel builder for ancestry, phasing and imputation.

WHY THIS FILE EXISTS
--------------------
``backend/external.py`` declares the panel ``onekg_sgdp`` and states what it may
legally contain. Nothing built it. ``backend/ancestry.py`` and
``backend/imputation.py`` both read ``external.panel_root()/onekg_sgdp/`` and
degrade to a documented "not built" payload when the files are absent, so the
application works without this script and gets meaningfully better with it.

This is the script that produces those files, from sources whose terms are
actually known, with every refusal written down rather than left as an absence.

SOURCES AND LICENCES
--------------------
    1000 Genomes phase 3, via IGSR at EMBL-EBI
        Open with no restriction on use, per the project's data reuse
        statement. Distributed under EMBL-EBI terms; citation requested.
        Genome build GRCh37. Retrieval URLs are listed in ONEKG_VCF_TEMPLATE
        and ONEKG_PANEL_URL and were verified to return HTTP 200 on 2026-08-04.

        TRAP, VERIFIED AGAINST THE LIVE FILES ON 2026-08-04. The IGSR release
        genotype VCFs for the 20130502 release publish the ID column as "."
        for every record. All 1,103,547 chr22 records and the first million
        records of the genome-wide sites file were checked and not one carries
        an rsID. A builder that assumes otherwise silently keeps nothing,
        because --array-file joins on rsID. The Browning lab republishes the
        SAME phase 3 v5a callset, GT-only, WITH the ID column populated, as
        the standard Beagle reference panel, and that is what
        ONEKG_BEAGLE_TEMPLATE points at and what --onekg-source selects by
        default. It is also about 35 percent smaller because it carries no
        INFO block. The IGSR URLs stay recorded and selectable so the
        canonical location remains visible, exactly as build_full_reference.py
        keeps the retired GWAS endpoint in front of its live mirror.

    Simons Genome Diversity Project, PUBLIC tier only, 279 samples
        Unrestricted. The public tier is the part of SGDP released without a
        data access agreement. Genome build GRCh37.

    PLINK-format genetic maps for GRCh37, Browning lab mirror
        Public, redistributable, used only to interpolate cM for panel.map.

Only these three feed the panel. Everything else on the shortlist was refused,
and the refusals are enforced in code below rather than left to good intentions.

HARD EXCLUSIONS
---------------
Three, each with a reason the builder prints before it does anything:

  1. THE 21-SAMPLE RESTRICTED SGDP TIER. Its signed agreement includes "I will
     not use the data for any commercial purposes". DNAInsight is MIT and makes
     a commercial-use grant to every downstream user, so a panel built from that
     tier could not be shipped, described or reproduced under this project's
     terms. There is no flag that turns it on. See REFUSAL_SGDP_RESTRICTED.

  2. HGDP, UNLESS BOTH ``--include-hgdp`` AND ``--accept-consent-caveat``.
     Two flags, not one, and the second is not a licence acceptance. HGDP is
     legally open under the Fort Lauderdale Principles; the objection is a
     CONSENT objection. Nature Genetics, 24 November 2025, concluded that broad
     reuse may diverge from what participants consented to, and the PRIMED
     Consortium voted on 21 August 2024 to keep permitting its use while
     acknowledging "failure to obtain informed consent consistent with current
     standards from many participants". Accepting a licence does not answer a
     consent objection, which is exactly why ``--accept-terms`` alone is not
     enough. See REFUSAL_HGDP and :func:`hgdp_gate`.

  3. THE ALLEN ANCIENT DNA RESOURCE, ENTIRELY. Its terms were not readable at
     review time and the compendium aggregates datasets that each carry their
     own upstream terms. Unread is not permissive. There is no flag. See
     REFUSAL_AADR.

ANCESTRY-INFORMATIVE MARKER SELECTION, AND WHAT IT IS NOT
---------------------------------------------------------
``informative_markers.tsv`` needs a defensible per-population marker ranking.
Two statistics are implemented in pure Python and either can be selected with
``--statistic``. Both are computed from per-population allele counts read
straight out of the VCF, so there is nothing to trust but arithmetic.

    fst (default)
        Wright's fixation index in its ratio-of-heterozygosities form:

            Hs = sum_i w_i * 2 * p_i * (1 - p_i)      within-population
            Ht = 2 * pbar * (1 - pbar)                total
            Fst = (Ht - Hs) / Ht,  and 0 when Ht == 0

        where p_i is the alternate-allele frequency in population i, w_i is
        that population's share of the called allele copies, and pbar is the
        copy-weighted mean frequency. Per-population ranking uses the same
        function on a two-group partition, that population against every other
        population pooled, which is the usual "one versus rest" AIM criterion.

    informativeness
        Rosenberg et al. (2003) informativeness for assignment, In:

            In = sum_j [ -pbar_j ln pbar_j + (1/K) sum_i p_ij ln p_ij ]

        over alleles j and the K populations i, with pbar_j the UNWEIGHTED mean
        across populations. Measured in nats. It answers "how much does knowing
        this genotype tell me about which population this person came from",
        which is closer to the actual question than Fst is.

LIMITATIONS OF BOTH, STATED PLAINLY
    * Neither carries the Weir and Cockerham (1984) small-sample correction.
      With unequal or small population samples both are biased upward. The
      1000 Genomes populations run from 61 to 113 samples, which is small
      enough for that to matter at the margin, so a marker near the cut is not
      meaningfully better than the one just below it.
    * Both are computed per marker, independently. Linkage disequilibrium is
      not accounted for, so the top of the ranking will contain correlated
      markers from the same haplotype block and the effective information is
      lower than the marker count suggests.
    * Both treat every site as biallelic. Multiallelic records are skipped
      rather than collapsed, because collapsing them would silently change what
      the frequency means.
    * Fst is a property of the panel's population labels. Relabelling the
      populations changes every number here. A marker that separates CEU from
      YRI says nothing about separating two European populations.
    * The whole ranking inherits the panel's ancestry bias. Under-represented
      ancestries get fewer informative markers because there are fewer samples
      to estimate their frequencies from, not because fewer markers exist.

These limitations are why the AIM set is written as a LABELLING of the panel
rather than as a filter on it. Every marker that survives the array and MAF
filters goes into panel.vcf.gz and panel.map; informative_markers.tsv names the
subset that ranked highest per population. Nothing is discarded on the strength
of a statistic this rough.

MEMORY DISCIPLINE
-----------------
A single 1000 Genomes phase 3 chromosome is over a gigabyte compressed and
several times that decompressed. Nothing is ever read whole. Every VCF is
streamed through ``gzip.GzipFile`` on the socket, one line at a time, and the
per-population counts for a record are discarded as soon as the record is
written. The marker ranking uses a bounded heap (:class:`TopN`), so peak memory
is proportional to the number of markers kept, not to the number seen.

CONSENT AND OFFLINE CONTRACT
----------------------------
Nothing here runs at import time and nothing here is reachable from the running
application. ``python data/build_panel.py --help`` works with no network.
Without ``--accept-terms`` the script performs a dry run: it prints every URL,
every licence and the estimated download size, and fetches nothing. This is the
same gate ``backend/snpedia.py`` applies to its harvest, which returns HTTP 403
until ``accept_license`` is true, and reusing the shape means there is one rule
to understand rather than two.

Output goes to ``backend.external.panel_root()``, which is ``~/.dnainsight/
panels/``, deliberately OUTSIDE the repository tree. Panels are tens of
gigabytes of source data reduced to hundreds of megabytes of artefact; neither
belongs in git, and the SNPedia precedent already established that non-bundled
data lives in the user's home directory.

USAGE
-----
    python data/build_panel.py                       # dry run, the default
    python data/build_panel.py --accept-terms --array-file uploads/mine.txt
    python data/build_panel.py --accept-terms --chromosomes 21,22 --limit 50000
    python data/build_panel.py --accept-terms --include-hgdp --accept-consent-caveat
    python data/build_panel.py --statistic informativeness --dry-run
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import heapq
import math
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend import external                                      # noqa: E402
# Reused rather than reimplemented, deliberately. build_full_reference.py
# already decided what counts as an array file and what counts as an rsID in
# one; two answers to that question would drift inside a release, and the whole
# point of --array-file is that both builders restrict to the SAME positions.
from data.build_full_reference import (                           # noqa: E402
    load_array_rsids,
    resolve_columns,
    stream_gzip_lines,
    stream_text_lines,
)

try:  # pragma: no cover - exercised only when requests is absent
    import requests
except ImportError:  # pragma: no cover
    requests = None  # type: ignore[assignment]

__all__ = [
    "PANEL_ID", "TARGET_BUILD", "PANEL_OUTPUTS",
    "ONEKG_PANEL_URL", "ONEKG_VCF_TEMPLATE", "ONEKG_X_VCF_URL",
    "ONEKG_RELEASE", "ONEKG_VCF_BYTES", "ONEKG_HOMEPAGE",
    "ONEKG_BEAGLE_TEMPLATE", "ONEKG_BEAGLE_BYTES", "ONEKG_SOURCES",
    "SGDP_BASE", "SGDP_METADATA_URL", "SGDP_VCF_TEMPLATE", "SGDP_ENA_PROJECT",
    "SGDP_PUBLIC_SAMPLES", "SGDP_RESTRICTED_SAMPLES",
    "GENETIC_MAP_URL", "HGDP_BASE", "HGDP_REUSE_STATEMENT_URL", "AADR_URL",
    "LICENCES", "REFUSAL_SGDP_RESTRICTED", "REFUSAL_HGDP", "REFUSAL_AADR",
    "REFUSALS", "STATISTICS", "AUTOSOMES",
    "PanelRefused", "TopN",
    "build_parser", "hgdp_gate", "refuse_restricted_sgdp", "refuse_aadr",
    "allele_frequency", "wright_fst", "population_fst", "informativeness",
    "marker_statistic", "parse_panel_file", "parse_sgdp_metadata",
    "vcf_sample_names", "parse_vcf_record", "allele_counts",
    "interpolate_cm", "sha256_of", "download_plan", "print_plan",
    "write_populations_tsv", "write_q_columns", "write_informative_markers",
    "write_build_txt", "panel_dir", "main",
]


# ---------------------------------------------------------------------------
# Identity and output contract
#
# The file names are not free choices. backend/ancestry.py and
# backend/imputation.py open these exact paths, and external.PANELS declares
# the first three as the completeness test for panel_status(). Renaming one
# here silently turns a built panel into a "partial" one.
# ---------------------------------------------------------------------------

PANEL_ID = "onekg_sgdp"
TARGET_BUILD = "GRCh37"

PANEL_OUTPUTS: tuple[str, ...] = (
    "panel.vcf.gz",
    "panel.map",
    "populations.tsv",
    "informative_markers.tsv",
    "q_columns.tsv",
    "BUILD.txt",
)


# ---------------------------------------------------------------------------
# Source endpoints
#
# Every URL below carries the date it was checked, for the same reason
# backend/external.py dates its licence strings: an endpoint nobody re-reads
# drifts, and the GWAS Catalog in build_full_reference.py is the standing
# example of a documented URL that quietly started returning 404.
# ---------------------------------------------------------------------------

ONEKG_HOMEPAGE = "https://www.internationalgenome.org/"
ONEKG_RELEASE = "20130502"
_ONEKG_BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/"
               f"{ONEKG_RELEASE}/")

# VERIFIED HTTP 200 ON 2026-08-04.
ONEKG_PANEL_URL = _ONEKG_BASE + "integrated_call_samples_v3.20130502.ALL.panel"
ONEKG_VCF_TEMPLATE = (_ONEKG_BASE + "ALL.chr{chrom}."
                      "phase3_shapeit2_mvncall_integrated_v5b.20130502."
                      "genotypes.vcf.gz")
# chrX is NOT v5b. It is v1c, and a template that assumes otherwise 404s on the
# one chromosome most likely to be requested by accident.
ONEKG_X_VCF_URL = (_ONEKG_BASE + "ALL.chrX."
                   "phase3_shapeit2_mvncall_integrated_v1c.20130502."
                   "genotypes.vcf.gz")

# The same phase 3 v5a callset with the ID column populated, republished by the
# Browning lab as the standard Beagle reference panel. GT-only, so it is
# smaller than the IGSR release as well as more useful. VERIFIED HTTP 200 AND
# rsIDs PRESENT ON 2026-08-04.
ONEKG_BEAGLE_TEMPLATE = ("https://bochet.gcc.biostat.washington.edu/beagle/"
                         "1000_Genomes_phase3_v5a/b37.vcf/"
                         "chr{chrom}.1kg.phase3.v5a.vcf.gz")

ONEKG_SOURCES: tuple[str, ...] = ("beagle", "igsr")

AUTOSOMES: tuple[str, ...] = tuple(str(i) for i in range(1, 23))

# Compressed sizes in bytes, measured by HTTP range probe on 2026-08-04. These
# exist so --dry-run can state the real cost of a build before anybody starts
# one, rather than a guess. A user on a metered connection is entitled to know
# that the unfiltered autosomal set is 13.8 GiB before it starts arriving.
ONEKG_VCF_BYTES: dict[str, int] = {
    "1": 1165011543, "2": 1255861869, "3": 1058895733, "4": 1071020149,
    "5": 946995712, "6": 959144250, "7": 870098463, "8": 824304010,
    "9": 644258549, "10": 741683063, "11": 734583405, "12": 710187972,
    "13": 533918831, "14": 485261255, "15": 438397058, "16": 473472094,
    "17": 415901551, "18": 418186091, "19": 344793450, "20": 327135853,
    "21": 209774472, "22": 205612353, "X": 1908425030,
}

# Same measurement, same date, for the GT-only Browning mirror.
ONEKG_BEAGLE_BYTES: dict[str, int] = {
    "1": 754763598, "2": 807495004, "3": 688092385, "4": 706626850,
    "5": 614138429, "6": 639545208, "7": 570868069, "8": 533370343,
    "9": 418568075, "10": 487361943, "11": 477505002, "12": 463932679,
    "13": 351679176, "14": 316379099, "15": 284598111, "16": 302859853,
    "17": 268106491, "18": 273973424, "19": 228198264, "20": 212189048,
    "21": 139082504, "22": 135429468, "X": 329578525,
}

# Simons Genome Diversity Project, public tier.
#
# HONEST NOTE ON VERIFICATION: sharehost.hms.harvard.edu presented a
# certificate chain this environment could not validate on 2026-08-04, so the
# per-sample VCF filename template below could NOT be confirmed against the
# live listing. It is recorded as the documented shape and both the directory
# and the template are overridable from the command line, so a user whose copy
# is laid out differently is not stuck. The dry run prints this caveat.
SGDP_BASE = "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
SGDP_METADATA_URL = SGDP_BASE + "SGDP_metadata.279public.21signed.csv"
SGDP_VCF_TEMPLATE = SGDP_BASE + "vcf_variants/{sample}.annotated.nh2.variants.vcf.gz"
# The same study is indexed at EMBL-EBI, which is reachable when the Harvard
# mirror is not. Recorded so the panel can state a second provenance.
SGDP_ENA_PROJECT = "https://www.ebi.ac.uk/ena/browser/view/PRJEB9586"
SGDP_PUBLIC_SAMPLES = 279
SGDP_RESTRICTED_SAMPLES = 21

# PLINK-format cM maps on GRCh37. panel.map is a genetic map, not a marker
# list: Beagle and FLARE both need centimorgans, and interpolating them from a
# published map is the only honest way to get them from a VCF that has none.
# VERIFIED HTTP 200 ON 2026-08-04.
GENETIC_MAP_URL = ("https://bochet.gcc.biostat.washington.edu/beagle/"
                   "genetic_maps/plink.GRCh37.map.zip")
GENETIC_MAP_MEMBER = "plink.chr{chrom}.GRCh37.map"

# Named so the exclusions are recorded decisions rather than oversights, in the
# same spirit as BLOCKED in backend/external.py. Neither URL is ever fetched by
# default.
HGDP_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/"
HGDP_REUSE_STATEMENT_URL = HGDP_BASE + "README_HGDP_datareuse_statement.md"
# The gnomAD HGDP plus 1000 Genomes joint callset is the only openly published
# HGDP genotype set with usable population labels. VERIFIED HTTP 206 on
# 2026-08-04. It is GRCh38, which is a second reason it is not mixed into a
# GRCh37 panel without the user saying so twice.
HGDP_VCF_TEMPLATE = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
                     "release/3.1.2/vcf/genomes/"
                     "gnomad.genomes.v3.1.2.hgdp_tgp.chr{chrom}.vcf.bgz")
AADR_URL = ("https://reich.hms.harvard.edu/"
            "allen-ancient-dna-resource-aadr-downloadable-genotypes-"
            "present-day-and-ancient-dna-data")

USER_AGENT = "DNAInsight/3.0 (+https://github.com/dnainsight) panel builder"
HTTP_TIMEOUT = 600
CHUNK = 1 << 20


LICENCES: dict[str, str] = {
    "onekg": (
        "1000 Genomes Project phase 3, originating from IGSR at EMBL-EBI. "
        "Open with no restriction on use, per the project's data reuse "
        "statement. Distributed under EMBL-EBI terms of use, which permit "
        "research and commercial reuse and redistribution. Citation of the "
        "project is requested, not required. No SPDX identifier applies. The "
        "Browning lab mirror carries the identical phase 3 v5a callset under "
        "the same open terms and adds no terms of its own."
    ),
    "sgdp": (
        "Simons Genome Diversity Project, PUBLIC tier only "
        f"({SGDP_PUBLIC_SAMPLES} samples). Released without a data access "
        "agreement and unrestricted in use. Citation of Mallick et al. 2016 "
        "is requested. No SPDX identifier applies. The separate "
        f"{SGDP_RESTRICTED_SAMPLES}-sample restricted tier is NOT used and "
        "cannot be enabled from this builder."
    ),
    "genetic_map": (
        "PLINK-format GRCh37 genetic maps published by the Browning lab, "
        "University of Washington. Freely redistributable and used here only "
        "to interpolate centimorgan positions for panel.map."
    ),
    "hgdp": (
        "HGDP-CEPH via IGSR. Open access under the Fort Lauderdale Principles. "
        "The licence is not the problem. See REFUSAL_HGDP."
    ),
}


# ---------------------------------------------------------------------------
# Refusals
#
# Each refusal is a constant with a reason, so the text a user sees, the text
# in backend/external.py PANELS['onekg_sgdp']['excluded'] and the text in the
# manifest all say the same thing. A refusal explained in three places and
# worded three ways is a refusal nobody trusts.
# ---------------------------------------------------------------------------

REFUSAL_SGDP_RESTRICTED = (
    "REFUSED: the 21-sample restricted tier of the Simons Genome Diversity "
    "Project.\n"
    "  Why: access requires a signed agreement whose terms include 'I will not "
    "use the data for any commercial purposes'.\n"
    "  Consequence: DNAInsight is MIT and grants every downstream user the "
    "right to sell what they build. A panel containing that tier could not be "
    "described, shipped or reproduced under those terms, and the restriction "
    "would travel silently to anyone who rebuilt it.\n"
    "  There is no flag that enables this. It is excluded by construction.\n"
    "  Only the " + str(SGDP_PUBLIC_SAMPLES) + "-sample public tier is used."
)

REFUSAL_HGDP = (
    "REFUSED BY DEFAULT: the Human Genome Diversity Project (HGDP-CEPH).\n"
    "  This is NOT a licence objection. HGDP is openly accessible through IGSR "
    "under the Fort Lauderdale Principles, and --accept-terms does not answer "
    "the objection, which is why a second and different flag is required.\n"
    "  Why: consent. Nature Genetics, 24 November 2025, concluded that broad "
    "reuse of HGDP may diverge from what participants consented to. The PRIMED "
    "Consortium voted on 21 August 2024 to keep permitting its use while "
    "acknowledging 'failure to obtain informed consent consistent with current "
    "standards from many participants'.\n"
    "  To include it anyway you must pass BOTH --include-hgdp AND "
    "--accept-consent-caveat. Accepting a licence is not accepting this; that "
    "is the whole point of the second flag.\n"
    "  Note also that the only openly published HGDP genotype callset is on "
    "GRCh38 while this panel is " + TARGET_BUILD + ", so it is written to the "
    "separate hgdp_optional panel and never merged into panel.vcf.gz."
)

REFUSAL_AADR = (
    "REFUSED: the Allen Ancient DNA Resource (AADR).\n"
    "  Why: its terms were not readable at review time, and the compendium "
    "aggregates many datasets that each carry their own upstream terms.\n"
    "  Unread is not permissive. A compendium cannot grant rights its "
    "components did not grant, and this project does not bundle data whose "
    "terms nobody has read.\n"
    "  There is no flag that enables this. Excluded entirely until the terms "
    "are read and recorded in data/DATA_SOURCES.md.\n"
    "  Recorded location: " + AADR_URL
)

REFUSALS: dict[str, str] = {
    "sgdp_restricted": REFUSAL_SGDP_RESTRICTED,
    "hgdp": REFUSAL_HGDP,
    "aadr": REFUSAL_AADR,
}

STATISTICS: tuple[str, ...] = ("fst", "informativeness")


class PanelRefused(RuntimeError):
    """A source the builder will not fetch. Carries the printable reason."""

    def __init__(self, source: str, reason: str) -> None:
        super().__init__(reason)
        self.source = source
        self.reason = reason


# ---------------------------------------------------------------------------
# Gates
# ---------------------------------------------------------------------------

def refuse_restricted_sgdp() -> str:
    """Always refuse the restricted SGDP tier and return the reason.

    A function rather than a bare constant so that any future caller that wants
    the tier has to go through something that cannot return success.
    """
    return REFUSAL_SGDP_RESTRICTED


def refuse_aadr() -> str:
    """Always refuse the Allen Ancient DNA Resource and return the reason."""
    return REFUSAL_AADR


def hgdp_gate(include_hgdp: bool, accept_consent_caveat: bool) -> tuple[bool, str]:
    """Decide whether HGDP may be fetched. Returns ``(allowed, reason)``.

    Both flags are required and they are deliberately not interchangeable.
    ``--include-hgdp`` says "I want this data". ``--accept-consent-caveat`` says
    "I have read the consent objection and I am choosing to proceed anyway".
    A single flag would collapse those into one click, and a licence acceptance
    would answer neither, because the objection is not about the licence.
    """
    if not include_hgdp and not accept_consent_caveat:
        return False, REFUSAL_HGDP
    if include_hgdp and not accept_consent_caveat:
        return False, (
            REFUSAL_HGDP + "\n\n  --include-hgdp was passed without "
            "--accept-consent-caveat. Wanting the data is not the same as "
            "having read why it is gated, so the request is refused."
        )
    if accept_consent_caveat and not include_hgdp:
        return False, (
            REFUSAL_HGDP + "\n\n  --accept-consent-caveat was passed without "
            "--include-hgdp. Nothing was requested, so nothing is included."
        )
    return True, (
        "HGDP included on an explicit double opt-in: --include-hgdp and "
        "--accept-consent-caveat were both passed. The consent caveat above "
        "still applies and is recorded verbatim in BUILD.txt, so anybody who "
        "later finds this panel on disk can see what was accepted."
    )


# ---------------------------------------------------------------------------
# Population genetics, in pure Python
#
# No numpy, on purpose: rebuilding a panel must not require a scientific stack
# or a compiler toolchain. Every function here is small enough to check by hand
# against the definitions in the module docstring, and the test suite does
# exactly that.
# ---------------------------------------------------------------------------

def allele_frequency(counts: tuple[int, int]) -> float:
    """Alternate-allele frequency from ``(alt_copies, called_copies)``.

    Returns 0.0 when nothing was called, which keeps a marker with no data out
    of the ranking instead of crashing on a division. A no-call site is not an
    invariant site, but for ranking purposes both are uninformative.
    """
    alt, total = int(counts[0]), int(counts[1])
    return (alt / total) if total > 0 else 0.0


def wright_fst(counts: dict[str, tuple[int, int]]) -> float:
    """Wright's Fst across populations, as the ratio of heterozygosities.

    ``counts`` maps a population code to ``(alt_copies, called_copies)``.

        Hs  = sum_i w_i * 2 * p_i * (1 - p_i),   w_i = n_i / sum(n)
        Ht  = 2 * pbar * (1 - pbar),             pbar = sum(alt) / sum(n)
        Fst = (Ht - Hs) / Ht

    pbar is copy-weighted rather than a plain mean of the p_i, because the
    populations in this panel differ in size by nearly a factor of two and an
    unweighted mean would let a small population move the total heterozygosity
    as much as a large one.

    Returns 0.0 when the site is invariant across the whole panel (Ht == 0),
    which is the correct answer rather than a NaN: an invariant site
    differentiates nothing. Populations with no called copies are ignored.
    Fewer than two populations with data also gives 0.0, because Fst is a
    statement about differentiation BETWEEN groups and there are none.
    """
    usable = {code: (int(a), int(n)) for code, (a, n) in counts.items() if int(n) > 0}
    if len(usable) < 2:
        return 0.0

    total_copies = sum(n for _a, n in usable.values())
    total_alt = sum(a for a, _n in usable.values())
    if total_copies <= 0:
        return 0.0

    pbar = total_alt / total_copies
    ht = 2.0 * pbar * (1.0 - pbar)
    if ht <= 0.0:
        return 0.0

    hs = 0.0
    for alt, n in usable.values():
        p = alt / n
        hs += (n / total_copies) * 2.0 * p * (1.0 - p)

    fst = (ht - hs) / ht
    # Clamp only the floating-point dust. A genuinely negative Fst is not
    # possible in this estimator, so anything below zero is rounding.
    if fst < 0.0:
        return 0.0
    return fst


def population_fst(counts: dict[str, tuple[int, int]], code: str) -> float:
    """One-versus-rest Fst for a single population.

    The standard AIM criterion: collapse every other population into one pooled
    group and ask how far this population sits from it. Implemented by calling
    :func:`wright_fst` on a two-group partition, so there is one estimator in
    this module and not two that could disagree.

    Returns 0.0 when ``code`` has no called copies or when there is no other
    population to compare it against.
    """
    mine = counts.get(code)
    if not mine or int(mine[1]) <= 0:
        return 0.0
    rest_alt = 0
    rest_total = 0
    for other, (alt, n) in counts.items():
        if other == code:
            continue
        rest_alt += int(alt)
        rest_total += int(n)
    if rest_total <= 0:
        return 0.0
    return wright_fst({code: (int(mine[0]), int(mine[1])),
                       "REST": (rest_alt, rest_total)})


def _x_log_x(value: float) -> float:
    """``x * ln(x)``, defined as 0 at x == 0 by the usual limit."""
    if value <= 0.0:
        return 0.0
    return value * math.log(value)


def informativeness(counts: dict[str, tuple[int, int]]) -> float:
    """Rosenberg et al. (2003) informativeness for assignment, in nats.

        In = sum_j [ -pbar_j ln pbar_j + (1/K) sum_i p_ij ln p_ij ]

    over the two alleles j of a biallelic site and the K populations i, with
    pbar_j the UNWEIGHTED mean frequency across populations. Unweighted is the
    published definition and it is the right one here: In asks how much a
    genotype tells you about which population somebody came from, and that
    question does not become less interesting because a population is small.
    That is the deliberate difference from :func:`wright_fst`, which weights.

    The maximum for a biallelic marker with two populations fixed for opposite
    alleles is ln(2), about 0.6931. Returns 0.0 for fewer than two populations
    with data.
    """
    freqs = [allele_frequency(c) for c in counts.values() if int(c[1]) > 0]
    k = len(freqs)
    if k < 2:
        return 0.0

    total = 0.0
    for alleles in (freqs, [1.0 - p for p in freqs]):
        pbar = sum(alleles) / k
        total += -_x_log_x(pbar) + sum(_x_log_x(p) for p in alleles) / k
    # Same dust clamp as Fst. In is non-negative by construction.
    return total if total > 0.0 else 0.0


def marker_statistic(counts: dict[str, tuple[int, int]],
                     statistic: str = "fst") -> float:
    """Dispatch to the global marker statistic named by ``--statistic``."""
    if statistic == "informativeness":
        return informativeness(counts)
    if statistic == "fst":
        return wright_fst(counts)
    raise ValueError(f"unknown statistic {statistic!r}; expected one of "
                     f"{', '.join(STATISTICS)}")


def _population_statistic(counts: dict[str, tuple[int, int]], code: str,
                          statistic: str) -> float:
    """Per-population ranking value under the selected statistic.

    For informativeness the one-versus-rest partition is used too, so the two
    statistics rank the same kind of thing and a user switching between them
    gets a comparable file rather than a differently shaped one.
    """
    if statistic == "fst":
        return population_fst(counts, code)
    mine = counts.get(code)
    if not mine or int(mine[1]) <= 0:
        return 0.0
    rest_alt = sum(int(a) for c, (a, _n) in counts.items() if c != code)
    rest_total = sum(int(n) for c, (_a, n) in counts.items() if c != code)
    if rest_total <= 0:
        return 0.0
    return informativeness({code: (int(mine[0]), int(mine[1])),
                            "REST": (rest_alt, rest_total)})


class TopN:
    """A bounded max-selection heap. Keeps the N highest-scoring items.

    The panel has tens of millions of candidate markers and we want a few
    thousand per population. Sorting the lot would cost memory proportional to
    the source; a heap of fixed size N costs memory proportional to the answer.
    Ties are broken by insertion order, which is genomic order here, so a
    rebuild from the same inputs produces the same file.
    """

    def __init__(self, size: int) -> None:
        self.size = max(0, int(size))
        self._heap: list[tuple[float, int, Any]] = []
        self._seq = 0

    def add(self, score: float, item: Any) -> None:
        """Offer one scored item to the selection."""
        if self.size == 0:
            return
        self._seq += 1
        # Negate the sequence so that, among equal scores, the EARLIER item
        # compares larger and therefore survives. Without this the survivor of
        # a tie depends on heap internals and the output stops being stable.
        entry = (float(score), -self._seq, item)
        if len(self._heap) < self.size:
            heapq.heappush(self._heap, entry)
        elif entry > self._heap[0]:
            heapq.heapreplace(self._heap, entry)

    def items(self) -> list[tuple[float, Any]]:
        """Selected items, highest score first."""
        ordered = sorted(self._heap, key=lambda e: (-e[0], -e[1]))
        return [(score, item) for score, _seq, item in ordered]

    def __len__(self) -> int:
        return len(self._heap)


# ---------------------------------------------------------------------------
# Source parsing
# ---------------------------------------------------------------------------

def parse_panel_file(lines: Iterable[str]) -> dict[str, tuple[str, str]]:
    """Parse the 1000 Genomes ``.panel`` file into sample -> (pop, superpop).

    Four whitespace-separated columns: sample, pop, super_pop, gender. Columns
    are resolved by NAME through the shared helper rather than by position, for
    the same reason build_full_reference.py resolves ClinVar by name: the file
    has gained a trailing column before and positions are not a contract.
    """
    wanted = {
        "sample": ("sample", "sample_id", "id"),
        "pop": ("pop", "population"),
        "superpop": ("super_pop", "superpop", "super_population"),
    }
    columns: dict[str, int] = {}
    out: dict[str, tuple[str, str]] = {}
    for line in lines:
        text = str(line or "").strip()
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        if not columns:
            columns = resolve_columns(fields, wanted, ("sample", "pop"))
            continue
        if len(fields) <= columns["pop"]:
            continue
        sample = fields[columns["sample"]].strip()
        pop = fields[columns["pop"]].strip().upper()
        idx = columns.get("superpop")
        superpop = (fields[idx].strip().upper()
                    if idx is not None and idx < len(fields) else "")
        if sample and pop:
            out[sample] = (pop, superpop or "UNKNOWN")
    return out


def parse_sgdp_metadata(lines: Iterable[str]) -> dict[str, tuple[str, str]]:
    """Parse the SGDP public metadata CSV into sample -> (population, region).

    SGDP publishes a comma-separated metadata table whose header has been
    reworded between releases, so the columns are resolved by name against a
    candidate list. A sample whose row does not carry both a usable identifier
    and a population is skipped rather than guessed at: an unlabelled sample in
    a supervised ancestry panel is worse than a missing one, because it becomes
    a population of one with a made-up name.
    """
    wanted = {
        "sample": ("Sample ID (Illumina)", "Sample_ID", "SampleID",
                   "Sample ID", "Illumina ID", "sample"),
        "population": ("Population ID", "Population_ID", "Population",
                       "population"),
        "region": ("Region", "region", "Continent"),
    }
    columns: dict[str, int] = {}
    out: dict[str, tuple[str, str]] = {}
    for line in lines:
        text = str(line or "").rstrip("\r\n")
        if not text.strip() or text.startswith("#"):
            continue
        fields = [f.strip().strip('"') for f in text.split(",")]
        if not columns:
            try:
                columns = resolve_columns(fields, wanted, ("sample", "population"))
            except KeyError:
                # Not the header yet, or a header this builder does not know.
                # Keep looking rather than crashing on a preamble line.
                continue
            continue
        if len(fields) <= columns["population"]:
            continue
        sample = fields[columns["sample"]]
        pop = fields[columns["population"]].strip().upper().replace(" ", "_")
        idx = columns.get("region")
        region = (fields[idx].strip().upper().replace(" ", "_")
                  if idx is not None and idx < len(fields) else "")
        if sample and pop:
            out[sample] = (pop, region or "SGDP")
    return out


def vcf_sample_names(header_line: str) -> list[str]:
    """Sample names from a VCF ``#CHROM`` header line.

    The first nine columns are fixed VCF fields; everything after FORMAT is a
    sample. A header with no FORMAT column carries no genotypes at all, and
    that returns an empty list rather than slicing something that is not there.
    """
    fields = str(header_line or "").rstrip("\n").split("\t")
    if len(fields) < 10 or fields[8].upper() != "FORMAT":
        return []
    return [f.strip() for f in fields[9:]]


def parse_vcf_record(line: str) -> dict | None:
    """Parse one VCF data line into a record dict, or None when unusable.

    Returns None for header lines, for multiallelic records and for anything
    that is not a single-nucleotide substitution. Both exclusions are
    deliberate:

      * a multiallelic record has no single alternate frequency, so every
        statistic in this module would be computed on a collapsed allele that
        does not exist;
      * an indel cannot be matched against a consumer array export, which
        reports substitutions only, so it could never be used by the
        --array-file path anyway.

    ``rsid`` is lower-cased so it joins directly against the array coverage set
    produced by build_full_reference.load_array_rsids.
    """
    text = str(line or "")
    if not text or text.startswith("#"):
        return None
    fields = text.rstrip("\n").split("\t")
    if len(fields) < 10:
        return None
    ref = fields[3].strip().upper()
    alt = fields[4].strip().upper()
    if len(ref) != 1 or len(alt) != 1:
        return None
    if ref not in "ACGT" or alt not in "ACGT":
        return None
    rsid = fields[2].strip().lower()
    try:
        pos = int(fields[1])
    except (TypeError, ValueError):
        return None
    return {
        "chrom": fields[0].strip(),
        "pos": pos,
        "rsid": rsid if rsid.startswith("rs") else "",
        "ref": ref,
        "alt": alt,
        "genotypes": fields[9:],
        "fixed": fields[:9],
    }


def allele_counts(genotypes: Sequence[str],
                  samples: Sequence[str],
                  sample_pops: dict[str, tuple[str, str]],
                  ) -> dict[str, tuple[int, int]]:
    """Per-population ``(alt_copies, called_copies)`` for one VCF record.

    Handles phased ``0|1`` and unphased ``0/1`` alike, and haploid calls such
    as the male X, which contribute one copy rather than two. Missing calls
    contribute nothing to either total, so a population's frequency is over the
    copies actually called and not over the copies that could have been.

    Only the first colon-delimited subfield is read, which is GT by VCF
    definition, so a record carrying DS or GP costs nothing extra.
    """
    out: dict[str, list[int]] = {}
    limit = min(len(genotypes), len(samples))
    for i in range(limit):
        pop = sample_pops.get(samples[i])
        if pop is None:
            continue
        cell = genotypes[i]
        gt = cell.split(":", 1)[0] if ":" in cell else cell
        gt = gt.strip()
        if not gt or gt[0] == ".":
            continue
        bucket = out.setdefault(pop[0], [0, 0])
        for token in gt.replace("|", "/").split("/"):
            if not token or token == ".":
                continue
            if token == "0":
                bucket[1] += 1
            elif token == "1":
                bucket[0] += 1
                bucket[1] += 1
            # Any other index means a third allele at a site this builder
            # already restricted to biallelic. Counting it as alt would be a
            # quiet lie, so it is dropped from both totals.
    return {code: (vals[0], vals[1]) for code, vals in out.items()}


def interpolate_cm(position: int,
                   positions: Sequence[int],
                   centimorgans: Sequence[float]) -> float:
    """Linearly interpolate a centimorgan value for one base-pair position.

    ``positions`` must be sorted ascending and the two sequences parallel.
    Outside the mapped range the nearest end value is returned rather than an
    extrapolation: recombination rate at a telomere is not linear and inventing
    a value there would put fake structure into the map that Beagle would then
    trust.

    Returns 0.0 when there is no map at all, which Beagle reads as "no genetic
    distance information", and which is honest.
    """
    if not positions:
        return 0.0
    if position <= positions[0]:
        return float(centimorgans[0])
    if position >= positions[-1]:
        return float(centimorgans[-1])
    i = bisect.bisect_left(positions, position)
    if i < len(positions) and positions[i] == position:
        return float(centimorgans[i])
    left_bp, right_bp = positions[i - 1], positions[i]
    left_cm, right_cm = float(centimorgans[i - 1]), float(centimorgans[i])
    span = right_bp - left_bp
    if span <= 0:
        return left_cm
    return left_cm + (right_cm - left_cm) * ((position - left_bp) / span)


def sha256_of(path: Path, chunk: int = CHUNK) -> str:
    """SHA-256 of a file, read in chunks so a multi-gigabyte panel is fine."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# The plan, and the consent gate
# ---------------------------------------------------------------------------

def panel_dir(panel_id: str = PANEL_ID) -> Path:
    """Target directory for a panel, under external.panel_root()."""
    return external.panel_root() / panel_id


def _human_bytes(count: int) -> str:
    """Format a byte count the way a person reading a download prompt reads."""
    value = float(count)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:,.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0
    return f"{value:,.1f} TiB"


def download_plan(args: argparse.Namespace) -> dict:
    """Everything this invocation would download, with sizes and licences.

    Built before any network call and printed verbatim, so the user sees the
    plan whether or not they then approve it. The dry run and the real run
    share this function on purpose: a plan that is only produced in dry-run
    mode is a plan that can stop describing the real run.
    """
    chromosomes = _selected_chromosomes(args)
    source = getattr(args, "onekg_source", "beagle")
    items: list[dict] = []

    onekg_bytes = sum(_onekg_vcf_bytes(c, source) for c in chromosomes)
    onekg = {
        "source": f"1000 Genomes phase 3 ({source} distribution)",
        "what": f"sample panel plus {len(chromosomes)} genotype VCF(s)",
        "urls": ([ONEKG_PANEL_URL]
                 + [_onekg_vcf_url(c, source) for c in chromosomes]),
        "bytes": onekg_bytes,
        "licence": LICENCES["onekg"],
        "homepage": ONEKG_HOMEPAGE,
    }
    if source == "igsr":
        onekg["caveat"] = (
            "The IGSR release VCFs publish the ID column as '.' for every "
            "record, verified against the live files on 2026-08-04. With this "
            "source --array-file matches nothing and informative_markers.tsv "
            "comes out empty. Use --onekg-source beagle unless you have a "
            "specific reason not to."
        )
    items.append(onekg)

    if not args.skip_sgdp:
        items.append({
            "source": f"SGDP public tier ({SGDP_PUBLIC_SAMPLES} samples)",
            "what": ("metadata CSV plus one per-sample VCF each"
                     if not args.sgdp_dir else
                     f"metadata CSV; genotypes read from {args.sgdp_dir}"),
            "urls": [(args.sgdp_metadata or SGDP_METADATA_URL),
                     (args.sgdp_vcf_template or SGDP_VCF_TEMPLATE)],
            # Deliberately None rather than a made-up number. The per-sample
            # VCF sizes could not be measured because the host's certificate
            # chain did not validate at review time, and a fabricated estimate
            # in a consent prompt is worse than an honest "unknown".
            "bytes": None,
            "licence": LICENCES["sgdp"],
            "homepage": SGDP_ENA_PROJECT,
            "caveat": (
                "The per-sample VCF filename template could NOT be verified "
                "against the live listing on 2026-08-04: the host presented a "
                "certificate chain that did not validate. Override it with "
                "--sgdp-vcf-template, or point --sgdp-dir at a local copy."
            ),
        })

    items.append({
        "source": "PLINK GRCh37 genetic maps (Browning lab)",
        "what": "one zip, cached then read per chromosome",
        "urls": [GENETIC_MAP_URL],
        "bytes": None,
        "licence": LICENCES["genetic_map"],
        "homepage": "https://faculty.washington.edu/browning/beagle/beagle.html",
    })

    hgdp_allowed, hgdp_reason = hgdp_gate(args.include_hgdp,
                                          args.accept_consent_caveat)
    if hgdp_allowed:
        items.append({
            "source": "HGDP-CEPH (double opt-in)",
            "what": f"{len(chromosomes)} joint-callset VCF(s), GRCh38",
            "urls": [HGDP_REUSE_STATEMENT_URL]
                    + [HGDP_VCF_TEMPLATE.format(chrom=c) for c in chromosomes],
            "bytes": None,
            "licence": LICENCES["hgdp"],
            "homepage": HGDP_BASE,
            "caveat": hgdp_reason,
        })

    known = [i["bytes"] for i in items if i["bytes"]]
    return {
        "panel": PANEL_ID,
        "build": TARGET_BUILD,
        "out_dir": str(panel_dir(args.panel_id)),
        "chromosomes": list(chromosomes),
        "onekg_source": source,
        "statistic": args.statistic,
        "array_file": args.array_file or "",
        "items": items,
        "known_bytes": sum(known),
        "outputs": list(PANEL_OUTPUTS),
        "refusals": {
            "sgdp_restricted": refuse_restricted_sgdp(),
            "aadr": refuse_aadr(),
            **({} if hgdp_allowed else {"hgdp": hgdp_reason}),
        },
        "hgdp_included": hgdp_allowed,
    }


def print_plan(plan: dict, *, dry_run: bool) -> None:
    """Print the download plan, the licences and every refusal.

    Printed on every run, not only on a dry run. A user who passes
    --accept-terms is agreeing to something specific and is entitled to see
    what it was, in the same terminal, before the first byte moves.
    """
    bar = "=" * 74
    print(bar)
    print("DNAInsight v3.0 reference panel builder")
    print(bar)
    print(f"panel:        {plan['panel']}")
    print(f"genome build: {plan['build']}")
    print(f"output dir:   {plan['out_dir']}")
    print(f"chromosomes:  {', '.join(plan['chromosomes'])}")
    print(f"1000G source: {plan['onekg_source']}")
    print(f"AIM statistic:{plan['statistic']:>13}")
    print(f"array filter: {plan['array_file'] or '(none: full panel, much larger)'}")
    print()

    print("WILL DOWNLOAD")
    print("-" * 74)
    for item in plan["items"]:
        size = (_human_bytes(item["bytes"]) if item["bytes"]
                else "size not published, not estimated")
        print(f"  {item['source']}")
        print(f"    what:     {item['what']}")
        print(f"    size:     {size}")
        print(f"    homepage: {item['homepage']}")
        for url in item["urls"][:3]:
            print(f"    url:      {url}")
        if len(item["urls"]) > 3:
            print(f"    url:      ... and {len(item['urls']) - 3} more of the "
                  f"same shape")
        print(f"    licence:  {item['licence']}")
        if item.get("caveat"):
            print(f"    caveat:   {item['caveat']}")
        print()
    print(f"  Total of the measurable parts: {_human_bytes(plan['known_bytes'])}")
    print("  Parts marked 'not estimated' are additional to that figure.")
    print()

    print("WILL NOT DOWNLOAD")
    print("-" * 74)
    for key in sorted(plan["refusals"]):
        for line in plan["refusals"][key].splitlines():
            print(f"  {line}")
        print()

    print("WILL WRITE")
    print("-" * 74)
    for name in plan["outputs"]:
        print(f"  {plan['out_dir']}/{name}")
    print()
    print("  Nothing is written inside the repository. This directory is")
    print("  outside the tree on purpose, exactly like the SNPedia cache.")
    print()

    if dry_run:
        print(bar)
        print("DRY RUN. Nothing was downloaded and nothing was written.")
        print("Re-run with --accept-terms to proceed. Passing --accept-terms")
        print("means you accept the licences printed above for your own copy;")
        print("DNAInsight ships none of this data.")
        print(bar)


# ---------------------------------------------------------------------------
# Transport
# ---------------------------------------------------------------------------

def _require_requests() -> None:
    if requests is None:  # pragma: no cover - only when the dep is missing
        raise RuntimeError("the requests package is not installed")


def _onekg_vcf_url(chrom: str, source: str = "beagle") -> str:
    """Genotype VCF URL for one 1000 Genomes phase 3 chromosome.

    ``beagle`` is the default because the IGSR release strips the ID column,
    and without rsIDs the --array-file join and informative_markers.tsv are
    both empty. ``igsr`` returns the canonical release URL, and chrX there is
    v1c rather than v5b, which is a filename trap worth having in one place.
    """
    if source == "igsr":
        return (ONEKG_X_VCF_URL if chrom.upper() == "X"
                else ONEKG_VCF_TEMPLATE.format(chrom=chrom))
    return ONEKG_BEAGLE_TEMPLATE.format(chrom=chrom.upper() if chrom.upper() == "X"
                                        else chrom)


def _onekg_vcf_bytes(chrom: str, source: str = "beagle") -> int:
    """Measured compressed size of one chromosome VCF, or 0 when unknown."""
    table = ONEKG_VCF_BYTES if source == "igsr" else ONEKG_BEAGLE_BYTES
    return table.get(chrom.upper(), 0)


def download_file(url: str, dest: Path, timeout: int = HTTP_TIMEOUT) -> int:
    """Stream a URL to disk in chunks. Returns bytes written.

    Used only for the genetic map zip. A zip index sits at the END of the file,
    so unlike a gzip it cannot be read from a forward-only socket, which is the
    same reason build_full_reference.py caches the GWAS archive before reading
    it. The cache lives under the panel directory, never under data/.
    """
    _require_requests()
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with requests.get(url, stream=True, timeout=timeout,
                      headers={"User-Agent": USER_AGENT}) as resp:
        resp.raise_for_status()
        with open(dest, "wb") as fh:
            for block in resp.iter_content(chunk_size=CHUNK):
                if not block:
                    continue
                fh.write(block)
                written += len(block)
    return written


def _load_genetic_map(archive: Path, chrom: str) -> tuple[list[int], list[float]]:
    """Read one chromosome's PLINK map from the cached zip.

    Returns parallel ``(positions, centimorgans)`` lists. One chromosome at a
    time, so the whole genome-wide map is never resident. A missing member
    returns empty lists and every cM comes out 0.0, which is a usable panel
    with a stated limitation rather than a failed build.
    """
    member = GENETIC_MAP_MEMBER.format(chrom=chrom)
    positions: list[int] = []
    centimorgans: list[float] = []
    try:
        with zipfile.ZipFile(archive) as zf:
            names = {n.rsplit("/", 1)[-1]: n for n in zf.namelist()}
            if member not in names:
                return positions, centimorgans
            with zf.open(names[member]) as handle:
                for raw in handle:
                    fields = raw.decode("utf-8", "replace").split()
                    if len(fields) < 4:
                        continue
                    try:
                        centimorgans.append(float(fields[2]))
                        positions.append(int(fields[3]))
                    except ValueError:
                        continue
    except (OSError, zipfile.BadZipFile):
        return [], []
    return positions, centimorgans


# ---------------------------------------------------------------------------
# The build itself
# ---------------------------------------------------------------------------

def _selected_chromosomes(args: argparse.Namespace) -> tuple[str, ...]:
    """Resolve --chromosomes into an ordered tuple of chromosome names."""
    raw = str(getattr(args, "chromosomes", "") or "").strip()
    if not raw or raw.lower() == "auto":
        return AUTOSOMES
    out: list[str] = []
    for token in raw.replace(";", ",").split(","):
        token = token.strip().upper().removeprefix("CHR")
        if token and token not in out and token in ONEKG_VCF_BYTES:
            out.append(token)
    return tuple(out) or AUTOSOMES


def _scratch_dir(base: Path) -> Path:
    """Working directory for caches and per-sample columns, under the panel."""
    return base / "_scratch"


def _stream_chromosome(url: str) -> Iterator[str]:
    """Yield decompressed VCF lines from a remote bgzipped file.

    bgzip output is a sequence of concatenated gzip members, which
    gzip.GzipFile reads transparently, so the same streaming helper that
    build_full_reference.py uses for ClinVar works here unchanged. That is the
    entire reason nothing needs a compiled bgzf reader.
    """
    return stream_gzip_lines(url, timeout=HTTP_TIMEOUT)


def _build_onekg_pass(chromosomes: Sequence[str],
                      sample_pops: dict[str, tuple[str, str]],
                      coverage: set[str],
                      out_vcf: Path,
                      out_map: Path,
                      map_archive: Path | None,
                      statistic: str,
                      min_maf: float,
                      per_population: int,
                      limit: int,
                      source: str = "beagle") -> dict:
    """Single streaming pass over the 1000 Genomes VCFs.

    Writes panel.vcf.gz and panel.map as it goes and accumulates the per
    population marker ranking in bounded heaps. One pass, because a second pass
    would mean downloading 13.8 GiB twice, and the ranking does not need to
    change which records are written: the AIM list labels the panel, it does
    not filter it.

    Returns a summary dict for BUILD.txt.
    """
    populations = sorted({pop for pop, _sup in sample_pops.values()})
    heaps = {code: TopN(per_population) for code in populations}
    index: dict[tuple[str, int, str, str], int] = {}
    order: list[str] = []
    seen = 0
    kept = 0
    unnamed = 0
    header_written = False

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_vcf, "wt", encoding="utf-8") as vcf_out, \
            open(out_map, "w", encoding="utf-8") as map_out:
        for chrom in chromosomes:
            url = _onekg_vcf_url(chrom, source)
            print(f"[1000g] streaming chr{chrom}: {url}", flush=True)
            positions, centimorgans = ([], [])
            if map_archive is not None:
                positions, centimorgans = _load_genetic_map(map_archive, chrom)
                print(f"[1000g] chr{chrom} genetic map rows: {len(positions):,}")
            samples: list[str] = []
            chrom_seen = 0

            for line in _stream_chromosome(url):
                if line.startswith("##"):
                    if not header_written:
                        vcf_out.write(line + "\n")
                    continue
                if line.startswith("#CHROM"):
                    samples = vcf_sample_names(line)
                    if not header_written:
                        # Restrict the written header to samples the panel
                        # actually labels. A column with no population is a
                        # column a supervised estimator cannot use.
                        keep = [s for s in samples if s in sample_pops]
                        vcf_out.write("\t".join(line.split("\t")[:9] + keep) + "\n")
                        order.extend(keep)
                        header_written = True
                    continue

                record = parse_vcf_record(line)
                if record is None:
                    continue
                seen += 1
                chrom_seen += 1
                if seen % 200000 == 0:
                    print(f"[1000g] {seen:>12,} records read, {kept:>9,} kept",
                          flush=True)
                if limit and chrom_seen > limit:
                    break
                if not record["rsid"]:
                    # No rsID means the marker cannot be joined against an
                    # array export or named in informative_markers.tsv. It is
                    # still a valid imputation marker, so it is counted and
                    # kept when no coverage filter is in force, but a coverage
                    # filter has nothing to test it against and drops it.
                    unnamed += 1
                    if coverage:
                        continue
                elif coverage and record["rsid"] not in coverage:
                    continue

                counts = allele_counts(record["genotypes"], samples, sample_pops)
                total_alt = sum(a for a, _n in counts.values())
                total_n = sum(n for _a, n in counts.values())
                if total_n <= 0:
                    continue
                freq = total_alt / total_n
                maf = min(freq, 1.0 - freq)
                if maf < min_maf:
                    continue

                key = (record["chrom"], record["pos"], record["ref"], record["alt"])
                if key in index:
                    continue
                index[key] = kept

                columns = {s: i for i, s in enumerate(samples)}
                kept_calls = [record["genotypes"][columns[s]]
                              for s in order if s in columns]
                vcf_out.write("\t".join(record["fixed"] + kept_calls) + "\n")

                cm = interpolate_cm(record["pos"], positions, centimorgans)
                # PLINK map marker names must be unique. When the source gave
                # no rsID the position is used, which is unique by definition
                # and is honest about being a coordinate rather than an
                # identifier.
                name = record["rsid"] or f"{record['chrom']}:{record['pos']}"
                map_out.write(f"{record['chrom']}\t{name}\t"
                              f"{cm:.6f}\t{record['pos']}\n")
                kept += 1

                if not record["rsid"]:
                    continue
                for code in populations:
                    value = _population_statistic(counts, code, statistic)
                    if value > 0.0:
                        heaps[code].add(value, record["rsid"])

    return {
        "records_seen": seen,
        "markers_kept": kept,
        "markers_without_rsid": unnamed,
        "samples": order,
        "populations": populations,
        "heaps": heaps,
        "index": index,
    }


def _sgdp_column(url_or_path: str,
                 index: dict[tuple[str, int, str, str], int],
                 marker_count: int,
                 local: bool) -> bytearray:
    """Genotype column for one SGDP sample, aligned to the panel marker index.

    One byte per panel marker: 0, 1 or 2 alternate copies, 255 for no call.
    Written to disk by the caller rather than held, so merging 279 samples
    costs 279 small files and a constant amount of memory instead of a matrix.

    The sample's own VCF is streamed and every record is looked up by
    (chrom, pos, ref, alt). Matching on position alone would silently merge a
    different alternate allele into the panel, which inverts a genotype.
    """
    column = bytearray([255]) * marker_count
    if local:
        opener = gzip.open(url_or_path, "rt", encoding="utf-8", errors="replace")
        lines: Iterable[str] = opener
    else:
        lines = _stream_chromosome(url_or_path)
    for line in lines:
        record = parse_vcf_record(line)
        if record is None:
            continue
        slot = index.get((record["chrom"], record["pos"],
                          record["ref"], record["alt"]))
        if slot is None:
            continue
        counts = allele_counts(record["genotypes"], ["S"], {"S": ("X", "X")})
        alt, called = counts.get("X", (0, 0))
        column[slot] = alt if called else 255
    return column


def _merge_columns(panel_vcf: Path,
                   merged_vcf: Path,
                   columns: Sequence[tuple[str, Path]]) -> int:
    """Rewrite panel.vcf.gz with extra per-sample genotype columns appended.

    The panel is streamed in and out one record at a time while one byte is
    read from each column file per record, so peak memory is one VCF line plus
    one byte per extra sample. That is what makes 279 extra samples affordable
    on a personal machine.
    """
    if not columns:
        return 0
    handles = [open(path, "rb") for _name, path in columns]
    names = [name for name, _path in columns]
    written = 0
    try:
        with gzip.open(panel_vcf, "rt", encoding="utf-8") as src, \
                gzip.open(merged_vcf, "wt", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("##"):
                    dst.write(line)
                    continue
                if line.startswith("#CHROM"):
                    dst.write(line.rstrip("\n") + "\t" + "\t".join(names) + "\n")
                    continue
                calls = []
                for handle in handles:
                    byte = handle.read(1)
                    value = byte[0] if byte else 255
                    calls.append({0: "0|0", 1: "0|1", 2: "1|1"}.get(value, ".|."))
                dst.write(line.rstrip("\n") + "\t" + "\t".join(calls) + "\n")
                written += 1
    finally:
        for handle in handles:
            handle.close()
    return written


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_populations_tsv(path: Path,
                          sample_pops: dict[str, tuple[str, str]],
                          order: Sequence[str] = ()) -> int:
    """Write populations.tsv as sample, population, superpopulation.

    The column layout is dictated by backend.ancestry.parse_population_map,
    which also treats a first field of 'sample' as a header and skips it.
    ``order`` fixes the row order to the VCF's own sample order so that a
    reader can zip this file against the VCF columns without sorting.
    """
    rows = [s for s in order if s in sample_pops] or sorted(sample_pops)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("sample\tpopulation\tsuperpop\n")
        for sample in rows:
            pop, superpop = sample_pops[sample]
            fh.write(f"{sample}\t{pop}\t{superpop}\n")
    return len(rows)


def write_q_columns(path: Path, populations: Sequence[str],
                    superpops: dict[str, str] | None = None) -> int:
    """Write q_columns.tsv, the .Q column to population-code mapping.

    A fastmixture .Q matrix carries no column names. In supervised and
    projection modes the columns follow the order of the population labels the
    model was fitted with, and that order is this file. Without it
    backend.ancestry falls back to calling the columns component_1,
    component_2 and so on and says out loud that they are unlabelled, which is
    correct but much less useful. Guessing the mapping would produce a
    confident, wrong ancestry report.
    """
    labels = superpops or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# fastmixture .Q column order for this panel.\n")
        fh.write("# One population code per line, column 1 first.\n")
        for code in populations:
            fh.write(f"{code}\t{labels.get(code, 'UNKNOWN')}\n")
    return len(populations)


def write_informative_markers(path: Path,
                              selections: dict[str, list[tuple[float, str]]],
                              statistic: str) -> int:
    """Write informative_markers.tsv as population code then rsID.

    Two columns, whitespace separated, which is exactly what
    backend.ancestry._read_informative_markers expects. The score is carried in
    a trailing comment column rather than a third field, because a third field
    would be silently ignored by that reader and a value nobody reads is a
    value that rots.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# Ancestry-informative markers, ranked by {statistic} "
                 f"(one versus rest).\n")
        fh.write("# Columns: population_code, rsid. See the module docstring "
                 "of data/build_panel.py\n")
        fh.write("# for the definition of the statistic and its limitations.\n")
        for code in sorted(selections):
            for _score, rsid in selections[code]:
                fh.write(f"{code}\t{rsid}\n")
                total += 1
    return total


def write_build_txt(path: Path, facts: dict) -> None:
    """Write BUILD.txt: build, versions, date, licences and output hashes.

    A panel found on disk months later must be able to state its own
    provenance, exactly like the meta table build_full_reference.py writes into
    reference.db. The hashes are of the other outputs; BUILD.txt cannot contain
    its own hash, and pretending otherwise would be a checksum nobody can
    verify.
    """
    lines: list[str] = []
    lines.append("DNAInsight reference panel")
    lines.append("=" * 74)
    lines.append(f"panel_id:        {facts['panel_id']}")
    lines.append(f"genome_build:    {facts['build']}")
    lines.append("builder:         data/build_panel.py")
    lines.append(f"builder_version: {facts['builder_version']}")
    lines.append(f"retrieved_at:    {facts['retrieved_at']}")
    lines.append(f"chromosomes:     {', '.join(facts['chromosomes'])}")
    lines.append(f"samples:         {facts['sample_count']}")
    lines.append(f"populations:     {facts['population_count']}")
    lines.append(f"markers:         {facts['marker_count']}")
    lines.append(f"aim_statistic:   {facts['statistic']}")
    lines.append(f"min_maf:         {facts['min_maf']}")
    lines.append(f"array_filter:    {facts['array_filter'] or 'none'}")
    lines.append("")
    lines.append("SOURCE VERSIONS")
    lines.append("-" * 74)
    for name, version in sorted(facts["source_versions"].items()):
        lines.append(f"  {name:<16} {version}")
    lines.append("")
    lines.append("LICENCES, VERBATIM")
    lines.append("-" * 74)
    for name, text in sorted(facts["licences"].items()):
        lines.append(f"  [{name}]")
        for chunk in text.split(". "):
            if chunk.strip():
                lines.append(f"    {chunk.strip().rstrip('.')}.")
        lines.append("")
    lines.append("EXCLUDED, AND WHY")
    lines.append("-" * 74)
    for name in sorted(facts["refusals"]):
        for line in facts["refusals"][name].splitlines():
            lines.append(f"  {line}")
        lines.append("")
    lines.append("OUTPUT CHECKSUMS (sha256)")
    lines.append("-" * 74)
    for name in sorted(facts["hashes"]):
        lines.append(f"  {facts['hashes'][name]}  {name}")
    lines.append("")
    lines.append("PANEL BIAS, STATED RATHER THAN HIDDEN")
    lines.append("-" * 74)
    lines.append("  Every openly licensed reference panel under-represents "
                 "non-European")
    lines.append("  ancestries. Results derived from this panel carry a "
                 "coverage figure and")
    lines.append("  the interface states the limitation. See "
                 "backend/ancestry.py.")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

BUILDER_VERSION = "3.0.0"


def build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the panel builder.

    Exposed as a function, like build_full_reference.build_parser, so the test
    suite can assert every documented flag parses without running a build or
    touching the network.
    """
    parser = argparse.ArgumentParser(
        prog="build_panel.py",
        description=("Build the 1000 Genomes plus public-tier SGDP reference "
                     "panel into ~/.dnainsight/panels/. Prints its download "
                     "plan and licences and refuses to fetch anything without "
                     "--accept-terms."),
        epilog=("Without --accept-terms this is a dry run: it prints the plan, "
                "the URLs, the licences and the estimated size, and downloads "
                "nothing."),
    )
    parser.add_argument("--accept-terms", action="store_true",
                        help="accept the printed licences and allow downloads. "
                             "Without this flag the run is a dry run.")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan and exit, even with --accept-terms")
    parser.add_argument("--array-file", default=None,
                        help="path to a raw array export, or a directory of "
                             "them, whose rsIDs restrict the panel to the "
                             "positions your chip actually reads")
    parser.add_argument("--chromosomes", default="auto",
                        help="comma-separated chromosomes, default all "
                             "autosomes. X is available and uses its own "
                             "release filename.")
    parser.add_argument("--onekg-source", default="beagle",
                        choices=list(ONEKG_SOURCES),
                        help="which distribution of phase 3 to stream. "
                             "'beagle' (default) is the same v5a callset with "
                             "rsIDs present; 'igsr' is the canonical release, "
                             "whose ID column is empty, which makes "
                             "--array-file match nothing.")
    parser.add_argument("--statistic", default="fst", choices=list(STATISTICS),
                        help="ancestry-informative marker statistic "
                             "(default fst)")
    parser.add_argument("--markers-per-population", type=int, default=5000,
                        metavar="N",
                        help="how many top-ranked markers to name per "
                             "population in informative_markers.tsv")
    parser.add_argument("--min-maf", type=float, default=0.01, metavar="F",
                        help="drop markers whose panel-wide minor allele "
                             "frequency is below F (default 0.01)")
    parser.add_argument("--skip-sgdp", action="store_true",
                        help="build from 1000 Genomes alone")
    parser.add_argument("--sgdp-dir", default=None,
                        help="read SGDP per-sample VCFs from this local "
                             "directory instead of downloading them")
    parser.add_argument("--sgdp-vcf-template", default=None, metavar="URL",
                        help="override the SGDP per-sample VCF URL template; "
                             "must contain {sample}")
    parser.add_argument("--sgdp-metadata", default=None, metavar="PATH_OR_URL",
                        help="read the SGDP public-tier metadata CSV from a "
                             "local path or an alternative URL. Needed with "
                             "--sgdp-dir when the Harvard host is not "
                             "reachable, which it was not at review time.")
    parser.add_argument("--include-hgdp", action="store_true",
                        help="request HGDP. Refused unless "
                             "--accept-consent-caveat is also passed.")
    parser.add_argument("--accept-consent-caveat", action="store_true",
                        help="acknowledge the HGDP consent objection. This is "
                             "not a licence acceptance and --accept-terms "
                             "does not imply it.")
    parser.add_argument("--panel-id", default=PANEL_ID,
                        help=f"panel directory name (default {PANEL_ID})")
    parser.add_argument("--limit", type=int, default=0, metavar="N",
                        help="stop after N records per chromosome, for smoke "
                             "tests. 0 means no limit.")
    parser.add_argument("--keep-scratch", action="store_true",
                        help="keep the cache and per-sample columns after the "
                             "build, for debugging")
    return parser


def _sgdp_sources(args: argparse.Namespace,
                  metadata: dict[str, tuple[str, str]]) -> list[tuple[str, str, bool]]:
    """Resolve SGDP samples to ``(sample, url_or_path, is_local)`` triples."""
    out: list[tuple[str, str, bool]] = []
    if args.sgdp_dir:
        base = Path(args.sgdp_dir)
        for sample in sorted(metadata):
            for pattern in (f"{sample}*.vcf.gz", f"{sample}*.vcf.bgz"):
                found = sorted(base.glob(pattern))
                if found:
                    out.append((sample, str(found[0]), True))
                    break
        return out
    template = args.sgdp_vcf_template or SGDP_VCF_TEMPLATE
    for sample in sorted(metadata):
        out.append((sample, template.format(sample=sample), False))
    return out


def main(argv: list[str] | None = None) -> int:
    """Build the panel. Returns a process exit code.

    Exit codes: 0 for a completed build or a dry run, 2 for a refused request
    such as --include-hgdp without --accept-consent-caveat, 1 for a build that
    started and failed.
    """
    args = build_parser().parse_args(argv)
    plan = download_plan(args)

    # THE GATE. A dry run is the default, and --accept-terms is the only thing
    # that turns it off. This mirrors backend/snpedia.py, which returns HTTP
    # 403 with the licence notice until accept_license is true.
    dry_run = bool(args.dry_run or not args.accept_terms)
    print_plan(plan, dry_run=dry_run)
    if dry_run:
        return 0

    if args.include_hgdp != args.accept_consent_caveat:
        # One flag without the other. Refuse loudly rather than quietly
        # dropping the request, so a scripted build cannot end up believing it
        # got HGDP when it did not.
        allowed, reason = hgdp_gate(args.include_hgdp, args.accept_consent_caveat)
        if not allowed:
            print(reason)
            print("\nRefusing to continue with a half-answered consent gate.")
            return 2

    _require_requests()
    base = panel_dir(args.panel_id)
    scratch = _scratch_dir(base)
    base.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    chromosomes = _selected_chromosomes(args)
    started = datetime.now(timezone.utc)

    print("array coverage set:")
    coverage = load_array_rsids(args.array_file)
    if coverage:
        print(f"  coverage set size: {len(coverage):,} distinct rsIDs")
    else:
        bang = "!" * 66
        print(f"  {bang}")
        print("  !! NO ARRAY FILE. Building the FULL panel, which is tens of")
        print("  !! gigabytes of download and far larger on disk than any")
        print("  !! consumer array needs. Pass --array-file <your raw export>")
        print("  !! to restrict the panel to positions your chip actually reads.")
        print(f"  {bang}")

    print(f"[1000g] fetching sample panel: {ONEKG_PANEL_URL}")
    sample_pops = parse_panel_file(stream_text_lines(ONEKG_PANEL_URL,
                                                     timeout=HTTP_TIMEOUT))
    print(f"[1000g] {len(sample_pops):,} samples, "
          f"{len({p for p, _s in sample_pops.values()})} populations")

    sgdp_meta: dict[str, tuple[str, str]] = {}
    if not args.skip_sgdp:
        try:
            source = args.sgdp_metadata or SGDP_METADATA_URL
            local = Path(source)
            if local.exists():
                print(f"[sgdp] reading metadata from {local}")
                lines: Iterable[str] = local.read_text(
                    encoding="utf-8", errors="replace").splitlines()
            else:
                print(f"[sgdp] fetching metadata: {source}")
                lines = stream_text_lines(source, timeout=HTTP_TIMEOUT)
            sgdp_meta = parse_sgdp_metadata(lines)
            print(f"[sgdp] {len(sgdp_meta):,} public-tier samples")
        except Exception as exc:
            # A missing SGDP tier degrades the panel; it does not fail the
            # build. The 1000 Genomes part is still a usable panel and the
            # omission is recorded in BUILD.txt rather than hidden.
            print(f"  WARNING: SGDP metadata unavailable "
                  f"({type(exc).__name__}: {exc}). Continuing without SGDP.")
            sgdp_meta = {}

    map_archive: Path | None = scratch / "plink.GRCh37.map.zip"
    try:
        print(f"[map] fetching {GENETIC_MAP_URL}")
        size = download_file(GENETIC_MAP_URL, map_archive)
        print(f"[map] cached {_human_bytes(size)} to {map_archive}")
    except Exception as exc:
        print(f"  WARNING: genetic map unavailable ({type(exc).__name__}: "
              f"{exc}). panel.map will carry 0.0 cM for every marker, which "
              f"Beagle reads as 'no genetic distance information'.")
        map_archive = None

    raw_vcf = scratch / "panel.1000g.vcf.gz"
    out_map = base / "panel.map"
    summary = _build_onekg_pass(
        chromosomes, sample_pops, coverage, raw_vcf, out_map, map_archive,
        args.statistic, args.min_maf, args.markers_per_population, args.limit,
        args.onekg_source)
    print(f"[1000g] {summary['records_seen']:,} records read, "
          f"{summary['markers_kept']:,} markers kept")
    if summary["markers_without_rsid"]:
        print(f"  NOTE: {summary['markers_without_rsid']:,} records carried no "
              f"rsID and are named by position in panel.map. They cannot "
              f"appear in informative_markers.tsv, which joins on rsID.")

    merged_samples = list(summary["samples"])
    all_pops = dict(sample_pops)
    columns: list[tuple[str, Path]] = []
    if sgdp_meta and summary["markers_kept"]:
        print(f"[sgdp] extracting {len(sgdp_meta)} genotype columns")
        for sample, source, is_local in _sgdp_sources(args, sgdp_meta):
            try:
                column = _sgdp_column(source, summary["index"],
                                      summary["markers_kept"], is_local)
            except Exception as exc:
                print(f"  WARNING: SGDP sample {sample} skipped "
                      f"({type(exc).__name__}: {exc})")
                continue
            path = scratch / f"sgdp.{sample}.gt"
            path.write_bytes(bytes(column))
            columns.append((sample, path))
            all_pops[sample] = sgdp_meta[sample]
        print(f"[sgdp] {len(columns)} sample column(s) extracted")

    out_vcf = base / "panel.vcf.gz"
    if columns:
        merged = _merge_columns(raw_vcf, scratch / "panel.merged.vcf.gz", columns)
        shutil.move(str(scratch / "panel.merged.vcf.gz"), str(out_vcf))
        merged_samples += [name for name, _p in columns]
        print(f"[merge] rewrote {merged:,} records with "
              f"{len(columns)} SGDP column(s)")
    else:
        shutil.move(str(raw_vcf), str(out_vcf))

    present = {s: all_pops[s] for s in merged_samples if s in all_pops}
    write_populations_tsv(base / "populations.tsv", present, merged_samples)

    pop_order: list[str] = []
    superpops: dict[str, str] = {}
    for sample in merged_samples:
        pop, superpop = present.get(sample, ("", ""))
        if pop and pop not in pop_order:
            pop_order.append(pop)
            superpops[pop] = superpop
    write_q_columns(base / "q_columns.tsv", pop_order, superpops)

    selections = {code: heap.items()
                  for code, heap in summary["heaps"].items() if len(heap)}
    written = write_informative_markers(base / "informative_markers.tsv",
                                        selections, args.statistic)
    print(f"[aim] {written:,} marker rows across {len(selections)} population(s)")

    hashes = {name: sha256_of(base / name)
              for name in PANEL_OUTPUTS
              if name != "BUILD.txt" and (base / name).exists()}
    licences = {"onekg": LICENCES["onekg"], "genetic_map": LICENCES["genetic_map"]}
    if columns:
        licences["sgdp"] = LICENCES["sgdp"]
    write_build_txt(base / "BUILD.txt", {
        "panel_id": args.panel_id,
        "build": TARGET_BUILD,
        "builder_version": BUILDER_VERSION,
        "retrieved_at": started.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "chromosomes": list(chromosomes),
        "sample_count": len(present),
        "population_count": len(pop_order),
        "marker_count": summary["markers_kept"],
        "statistic": args.statistic,
        "min_maf": args.min_maf,
        "array_filter": args.array_file or "",
        "source_versions": {
            "1000_genomes": (f"phase 3 v5a, release {ONEKG_RELEASE} "
                             f"({TARGET_BUILD}), {args.onekg_source} "
                             f"distribution"),
            "sgdp": (f"public tier, {len(columns)} samples merged"
                     if columns else "not included in this build"),
            "genetic_map": ("plink.GRCh37.map.zip" if map_archive
                            else "unavailable; cM recorded as 0.0"),
        },
        "licences": licences,
        "refusals": plan["refusals"],
        "hashes": hashes,
    })

    if not args.keep_scratch:
        shutil.rmtree(scratch, ignore_errors=True)

    print("-" * 74)
    for name in PANEL_OUTPUTS:
        target = base / name
        state = _human_bytes(target.stat().st_size) if target.exists() else "MISSING"
        print(f"  {name:<26} {state}")
    print(f"WROTE {base}")
    print("Reminder: this directory is outside the repository on purpose and "
          "must never be committed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
