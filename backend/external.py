"""
external.py -- out-of-tree external tool registry, licence gate and runner.

WHY THIS FILE EXISTS
--------------------
DNAInsight is MIT, and `data/DATA_SOURCES.md` records a deliberate rule: only
CC0 and US public domain data is bundled, so the whole repository stays
redistributable and the MIT grant survives for every downstream user.

Several v3.0 capabilities have no acceptable MIT or Apache-2.0 implementation.
The best available tools for global ancestry, phasing, imputation and Y
haplogroup calling are GPL-3.0. Vendoring one into this tree would relicense the
project by copyleft. Reimplementing four published algorithms badly would be
worse than not shipping them at all.

So the third option, and the one this module implements: DNAInsight ships ONLY
the adapter, which is MIT. The tool itself is installed by the user, on explicit
consent, into ``~/.dnainsight/tools/``, deliberately outside this repository
tree. That is precisely the pattern ``backend/snpedia.py`` already uses for
CC-BY-NC-SA content, and reusing it means there is one rule to understand rather
than two.

Three consequences, enforced here in code rather than left to good intentions:

  1. Every capability that needs an external tool degrades to a documented
     "unavailable" payload instead of raising. The UI hides the control; it
     never renders a dead button. This mirrors the existing
     ``/api/capabilities`` contract.

  2. No tool runs until its licence has been accepted. The gate raises
     ``LicenceRequired``, which the route layer turns into HTTP 403 carrying the
     licence notice, the same shape as the SNPedia harvest gate.

  3. Tools whose licence forbids redistribution or commercial use are BLOCKED
     and cannot be installed at all, even on explicit consent. ADMIXTURE,
     RFMix, yhaplo and yallHap are in that list. Naming them here makes the
     exclusion a recorded decision rather than an oversight, exactly as
     ``DATA_SOURCES.md`` section 9 does for PharmGKB.

OFFLINE CONTRACT
----------------
Nothing in this module touches the network at import time, during a scan, or
during any read path. ``install_hint()`` returns instructions; it does not
download. Fetching is the user's own action, performed outside the running
application. The app therefore continues to satisfy "runs fully offline" in the
sense the project has always meant it: the application makes no network calls,
while builders and installers run on explicit user action.

LICENCE STATEMENTS
------------------
The licence strings below were verified at source in August 2026 and each
carries the date. They are recorded here because a licence that nobody re-reads
drifts, and because ``spdx`` being absent is itself information: it means the
project publishes no machine-readable licence and the terms had to be read by
hand.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

__all__ = [
    "ExternalError", "LicenceRequired", "ToolUnavailable", "ToolBlocked",
    "TOOLS", "PANELS", "BLOCKED",
    "tools_root", "registry", "get", "is_blocked", "resolve", "is_available",
    "licence_accepted", "accept_licence", "revoke_licence", "licence_notice",
    "status", "status_all", "capability_report", "unavailable",
    "run", "install_hint", "panel_root", "panel_status", "reset_cache",
]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class ExternalError(Exception):
    """Base class for every failure raised by this module."""


class LicenceRequired(ExternalError):
    """The tool exists but the user has not accepted its licence.

    Carries ``tool_id`` and ``notice`` so the route layer can return 403 with
    the full text rather than a bare message.
    """

    def __init__(self, tool_id: str, notice: str) -> None:
        super().__init__(notice)
        self.tool_id = tool_id
        self.notice = notice


class ToolUnavailable(ExternalError):
    """The tool is not installed. Callers degrade; they do not crash."""

    def __init__(self, tool_id: str, message: str = "") -> None:
        super().__init__(message or f"External tool '{tool_id}' is not installed.")
        self.tool_id = tool_id


class ToolBlocked(ExternalError):
    """The tool is permanently excluded on licence grounds.

    This is not a "not yet installed" state and no amount of consent clears it.
    """

    def __init__(self, tool_id: str, reason: str) -> None:
        super().__init__(reason)
        self.tool_id = tool_id
        self.reason = reason


# ---------------------------------------------------------------------------
# Registry
#
# `composable` describes whether the tool's code could legally live inside this
# MIT tree. Everything marked False is invoked as a separate process across a
# licence boundary and is never imported, linked or vendored.
# ---------------------------------------------------------------------------

TOOLS: dict[str, dict] = {
    # -- Apache-2.0 and MIT: composable in principle, still out-of-tree in
    #    practice because they are large compiled artefacts, not because of
    #    licence. Recorded honestly so the distinction stays visible.
    "flare": {
        "id": "flare",
        "name": "FLARE",
        "purpose": "Local ancestry inference (per-segment ancestry assignment)",
        "capability": "ancestry_local",
        "licence": "Apache License 2.0",
        "spdx": "Apache-2.0",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/browning-lab/flare",
        "binaries": ["flare.jar"],
        "kind": "jar",
        "verified": "2026-08-04",
        "requires": ["java"],
        "notes": (
            "Requires phased study VCF, phased reference VCF, a population map "
            "file and a PLINK-format cM genetic map. Preprocess with Beagle."
        ),
    },
    "hap_ibd": {
        "id": "hap_ibd",
        "name": "hap-ibd",
        "purpose": "IBD segment detection on phased data",
        "capability": "ibd_phased",
        "licence": "Apache License 2.0",
        "spdx": "Apache-2.0",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/browning-lab/hap-ibd",
        "binaries": ["hap-ibd.jar"],
        "kind": "jar",
        "verified": "2026-08-04",
        "requires": ["java"],
        "notes": "Requires phased VCF with no missing alleles plus a cM map.",
    },
    "haplogrep": {
        "id": "haplogrep",
        "name": "HaploGrep 3",
        "purpose": "mtDNA haplogroup classification against PhyloTree",
        "capability": "haplogroup_mt",
        "licence": "MIT License",
        "spdx": "MIT",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/genepi/haplogrep3",
        "binaries": ["haplogrep3", "haplogrep3.exe", "haplogrep3.jar"],
        "kind": "binary",
        "verified": "2026-08-04",
        "requires": ["java"],
        "notes": (
            "Tree catalogue is separate from the binary. Record which tree "
            "version produced a call: a haplogroup is meaningless without it."
        ),
    },
    "cladefinder": {
        "id": "cladefinder",
        "name": "Clade Finder",
        "purpose": "Y clade determination from positive and negative SNP calls",
        "capability": "haplogroup_y_second_opinion",
        "licence": "MIT License",
        "spdx": "MIT",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/hprovyn/clade-finder",
        "binaries": ["cladefinder", "cladefinder.py"],
        "kind": "python",
        "verified": "2026-08-04",
        "requires": [],
        "notes": "Used as an independent second opinion against Yleaf.",
    },
    "shapeit": {
        "id": "shapeit",
        "name": "SHAPEIT5",
        "purpose": "Haplotype phasing",
        "capability": "phasing",
        "licence": "MIT License",
        "spdx": "MIT",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        # The original odelaneau/shapeit5 repository is disabled by GitHub
        # staff. Pinning to it is a live supply-chain hazard, so the live
        # location is recorded here and the dead one is named so nobody
        # "restores" it.
        "homepage": "https://github.com/odelaneau/shapeit",
        "superseded_url": "https://github.com/odelaneau/shapeit5 (disabled by GitHub staff)",
        "binaries": ["phase_common", "phase_rare", "shapeit5"],
        "kind": "binary",
        "verified": "2026-08-04",
        "requires": [],
        "notes": "Alternative to Beagle for the phasing step, and MIT rather than GPL.",
    },

    # -- GPL-3.0: NOT composable. Separate process only, never imported.
    "fastmixture": {
        "id": "fastmixture",
        "name": "fastmixture",
        "purpose": "Global ancestry estimation, supervised and projection modes",
        "capability": "ancestry_global",
        "licence": "GNU General Public License v3.0",
        "spdx": "GPL-3.0-only",
        "composable": False,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/Rosemeis/fastmixture",
        "binaries": ["fastmixture"],
        "kind": "python",
        "verified": "2026-08-04",
        "requires": [],
        "notes": (
            "Run in --projection mode against a fixed reference panel, which is "
            "the same shape as a DIY admixture calculator without inheriting an "
            "unlicensed community model file."
        ),
    },
    "beagle": {
        "id": "beagle",
        "name": "Beagle 5.5",
        "purpose": "Phasing and genotype imputation",
        "capability": "imputation",
        "licence": "GNU General Public License v3.0 or later",
        "spdx": "GPL-3.0-or-later",
        "composable": False,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://faculty.washington.edu/browning/beagle/beagle.html",
        "binaries": ["beagle.jar", "beagle.27Feb25.75f.jar"],
        "kind": "jar",
        "verified": "2026-08-04",
        "requires": ["java"],
        "notes": "Emits DR2 per variant. DR2 is carried through as a first-class field.",
    },
    "yleaf": {
        "id": "yleaf",
        "name": "Yleaf",
        "purpose": "Y chromosome haplogroup calling",
        "capability": "haplogroup_y",
        "licence": "GNU General Public License v3.0",
        "spdx": "GPL-3.0-only",
        "composable": False,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/genid/Yleaf",
        "binaries": ["Yleaf", "yleaf"],
        "kind": "python",
        "verified": "2026-08-04",
        "requires": [],
        "notes": (
            "Chosen over every alternative for one reason: it accepts PLINK and "
            "SNP-array input directly, which is the only format a consumer "
            "array user actually has. --tree selects YFull, ISOGG or FTDNA."
        ),
    },
    "ibis": {
        "id": "ibis",
        "name": "IBIS",
        "purpose": "IBD segment detection without phasing",
        "capability": "ibd_unphased",
        "licence": "GNU General Public License v3.0",
        "spdx": "GPL-3.0-only",
        "composable": False,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://github.com/williamslab/ibis",
        "binaries": ["ibis", "ibis.exe"],
        "kind": "binary",
        "verified": "2026-08-04",
        "requires": [],
        "notes": (
            "Phase-free, so it runs on raw array data with no imputation step. "
            "That is what makes household IBD possible in Wave 3 rather than "
            "Wave 4."
        ),
    },

    "samtools": {
        "id": "samtools",
        "name": "SAMtools",
        "purpose": "Targeted pileup extraction from BAM and CRAM alignments",
        "capability": "sequencing_pileup",
        "licence": "MIT License",
        "spdx": "MIT",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://www.htslib.org/",
        "binaries": ["samtools", "samtools.exe"],
        "kind": "binary",
        "verified": "2026-08-04",
        "requires": [],
        "notes": (
            "Used for targeted extraction at the reference positions DNAInsight "
            "already cares about, never for full variant calling. A 30x BAM is "
            "100 to 200 GB and re-calling it locally would be a different "
            "product."
        ),
    },

    # -- Local model runtime for the grounded assistant.
    "ollama": {
        "id": "ollama",
        "name": "Ollama",
        "purpose": "Local language model runtime for the grounded assistant",
        "capability": "assistant",
        "licence": "MIT License",
        "spdx": "MIT",
        "composable": True,
        "commercial_ok": True,
        "redistributable": True,
        "homepage": "https://ollama.com/",
        "binaries": ["ollama", "ollama.exe"],
        "kind": "binary",
        "verified": "2026-08-04",
        "requires": [],
        "notes": (
            "Reached over loopback only. Genotypes never leave the process; the "
            "assistant sends finding text and citations, never raw calls."
        ),
    },
}


# ---------------------------------------------------------------------------
# Permanently excluded tools.
#
# These are named rather than omitted. An absent entry looks like an oversight;
# a BLOCKED entry with a reason is a decision somebody made on a date.
# ---------------------------------------------------------------------------

BLOCKED: dict[str, dict] = {
    "admixture": {
        "id": "admixture",
        "name": "ADMIXTURE",
        "purpose": "Global ancestry estimation",
        "licence": "Academic use only. No LICENSE file exists in the repository; "
                   "Bioconda labels it 'Free for Academic Use' and the v1.4 "
                   "manual contains no licence section at all.",
        "spdx": None,
        "reason": (
            "Binary-only distribution under academic terms. Shipping it, or "
            "instructing users to install it as part of a redistributable "
            "product, would attach terms this project cannot honour. Use "
            "fastmixture, which is GPL-3.0 and permits commercial use."
        ),
        "replacement": "fastmixture",
        "verified": "2026-08-04",
    },
    "rfmix": {
        "id": "rfmix",
        "name": "RFMix v2",
        "purpose": "Local ancestry inference",
        "licence": "Academic research use only. Commercial users must contact "
                   "Stanford's Office of Technology Licensing.",
        "spdx": None,
        "reason": (
            "Explicit academic-only grant. Use FLARE, which is Apache-2.0 and "
            "carries no such restriction."
        ),
        "replacement": "flare",
        "verified": "2026-08-04",
    },
    "yhaplo": {
        "id": "yhaplo",
        "name": "yhaplo",
        "purpose": "Y haplogroup calling",
        "licence": "Non-commercial use only, per LICENSE.txt in the 23andMe "
                   "repository. ISOGG independently describes it as free for "
                   "non-commercial use.",
        "spdx": None,
        "reason": (
            "A non-commercial term would strip the right to sell anything built "
            "on DNAInsight, which is exactly the outcome the SNPedia exclusion "
            "exists to prevent. Use Yleaf."
        ),
        "replacement": "yleaf",
        "verified": "2026-08-04",
    },
    "yallhap": {
        "id": "yallhap",
        "name": "yallHap",
        "purpose": "Y haplogroup calling",
        "licence": "PolyForm Noncommercial License 1.0.0",
        "spdx": "PolyForm-Noncommercial-1.0.0",
        "reason": "Non-commercial. Same reasoning as yhaplo. Use Yleaf.",
        "replacement": "yleaf",
        "verified": "2026-08-04",
    },
    "diydodecad": {
        "id": "diydodecad",
        "name": "DIYDodecad and the community admixture calculators",
        "purpose": "Admixture estimation using Eurogenes, Dodecad, MDLP and "
                   "HarappaWorld model files",
        "licence": "DIYDodecad is published free of charge for non-commercial "
                   "use. The Eurogenes, MDLP and HarappaWorld model files "
                   "publish no licence at all.",
        "spdx": None,
        "reason": (
            "The stevenliuyi/admix runner is GPL-3.0, but its author states "
            "plainly that the model files are the property of their authors and "
            "are not covered by that licence. Unlicensed is not permissive. "
            "Build the panel from 1000 Genomes and the public SGDP tier "
            "instead, where the terms are actually known."
        ),
        "replacement": "panel:onekg_sgdp",
        "verified": "2026-08-04",
    },
}


# ---------------------------------------------------------------------------
# Reference panels.
#
# Panels are data, not tools, but the same rule applies and for the same reason,
# so they live behind the same gate.
# ---------------------------------------------------------------------------

PANELS: dict[str, dict] = {
    "onekg_sgdp": {
        "id": "onekg_sgdp",
        "name": "1000 Genomes phase 3 plus the public SGDP tier",
        "purpose": "Reference panel for global ancestry, phasing and imputation",
        "licence": (
            "1000 Genomes is open with no restriction on use, distributed via "
            "IGSR under EMBL-EBI terms; citation requested. The 279-sample "
            "public tier of the Simons Genome Diversity Project is "
            "unrestricted."
        ),
        "spdx": None,
        "commercial_ok": True,
        "verified": "2026-08-04",
        "files": ["panel.vcf.gz", "panel.map", "populations.tsv"],
        "excluded": {
            "sgdp_restricted": (
                "The 21-sample restricted SGDP tier sits behind a signed "
                "agreement whose terms include 'I will not use the data for any "
                "commercial purposes'. It is excluded by construction and the "
                "builder refuses to fetch it."
            ),
            "hgdp": (
                "HGDP is legally open but ethically contested. Nature Genetics, "
                "24 November 2025, concluded broad reuse may diverge from what "
                "participants consented to. The PRIMED Consortium voted on "
                "21 August 2024 to keep permitting it while acknowledging "
                "'failure to obtain informed consent consistent with current "
                "standards from many participants'. Available only behind an "
                "explicit second opt-in, never by default."
            ),
            "aadr": (
                "Allen Ancient DNA Resource terms were not readable at review "
                "time and the compendium aggregates datasets each carrying "
                "their own upstream terms. Excluded until read."
            ),
        },
        "note": (
            "Reference panels are ancestry-biased. Non-European ancestries are "
            "under-represented in every openly licensed panel that exists, so "
            "every result derived from this panel carries a coverage figure and "
            "the UI states the limitation rather than hiding it."
        ),
    },
    "hgdp_optional": {
        "id": "hgdp_optional",
        "name": "HGDP-CEPH, optional second opt-in",
        "purpose": "Broader population coverage for global ancestry",
        "licence": "Open access via IGSR under Fort Lauderdale Principles.",
        "spdx": None,
        "commercial_ok": True,
        "verified": "2026-08-04",
        "files": ["hgdp.vcf.gz", "hgdp_populations.tsv"],
        "ethics_gate": True,
        "note": (
            "Requires a separate acknowledgement beyond the licence gate, "
            "because the objection to HGDP is a consent objection and not a "
            "licence objection. Accepting the licence does not answer it."
        ),
    },
}


# ---------------------------------------------------------------------------
# Filesystem layout
#
# Everything lives under ~/.dnainsight/, outside the repository tree, exactly
# like the SNPedia cache. .gitignore carries a second line of defence anyway,
# because one careless `git add -A` is all it takes.
# ---------------------------------------------------------------------------

def _home_root() -> Path:
    override = os.environ.get("DNAINSIGHT_HOME")
    if override:
        return Path(override)
    return Path.home() / ".dnainsight"


def tools_root() -> Path:
    """Directory holding user-installed external tools."""
    return _home_root() / "tools"


def panel_root() -> Path:
    """Directory holding user-built reference panels."""
    return _home_root() / "panels"


def _consent_path() -> Path:
    return _home_root() / "licences_accepted.json"


def _read_consent() -> dict:
    path = _consent_path()
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        # A corrupt consent file must mean "nothing accepted", never
        # "everything accepted". Fail closed.
        return {}


def _write_consent(data: dict) -> None:
    path = _consent_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# Registry access
# ---------------------------------------------------------------------------

def registry() -> list[dict]:
    """Every installable tool, in a stable order."""
    return [dict(TOOLS[k]) for k in sorted(TOOLS)]


def get(tool_id: str) -> dict | None:
    """Registry entry for ``tool_id``, or None if it is not a known tool."""
    return dict(TOOLS[tool_id]) if tool_id in TOOLS else None


def is_blocked(tool_id: str) -> bool:
    return str(tool_id or "").strip().lower() in BLOCKED


def _require_known(tool_id: str) -> dict:
    key = str(tool_id or "").strip().lower()
    if key in BLOCKED:
        entry = BLOCKED[key]
        raise ToolBlocked(key, entry["reason"])
    if key not in TOOLS:
        raise ExternalError(f"Unknown external tool: {tool_id!r}")
    return TOOLS[key]


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

_resolve_cache: dict[str, str | None] = {}


def reset_cache() -> None:
    """Clear the resolution cache. Tests and the install flow call this."""
    _resolve_cache.clear()


def resolve(tool_id: str) -> Path | None:
    """Locate an installed tool.

    Search order, most specific first:
      1. ``DNAINSIGHT_TOOL_<ID>`` environment override, so a user with the tool
         already installed somewhere sensible is not forced to duplicate it.
      2. ``~/.dnainsight/tools/<id>/`` and its ``bin/`` subdirectory.
      3. The system PATH.

    Returns None when the tool is absent. It does NOT raise, because absence is
    the normal state for most users and callers must degrade rather than fail.
    """
    key = str(tool_id or "").strip().lower()
    if key in _resolve_cache:
        cached = _resolve_cache[key]
        return Path(cached) if cached else None

    entry = TOOLS.get(key)
    if entry is None:
        _resolve_cache[key] = None
        return None

    found: Path | None = None

    override = os.environ.get(f"DNAINSIGHT_TOOL_{key.upper()}")
    if override:
        candidate = Path(override)
        if candidate.exists():
            found = candidate

    if found is None:
        base = tools_root() / key
        for name in entry.get("binaries", []):
            for candidate in (base / name, base / "bin" / name):
                if candidate.exists():
                    found = candidate
                    break
            if found is not None:
                break

    if found is None:
        for name in entry.get("binaries", []):
            if name.endswith(".jar") or name.endswith(".py"):
                continue
            which = shutil.which(name)
            if which:
                found = Path(which)
                break

    _resolve_cache[key] = str(found) if found else None
    return found


def _requirements_met(entry: dict) -> tuple[bool, list[str]]:
    """Check interpreter or runtime prerequisites such as Java."""
    missing = [req for req in entry.get("requires", []) if not shutil.which(req)]
    return (not missing), missing


def is_available(tool_id: str) -> bool:
    """True when the tool is installed, its runtime exists and its licence is accepted.

    All three conditions matter. An installed tool whose licence has not been
    accepted is deliberately reported unavailable, so the user is never in a
    state where DNAInsight silently used something they did not agree to.
    """
    key = str(tool_id or "").strip().lower()
    if key in BLOCKED or key not in TOOLS:
        return False
    if resolve(key) is None:
        return False
    ok, _missing = _requirements_met(TOOLS[key])
    if not ok:
        return False
    return licence_accepted(key)


# ---------------------------------------------------------------------------
# Licence gate
# ---------------------------------------------------------------------------

def licence_notice(tool_id: str) -> str:
    """Full text shown to the user before anything is installed or run."""
    key = str(tool_id or "").strip().lower()
    if key in BLOCKED:
        b = BLOCKED[key]
        return (
            f"{b['name']} is permanently excluded from DNAInsight.\n\n"
            f"Licence: {b['licence']}\n\n"
            f"Reason: {b['reason']}\n\n"
            f"Use instead: {b.get('replacement') or 'no replacement recorded'}."
        )
    entry = TOOLS.get(key)
    if entry is None:
        return f"Unknown external tool: {tool_id!r}"
    lines = [
        f"{entry['name']} is a separate program, not part of DNAInsight.",
        "",
        f"Purpose: {entry['purpose']}",
        f"Licence: {entry['licence']}"
        + (f" ({entry['spdx']})" if entry.get("spdx") else ""),
        f"Homepage: {entry['homepage']}",
        f"Licence verified: {entry['verified']}",
        "",
    ]
    if not entry.get("composable", True):
        lines += [
            "This tool is copyleft. DNAInsight is MIT and does not include, "
            "link to, or vendor it. It is installed into your own home "
            "directory and invoked as a separate process, so its licence "
            "attaches to your copy of that program and not to DNAInsight.",
            "",
        ]
    lines += [
        "Obligations under this licence fall on you, the person installing it.",
        "DNAInsight ships only the adapter code, which is MIT.",
    ]
    if entry.get("notes"):
        lines += ["", f"Note: {entry['notes']}"]
    if entry.get("superseded_url"):
        lines += ["", f"Do not use: {entry['superseded_url']}"]
    return "\n".join(lines)


def licence_accepted(tool_id: str) -> bool:
    key = str(tool_id or "").strip().lower()
    record = _read_consent().get(key)
    return bool(record and record.get("accepted"))


def accept_licence(tool_id: str, *, accept: bool = False) -> dict:
    """Record acceptance of a tool's licence.

    ``accept`` must be passed explicitly as True. A default-True parameter here
    would mean a stray call silently grants consent, and consent is the one
    thing in this module that must never happen by accident.
    """
    entry = _require_known(tool_id)          # raises ToolBlocked where relevant
    key = entry["id"]
    if not accept:
        raise LicenceRequired(key, licence_notice(key))
    data = _read_consent()
    data[key] = {
        "accepted": True,
        "licence": entry["licence"],
        "spdx": entry.get("spdx"),
        "accepted_at": datetime.now(timezone.utc).isoformat(),
        "dnainsight_recorded_version": entry.get("verified"),
    }
    _write_consent(data)
    reset_cache()
    return dict(data[key])


def revoke_licence(tool_id: str) -> bool:
    """Withdraw acceptance. The tool immediately reports unavailable again."""
    key = str(tool_id or "").strip().lower()
    data = _read_consent()
    if key not in data:
        return False
    del data[key]
    _write_consent(data)
    reset_cache()
    return True


# ---------------------------------------------------------------------------
# Status reporting
# ---------------------------------------------------------------------------

def status(tool_id: str) -> dict:
    """Full state of one tool. Never raises for a blocked or unknown id."""
    key = str(tool_id or "").strip().lower()

    if key in BLOCKED:
        b = BLOCKED[key]
        return {
            "id": key, "name": b["name"], "state": "blocked",
            "available": False, "installed": False, "licence_accepted": False,
            "licence": b["licence"], "spdx": b.get("spdx"),
            "reason": b["reason"], "replacement": b.get("replacement"),
            "capability": None,
        }

    entry = TOOLS.get(key)
    if entry is None:
        return {"id": key, "state": "unknown", "available": False,
                "installed": False, "licence_accepted": False}

    path = resolve(key)
    runtime_ok, missing = _requirements_met(entry)
    accepted = licence_accepted(key)

    if path is None:
        state = "not_installed"
    elif not runtime_ok:
        state = "runtime_missing"
    elif not accepted:
        state = "licence_not_accepted"
    else:
        state = "ready"

    return {
        "id": key,
        "name": entry["name"],
        "purpose": entry["purpose"],
        "capability": entry["capability"],
        "state": state,
        "available": state == "ready",
        "installed": path is not None,
        "path": str(path) if path else None,
        "licence": entry["licence"],
        "spdx": entry.get("spdx"),
        "composable_with_mit": entry.get("composable", True),
        "licence_accepted": accepted,
        "missing_runtime": missing,
        "homepage": entry["homepage"],
        "licence_verified": entry["verified"],
    }


def status_all() -> dict:
    """State of every tool and every permanently excluded tool."""
    return {
        "tools": [status(k) for k in sorted(TOOLS)],
        "blocked": [status(k) for k in sorted(BLOCKED)],
        "panels": [panel_status(k) for k in sorted(PANELS)],
        "tools_root": str(tools_root()),
        "panel_root": str(panel_root()),
        "policy": (
            "DNAInsight is MIT and bundles only CC0 and US public domain data. "
            "External tools are installed into your own home directory, outside "
            "the repository tree, on explicit consent. Tools whose licence "
            "forbids redistribution or commercial use are permanently excluded "
            "and cannot be installed."
        ),
        "offline": (
            "No tool is contacted over the network by DNAInsight. Installation "
            "is your own action, performed outside the running application."
        ),
    }


def panel_status(panel_id: str) -> dict:
    """Whether a reference panel has been built, and under what terms."""
    key = str(panel_id or "").strip().lower()
    entry = PANELS.get(key)
    if entry is None:
        return {"id": key, "state": "unknown", "available": False}
    base = panel_root() / key
    present = [f for f in entry.get("files", []) if (base / f).exists()]
    complete = len(present) == len(entry.get("files", []))
    return {
        "id": key,
        "name": entry["name"],
        "purpose": entry["purpose"],
        "state": "ready" if complete else ("partial" if present else "not_built"),
        "available": complete,
        "path": str(base),
        "files_present": present,
        "files_expected": list(entry.get("files", [])),
        "licence": entry["licence"],
        "commercial_ok": entry.get("commercial_ok"),
        "ethics_gate": entry.get("ethics_gate", False),
        "excluded": entry.get("excluded", {}),
        "note": entry.get("note", ""),
        "licence_verified": entry["verified"],
    }


def capability_report() -> dict:
    """Capability flags for ``/api/capabilities``.

    The UI hides a control whose capability is False. That is the whole point:
    a user without Beagle installed should not see an imputation slider that
    does nothing, and should not be told imputation "failed" when it was simply
    never possible.
    """
    out: dict[str, bool] = {}
    for key, entry in TOOLS.items():
        out[entry["capability"]] = is_available(key)
    out["panel_onekg_sgdp"] = panel_status("onekg_sgdp")["available"]
    out["panel_hgdp"] = panel_status("hgdp_optional")["available"]
    return out


def unavailable(tool_id: str, capability: str, *, detail: str = "") -> dict:
    """The standard degraded payload.

    Every adapter returns this instead of raising when its tool is absent. The
    shape is fixed so the frontend renders one honest empty state rather than
    ten different error strings, and so a caller can always distinguish "we
    looked and found nothing" from "we could not look at all". That distinction
    is the same one the genoset engine already draws between unmatched and
    not testable.
    """
    st = status(tool_id)
    reasons = {
        "not_installed": (
            f"{st.get('name', tool_id)} is not installed. This analysis was not "
            f"attempted, which is different from finding nothing."
        ),
        "runtime_missing": (
            f"{st.get('name', tool_id)} is installed but requires "
            f"{', '.join(st.get('missing_runtime') or ['a runtime'])}, which was "
            f"not found."
        ),
        "licence_not_accepted": (
            f"{st.get('name', tool_id)} is installed but its licence has not "
            f"been accepted, so it was not run."
        ),
        "blocked": st.get("reason", "This tool is permanently excluded."),
        "unknown": f"Unknown external tool: {tool_id!r}",
    }
    return {
        "available": False,
        "capability": capability,
        "tool": st.get("name", tool_id),
        "tool_id": st.get("id", tool_id),
        "state": st.get("state"),
        "reason": detail or reasons.get(st.get("state"), "Unavailable."),
        "not_attempted": True,
        "results": [],
        "how_to_enable": install_hint(tool_id) if st.get("state") != "blocked" else None,
    }


def install_hint(tool_id: str) -> dict | None:
    """Instructions for the user. This function does not download anything."""
    key = str(tool_id or "").strip().lower()
    if key in BLOCKED:
        return None
    entry = TOOLS.get(key)
    if entry is None:
        return None
    target = tools_root() / key
    return {
        "tool": entry["name"],
        "homepage": entry["homepage"],
        "install_to": str(target),
        "expected_files": list(entry.get("binaries", [])),
        "requires": list(entry.get("requires", [])),
        "steps": [
            f"Download {entry['name']} from {entry['homepage']}.",
            f"Place the executable in {target} (create the folder if needed).",
            "Accept the licence in DNAInsight under Settings, External Tools, "
            "or POST accept_license true to the tools endpoint.",
        ],
        "licence": entry["licence"],
        "note": entry.get("notes", ""),
        "alternative": (
            f"Or set DNAINSIGHT_TOOL_{key.upper()} to an existing install path."
        ),
    }


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

# A hard ceiling so a runaway external process cannot hang a scan forever.
DEFAULT_TIMEOUT = 900


def run(tool_id: str,
        args: Sequence[str],
        *,
        timeout: int = DEFAULT_TIMEOUT,
        cwd: str | Path | None = None,
        env: dict | None = None,
        check: bool = True) -> subprocess.CompletedProcess:
    """Invoke an external tool as a separate process.

    The subprocess boundary is not a convenience. It is the licence boundary:
    DNAInsight never imports, links or vendors a GPL-3.0 tool, it executes one
    the user installed themselves. Keeping every invocation in this one function
    means that boundary is auditable in a single place.

    Raises ToolBlocked, LicenceRequired or ToolUnavailable. Callers that want to
    degrade rather than fail should check ``is_available`` first and return
    ``unavailable(...)``.
    """
    entry = _require_known(tool_id)
    key = entry["id"]

    if not licence_accepted(key):
        raise LicenceRequired(key, licence_notice(key))

    path = resolve(key)
    if path is None:
        raise ToolUnavailable(key)

    runtime_ok, missing = _requirements_met(entry)
    if not runtime_ok:
        raise ToolUnavailable(
            key, f"{entry['name']} requires {', '.join(missing)}, which was not found."
        )

    kind = entry.get("kind", "binary")
    if kind == "jar":
        cmd = ["java", "-jar", str(path), *[str(a) for a in args]]
    elif kind == "python":
        cmd = [sys.executable, str(path), *[str(a) for a in args]]
    else:
        cmd = [str(path), *[str(a) for a in args]]

    run_env = dict(os.environ)
    if env:
        run_env.update(env)

    try:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(cwd) if cwd else None,
            env=run_env,
            check=check,
        )
    except subprocess.TimeoutExpired as exc:
        raise ExternalError(
            f"{entry['name']} exceeded the {timeout}s timeout and was stopped."
        ) from exc
    except subprocess.CalledProcessError as exc:
        tail = (exc.stderr or "").strip().splitlines()[-5:]
        raise ExternalError(
            f"{entry['name']} exited with code {exc.returncode}. "
            + (" | ".join(tail) if tail else "No stderr output.")
        ) from exc


def guard(tool_id: str, capability: str):
    """Return None when the tool is ready, or the degraded payload when it is not.

    Adapters open with::

        blocked = external.guard("beagle", "imputation")
        if blocked is not None:
            return blocked

    which keeps the degradation contract identical across all ten modules.
    """
    if is_available(tool_id):
        return None
    return unavailable(tool_id, capability)
