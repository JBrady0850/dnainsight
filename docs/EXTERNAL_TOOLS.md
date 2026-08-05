---
created: 2026-08-04
modified: 2026-08-04
tags: [archivist, dnainsight, reference, decision]
aliases: [DNAInsight External Tools, EXTERNAL_TOOLS]
---

# DNAInsight External Tools

The licence architecture, and how to install a tool.

Referenced from `README.md`, `CHANGELOG.md` and `data/DATA_SOURCES.md`. The
registry this document describes is `backend/external.py`, and
`data/tools_manifest.json` is its published mirror, checked against it by a test
so the two cannot drift.

Last verified at source: 2026-08-04. Every licence line below carries that date
in the registry, because a licence nobody re-reads drifts.

## 1. The problem

DNAInsight is MIT, and `data/DATA_SOURCES.md` records the rule that follows from
that: only CC0 and US public domain data is bundled, so the repository stays
redistributable and the MIT grant survives for every downstream user, including
anyone who sells a service built on it.

Six of the ten v3.0 capabilities have no acceptable MIT or Apache-2.0
implementation. **The best available tools for global ancestry, imputation,
phasing and Y haplogroup calling are GPL-3.0**, and several of the obvious
alternatives are academic-only or non-commercial.

Three options existed. Two of them are wrong.

| Option | Consequence |
|---|---|
| Vendor the GPL-3.0 binary | Copyleft relicenses DNAInsight. Every downstream user loses the MIT grant, silently |
| Reimplement four published algorithms | A bad reimplementation of fastmixture is worse than not shipping ancestry at all |
| Ship only the adapter | The tool is a separate program the user installed. Its licence attaches to their copy of it |

## 2. The rule

**DNAInsight ships only the adapter, which is MIT. The tool is installed by the
user, on explicit consent, into `~/.dnainsight/tools/`, deliberately outside
this repository tree.**

**The subprocess boundary IS the licence boundary.** No external tool is
imported, linked or vendored. Every invocation goes through the single
`external.run` function, so the boundary is auditable in one place rather than
scattered across ten adapters.

This is not a new pattern. It is exactly what `backend/snpedia.py` already does
for CC-BY-NC-SA content: the data lives in `~/.dnainsight/`, the gate refuses
until the user accepts, and the repository never holds the thing it cannot hold.
Reusing it means there is one rule to understand rather than two.

Three consequences, enforced in code rather than left to good intentions.

1. **A capability that needs a missing tool degrades, it never raises.** It
   returns `available: false` with `not_attempted: true` and a plain-English
   reason. The UI hides the control; it never renders a dead button. See
   `docs/API_V3.md` section 1.
2. **No tool runs until its licence has been accepted.** The gate raises
   `LicenceRequired`, which the route layer turns into HTTP 403 carrying the
   full notice. Same shape as the SNPedia harvest gate.
3. **A tool whose licence forbids redistribution or commercial use is BLOCKED
   and cannot be installed at all, even on explicit consent.** Section 5.

## 3. Two tiers, recorded honestly

`composable` in the registry records whether a tool's code could **legally** live
inside this MIT tree. It is separate from whether it actually does.

### Composable in principle, Apache-2.0 or MIT

These carry no copyleft. They are still out-of-tree because they are large
compiled artefacts, and because one install rule is better than two. Recording
the distinction keeps it visible rather than letting "we did not vendor it" turn
into "we could not have".

| Tool | Licence | Capability | Runtime |
|---|---|---|---|
| FLARE | Apache-2.0 | local ancestry, `ancestry_local` | Java |
| hap-ibd | Apache-2.0 | phased IBD, `ibd_phased` | Java |
| HaploGrep 3 | MIT | mtDNA haplogroup, `haplogroup_mt` | Java |
| Clade Finder | MIT | Y second opinion, `haplogroup_y_second_opinion` | Python |
| SHAPEIT5 | MIT | phasing, `phasing` | none |
| SAMtools | MIT | BAM and CRAM pileup, `sequencing_pileup` | none |
| Ollama | MIT | assistant runtime, `assistant` | none |

**SHAPEIT5 supply-chain note.** The original repository
`odelaneau/shapeit5` is **disabled by GitHub staff**. Pinning to it is a live
hazard. The live location is `https://github.com/odelaneau/shapeit`. The dead URL
is recorded in the registry under `superseded_url` and printed in the licence
notice, so nobody "restores" it.

### Subprocess only, GPL-3.0

Never imported, never linked, never vendored.

| Tool | Licence | Capability | Runtime |
|---|---|---|---|
| fastmixture | GPL-3.0-only | global ancestry, `ancestry_global` | Python |
| Beagle 5.5 | GPL-3.0-or-later | imputation, `imputation` | Java |
| Yleaf | GPL-3.0-only | Y haplogroup, `haplogroup_y` | Python |
| IBIS | GPL-3.0-only | unphased IBD, `ibd_unphased` | none |

All four permit commercial use and redistribution. Copyleft is the only reason
they are out of tree, and the licence notice says so before the user accepts.

Why these four specifically:

- **fastmixture** replaces ADMIXTURE. Run in `--projection` mode against a fixed
  reference panel it is the same shape as a DIY admixture calculator, without
  inheriting an unlicensed community model file.
- **Yleaf** was chosen over every alternative for one reason: **it accepts PLINK
  and SNP-array input directly**, which is the only format a consumer array user
  actually has.
- **IBIS** is phase-free, so it runs on raw array data with no imputation step.
  That is what makes household IBD a Wave 3 feature rather than a Wave 4 one.
- **Beagle** emits DR2 per variant, which DNAInsight carries through as a
  first-class field rather than a footnote.

## 4. What breaks without each tool

Nothing below is claimed to work on a machine that does not have the tool
installed. This table is the honest state of a fresh install, which is: none of
them.

| Capability | Needs | Without it |
|---|---|---|
| Global ancestry | fastmixture **and** a built panel | not attempted, `problem` distinguishes which is missing |
| Local ancestry, chromosome painting | FLARE **and** a built panel **and** phased input | not attempted; unphased input is caught and named, not silently painted |
| Imputation | Beagle **and** a built panel | not attempted |
| Phasing | SHAPEIT5 or Beagle | not attempted |
| Y haplogroup depth | Yleaf, optionally Clade Finder for a second opinion | the bundled 49-marker backbone call only, flagged provisional |
| mtDNA haplogroup depth | HaploGrep 3 | the bundled 28-node backbone call only, flagged provisional |
| Phased IBD | hap-ibd | unphased IBD still works |
| Unphased IBD | nothing. Pure Python. IBIS optional | works |
| BAM and CRAM ingest | samtools | VCF and gVCF ingest still works |
| Assistant | Ollama and a local model | refuses |

Everything v1 and v2 did still works with no external tool at all.

## 5. Permanently blocked

These are named rather than omitted. **An absent entry looks like an oversight; a
BLOCKED entry with a reason is a decision somebody made on a date.** No consent
clears them, and `POST /api/v3/tools/<id>/licence` returns 409 rather than 403,
because 403 would imply a body exists that would work.

| Tool | Licence problem | Use instead |
|---|---|---|
| **ADMIXTURE** | Academic use only. **No LICENSE file exists in the repository at all.** Bioconda labels it "Free for Academic Use" and the v1.4 manual has no licence section | fastmixture |
| **RFMix v2** | Academic research use only. Commercial users must contact Stanford's Office of Technology Licensing | FLARE |
| **yhaplo** | Non-commercial only, per LICENSE.txt in the 23andMe repository. ISOGG independently describes it as free for non-commercial use | Yleaf |
| **yallHap** | PolyForm Noncommercial 1.0.0 | Yleaf |
| **DIYDodecad** and the Eurogenes, Dodecad, MDLP and HarappaWorld model files | DIYDodecad is free for non-commercial use. **The model files publish no licence at all** | build the panel from 1000 Genomes and the public SGDP tier |

Two of these deserve their reasoning spelled out.

**yhaplo and yallHap.** A non-commercial term would strip the right to sell
anything built on DNAInsight. That is exactly the outcome the SNPedia exclusion
exists to prevent, and it does not become acceptable because the artefact is a
program rather than a dataset.

**The community admixture calculators.** The `stevenliuyi/admix` runner is
GPL-3.0, which looks like it settles the question. It does not: the author states
plainly that the model files are the property of their authors and are **not**
covered by that licence. **Unlicensed is not permissive.** The alternative is a
panel built from 1000 Genomes and the public SGDP tier, where the terms are
actually known.

## 6. Reference panels

Panels are data, not tools, but the same rule applies for the same reason, so
they live behind the same gate. Built by `data/build_panel.py` into
`~/.dnainsight/panels/`.

### The clean panel: `onekg_sgdp`

**1000 Genomes phase 3 plus the 279-sample PUBLIC tier of SGDP.** 1000 Genomes is
open with no restriction on use. The SGDP public tier is unrestricted. Both
permit commercial use.

Files: `panel.vcf.gz`, `panel.map`, `populations.tsv`, plus
`informative_markers.tsv` and `q_columns.tsv`.

### Excluded, with reasons

**The 21-sample restricted SGDP tier** sits behind a signed agreement whose terms
include "I will not use the data for any commercial purposes". **It is excluded
by construction.** There is no flag that fetches it and no consent path that
enables it. A non-commercial term inside a bundled panel is the same failure as
a non-commercial term inside `data/`.

**HGDP is legally open and gated anyway.** The licence is not the problem. Two
published events are, and both are recorded with dates so the decision can be
re-examined rather than inherited:

- **Nature Genetics, 24 November 2025**, concluded that broad reuse may diverge
  from what participants consented to.
- **The PRIMED Consortium voted on 21 August 2024** to keep permitting its use,
  while acknowledging "failure to obtain informed consent consistent with current
  standards from many participants".

They point in opposite directions and a reader is entitled to both. HGDP sits
behind a **second opt-in** that is separate from the licence gate:
`--include-hgdp` is refused unless `--accept-consent-caveat` is also passed, and
`--accept-terms` does not imply it. **A consent objection is not answered by
accepting a licence**, so the two acknowledgements are not allowed to be the same
click.

**The Allen Ancient DNA Resource is excluded entirely.** Its terms were not
readable at review time and the compendium aggregates datasets each carrying
their own upstream terms. An unread licence is not a permissive licence.

### The panel limitation, stated rather than hidden

**Reference panels are ancestry-biased.** Non-European ancestries are
under-represented in every openly licensed panel that exists. Every result
derived from this panel therefore carries a per-population coverage figure, and
a population the panel cannot resolve is reported NOT RESOLVABLE rather than as
zero percent.

## 7. The offline contract

**The running application makes zero network calls.** Nothing in
`backend/external.py` touches the network at import time, during a scan, or on
any read path. `install_hint()` returns instructions; it does not download.

Builders and tool installers may download, but only when the user runs them, and
**both builders dry-run by default**:

```
python data/build_panel.py                 # prints the plan and licences, downloads nothing
python data/build_pgx_alleles.py           # prints the endpoints and the columns it refuses to transfer
```

Neither fetches anything until `--accept-terms` is passed, and `--dry-run`
overrides `--accept-terms` so a plan can always be re-printed safely.

## 8. Installing a tool

### 8.1 Where things go

```
~/.dnainsight/
├── tools/
│   ├── beagle/beagle.jar
│   ├── fastmixture/fastmixture
│   ├── yleaf/Yleaf
│   └── ...
├── panels/
│   ├── onekg_sgdp/panel.vcf.gz
│   └── chains/hg38ToHg19.over.chain.gz
└── licences_accepted.json
```

Nothing goes inside the repository. `.gitignore` carries a second line of
defence anyway, because one careless `git add -A` is all it takes.

### 8.2 Resolution order

`external.resolve()` looks in three places, most specific first:

1. **`DNAINSIGHT_TOOL_<ID>`**, an environment variable pointing at an existing
   install, so a user who already has samtools somewhere sensible is not forced
   to duplicate it. For example `DNAINSIGHT_TOOL_SAMTOOLS=/usr/bin/samtools`.
2. `~/.dnainsight/tools/<id>/` and its `bin/` subdirectory.
3. The system PATH, for executables only. `.jar` and `.py` artefacts are never
   resolved from PATH, because a jar on PATH is not a runnable command.

Absence returns None rather than raising. Absence is the normal state for most
users and every caller degrades rather than fails.

`DNAINSIGHT_HOME` relocates the whole tree.

### 8.3 The three steps

```
# 1. See what a tool is, what it needs, and its full licence notice
curl http://127.0.0.1:5050/api/v3/tools/beagle

# 2. Download it yourself from the homepage in that response, and place it at
#    the "install_to" path, creating the folder if needed.

# 3. Accept the licence
curl -X POST http://127.0.0.1:5050/api/v3/tools/beagle/licence \
     -H 'Content-Type: application/json' \
     -d '{"accept_license": true}'
```

Step 3 is refused with **403 and the full notice** until `accept_license` is
literally true. In the Python API, `accept_licence(tool_id, accept=True)` takes
its flag as a keyword defaulting to False: a default-True parameter would mean a
stray call silently grants consent, and consent is the one thing here that must
never happen by accident.

`DELETE` on the same path revokes. The capability reports unavailable
immediately.

Acceptance is recorded in `~/.dnainsight/licences_accepted.json` with the
licence text, the SPDX identifier and a UTC timestamp. **A corrupt consent file
means "nothing accepted", never "everything accepted".** It fails closed.

### 8.4 Per-tool notes

| Tool | Install target | Expect | Note |
|---|---|---|---|
| Beagle 5.5 | `tools/beagle/` | `beagle.jar` | needs Java. Emits DR2 per variant |
| FLARE | `tools/flare/` | `flare.jar` | needs Java, phased study VCF, phased reference VCF, a population map and a PLINK cM map. Preprocess with Beagle |
| hap-ibd | `tools/hap_ibd/` | `hap-ibd.jar` | needs Java. Requires phased VCF with no missing alleles plus a cM map |
| HaploGrep 3 | `tools/haplogrep/` | `haplogrep3` or `haplogrep3.jar` | needs Java. **The tree catalogue is separate from the binary.** Record which tree version produced a call |
| Clade Finder | `tools/cladefinder/` | `cladefinder.py` | used as an independent second opinion against Yleaf |
| SHAPEIT5 | `tools/shapeit/` | `phase_common`, `phase_rare` | use `odelaneau/shapeit`, **not** the disabled `shapeit5` repository |
| fastmixture | `tools/fastmixture/` | `fastmixture` | run in `--projection` mode against a fixed panel |
| Yleaf | `tools/yleaf/` | `Yleaf` | `--tree` selects YFull, ISOGG or FTDNA |
| IBIS | `tools/ibis/` | `ibis` | see the known gap in section 10 |
| SAMtools | `tools/samtools/` | `samtools` | usually already on PATH |
| Ollama | `tools/ollama/` | `ollama` | usually already on PATH. Reached over loopback only |

Every external process runs under a **900 second timeout**, so a runaway tool
cannot hang a scan forever. A timeout and a non-zero exit are both reported with
the stderr tail rather than as a bare failure.

### 8.5 Building the panel

```
python data/build_panel.py                                 # dry run, prints the plan
python data/build_panel.py --accept-terms --array-file my_raw_dna.txt
```

`--array-file` restricts the panel to the positions your chip actually reads,
which is the difference between a panel you can build on a laptop and one you
cannot.

`--onekg-source` defaults to `beagle`, and that default exists for a verified
reason: the IGSR 20130502 release genotype VCFs publish `ID` as `.` for every
record, so an rsID-keyed join against them matches nothing. `igsr` remains
selectable and prints the caveat.

HGDP needs both flags:

```
python data/build_panel.py --accept-terms --include-hgdp --accept-consent-caveat
```

## 9. Verifying the position

```
GET /api/v3/tools           every tool, every blocked entry, every panel
GET /api/v3/capabilities    the three axes: subsystems, external, panels
GET /api/v3/licence-audit   the bundling rule, checked at runtime
```

`/api/v3/licence-audit` is the runtime half of `data/DATA_SOURCES.md`. A
non-empty `violations` list means something non-redistributable reached a
bundled artefact. It is what stops that document from being prose somebody has
to remember.

## 10. Known gaps in this layer

**No tool's command-line arguments have been executed against an installed
binary.** Everything below is an assumption from documentation. Full detail is
in `docs/KNOWN_GAPS.md`; the summary is here so nobody reads this document alone
and takes it for tested.

- Beagle's `gt= ref= map= out= chrom= nthreads=` argument form, its
  `<out>.vcf.gz` output path, and `DR2=`, `AF=` and `IMP` in INFO.
- fastmixture's `--bfile --out --seed --projection --reference` flags and its
  `ancestry*.Q` output. **Without `q_columns.tsv` the columns are reported as
  `component_N` with `components_labelled: false`, because a numbered component
  is not a population.**
- **FLARE's per-haplotype ancestry is assumed to arrive as `AN1`/`AN2` FORMAT
  fields holding integer indices into the reference panel population order.
  This is the weakest assumption in the release.**
- CLI flags for Yleaf, HaploGrep, Clade Finder and IBIS are isolated in
  `YLEAF_ARGS`, `HAPLOGREP_ARGS`, `CLADEFINDER_ARGS` and `IBIS_ARGS`, so
  correcting one is a one-line change in a known place.
- **IBIS writes `.ped`/`.map`, but `IBIS_ARGS` references bed/bim/fam.** A PLINK
  conversion step, or a switch to IBIS text-input mode, is still needed.

## 11. Adding a tool

Two rules, and they are not negotiable.

1. **It must permit commercial use and redistribution.** If either answer is no,
   it goes in `BLOCKED` with a reason, a verification date and a named
   replacement. Not omitted. An omission reads as an oversight.
2. **It is invoked through `external.run` and nowhere else.** That one function
   is the licence boundary and it has to stay auditable in one place.

Then:

- Add the registry entry with `licence`, `spdx` (or `None`, which is itself
  information: it means the project publishes no machine-readable licence and
  the terms had to be read by hand), `composable`, `commercial_ok`,
  `redistributable`, `homepage`, `binaries`, `kind`, `requires` and `verified`.
- Give it a **capability name that does not collide with a DNAInsight module
  name**. Capabilities are merged into `available_subsystems()` under a `tool_`
  prefix precisely because `imputation` and `assistant` already collided once,
  and a user without Beagle was told DNAInsight's own imputation module was
  missing.
- Open the adapter with the standard guard, so degradation is identical across
  all ten modules:

```python
blocked = external.guard("beagle", "imputation")
if blocked is not None:
    return blocked
```

- Regenerate `data/tools_manifest.json`. A test compares it against
  `backend/external.py` field by field and will fail if you do not.
