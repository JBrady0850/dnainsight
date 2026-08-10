DNAInsight v3.3.0

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Bundled SNPs](https://img.shields.io/badge/Bundled_SNPs-122_curated-orange)
![Genosets](https://img.shields.io/badge/Genosets-65-orange)
![Tests](https://img.shields.io/badge/Tests-3425_passing-brightgreen)
[![Buy Me a Coffee](https://img.shields.io/badge/Buy_Me_a_Coffee-support_this_project-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/jbrady2852)

![DNAInsight dashboard showing a sample profile with findings, risk categories and next steps](DNAInsight.png)

Personal DNA analysis that runs entirely on your own computer. Read your raw DNA
file from any major provider, annotate it against curated clinical evidence,
explore the results with real filters, and generate reports you can print, keep
or hand to a clinician.

Nothing leaves your machine unless you explicitly ask it to.

---

What is new in v3.0

v2.0 told you what a variant means for you and ranked it. v3.0 tells you when
that meaning CHANGES, where your data runs out, and what it could not check.

| Capability | v2.0 | v3.0 |
|---|---|---|
| Reclassification | none, every scan stood alone | change ledger, 15 change kinds, dated additive addenda |
| Provenance | none | signed reproducible manifests, plus a runtime audit of the bundling rule |
| Input formats | consumer array exports | plus VCF, gVCF, BAM and CRAM, with the genome build detected and mismatches refused |
| Ancestry | none | global proportions, local ancestry, chromosome painting |
| Haplogroups | none | bundled Y and mtDNA backbones, both naming systems reported, with optional tools for depth |
| Relatives | trio checks on loaded files | IBD across loaded kits, relationship ranges, parental phasing, chromosome browser |
| Imputation | none | DR2 as a first-class field, imputed calls structurally capped below typed ones |
| Pharmacogenomics | CPIC level per variant | star-allele diplotypes for 9 genes, Indeterminate by default, prescription guard |
| Carrier screening | risk allele per variant | 11-gene panel with residual risk arithmetic |
| Assistant | none | grounded local model, refusal-first, genotypes never leave the process |
| Cross-vendor comparison (v3.2) | none | which two of your own kits disagree, out of how many shared positions, with strand artifacts and palindromic sites counted apart from real disagreement |
| Endpoints | 20 | plus 32 v3 paths |
| Tests | 1929 | 3425 |

`CHANGELOG.md` has the full history, including the v3.1 installer and interface
fixes and the v3.2 concordance endpoint.

---

The honesty features, which are the point

Most consumer DNA tools fail in the same four ways. DNAInsight addresses each one
directly, and these behaviours are covered by tests so they cannot regress.

**1. It will not alarm you about a variant you do not carry.**
A ClinVar classification describes an ALLELE, not a position. Showing
"pathogenic" to someone carrying two reference copies is simply wrong, and it is
the most common way these tools frighten people for no reason. Non-carriers are
scored down to a quarter, their repute is cleared, and the card says plainly that
you do not carry the reported variant.

**2. It admits when the strand cannot be verified.**
Testing companies report the plus strand of GRCh37. Reference databases sometimes
store the other one. For an A/T or C/G genotype, complementing gives you the
other allele you already have, so no metadata can settle which reading is right.
Those sites are capped at magnitude 2, badged "strand ambiguous", and explained.
13 of the 122 bundled variants are affected.

**3. It distinguishes "not present" from "never checked".**
A genoset is a rule over several positions. If your array did not read one of
them, the rule cannot be evaluated at all. Those appear in a separate section
headed "not testable on your array", never mixed in with rules that were checked
and found absent.

**4. It never labels a trait good or bad.**
Traits and polygenic scores always render neutral grey. Eye colour is not a
verdict. A no-call scores exactly zero, because a failed probe is not a finding.

v3.0 applies the same rule to four new places. An imputed call can never outrank
a measured one, so it carries its DR2 and ceilings below the typed maximum. A
pharmacogenomic result defaults to Indeterminate rather than Normal, because an
unread position has not shown a variant is absent. The phrase "not a carrier" is
forbidden in code; the answer is always "not a carrier for the N variants
tested", with the residual risk. And a population the panel cannot resolve is
reported NOT RESOLVABLE, never zero percent, because zero percent is a
measurement and unmeasurable is not.

---

Why the magnitude is not the SNPedia magnitude

SNPedia's Magnitude and Repute are hand-curated by wiki editors and licensed
CC-BY-NC-SA-3.0-US. Redistributing them would force this repository to
non-commercial share-alike. So DNAInsight computes its own from CC0 and public
domain evidence: CPIC guideline level, ClinVar review status, replicated GWAS
support, publication counts, population frequency and your carrier status.

The 0 to 10 shape is intentional, so anyone who has read a Promethease report can
read this one. The numbers are not the same, and the interface says so.

Every card can expand to show exactly how its number was produced, step by step.
An opaque interest score would be worse than none.

If you want the real SNPedia values, the Database view has an opt-in fetch that
runs on your machine and writes to `~/.dnainsight/`, outside this project folder.
It is licence-gated and refuses to run until you accept the terms. Nothing it
downloads is ever committed. See `data/DATA_SOURCES.md`.

---

Supported DNA providers

| Provider | Format | Notes |
|---|---|---|
| AncestryDNA | .txt tab-delimited | V1 and V2 arrays |
| 23andMe | .txt tab-delimited | all array versions |
| MyHeritage | .csv or .txt | |
| FamilyTreeDNA | .csv | the SNP test, not the STR test |
| LivingDNA | .txt | |
| Generic TSV | .txt / .csv | auto-detected column layout |
| VCF, gVCF | plain or gzipped | new in v3.0, read by streaming |
| BAM, CRAM | alignment | new in v3.0, targeted pileup only, needs samtools |

Upload the uncompressed file. Extract the zip from your provider first.

Genome build is detected from contig LENGTHS, the one header field that cannot
quietly disagree with the coordinates in the file. Anything that is not GRCh37 is
refused rather than annotated, because mixing builds is the most common way this
class of tool produces confidently wrong answers.

---

Requirements

- Python 3.10 or newer
- About 60 MB of disk space
- Internet optional. Everything core works offline.

Two dependencies, Flask and requests. No compiler, no scientific stack, no
database server.

**Optional, for six v3 features only.** Ancestry, local ancestry, imputation, Y
and mtDNA haplogroup depth, and the assistant each need a third-party program.
DNAInsight cannot bundle them, because the best available tools are GPL-3.0 and
vendoring one would relicense this repository by copyleft. So DNAInsight ships
only the adapter, which is MIT, and you install the tool yourself into
`~/.dnainsight/tools/` on explicit consent. Until you do, the capability reports
`not_attempted` and the interface hides the control rather than pretending to
have run. `docs/EXTERNAL_TOOLS.md` is the install guide and the licence record.

---

Download

If you have never used GitHub, this is the only part that is unfamiliar. It
takes about a minute.

1. Open **https://github.com/JBrady0850/dnainsight** in a browser.
2. Click the green **Code** button near the top right.
3. Choose **Download ZIP** from the menu that opens.
   The direct link is
   https://github.com/JBrady0850/dnainsight/archive/refs/heads/main.zip
4. Extract the ZIP somewhere you can find again, such as your Desktop or
   Documents folder. You will get a folder named `dnainsight-main`.

**On Windows, do this before extracting.** Right-click the downloaded ZIP,
choose **Properties**, tick **Unblock** at the bottom, then **OK**. Windows
marks files that came from the internet, and without this step the installer
may be blocked with no useful explanation. Then right-click the ZIP and choose
**Extract All**.

Do not double-click `install.bat` while you are still looking inside the ZIP.
Windows will preview a ZIP like a folder, but scripts run from that preview
cannot find the files they need. Extract first, then open the extracted folder.

If Windows shows a blue **"Windows protected your PC"** box, that is SmartScreen
reacting to a script it has not seen before, not a virus warning. Click **More
info**, then **Run anyway**.

**On macOS**, double-click the ZIP to extract it. You will need Terminal for the
next step: open Terminal, type `cd ` with a trailing space, then drag the
extracted folder onto the Terminal window and press Enter. That puts you in the
right place without typing a path.

**If you already use git**, skip all of the above:

```
git clone https://github.com/JBrady0850/dnainsight.git
cd dnainsight
```

---

Install

Open the extracted folder, then:

Windows: double-click `install.bat`
macOS and Linux: `bash install.sh`

The installer checks for Python and installs it if it is missing, fetches the
two dependencies, builds the bundled reference, and creates a launcher. It
offers to start DNAInsight when it finishes.

Or manually:

```
pip install -r requirements.txt
python app.py
```

DNAInsight opens at http://127.0.0.1:5050 . Press Ctrl+C to stop.

Both installers finish by importing the application and building it, so a
successful banner means it actually starts. If verification fails they say so and
stop rather than reporting success.

To run the test suite as well:

```
pip install -r requirements.txt -r requirements-dev.txt
python -m pytest tests -q
```

---

Walkthrough

**1. Get your raw data.** AncestryDNA: your name, then DNA, Settings, Download
DNA Data. 23andMe: Tools, Browse Raw Data, Download. MyHeritage: DNA, Manage DNA
Kits, Download. FamilyTreeDNA: myDNA, Chromosome Browser, Download Raw Data.

**2. Add a profile.** Click "+ Add profile", enter your details and upload the
extracted .txt file. DNAInsight reports how many variants it read, typically
600,000 to 700,000.

**3. Run a scan.** Choose which subsystems to include. Everything except the
online API option works with no internet. A scan of the bundled reference is
effectively instant; the optional API pass can take 30 to 60 minutes because it
batches every rsID.

**4. Explore the findings.** Open the Filters panel. Drag the magnitude slider,
untick reputes, pick a gene, change the reference population, or type a query:

```
chr7                 everything on chromosome 7
chr7:1000-2000       a position range
/MAG>=3              magnitude 3 and above
/STARS>=2            ClinVar review stars
/CLNSIG=5,4          pathogenic and likely pathogenic only
/flipped             calls whose alleles were complemented
/ambiguous           calls whose strand cannot be verified
/carrier             variants you actually carry
/imputed             predicted calls, never measured on your array
/typed               measured calls only
/dr2>=0.9            imputation quality threshold
```

Escape resets every filter. Ctrl+H opens help. F opens the panel.

**5. Generate reports.**

- **Genetic health report.** For you. Plain language, a glossary, and a clear
  split between what needs a prescriber, what needs a lifestyle change, and what
  is background.
- **Doctor discussion report.** For a clinician. Leads with a prescription
  critical table sorted by CPIC level, groups interactions by drug, states its
  own limitations, and includes an AI analysis prompt block.
- **Interactive offline report.** One HTML file containing the findings AND a
  working filter engine. No server, no internet, no DNAInsight install. It makes
  zero network requests. This is the copy worth keeping.
- **Report addendum.** New in v3.0. Dated and additive, and it never rewrites the
  original, because a clinician who acted on the January report has to be able to
  see exactly what the January report said.

**6. Scan again later.** Every scan writes a ledger snapshot, so the next one
tells you what changed FOR YOU rather than what changed in the databases: a
variant of uncertain significance that became pathogenic and that you carry, a
CPIC level that moved, a finding that became evaluable. `POST
/api/profiles/<id>/manifest` emits a signed record of exactly which database
versions produced your report, so it can be reproduced or shown to have drifted.

---

Pooling several DNA files

Different testing chips read different positions, so two files from the same
person cover more together than either alone.

Go to DNA files, add another file, and set the relationship to "Yourself". The
files are pooled.

**Where the two files disagree, both readings are kept and shown side by side.**
Nothing is voted on and no winner is picked, because a disagreement between two
arrays is information about reliability that a silent merge would destroy.

Set any other relationship, such as Mother, and that file is used for comparison
only. It never changes your own calls. With both parents loaded you also get
offspring transmission probability and Mendelian consistency checking.

New in v3.0, loaded kits are also compared for shared DNA segments, with a shared
cM total and a relationship range. This compares only the files you loaded here.
There is no matching database and nothing to opt out of.

---

What is in the bundled reference

`data/snp_reference.json`, 122 curated variants, chosen for actionability rather
than count. Each carries a plain-English interpretation plus a risk allele, CPIC
level, ClinVar review stars, publication count, topics and medicines. That
evidence layer is what makes ranking and offline carrier awareness possible.

| Category | Focus | Key genes |
|---|---|---|
| Pharmacogenomics | CPIC-level drug response: warfarin, statins, antidepressants, opioids, thiopurines, fluoropyrimidines | CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, NUDT15, DPYD, UGT1A1, G6PD |
| Metabolic | obesity, type 2 diabetes, nutrient processing, iron overload | FTO, TCF7L2, PPARG, MTHFR, HFE, SLC30A8, MTNR1B |
| Cardiovascular | clotting, lipids, coronary risk | F5, F2, APOE, LPA, 9p21, APOA5 |
| Inflammation | chronic inflammation, autoimmune susceptibility | IL6, TNFA, IL10, CTLA4, IL6R, ERAP1, HLA-DQA1 |
| Neurological | mood, folate cycle, stress response, neurotransmitters | MTHFR, COMT, BDNF, MAOA, SLC6A4, FKBP5, OXTR |
| Detox | oxidative stress, alcohol, nicotine | GSTP1, SOD2, NQO1, ALDH2, CYP1A2, CHRNA3 |

Coverage: 115 of 122 have a risk allele, 46 have a CPIC assignment of which 24
are Level A, and 118 have population frequencies across 16 populations.

Annotation draws on CPIC, ClinVar, gnomAD, 1000 Genomes and the GWAS Catalog,
with the PGS Catalog for polygenic scores. PharmGKB and ClinPGx are deliberately
used for nothing; `data/DATA_SOURCES.md` section 9 records why.

**Want more?** `python data/build_full_reference.py --array-file <your raw file>`
builds a much larger local database from ClinVar, the GWAS Catalog and CPIC,
filtered to the positions your array actually reads. It is gitignored because it
is large and fully reproducible.

**What is not verified.** `docs/KNOWN_GAPS.md` lists every figure in this release
that was not machine-checked at source: the bundled Y and mtDNA markers, 17 star
alleles, and the carrier frequencies. Everything in it ships and works. That is
exactly why it is written down.

---

Privacy

- Your DNA is stored only in `dnainsight.db` on your computer.
- The optional MyVariant.info pass sends rsIDs only. Never genotypes.
- The optional SNPedia harvest sends page titles only, and writes outside this
  folder.
- The interactive report makes no network requests at all.
- Delete a profile and its data goes with it.
- The local assistant talks to a model on loopback or it does not answer.
  Genotypes are stripped before anything leaves the process, and the prompt is
  re-scanned afterward. If any survived, nothing is sent.
- The running application makes no network calls at all. Builders and external
  tools only run when you run them.

`.gitignore` blocks `uploads/`, every `.db`, and anything with `snpedia` in its
name, so your genetic data and any non-commercial cache cannot be committed by
accident.

---

Disclaimer

DNAInsight is not a medical device and does not provide medical advice. Consumer
DNA arrays are not clinical-grade tests and cover far less than clinical exome or
genome sequencing, so a negative result here does not rule anything out.

Do not start, stop or change any medication based on this software. Confirm any
significant finding with a clinically validated test and discuss it with a
licensed clinician, pharmacist or genetic counsellor. The American Board of
Genetic Counseling maintains a directory at findageneticcounselor.com .

Three v3.0 features need naming. Carrier screening here is not clinical carrier
screening, and anyone making a reproductive decision needs a test ordered through
a clinician. Pharmacogenomic diplotypes are not clinical pharmacogenomic testing,
and CYP2D6 is always provisional because an array cannot see copy number.
Imputed calls are never confirmatory, whatever their DR2.

For educational and research use only.

---

Project layout

```
dnainsight/
├── app.py                     Flask entry point, initialises the schema
├── requirements.txt           two runtime dependencies
├── requirements-dev.txt       adds pytest and flake8
├── backend/
│   ├── parsers.py             provider detection and raw file parsing
│   ├── sequencing.py          VCF, gVCF, BAM, CRAM, build detection, liftover
│   ├── merge.py               multi-file pooling, conflicts, trio
│   ├── orientation.py         strand reconciliation and ambiguity
│   ├── scanner.py             offline annotation and the API pass
│   ├── frequency.py           population frequency, strand tolerant
│   ├── genosets.py            criteria parser and evaluator
│   ├── traits.py              traits, ABO and RhD
│   ├── prs.py                 polygenic scores
│   ├── snpedia.py             opt-in local cache, licence gated
│   ├── scoring.py             magnitude, repute, confidence
│   ├── pipeline.py            scan orchestrator, correct stage order
│   ├── filters.py             filtering, sorting, faceting
│   ├── database.py            SQLite access
│   ├── external.py            tool registry, licence gate, subprocess runner
│   ├── ledger.py              reclassification snapshots and addenda
│   ├── provenance.py          source graph, signed manifests, licence audit
│   ├── haplogroups.py         Y and mtDNA backbones plus adapters
│   ├── relatedness.py         IBD, cM estimation, relationship ranges
│   ├── imputation.py          Beagle adapter, DR2, the magnitude cap
│   ├── ancestry.py            global and local ancestry, panel manifest
│   ├── diplotype.py           CPIC star alleles, prescription guard
│   ├── carrier.py             carrier panel, residual risk, ACMG coverage
│   ├── assistant.py           grounded local model, refusal-first
│   ├── genetic_report.py      report for you
│   ├── doctor_report.py       report for a clinician
│   ├── interactive_report.py  self-contained offline report
│   └── routes.py, routes_v2.py, routes_v3.py
├── data/
│   ├── build_reference.py     curated table plus evidence overlay
│   ├── evidence_overlay.py    risk alleles, CPIC levels, stars, publications
│   ├── build_genosets.py      the 65-rule corpus
│   ├── build_frequencies.py   population frequencies from Ensembl
│   ├── build_prs.py           polygenic models, with a licence gate
│   ├── build_full_reference.py  optional large local database
│   ├── build_panel.py         1000 Genomes plus public-tier SGDP panel
│   ├── build_pgx_alleles.py   CPIC allele tables
│   └── DATA_SOURCES.md        licence record for every source
├── frontend/index.html        the single page app
├── docs/
│   ├── API_V2.md, API_V3.md   the API contracts
│   ├── EXTERNAL_TOOLS.md      licence architecture and install guide
│   └── KNOWN_GAPS.md          every figure this release did not verify
├── tests/                     3280 tests
└── tools/                     release gate and verification harnesses (dev-only)
```

---

Contributing

Evidence bar for a new variant: CPIC Level A or B, or ClinVar pathogenic at 2
review stars or better, or a GWAS association replicated in two or more
independent studies. It must plausibly be on 23andMe v4/v5 or AncestryDNA v2, and
it needs a plain-English interpretation aimed at a non-expert.

Add the row to `data/build_reference.py`, add its evidence to
`data/evidence_overlay.py` including the risk allele on the GRCh37 plus strand,
then run:

```
python data/build_reference.py
python -m pytest tests -q
```

Before opening a pull request, run the full release gate:

```
python tools/golive.py
```

It rebuilds every derived artifact, audits every source file for duplication,
runs the module smoke tests, the strand regression, the pipeline contract check,
the filter engine check, a 42-endpoint API sweep, the interactive report
verification, the GitHub Actions runtime check, the lint gate, a clean-clone CI
simulation, the harness isolation guard and the full test suite.

`CONTRIBUTING.md` documents the design decisions that should not be quietly
reversed and the data-source parsing traps worth reading before you touch the
reference builders.

**The most valuable contribution right now is verification.**
`docs/KNOWN_GAPS.md` lists every unverified figure with the file and field it
lives in. Confirming one Y marker against ISOGG, one mtDNA position against
PhyloTree, one star allele against the CPIC tables, or one carrier frequency
against a citable source, and flipping its `verified` flag, is worth more than a
new feature. Bring the source with the pull request.

Any new external tool must permit commercial use and redistribution, or it goes
in `BLOCKED` with a reason and a named replacement. It is invoked through
`external.run` and nowhere else, because that one function is the licence
boundary and it has to stay auditable in one place.

---

License

MIT. The code is freely reusable. Bundled data is CC0 or US public domain, so the
whole repository is redistributable. Per-source licences are recorded in
`data/DATA_SOURCES.md`.

External tools keep their own licences, which attach to your copy of those
programs and not to DNAInsight, because you install them and they run as
separate processes. `GET /api/v3/licence-audit` checks the bundling rule at
runtime, so the document and the code cannot drift apart in silence.
