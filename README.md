DNAInsight v2.0.0

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Bundled SNPs](https://img.shields.io/badge/Bundled_SNPs-122_curated-orange)
![Genosets](https://img.shields.io/badge/Genosets-65-orange)
![Tests](https://img.shields.io/badge/Tests-1929_passing-brightgreen)

Personal DNA analysis that runs entirely on your own computer. Read your raw DNA
file from any major provider, annotate it against curated clinical evidence,
explore the results with real filters, and generate reports you can print, keep
or hand to a clinician.

Nothing leaves your machine unless you explicitly ask it to.

---

What is new in v2.0

v1.x told you what a position meant. v2.0 tells you what it means FOR YOU, ranks
it, and is honest about what it cannot establish.

| Capability | v1.2 | v2.0 |
|---|---|---|
| Interest ranking | none, every finding looked equal | DNAInsight magnitude, 0 to 10, with a per-finding audit trail |
| Direction of effect | none | repute: good, bad, or deliberately unset |
| Carrier awareness | online API only | offline too, from a bundled risk allele per variant |
| Multi-SNP rules | none | 65 genosets with a full boolean criteria engine |
| Population frequency | none | 16 populations, your exact genotype, observed or Hardy-Weinberg |
| Strand handling | none | full reconciliation, with unverifiable sites flagged not guessed |
| Multiple DNA files | one per profile | pool several, keep disagreements, compare relatives |
| Polygenic scores | none | 7 models with mandatory caveats and coverage honesty |
| Traits and blood type | none | 18 traits plus ABO and RhD, refusing to guess when uncertain |
| Filtering | search box and a gene dropdown | server-side engine, sliders, facets, 20 sort orders, query grammar |
| Reports | 2 static | 3, including a self-contained offline interactive report |
| Tests | 138 | 1929 |

---

The honesty features, which are the point

Most consumer DNA tools fail in the same four ways. v2.0 addresses each one
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

Upload the uncompressed file. Extract the zip from your provider first.

---

Requirements

- Python 3.10 or newer
- About 60 MB of disk space
- Internet optional. Everything core works offline.

Two dependencies, Flask and requests. No compiler, no scientific stack, no
database server.

---

Install

Windows: double-click `install.bat`
macOS and Linux: `bash install.sh`

Or manually:

```
pip install -r requirements.txt
python app.py
```

DNAInsight opens at http://127.0.0.1:5050 . Press Ctrl+C to stop.

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

---

What is in the bundled reference

`data/snp_reference.json`, 122 curated variants, chosen for actionability rather
than count. Each carries a plain-English interpretation plus, new in v2.0, a risk
allele, CPIC level, ClinVar review stars, publication count, topics and
medicines. That evidence layer is what makes ranking and offline carrier
awareness possible.

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

**Want more?** `python data/build_full_reference.py --array-file <your raw file>`
builds a much larger local database from ClinVar, the GWAS Catalog and CPIC,
filtered to the positions your array actually reads. It is gitignored because it
is large and fully reproducible.

---

Privacy

- Your DNA is stored only in `dnainsight.db` on your computer.
- The optional MyVariant.info pass sends rsIDs only. Never genotypes.
- The optional SNPedia harvest sends page titles only, and writes outside this
  folder.
- The interactive report makes no network requests at all.
- Delete a profile and its data goes with it.

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

For educational and research use only.

---

Project layout

```
dnainsight/
├── app.py                     Flask entry point, initialises the schema
├── requirements.txt           two runtime dependencies
├── requirements-dev.txt       adds pytest
├── backend/
│   ├── parsers.py             provider detection and raw file parsing
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
│   ├── routes.py              v1 endpoints, unchanged behaviour
│   ├── routes_v2.py           v2 endpoints
│   ├── genetic_report.py      report for you
│   ├── doctor_report.py       report for a clinician
│   └── interactive_report.py  self-contained offline report
├── data/
│   ├── build_reference.py     curated table plus evidence overlay
│   ├── evidence_overlay.py    risk alleles, CPIC levels, stars, publications
│   ├── build_genosets.py      the 65-rule corpus
│   ├── build_frequencies.py   population frequencies from Ensembl
│   ├── build_prs.py           polygenic models, with a licence gate
│   ├── build_full_reference.py  optional large local database
│   └── DATA_SOURCES.md        licence record for every source
├── frontend/index.html        the single page app
├── docs/
│   └── API_V2.md              the authoritative API contract
├── tests/                     1929 tests
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

---

License

MIT. The code is freely reusable. Bundled data is CC0 or US public domain, so the
whole repository is redistributable. Per-source licences are recorded in
`data/DATA_SOURCES.md`.
