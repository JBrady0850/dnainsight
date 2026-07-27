---
created: 2026-07-26
modified: 2026-07-27
tags: [archivist, dnainsight, project, decision, reference]
aliases: [DNAInsight v2 handoff, v2 build record, v2 remaining work]
---

# DNAInsight v2.0 handoff

State of the build, what is verified, and what a later session needs to know so
it does not re-derive anything or reverse a deliberate decision.

## Verification status, as of 2026-07-27

| Check | Command | Result |
|---|---|---|
| Unit tests | `python -m pytest tests -q` | **1885 passed**, 0 failed, 0 xfail |
| Source integrity | `python _build\audit.py` | all clean except one known false positive |
| Module smoke | `python _build\smoke.py` | 8 of 8 modules ok |
| Strand regression | `python _build\vfreq.py` | fixed and verified |
| Pipeline contract | `python _build\vpipe.py` | every contract key present, all invariants hold |
| Filter engine | `python _build\vfilters.py` | every documented rule verified |
| Interactive report | `python _build\vreport2.py` | self-contained, zero network requests |
| Static reports | `python _build\vreports.py` | genetic and doctor render every v2 field |
| Frontend contract | `python _build\vfrontend.py` | all 13 checklist items present |
| Database safety | `python _build\vdbloss.py` | user data survives import and `create_app()` |
| API sweep | `pwsh _build\gate.ps1` | **42 of 42 endpoint checks passed** |

One standing gate runs every check above and prints a single verdict:

```
python _build\golive.py
```

`backend\parsers.py` is reported DUPLICATED by the auditor. That is a **false
positive**: its two TSV parsers are legitimately similar, and the file is byte
identical to the v1.2 original at 234 lines. Do not "repair" it.

## What is done and working

Every item scoped for v2.0 is built, tested and verified. The list below is the
delivered surface, not a plan.

- **Eleven new backend modules**: `scoring.py`, `pipeline.py`, `filters.py`,
  `routes_v2.py`, `orientation.py`, `genosets.py`, `prs.py`, `merge.py`,
  `snpedia.py`, `frequency.py` extensions, `interactive_report.py`.
- **20 new v2 endpoints**, 42 of 42 endpoint checks green.
- **65 genosets** with the full criteria grammar: `and`, `or`, `not`,
  `atleast(N, ...)`, exact `rs1234(A;T)` and at-least-one `rs1234(T)`, and
  cross-genoset references.
- **7 polygenic models** over 42 variants.
- **18 traits**, including the ABO indel logic.
- **16 populations** covering 118 of 122 bundled rsIDs.
- **An evidence overlay** of 121 rows giving 115 risk alleles and 46 CPIC
  assignments, 24 of them Level A.
- **A computed magnitude, repute and confidence system** with a per-finding
  audit trail in `magnitude_factors`.
- **A rebuilt frontend**, `frontend/index.html` at 2421 lines, covering all 13
  items of the `docs/API_V2.md` section 5 checklist. The v1.2 file is preserved
  at `frontend/index.html.v1-backup`.
- **A self-contained offline interactive report** that opens with no server and
  makes zero network requests.
- **Both static report generators** rewritten to render magnitude, repute,
  confidence, review stars, CPIC level, population frequency, the genoset
  section and the false-positive warning at magnitude 6 or above.
- **Tier 2 full-array reference builder**, `data/build_full_reference.py`,
  producing a gitignored `data/reference.db` from ClinVar, GWAS and CPIC.
- **Fourteen pytest suites**, 1885 tests.
- **`data/DATA_SOURCES.md`** recording the licence for every bundled source.
- **`requirements.txt` reconciled** to version ranges rather than exact pins.
- **`README.md`** rewritten for v2, including the licensing position.

`docs/API_V2.md` is the authoritative contract. Code that disagrees with it is
wrong.

## Defects found and fixed during this build, worth not reintroducing

1. **`_resolve_db_path()` destroyed the user's database at import time.** The old
   write probe connected to the real database, created a test table inside it,
   then called `unlink()` on it. Every launch wiped every profile. It was
   intermittent on Windows only because WAL locks sometimes made the unlink fail.
   The fix opens an existing database read-only with a `PRAGMA schema_version`
   round trip and probes writability with a separate pid-suffixed throwaway file
   that is the only thing ever unlinked. Guarded permanently by
   `tests/test_database.py` and `_build/vdbloss.py`.
2. **Report filenames collided.** Four reports generated inside the same second
   shared one path, so report id 1 served report 4's content. Fixed with
   `_unique_report_name()`, applied to both the v1 and the v2 report routes.
3. **`clinvar_sig_code` matched "pathogenic" inside "pathogenicity".** Records
   classified "conflicting classifications of pathogenicity" scored 5 and were
   coloured Bad. The conflicting guard now runs before the compound fallback.
4. **`_ensure_contract_keys` omitted 13 keys** on non-SNP entities, and its None
   coercion used `setdefault`, which cannot replace an existing None, so the
   coercion never ran at all.
5. **`build_facets` omitted the documented `clinvar_diseases` bucket.**
6. **rs8176719 `--` was read as a homozygous deletion**, producing a confident
   blood group O call from a failed 23andMe probe. `-` was removed from
   `_DELETION_TOKENS` and an explicit `_INDEL_NOCALL_TOKENS` guard was added.
7. **Frequency lookup was strand-naive.** rs1801133 returned `unavailable` for
   CEU although the data was present, because Ensembl stores G/A and 23andMe
   reports C/T. `resolve_strand()` now returns 41.84 percent with
   `flipped: True`. `aggregate_frequency` in GLOBAL mode had the same fault
   while MAX, AVG and MIN did not.
8. **`load_frequencies(path)` did not persist**, so a caller-supplied source was
   silently discarded. Fixed with `reset_source()`.
9. **`create_app()` never initialised the schema**, so any WSGI host or test
   client hit `no such table: profiles`.
10. **Two mislabelled genes**: rs1800544 was ADRB3 and is ADRA2A, rs30187 was
    CRP and is ERAP1.

## Load-bearing details

### Pipeline stage order

`backend/pipeline.py` runs merge, then bundled and API annotation, then
**orientation, frequency, snpedia, scoring**, then genosets, traits and
polygenic scores. Two of those adjacencies are not stylistic:

- Strand resolution must precede frequency lookup, or minus-strand frequencies
  silently return unavailable.
- Frequency must precede scoring, or the rarity adjustment has nothing to read
  and every magnitude flattens toward the default of 1.0.

### Tier 2 reference builder

`data/reference.db` must be gitignored because GitHub rejects files over 100 MB.
These entries are already present in `.gitignore`: `data/reference.db`,
`data/reference.db-wal`, `data/reference.db-shm`, `data/*.txt.gz`,
`data/_cache/`, `*.bak`.

Parsing traps already verified and worth honouring:

- The ClinVar rsID column is literally `RS# (dbSNP)`, with a space and a hash.
- The gene column differs BETWEEN ClinVar tables. `variant_summary.txt.gz` calls
  it `GeneSymbol`; `gene_specific_summary.txt` calls it `#Symbol`. Verified
  against the live files on 2026-07-27. The builder accepts both spellings and
  resolves every column by name, never by position, precisely because two
  ClinVar tables describing the same genes disagree about the column name.
  `PhenotypeIDS` ends in a capital S.
- Use `PositionVCF`, not `Start` and `Stop`, which are right-shifted.
- Filter to `Assembly == "GRCh37"`. `na` and NCBI36 rows exist. A `-1` rsID means
  no dbSNP mapping.
- Review status: the docs publish "no classification for the individual variant"
  but the DATA emits "single". Code against the data string. The full mapping is
  in `backend/scoring.py` as `REVIEW_STATUS_STARS`.
- GWAS `P-VALUE` underflows to 0 for thousands of rows, so use `PVALUE_MLOG`. The
  sign of beta lives in the `95% CI (TEXT)` column, not the numeric one.
- The GWAS endpoint `https://www.ebi.ac.uk/gwas/api/search/downloads/alternative`
  RETURNS 404. Verified on 2026-07-27. EBI retired that route and now publishes
  the ontology-annotated association table as
  `gwas-catalog-associations_ontology-annotated-full.zip` under
  `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/`. The builder tries
  the retired URL first so a restored service is picked up automatically, then
  falls through to the FTP mirror. GWAS is the only moving target among these
  sources, so it carries a candidate list rather than a single address.
- CPIC is CC0 and safe, but do NOT copy the `clinpgxlevel` or `pgxtesting`
  columns out of the CPIC dump. They arrive inside a CC0 file, which makes them
  look safe, but they are PharmGKB-sourced and carry a no-commercial-sale
  clause. They are excluded at the wire level by `CPIC_SELECT`, not filtered
  afterwards.

## The build hazard, restated because it will bite again

The MCP file bridge **double-applies append-mode writes AND `edit_block`
replacements**. `rewrite(A)` then `append(B)` lands on disk as `A+B+B`. An
`edit_block` can insert its replacement twice. This corrupted seven Python files
and one Markdown file during this build.

Duplicated Python **still imports cleanly**, because later definitions shadow
earlier ones. The corruption is invisible at runtime and cannot be caught by
"does it run". Detect it structurally:

```
python _build\audit.py                    parse AST, find repeated defs and blocks
python _build\repair.py --dry-run         show what would be collapsed
python _build\repair.py                   collapse adjacent duplicate runs, with .bak
```

Rules for any future session using that bridge:

1. Never pass mode `append`.
2. Write each file in one `rewrite` call. If it is too large, write numbered
   `_build\<name>.partNN` files with `rewrite` and concatenate with
   `python _build\assemble.py <relative/target>`.
3. Run `_build\audit.py` after **every** `edit_block`, and re-read the region of
   any Markdown file you edit, since the auditor only parses Python.
4. Never use inline `python -c` or inline pwsh expressions. The bridge wraps
   commands in PowerShell and they fail to parse. Put code in a file, run the
   file.
5. Never use shell redirects or pipes inside `cmd /c`. The bridge re-executes
   concurrently and the two runs fight over the output file. Capture inside a
   single Python parent process instead, the way `_build/golive.py` does.

## Design decisions that should not be quietly reversed

1. **No SNPedia data in the repository.** It is CC-BY-NC-SA-3.0-US. Committing a
   harvested cache would relicense the whole repository and foreclose commercial
   use permanently. The harvester is opt-in, licence-gated, and writes to
   `~/.dnainsight/`, outside the repo.
2. **No PharmGKB or ClinPGx bulk data**, for the no-commercial-sale clause.
3. **Conflicting pooled calls are both retained.** No voting, no winner.
4. **A missing SNP evaluates false in a genoset**, and the genoset is reported as
   not testable rather than absent.
5. **Traits and polygenic scores never get a repute.** A trait is not good or bad.
6. **Non-carriers are down-weighted and their repute cleared.** A classification
   describes an allele, not a position.
7. **Palindromic A/T and C/G sites are flagged, not guessed**, and capped at
   magnitude 2.
8. **`0.0` frequency and `null` frequency are different facts** and must render
   differently. Not observed is not the same as not known.
9. **A null magnitude sorts and filters as 1**, not 0.
10. **No `localStorage`, `sessionStorage` or any browser storage API** in the
    frontend or the interactive report. UI state lives in a plain object.
11. **No personal genetic data in the repository.** `uploads/` content and any
    `.db` file is somebody's genome and cannot be recalled once pushed.

## Remaining optional work, none of it blocking

1. **Run the Tier 2 builder.** The code and its tests are done, but the database
   itself is a local artifact by design and has to be generated on the machine
   that will use it: `python data\build_full_reference.py`. It downloads roughly
   400 MB and takes a while. The app runs correctly on Tier 1 alone.
2. **Run the SNPedia harvester** if the user wants crowd Magnitude and Repute
   alongside the computed values: `python -m backend.snpedia --harvest`. Opt-in,
   licence-gated, writes to `~/.dnainsight/`.
3. **Decide the fate of `_build/`.** It is developer scratch and is gitignored.
   Four files encode hard-won knowledge and would be worth promoting to a
   committed `tools/` directory: `audit.py`, `repair.py`, `assemble.py` and
   `golive.py`.
4. **Populations for the last 4 of 122 bundled rsIDs.** No frequency source
   carries them; they correctly render as null rather than 0.0.
