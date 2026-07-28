# Changelog

All notable changes to DNAInsight are documented here. This project follows
[Semantic Versioning](https://semver.org/).

## [2.0.0] - 2026-07-26

Feature parity with the Promethease and SNPedia model, built entirely from CC0
and public-domain data so the repository stays redistributable.

### Licensing position, decided first because it shapes everything else

SNPedia is licensed `CC-BY-NC-SA-3.0-US`. Its Magnitude and Repute fields are
hand-curated by wiki editors and have no open-data equivalent. Bundling them
would force the derived database to non-commercial share-alike and permanently
foreclose commercial use of DNAInsight. So:

- **Nothing SNPedia-derived ships in this repository.** DNAInsight computes its
  own magnitude and repute from CPIC (CC0), ClinVar (US public domain), gnomAD
  and 1000 Genomes (CC0) and the GWAS Catalog.
- An **opt-in harvester** lets a user build a local SNPedia cache on their own
  machine for their own personal, non-commercial use. It is licence-gated, and
  it writes to `~/.dnainsight/snpedia_cache.db`, deliberately outside the repo.
- `data/DATA_SOURCES.md` records the licence for every source.

### Added: interest scoring that can be audited

- **`backend/scoring.py`** computes a 0 to 10 DNAInsight magnitude, a
  Good / Bad / unset repute and a four-level confidence. The scale intentionally
  mirrors the shape people already know from Promethease reports, but the numbers
  are ours and the UI labels them as such.
- Every score carries a **`magnitude_factors`** audit trail listing each step
  that fired, so a number can be explained rather than taken on trust.
- **Carrier awareness.** A ClinVar classification describes an allele, not a
  position. A non-carrier is multiplied down to a quarter and has its repute
  cleared, instead of being shown "pathogenic" for a variant they do not have.
  This is the single largest honesty improvement in the release.
- **No-calls score zero** and **palindromic sites are capped at 2.0** and flagged
  dubious, because an unverifiable call must not outrank a verifiable one.

### Added: the genoset engine

- **`backend/genosets.py`** implements the full criteria grammar: `and`, `or`,
  `not`, `atleast(N, ...)`, exact `rs1234(A;T)`, homozygous `rs1234(T;T)`,
  at-least-one `rs1234(T)`, cross-genoset references and arbitrary nesting.
- A **missing SNP evaluates false**, never null and never imputed to the
  population-major allele.
- **65 authored genosets** in `data/genosets.json`, including all six APOE
  diplotypes, warfarin metabolizer classes and SLCO1B1 statin risk. Every rule
  references only rsIDs that exist in the bundled reference.
- A genoset whose required positions were not genotyped is reported as **not
  testable**, which is a different fact from absent. Conflating those two is how
  a report tells someone they do not have something it never checked.

### Added: strand reconciliation, the dominant correctness risk

- **`backend/orientation.py`** handles plus and minus orientation and
  `StabilizedOrientation` semantics.
- Consumer arrays report the GRCh37 plus strand; Ensembl and dbSNP store many
  variants on the minus strand. rs1801133 is the canonical case: 23andMe calls it
  C/T, Ensembl stores G/A. Before this release the frequency lookup silently
  returned nothing for every such variant, which reads as missing data.
- **A/T and C/G heterozygotes are irreducibly ambiguous** and are now flagged
  rather than guessed. 13 of the 122 bundled SNPs are affected.

### Added: population frequency

- **`backend/frequency.py`** with 16 populations from 1000 Genomes Phase 3 via
  Ensembl, covering 118 of the 122 bundled rsIDs.
- Genotype frequency per population, observed where available and derived under
  Hardy-Weinberg otherwise, with the method always reported.
- A **rarity colour ramp** and coarse bands, plus a strict distinction between
  `0.0` (not observed in that panel) and `null` (unknown). These are different
  facts and are rendered differently.
- A population selector is a **frequency denominator, not ancestry inference**,
  and says so.

### Added: multi-file pooling, comparison and trio

- **`backend/merge.py`** pools any number of `self` files into one genotype set.
- **Conflicting calls are both retained and surfaced.** There is no voting, no
  confidence weighting and no automatic winner, because a disagreement between
  two arrays is information about reliability that a silent merge destroys.
- Eight relationship roles, comparison rows per relative, Mendelian violation
  detection and offspring transmission probability.

### Added: polygenic scores

- **`backend/prs.py`** with seven authored additive models (type 2 diabetes,
  coronary artery disease, BMI, venous thromboembolism, LDL, homocysteine,
  inflammation), 42 variant rows, weights as natural-log odds ratios with the
  source odds ratio recorded per row.
- Reference means and standard deviations are derived analytically under
  Hardy-Weinberg and validated to within 1e-6 by `--validate`.
- A `--from-pgs` mode fetches PGS Catalog scores with a **hard licence gate**
  that refuses NonCommercial, NoDerivatives and academic-only scores.
- Every result carries mandatory caveats and is marked unreliable below 90
  percent variant coverage.

### Added: traits and blood type

- **`backend/traits.py`** with 18 traits and ABO / Rh prediction.
- Traits are **never** assigned a repute. A trait is not good or bad.
- Blood type degrades honestly: when the decisive rs8176719 tag is missing or
  uncalled the answer is "not determinable", never a guess.

### Added: the filter engine and 20 new endpoints

- **`backend/filters.py`** implements server-side filtering, 10 sort keys in both
  directions, faceting with counts, and the free-text grammar including
  `chr7:1000-2000` region search and `/CLNSIG=`, `/STARS>=`, `/MAG>=`,
  `/COUNT>=`, `/flipped`, `/ambiguous` operators.
- Frequency and publication filters are **exempt for genosets, traits and
  scores**, which have no single position to have a frequency at.
- A **null magnitude sorts and filters as 1**, not 0.
- **`backend/routes_v2.py`** adds 20 endpoints: capabilities, populations,
  sources CRUD, the v2 scan, filtered findings, facets, genosets, traits, prs,
  pgx, conflicts, trio, qc, filtered export in three formats and the
  licence-gated SNPedia admin surface.
- **`backend/pipeline.py`** orchestrates the scan in the one correct order:
  merge, annotate, resolve strand, then frequency, then scoring. Strand must
  precede frequency, and frequency must precede scoring, or the score is
  computed from the wrong numbers.
- `docs/API_V2.md` is the authoritative contract for both backend and frontend.

### Added: evidence layer

- **`data/evidence_overlay.py`** adds risk allele, CPIC level, ClinVar review
  stars, publication count, topics and medicines for 121 of the 122 bundled
  rsIDs: 115 risk alleles, 46 CPIC assignments of which 24 are Level A.
- This is what makes the offline reference **carrier-aware** rather than
  allele-general, closing the limitation the v1.2 README acknowledged.
- Before it, every finding scored an identical base of 1.0. After it, DPYD \*2A
  homozygous scores 9.3 and a background trait scores 1.0.

### Fixed

- **rs8176719 `--` read as a homozygous deletion.** `--` is the token 23andMe and
  AncestryDNA write for a failed probe, and `-` was in the deletion token set.
  A failed read of the one decisive ABO tag was reported as a confident group O.
  Now every no-call spelling degrades to unknown with confidence none.
- **Strand-naive frequency lookup.** Every minus-strand variant returned
  `unavailable` although the data was present. rs1801133 CEU now returns 41.84
  percent with `flipped` recorded.
- **`aggregate_frequency` GLOBAL mode was strand-naive** while MAX, AVG and MIN
  were not, so GLOBAL returned null exactly where the others returned numbers.
- **`load_frequencies(path)` did not persist.** Every accessor re-resolved to the
  module default, so a caller pointing the module at a fixture was silently read
  back off the bundled file. Added `reset_source()`.
- **`create_app()` did not initialise the database schema**, so the app only
  worked when launched through `main()`. Any WSGI host or test client failed with
  "no such table: profiles".
- **Two mislabelled genes** in the curated reference: rs1800544 was labelled
  ADRB3 but is ADRA2A (ADRB3 Trp64Arg is rs4994), and rs30187 was labelled CRP
  but is an ERAP1 variant. Both now carry a note about the common mislabelling.

- **DATA LOSS: `_resolve_db_path()` deleted the user's database at import time.**
  The write probe connected to the real database file, created a table inside it,
  then called `unlink()` on it. Every launch destroyed every stored profile,
  finding and report. It was intermittent on Windows only because WAL locks
  sometimes made the unlink fail, which is why it survived to this point. The
  probe now opens an existing database read-only and verifies it with a
  `PRAGMA schema_version` round trip, touching nothing; writability is proved
  against a separate pid-suffixed throwaway file which is the only path ever
  unlinked. Guarded by `tests/test_database.py` and `tools/vdbloss.py`, which
  writes a canary profile and asserts it survives both a bare import and a full
  `create_app()`.
- **Report filenames collided.** The generated name carried a one-second
  timestamp, so four reports produced inside the same second shared a single
  path and report id 1 served report 4's content. `_unique_report_name()` now
  appends a counter until the path is free, on both the v1 and v2 report routes.
- **`clinvar_sig_code` matched "pathogenic" inside "pathogenicity".** Records
  classified "conflicting classifications of pathogenicity" scored as pathogenic
  and were coloured Bad, sweeping genuinely disputed variants into the default
  whitelist. The conflicting-record guard now runs before the compound fallback.
- **`_ensure_contract_keys` omitted 13 contract keys** on genoset, trait and
  polygenic-score findings, and its None coercion used `setdefault`, which cannot
  replace an existing None, so the coercion never ran at all.
- **`build_facets` omitted the documented `clinvar_diseases` bucket.**

- **CI ran on the deprecated Node 20 action runtime.** `actions/checkout@v4` and
  `actions/setup-python@v5` both bundle Node 20, which GitHub deprecated on
  2025-09-19 and now force-runs on Node 24 while printing a warning on every job.
  Bumped to `actions/checkout@v5` and `actions/setup-python@v6`, the releases that
  declare Node 24 natively. `tools/vactions.py` now fails the gate on any Node 20
  action pin, so this cannot silently return.
- **The CI lint job would have failed on a clean clone.** The go-live gate ran
  pytest but never ran flake8, so thirteen findings were invisible locally while
  they would have turned the `lint` job red. Cleared all thirteen. Two were real
  defects rather than formatting, and both are listed below. Added
  `tools/vlint.py`, which runs the lint job's own command, and
  `tools/vci.py`, which copies exactly the files `git ls-files` would publish
  into a temp tree and runs the workflow's steps there. A green local run proves
  nothing about a fresh clone, because the working tree holds gitignored
  artifacts a clone will not have.
- **`harvest_genosets` silently ignored its own `rate_limit` argument.** It built
  a `_RateLimiter` from the caller's value and then never used it, so every
  request fell back to the module-level 2.0/s default. `harvest` had the same gap
  for its per-page fetches: it threaded the limiter into target enumeration only.
  `fetch_subject` and `fetch_wikitext` now accept and forward a limiter, so both
  harvesters honour the argument. flake8 saw this as an unused local; it was a
  parameter that lied.
- **DATA LOSS RISK: the upload destination was unbounded and overwrote silently.**
  The path was `f"{profile_name}_{filename}"` with no length cap and no collision
  check. Two consequences, both reproduced. First, re-uploading a file the app had
  itself named prepended the profile name again, so the component grew about
  fifteen characters per cycle; at 246 characters the next write crossed the
  filesystem's 255-byte limit for a single component and Windows raised OSError
  EINVAL, returning HTTP 500. Second, two profiles sharing a name and a filename
  resolved to one path, so the second upload replaced the first person's raw
  export with no error and no warning, and raw DNA cannot be recovered once
  replaced. `_bounded_upload_path()` in `backend/routes.py` now truncates the stem
  to 100 characters and appends a counter until the path is free;
  `backend/routes_v2.py` imports that one implementation rather than keeping a
  second copy of the rule. Guarded by `tests/test_uploads.py`, 31 tests.
- **The verification gate was not idempotent.** Four harnesses selected their
  sample as `sorted(uploads/*.txt)[0]`, the alphabetically first file, which after
  one run was one of their own artifacts. That is what drove the filename growth
  above, and it meant the gate's result depended on how many times it had been
  run: it passed while names were short and failed three stages at once when they
  were not. Selection now uses `tools/sample.py`, which picks the SHORTEST name.
  A derived artifact is strictly longer than its source, because derivation only
  prepends, so the shortest name can never be output derived from itself. No
  prefix blacklist is used, because a blacklist stops working the moment someone
  adds a harness and forgets to update it.

- **Verification harnesses wrote into the real database and uploads folder.**
  `backend/database.py`'s highest-priority DB_PATH candidate, and
  `backend/routes.py` / `routes_v2.py`'s `UPLOAD_DIR` and `REPORTS_DIR`, are the
  exact locations the installed app itself uses. Six `tools` harnesses booted
  the real Flask app to verify endpoints and reports, so they wrote real rows
  and real files there: `DURABILITY CANARY`, `Report Verify`, `Collision Test`,
  `Static Report Verify` and `V2 Verify` profiles, plus 99 files in `uploads/`.
  Nothing distinguishes a harness row from a real one in `list_profiles()`, so
  these appeared in the app's own UI indistinguishable from user data. Added an
  explicit `DNAINSIGHT_DB_PATH` / `DNAINSIGHT_UPLOAD_DIR` / `DNAINSIGHT_REPORTS_DIR`
  environment override, checked before each hardcoded default, and
  `tools/isolated_db.py`, which every harness now calls before importing `app`
  to redirect all three to a throwaway per-run directory. `tools/vdbisolation.py`
  proves it: runs all six harnesses and asserts the real database's row count and
  mtime, and the real `uploads/` file count, are unchanged. Cleared the harness
  rows and files already present.

### Testing

- 1929 tests passing, up from 138. New suites for scoring, filters,
  pipeline, orientation, genosets, frequency, prs, merge, traits, snpedia
  and database path resolution.
- Verification harnesses under `tools/`: a structural auditor, a Markdown
  duplicate detector, a repair pass, a pipeline contract check, a filter
  engine check, an interactive-report check, a static-report check, a
  frontend contract check, a data-loss canary and a 42-check API sweep.
- `python tools/golive.py` runs every check above and prints one verdict.

### Known build hazard, documented for future sessions

The MCP file bridge used to author this release **double-applies append-mode
writes and `edit_block` replacements**: `rewrite(A)` followed by `append(B)`
lands on disk as `A+B+B`, and a single `edit_block` can insert its replacement
twice. This silently duplicated seven source files and one Markdown file. Duplicated Python still imports cleanly
because later definitions shadow earlier ones, so the corruption is invisible at
runtime and must be detected structurally. `tools/audit.py` detects it and
`tools/repair.py` reverses it losslessly. Never use append mode against this
repository from that bridge.

### Repository cleanup for public release

- **`_build/` reorganized into a committed `tools/` directory.** Twenty-six
  scripts that are permanent regression guards or part of the release gate
  (`golive.py` and everything it calls, `audit.py` / `repair.py` / `assemble.py`,
  `isolated_db.py`, `sample.py`) are now tracked source instead of developer
  scratch. Roughly ninety one-off patch scripts, ad-hoc diagnostics, a separate
  project's leftover files, and stale logs and backups were deleted rather than
  moved. `_build/` remains gitignored for future one-off scratch, and `golive.py`
  gained a new stage that runs `tools/vdbisolation.py`, so a harness that starts
  writing into the real database or `uploads/` again fails the gate.
- **`docs/HANDOFF_V2.md` retired.** It was an internal engineering handoff note,
  not user or contributor documentation, and it pointed at `_build/` paths that
  no longer exist. Its durable content — pipeline stage ordering, the Tier 2
  reference builder's parsing traps, and the design decisions that should not be
  quietly reversed — is now in `CONTRIBUTING.md` under Architecture Notes.
- **`frontend/index.html.v1-backup` removed.** The v1.2 frontend is preserved in
  git history; keeping a redundant static copy in the working tree served no
  purpose once v2 shipped.
- **Working-directory test artifacts cleared.** `uploads/`, `reports_output/`,
  and the local `dnainsight.db` are gitignored but accumulate real files across
  development and gate runs. Cleared all three so the working tree matches what
  a fresh clone starts with.
- **The release gate depended on stray leftover files it never created.**
  Clearing `uploads/` above broke six end-to-end harnesses
  (`vpipe.py`, `vfilters.py`, `vreport.py`, `vreport2.py`, `vserver.py`,
  `vreports.py`): each needs a raw-DNA file to POST, and had always silently
  relied on whatever manually-uploaded file happened to be sitting in
  `uploads/` from an earlier session. A genuinely fresh clone would fail the
  same way on its first run. Added `tests/fixtures/sample_23andme.txt`, a
  synthetic file built from the bundled reference's own rsIDs (no real
  person's genetic data, so it is safe to track) with a deliberate subset of
  rsIDs omitted so some multi-SNP genosets stay genuinely not-testable. New
  `tools/sample.py:ensure_fixture()` seeds it into `uploads/` on first use;
  every harness that needs an upload now calls it instead of assuming one
  exists. `tools/vdbisolation.py` also crashed on a missing `dnainsight.db`
  (it read `list_profiles()` without initialising the schema first); it now
  calls `init_db()`, which is idempotent, before reading. `tools/vcomplete.py`
  still required the now-deleted `docs/HANDOFF_V2.md`; its required-files list
  is updated to match the current tree.

## [1.2.0] - 2026-07-03

### Security
- **Upload path traversal fixed.** Uploaded filenames are now sanitized with
  `werkzeug.utils.secure_filename` before being written to disk, preventing a
  crafted `filename` from escaping the `uploads/` directory.
- **Upload size cap.** Requests are limited to 64 MB via Flask
  `MAX_CONTENT_LENGTH`, plus an explicit per-file check, preventing
  memory-exhaustion from oversized uploads. Oversized requests return `413`.
- **Upload type validation.** Only `.txt`, `.csv`, and `.tsv` raw exports are
  accepted; compressed/binary uploads are rejected with a clear message.

### Added
- **Zygosity on every finding.** Each finding now reports `homozygous`,
  `heterozygous`, `hemizygous`, or `no_call`, computed directly from the two
  alleles. Shown in the findings table, both report types, and CSV/JSON exports.
- **Carrier-aware API annotation.** When the MyVariant.info API returns an
  authoritative alternate allele (`dbsnp.alt`), the scanner checks how many
  copies the person actually carries. ClinVar classifications for variants the
  person does **not** carry are downgraded to informational and annotated as a
  reference genotype, rather than raising a false alarm.
- **Single-SNP lookup.** New endpoint `GET /api/profiles/<id>/lookup/<rsid>`
  and a lookup box in the Findings view answer "what is my genotype at rsX?"
  without running a full scan.
- **Reference integrity guard.** `build_reference.py` now detects and warns on
  duplicate rsIDs (previously dropped silently), and a new test suite enforces
  uniqueness, row shape, and valid category values.

### Fixed
- **`_split_genotype` no-call ordering bug.** Two-character no-call tokens such
  as `--` and `00` were split into (`-`,`-`) / (`0`,`0`) before the no-call
  check ran; they now normalize to (`N`,`N`) first.
- **Version string inconsistencies.** The startup banner and report footer
  showed `v1.0` while the API reported `1.1.0`. Version is now single-sourced
  from `backend.__version__` and consistent everywhere.

### Changed
- `findings` table gains a `zygosity` column (auto-migrated on existing
  databases).
- Bumped bundled reference and application version to `1.2.0`.

## [1.1.0] - 2026-06-24

### Added
- Expanded bundled reference to 122 curated medical SNPs (CPIC Level A
  pharmacogenomics: DPYD, TPMT, NUDT15, UGT1A1, G6PD, and others).
- Background database refresh from MyVariant.info with progress polling.
- JSON and CSV export of findings.
- Staleness banner prompting monthly reference refresh.

## [1.0.0] - 2026-06-22

### Added
- Initial release: multi-provider raw DNA parsing (AncestryDNA, 23andMe,
  MyHeritage, FamilyTreeDNA, LivingDNA), offline-first SNP annotation,
  three-silo finding classification, and genetic + doctor HTML reports.
