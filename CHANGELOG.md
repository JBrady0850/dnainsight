# Changelog

All notable changes to DNAInsight are documented here. This project follows
[Semantic Versioning](https://semver.org/).

## [3.3.0] - 2026-08-10

### Changed
- **`Y_BACKBONE` now records what a marker IS, not only what its alleles are.**
  Entries gained `variant_type` ("snv" by default), `ancestral_seq` and
  `derived_seq`. A marker that is not a base substitution clears the
  single-base `ancestral` and `derived` fields and carries whole sequences
  instead.

  Widening the existing fields to hold "GGGG" was considered and rejected. Every
  reader of those fields, `marker_state` included, treats them as one base, so a
  widened value would have been compared against an array call and silently
  never matched. A field that cannot match is worse than a field that is
  absent, because absence is visible.

### Fixed
- **M17 is a deletion, not a G>A substitution.** Confirmed against dbSNP rs3908:
  `snp_class: delins`, SPDI `NC_000024.10:19571278:GGGG:GGG`, HGVS
  `g.19571282del`, `SEQ=[G/-]`. A single-base deletion inside a four-base G
  homopolymer. There is no A allele at that site and there never was.
- **M91 is a 9T to 8T length polymorphism, not an A>T substitution.** Stated as
  such in the Karafet et al. 2008 text. No rsID has been established, so the
  sequences come from the publication and the row stays unverified.
- **Both corrections change what the tree can honestly claim.** M91 defines BT,
  which sits on the path to every non-A haplogroup, so BT now reports as
  ASSUMED rather than confirmed on essentially every call. That is not a
  regression; it is the first accurate statement this project has made about
  that node. Four tests that encoded the old assumption were updated to assert
  the corrected behaviour, with the reason written into each.

### Added
- `untypeable_markers()` and `untypeable_markers` on the haplogroups payload,
  reported separately from `unverified_markers`. An unverified marker might be
  right and nobody checked. An untypeable one is an indel that no array base
  call can ever satisfy, so its ceiling is STRUCTURAL rather than a matter of
  coverage. Collapsing the two would tell a user their array fell short when the
  marker was never callable at all, which is the same "not present" against
  "never checked" distinction the rest of the project already holds.
- The dbSNP audit folded into the table as a visible provenance block:
  `assembly`, `ref_carries` and `dbsnp_checked` on 18 rows, and a `multi_allelic`
  flag on M45, M343, M269 and P312.

  **Nothing was marked `verified`.** dbSNP settles class, position and the
  reference/alternate pair; it cannot settle ancestral against derived, which is
  the assignment `verified` is about. The audit measured why that matters: the
  GRCh38 reference Y carries the DERIVED allele at 10 of the 17 determinable
  nodes, so a builder mapping `ref` onto `ancestral` would invert 59 percent of
  this table.
- The unresolved M20 conflict recorded in its note rather than quietly patched.
  dbSNP rs3911 carries A/G single-allelic; the table records ancestral A,
  derived C. One of the rsID assignment and the allele pair is wrong and the
  audit cannot say which, so neither was changed.
- 15 tests, 3425 total.

### Fixed, unrelated to the backbone
- **`datetime.utcnow()` removed from all 9 call sites** in `backend/database.py`
  and `backend/routes.py`. It is deprecated and scheduled for removal.

  The obvious replacement is not a drop-in. `datetime.now(timezone.utc).isoformat()`
  appends "+00:00", so every row written after the change would carry a
  different string shape from every row written before it, inside columns this
  application sorts and displays, and one call site appends a literal "Z" which
  would have produced "+00:00Z". `database._utc_now_iso()` converts back to a
  naive value before formatting, so the stored string is byte-identical to every
  release up to 3.2.2 and no migration is needed.

### Decided
- **`REFUSAL_AADR` stays in force.** The Harvard Dataverse record states CC0 1.0,
  which answers the first of the two grounds the exclusion rested on. The
  second, that a compendium cannot grant rights its components did not grant, is
  a question about the constituent studies and nothing read so far addresses it.
  A licence declaration by an aggregator is evidence about the aggregator's
  intent, not proof the underlying rights existed to be granted.
  `data/DATA_SOURCES.md` now records the CC0 finding, the citation, both grounds
  and what would settle the open one, so the next reader does not re-litigate
  ground 1 and mistake it for the whole question.
- **NCBI dbSNP recorded as a source**, US Government work and public domain,
  with an explicit note on the reference-against-ancestral distinction and on
  the 31 markers it cannot reach because it does not index marker names.

### Still open
- **Karafet 2008 Supplemental Table 1 remains unobtained.** The PMC article page
  now returns a CAPTCHA challenge rather than the article, and no route around
  that was attempted. Without it, 31 of 49 markers carry no rsID and cannot be
  audited, including M60, M175 and M267.

## [3.2.2] - 2026-08-10

### Fixed
- **`ledger._utc_now()` now issues strictly increasing timestamps.** A
  microsecond FORMAT is not microsecond RESOLUTION, and the difference was
  silently corrupting the reclassification ledger on Windows.

  `datetime.now()` rides on `time.time()`, which on Windows is served by
  `GetSystemTimeAsFileTime` and advances in steps of roughly 15.6 ms. Snapshots
  written inside one step all received the SAME `created_at` string, so the old
  docstring's promise that "two snapshots taken in the same second still order
  deterministically" was false on the platform this project is most used on.

  The consequence was not a cosmetic tie. `changes_for(since=...)` skips any
  comparison whose newer snapshot is at or before the cutoff, so two consecutive
  snapshots sharing a stamp made a real reclassification report as **no change
  at all**. A user asking what changed since their last scan could be told
  nothing had, while the ledger held the change.

  Stamps are now issued under a lock and never repeat or move backwards. A
  collision advances one microsecond past the last issued value, which also
  absorbs a clock stepped backwards by NTP. Every stamp stays fixed width, so
  lexical order still equals chronological order and `ORDER BY created_at` is
  untouched.

  **How this was found.** An intermittent failure of
  `test_since_filters_out_earlier_comparisons` on Windows, on a suite that was
  green minutes earlier and green 60 consecutive times on Linux. It was
  reproduced deliberately by quantising the clock to a 15.6 ms tick before it
  was treated as anything other than a flake.
- **`changes_for` declares an ambiguous cutoff instead of guessing.** Databases
  written before this release can still hold tied timestamps. A bare timestamp
  cannot say WHICH of several snapshots sharing it the caller meant, and no
  tie-break invented in `changes_for` recovers an identity the cutoff never
  carried: guessing the earliest silently drops a comparison, guessing the
  latest silently adds one. Both were prototyped and both were rejected.

  The filter is therefore unchanged, and the payload gained `since_ambiguous`
  and `since_note`. A caller told "nothing changed" is now also told the cutoff
  matched several snapshots and the answer may be incomplete. Absent and zero
  stay separate, which is the rule the rest of this module already follows.

### Added
- 9 tests in `tests/test_ledger.py`, 3410 total. They pin that a burst of stamps
  never repeats and never goes backwards, that a frozen coarse clock still
  yields distinct increasing stamps, that a backward clock step does not reorder
  snapshots, and that an empty result from a tied cutoff is never silent.

## [3.2.1] - 2026-08-10

### Added
- `tools/audit_y_dbsnp.py` and `docs/Y_BACKBONE_AUDIT.md`. The Y backbone audit
  the v3.1.1 notes said was needed, actually run, against NCBI dbSNP, which is
  a US Government work and public domain.

  **`Y_BACKBONE` is unchanged. Every row remains `verified: false`.** The audit
  reports; it does not correct. Two of its findings need a schema decision and a
  primary source respectively before anything can be written back.
- `tests/test_y_dbsnp_audit.py`, 25 tests, covering the classification logic
  offline against real dbSNP records captured on the run date.

### Findings
- **M17 is an indel, not a base substitution.** dbSNP rs3908 reports
  `snp_class: delins`, SPDI `NC_000024.10:19571278:GGGG:GGG`, HGVS
  `g.19571282del`, `SEQ=[G/-]`. `Y_BACKBONE["R1a1a"]` records ancestral G,
  derived A. There is no A allele at that site.

  M17 defines R1a1a, so a genotyping rule expecting a G or A base call cannot
  fire against a deletion and that node is currently unreachable by the logic
  meant to reach it. This is the first of the four indel leads corroborated
  against a citable accession rather than recall.
- **M20 conflicts with dbSNP on the allele pair.** rs3911 carries A/G,
  single-allelic. `Y_BACKBONE["L"]` records ancestral A, derived C. There is no
  C allele and complementing does not reconcile the two. One of them is wrong
  and the audit cannot say which.
- **The reference Y carries the derived allele at 10 of 17 determinable nodes.**
  The v3.1.1 entry argued that mapping `ref` onto `ancestral` "would have
  inverted roughly half the tree with the entire suite still green". That was an
  argument. It is now a measurement: 59 percent. `ref_carries` stays in force.
- **Coverage is 18 of 49 markers, and the other 31 are reported as unaudited.**
  They carry no rsID, and dbSNP does not index marker names: `esearch` for M91,
  M175 and M267 each returns zero hits. M91, M60, M175 and M267 all sit in that
  unreachable set, so the retracted audit's claims about them remain
  unsupported. Karafet 2008 Supplemental Table 1 is still the missing artifact.

### Notes on the implementation
- `parse_spdi` exists as its own tested function because the dbSNP `spdi` field
  is a comma-separated list when a site is multi-allelic. A first draft split on
  ":" and expected four parts, which blanked the allele pair for M45, M343, M269
  and P312 and reported all four as conflicts against named markers. Every one
  was fine, with the recorded pair a subset of a larger observed set. That draft
  was discarded before it was acted on and the parsing is now pinned by test.
- `CONFLICT` is the residue of the classification, never its default. Class
  errors, multi-allelic subsets and strand differences are each ruled out first,
  for the same reason `backend/concordance.py` orders its verdicts that way: the
  bucket that becomes a published claim about real data is the one nothing else
  accounted for.

## [3.2.0] - 2026-08-10

### Added
- `backend/concordance.py` and `GET /api/profiles/<pid>/concordance`.
  Cross-vendor agreement across the kits already loaded into a profile: which
  two files disagree, out of how many positions they could have disagreed over,
  and what kind of disagreement it is.

  `merge.py` has pooled kits and retained conflicts since v2.0, and it
  deliberately refuses to reconcile them. What it never said was WHICH two files
  disagreed or what the denominator was, which is the number a user with two
  kits actually wants.

  The module is derived, read-only and additive. It never mutates the merged set
  and nothing it returns changes which genotype the rest of the application uses.
- Conflict classification, which is the reason this is a module and not a field.
  Most apparent vendor disagreement is strand orientation, not error. One company
  reports a SNP on the plus strand and another on the minus strand, so the same
  person reads AA in one file and TT in the other and nothing is wrong with
  either. Publishing that as a vendor error rate would be a false accusation
  about a named company in the most credible form one can take, a statistic.

  So every disagreement is classified before it is counted:

  - `orientation_artifact`, one call is the exact complement of the other and
    neither is an irreducible heterozygote, so the complement explains it
  - `indeterminate`, at least one call sits on a palindromic A/T or C/G
    heterozygote, where a strand flip and a real difference cannot be told apart
  - `genuine`, the residue, and only once every other explanation is ruled out

  The indeterminate bucket is never folded into either neighbour. Folding it
  into genuine overstates vendor disagreement; folding it into artifact hides
  real disagreement. It is the same treatment "not testable on your array"
  already gets in `genosets.py` and "strand ambiguous" gets in `scoring.py`.
- `tests/test_concordance.py`, 44 tests, written before the module. Four lines
  are pinned there specifically: a palindromic disagreement is never resolved,
  an unparsed genotype can never become evidence against a company, the three
  buckets always sum to the conflict total, and AG against CT is an artifact
  while AA against GG is genuine.

### Changed
- `pipeline.available_subsystems()` reports `concordance`, so the capability map
  can gate the control the same way it gates every other v3 subsystem.
- `backend/routes_v3.py` gained section 6 and renumbered the sections after it.

### Notes on the implementation, recorded because they are easy to get wrong
- All strand logic is delegated to `backend/orientation.py`, which already owns
  the complement table and the ambiguity rule. A second copy would drift, and
  drift in strand handling produces results that are internally consistent and
  externally backwards, which is the same failure mode the v3.1.1 `ref_carries`
  guard exists to prevent.
- Genotypes are validated against ACGT specifically, **not** against the keys of
  `orientation.COMPLEMENT`. That table tolerates no-call and indel symbols on
  purpose, so that flipping a whole file leaves them intact. Inheriting that
  tolerance here would have let `NN` become evidence that two companies
  disagree.
- No rate is reported without its denominator. Every rate travels with `shared`,
  and a pair with nothing in common reports `comparable: false` with a rate of
  `null`, never 0.0 and never 100.0. "They never agreed" and "we never compared
  them" are different claims.
- `findings_covered` is `null` when no findings were supplied and `0` when
  findings were supplied and none were covered. Null means nobody asked.
- A single kit returns `available: false` with `not_attempted: true` and
  `totals: null`. Totals of zero would claim a comparison ran and found nothing,
  which did not happen. One kit is an absent comparison, not a failed one.
- Same-provider pairs are compared and flagged with `same_provider`, never
  dropped. Two kits from one vendor years apart ran on different chips and their
  agreement with each other is a real number.
- Kits with no declared provider each form their own coverage group. Pooling
  every undeclared kit into one bucket would invent a vendor that agreed with
  itself.
- A known property, stated rather than hidden: a palindromic heterozygote reads
  the same on either strand, so two kits always agree there. Those positions
  count towards agreement and the agreement rate is very slightly optimistic as
  a result. They are left in because excluding positions from a denominator for
  being too easy is its own distortion, and the indeterminate count sits next to
  the rate.

## [3.1.1] - 2026-08-09

### Added
- Both haplogroup naming systems on every call. `snp_name()` and
  `equivalent_names()` render a node in the SNP-based form as well as the
  letter-number one, so J1 also reports as **J-M267** and R1b1a1b as
  **R-M269**. Every Y payload gained `equivalent_names` and `also_written_as`.

  This came from a real misreading. A widely shared claim held that consumer
  vendors "misclassify J1 as generic J-M267". M267 is the SNP that defines J1,
  so the two strings name the same node and there was no misclassification to
  report. A DNAInsight result sat next to a FamilyTreeDNA result with no way to
  see they agreed. Low resolution and a wrong call are different failures and
  this project exists to keep them apart, so the fix belongs in the payload
  rather than in documentation.
- `assembly`, `ref_carries` and `dbsnp_checked` on every Y backbone entry, and
  the rule that no entry may claim `verified` without all three.

  dbSNP is the obvious source for confirming that table and it cannot answer the
  question that matters. It reports **reference over alternate**; the table
  records **ancestral over derived**. On the Y these routinely disagree, because
  the GRCh38 reference Y comes from a lineage carrying the derived allele at
  many backbone nodes. dbSNP gives rs2032595 (M168) as chrY:12702062 T>C
  forward while the table gives ancestral C, derived T, and both are correct.
  A builder mapping `ref` onto `ancestral` would have inverted roughly half the
  tree with the entire suite still green, since the data would be internally
  consistent and externally backwards. The same shape as the CPIC positive-strand
  conflicts in `docs/KNOWN_GAPS.md`.
- `tests/test_haplogroup_nomenclature.py`, 28 tests, and `tests/test_version_consistency.py`, 9 tests. J1 and J2 equivalence sets
  may never intersect, no name may be claimed by two nodes, and nothing may be
  marked verified without a recorded reference orientation.
- `NOTICE` at the repository root, for attribution-licensed data. Currently
  lists none.

### Changed
- `data/DATA_SOURCES.md` now permits **CC-BY** alongside CC0 and public domain,
  on the owner's explicit instruction of 2026-08-09. Share-alike and
  non-commercial terms remain refused. The rule existed to stop a licence
  stripping the MIT grant from downstream users, and attribution does not do
  that. No CC-BY data is bundled yet; `licence_audit()` is deliberately left
  strict until the first such file actually lands.
- `unverified_markers()` rows now carry `assembly`, `ref_carries` and
  `dbsnp_checked`, so the audit list states the remaining work rather than only
  naming the row.

### Fixed
- A stale competitive claim in the `haplogroups.py` docstring asserting that
  MyHeritage does not offer haplogroups. MyHeritage added a Y-DNA tier in 2026,
  capped at Intermediate by product decision rather than by the assay. The
  corrected passage carries the date it was checked, because claims about other
  vendors go stale silently.
- `tools/vversion.py` compared every version string in the repository against
  every other one, including the built data artefacts, whose versions are
  deliberately independent of the application's. It could therefore never pass,
  and had reported a mismatch on every run of the release gate since v2.0 while
  being carried as a standing "safe to ship" warning. A warning that is always
  on is a warning nobody reads, so a genuine version skew would have looked
  identical to every other run. Application and artefact declarations are now
  checked within their own groups. **The release gate reports zero blockers and
  zero warnings for the first time.**

## [3.1.0] - 2026-08-04

### Added
- Installer verification. `install.bat` and `install.sh` both gained a fifth
  step that imports the backend, builds the Flask app and parses the bundled
  reference before declaring success. Building the app is the load-bearing
  check: a missing data file or a broken blueprint surfaces there and nowhere
  earlier. Failure is fatal and names the command that shows the real error.
- `tests/test_install_scripts.py`, 47 tests locking both installers to reality.
  Referenced paths must exist, the two scripts must agree on builders and
  requirements, step counters must be contiguous, PEP 668 escape hatches must
  be present, and the installers' own three checks run in-process. All static,
  so a Linux runner can test the Windows batch file.
- Three states in the dashboard capability table instead of two: available,
  not built, and needs a separate tool. The table also lists all sixteen
  subsystems rather than the five it was hardcoded to in v2.

### Fixed
- `install.bat` contained `echo ... Settings > Database to update`. The `>` is
  a redirect, so the tip printed truncated and a file literally named
  `Database` was created in every Windows install directory. Shipped in v1 and
  v2, found by running the installer rather than reading it.
- The install banner lost its exclamation mark to `EnableDelayedExpansion`.
  Removed rather than escaped: the first fix used `^!`, which the parser still
  ate, and the test guarding it checked the input instead of the output.
- `install.bat` called bare `pip`. A machine can have `python.exe` on PATH
  without `pip.exe`, and it then failed two steps after python was proven to
  work. Now `python -m pip` throughout.
- `data/build_reference.py` announced "Built bundled reference v2.0.0" during a
  v3 install. It now prints the data version and the application version
  separately, because they are different facts.
- `requirements-dev.txt` pinned `pytest>=8.0,<9.0`, which fails on any current
  environment and had done since before v3.0, so a clean checkout could not run
  its own release gate. Raised to `>=8.4,<10.0` after verifying the full suite
  on 9.1.1.
- flake8 was commented out of the dev requirements as optional while CI runs it
  as a job and `tools/golive.py` treats it as a blocker. Now a real dependency.
- Duplicated blocks in `backend/ancestry.py`, `backend/haplogroups.py`,
  `backend/provenance.py` and `backend/routes_v3.py`. Three were genuine
  copy-paste and were factored into named helpers. The provenance case was the
  opposite, two parallel licence-table literals where the columns are the
  point, so it moved to a keyword-only constructor with no defaults: omitting a
  field is now a TypeError at import rather than an unstated licence claim.
- One flake8 E131 in `backend/imputation.py`.
- `tests/test_install_scripts.py` ran `bash -n` against a Windows path on
  Windows, where `shutil.which` finds the WSL launcher. Skipped on Windows,
  where the check would be testing the developer's WSL setup rather than this
  repository.

### Changed
- `tools/vfrontend.py` integrity pin updated for the capability table change,
  with a comment recording why. Any other movement in those three numbers is
  unreviewed drift and should keep failing.
- `DNAInsight.png` regenerated from a running instance against a 640,000 marker
  export, palettised from 425 KB to 156 KB, and de-identified: the previous
  capture carried a real name and date of birth into an image that ships in a
  public repository.
- Page title dropped its stale "v2".
- `tools/golive.py` closing suggestion no longer tells you to branch for v2.0.

### Known
- VERSION CONSISTENCY reports the bundled data artefacts at 2.0.0 against an
  application at 3.1.0. That is correct and deliberate. The data did not
  change, and `backend/provenance.py` documents why an artefact version must
  not track the application version.

## [3.0.0] - 2026-08-04

Ten capabilities, delivered in five waves: a reclassification ledger, a
provenance graph with signed manifests, sequencing ingest, haplogroups,
household IBD, imputation, ancestry, star-allele pharmacogenomics, carrier
screening with residual risk, and a grounded local assistant.

Six of the ten need a tool DNAInsight cannot legally bundle. That constraint,
not the algorithms, is what shaped this release.

### Licensing position, decided first because it shapes everything else

v2.0 established that nothing SNPedia-derived ships in this repository. v3.0
extends the same rule from data to executables.

The best available tools for global ancestry, imputation, phasing and Y
haplogroup calling are GPL-3.0, and several of the obvious alternatives are
academic-only or non-commercial. Vendoring a GPL-3.0 binary would relicense
DNAInsight by copyleft. Reimplementing four published algorithms badly would be
worse than not shipping them.

- **DNAInsight ships only the adapter, which is MIT.** The tool is installed by
  the user, on explicit consent, into `~/.dnainsight/tools/`, deliberately
  outside this repository tree. **The subprocess boundary IS the licence
  boundary.** No external tool is ever imported, linked or vendored, and every
  invocation goes through the single `external.run` function so that boundary is
  auditable in one place.
- This is the pattern `backend/snpedia.py` already used for CC-BY-NC-SA data, so
  there is one rule to understand rather than two.
- **Five tools are permanently blocked and cannot be installed even on
  consent**, each with a recorded reason and a named replacement: ADMIXTURE,
  RFMix v2, yhaplo, yallHap, and DIYDodecad with the Eurogenes, Dodecad, MDLP
  and HarappaWorld model files. Naming them makes the exclusion a decision
  somebody made on a date rather than an oversight, exactly as
  `data/DATA_SOURCES.md` section 9 does for PharmGKB.
- **The offline contract is unchanged.** The running application makes zero
  network calls. Builders and tool installers may download, but only when the
  user runs them, and both builders dry-run by default until `--accept-terms` is
  passed.
- `docs/EXTERNAL_TOOLS.md` is the full architecture and install guide.

### Added: the reclassification ledger

- **`backend/ledger.py`** records an insert-only snapshot of every finding's
  comparable clinical state after each scan, and diffs consecutive snapshots
  into **15 classified change kinds**. `vus_resolved_pathogenic` ranks highest
  in the system by design: a variant of uncertain significance that became
  pathogenic **and that this person carries** is the highest-value event in
  personal genomics, and no consumer product surfaces it as an event.
- Every change record names the field, both values, and whether the move was
  up, down or lateral. "Changed" on its own is useless to a reader deciding
  whether to worry.
- The fingerprint deliberately excludes prose. Interpretation text churns
  constantly, and a naive whole-record diff would report a hundred clinically
  identical findings and bury the one real reclassification.
- **Addenda are dated and additive.** `supersedes` is always None and no code
  path in the module writes to a prior snapshot or a prior report. A clinician
  who acted on the January report must be able to see exactly what the January
  report said, because their decision is only defensible against the evidence
  that existed at the time.
- `routes_v2`'s scan now records a ledger snapshot on completion, defensively,
  so a ledger failure can never turn a successful scan into a failed one.

### Added: provenance, the signed manifest and a runtime licence audit

- **`backend/provenance.py`** carries `SOURCES`, a machine-readable mirror of
  `data/DATA_SOURCES.md`, and `licence_audit()`, which checks every declared
  bundled artefact against the rules that document states. That turns
  `DATA_SOURCES.md` from prose somebody has to remember into an enforced runtime
  contract, exposed at `GET /api/v3/licence-audit`.
- **Signed reproducible report manifests.** HMAC-SHA256 over the canonical
  manifest, keyed by a secret generated on and stored only on this machine. The
  scope is stated in the payload: this proves the report was not altered after
  generation **on this machine**. It is not a public-key attestation and does
  not prove authorship to a third party. Overclaiming that would be worse than
  not signing.
- `verify_manifest` returns a structured verdict naming the field that drifted,
  never a bare boolean. "Verification failed" tells a clinician nothing they can
  act on.
- **Conflict detection displays disagreement and never resolves it.** Every
  conflict record carries both positions with `verdict` None and `resolved`
  False, mirroring how `backend/merge.py` already treats two pooled files that
  disagree at a position.

### Added: sequencing ingest

- **`backend/sequencing.py`** reads VCF, gVCF and gzipped input by streaming,
  so a whole-genome file does not have to fit in memory.
- **Genome build is detected from contig LENGTHS**, which is the only header
  field that cannot quietly disagree with the coordinates in the body of the
  file. An `assembly=` tag is a claim, a `##reference` line is frequently a
  stale path, and a `chr` prefix is evidence about the distributor rather than
  the assembly. Lengths from both builds produce `None` with confidence
  `conflict`, because such a file is either concatenated or hand-edited and
  neither is safe to annotate.
- **`BuildMismatch` refuses the file rather than annotating it.** Mixing builds
  is the most common way this class of tool produces confidently wrong answers,
  and it now gets the same loudness the project already gives strand ambiguity.
- A **UCSC chain-format liftover**. With no chain file present it returns the
  standard unavailable payload rather than passing coordinates through
  untranslated. Where a chain maps to the minus strand the alleles are
  complemented, because moving a coordinate while leaving a plus-strand allele
  alone produces a genotype that is wrong at a position that is right.
- **Targeted BAM and CRAM pileup extraction via samtools**, at the reference
  positions DNAInsight already cares about, never full variant calling. A 30x
  BAM is 100 to 200 GB and re-calling it locally would be a different product.
- Skipped records are reported with all five reasons, zero-valued when nothing
  was skipped. A user who uploads a four-million-record gVCF and gets 600,000
  calls is owed the arithmetic.

### Added: haplogroups

- **`backend/haplogroups.py`** ships bundled backbone trees: **49 Y markers and
  28 mtDNA nodes**, callable with no external tool installed.
- A **tri-state marker call**: derived, ancestral, or `NOT_ON_ARRAY`. The third
  state is the point. An array that never read a marker has not shown the marker
  is ancestral.
- A **computed resolution ceiling**, not a stored constant. It describes this
  array and this person, and it says in one plain sentence where their data runs
  out.
- Adapters for **Yleaf**, **HaploGrep 3** and **Clade Finder**. Where two tools
  disagree the disagreement is surfaced, not reconciled.
- **Tree versioning on every response.** A haplogroup is meaningless without the
  tree that produced it.
- The bundled backbone is provisional and every response says so. See
  `docs/KNOWN_GAPS.md`; `verified_only=true` refuses every unverified marker and
  yields an unresolved call, which is the honest state.

### Added: household genomics

- **`backend/relatedness.py`** implements IBS-based IBD detection in pure
  Python, with an **IBIS** adapter for users who install it.
- **Relationships are returned as RANGES, never one answer.** A half sibling, a
  grandparent and an aunt all share about a quarter of their DNA and no total
  separates them. Returning one of the three would be inventing certainty.
- **Every centimorgan figure is flagged `cm_estimated`.** The genetic map is
  approximate whole-chromosome averages, and a number derived from it must not
  read as a measurement.
- **Role-disagreement detection.** A declared sibling sharing about 1,700 cM is
  a half sibling, and the person deserves to be told that rather than shown a
  number and left to work it out.
- Parental phasing where the trio makes it determinable, ambiguous where it does
  not, plus chromosome browser data shaped for the offline report's SVG.
- **This is deliberately not a matching service.** GEDmatch's profile count is a
  network effect, not a software moat. Only kits loaded on this machine are
  compared, which needs no database and is immune to the opt-out circumvention
  documented at GEDmatch in 2023 and at MyHeritage in 2025. The response says so
  in `scope`, so the limitation is never mistaken for a bug.

### Added: imputation

- **`backend/imputation.py`** adapts Beagle 5.5 and carries **DR2 as a
  first-class, filterable field** rather than a footnote.
- **`IMPUTED_MAGNITUDE_CAP` is 3.0** below the DR2 threshold, and even a
  perfectly imputed call ceilings at 9.5 against a typed ceiling of 10.0. **An
  imputed call can never reach parity with a typed one**, and that property is
  a structural guarantee asserted in the tests, not a convention.
- **The cap is written into the magnitude audit trail as a named step**, and the
  step is appended even when the cap did not bind. "The ceiling was considered
  and did not bind" is information; a silent no-op would leave the reader unable
  to tell the rule ran at all.
- **`assert_no_imputed_pathogenic_without_quality`** is invariant 1 of this
  project, "do not alarm a non-carrier", restated for imputation. Imputation
  multiplies that risk by volume, so it is a checked guarantee with its own
  endpoint rather than a comment somebody wrote once.
- `filters.py` gained `/imputed`, `/typed`, `/provisional` and `/dr2>=N`. A
  finding with no DR2 never matches a DR2 comparison, and a missing `imputed`
  key means typed rather than unknown.

### Added: ancestry

- **`backend/ancestry.py`** adapts fastmixture in projection mode for global
  ancestry and FLARE for local ancestry, plus chromosome painting.
- **A population the array cannot resolve is reported NOT RESOLVABLE with a null
  proportion, never as zero percent.** Zero percent is a measurement: it says we
  looked and found none. Not resolvable says we could not look. Reporting one as
  the other is how every incumbent turns a model artefact into an apparent fact
  about a person.
- **Confidence intervals always carry the method that produced them.** With
  bootstrap replicates at the default of 0 the intervals are Wilson
  approximations from marker counts and are labelled as such, because presenting
  an approximation as a bootstrap is a quiet overclaim.
- **Per-population marker coverage** on every result.
- **`panel_manifest` publishes the model**: populations, per-population sample
  counts, source, licence, build, marker count and a content hash. An incumbent
  that says "4,500 regions" cannot produce that dict for its own product.
- FLARE requires phased input. Unphased input is caught and reported as
  `input_not_phased` with the fix, rather than silently producing a picture that
  looks like a painted chromosome without being one.

### Added: star-allele pharmacogenomics

- **`backend/diplotype.py`** calls CPIC star alleles for **9 genes**: CYP2C19,
  CYP2C9, VKORC1, SLCO1B1, TPMT, NUDT15, DPYD, UGT1A1 and CYP2D6.
- **Indeterminate is the default, not Normal.** It is returned whenever a
  chromosome was filled with the reference allele while some allele of that gene
  was untestable, whenever an activity value is unknown, and whenever nothing
  could be called. Defaulting to Normal is how this class of tool tells somebody
  they metabolise a drug fine when nobody checked. The cost of Indeterminate is
  a user who has to ask a pharmacist; the cost of a wrong Normal is a user who
  does not.
- Alleles carry a tri-state: present, absent, or **untestable**. An array that
  does not carry rs4244285 has not shown CYP2C19\*2 is absent, it has shown
  nothing.
- **CYP2D6 is always provisional.** Arrays cannot see copy number or hybrid
  alleles, so \*2xN duplications and \*5 whole-gene deletions are invisible.
- **The prescription guard names banned imperative language and a test
  enforces it.** No output string may instruct anyone to begin, continue or
  change a medicine. `audit_language` exists so that claim is testable rather
  than believed.

### Added: carrier screening with residual risk

- **`backend/carrier.py`** implements the residual risk arithmetic
  `f(1-DR)/(1-f*DR)`, derived from Bayes in the module docstring, for an
  11-gene panel.
- **The bare phrase "not a carrier" is forbidden in code**, along with
  "non-carrier", "no risk", "rules out" and five others. The answer is always
  "not a carrier for the N variants tested", where N is the number this file
  could actually read. CFTR alone has thousands of catalogued variants and a
  consumer array reads a few dozen.
- **When a detection rate is unknown, `residual_risk` returns None with a
  reason** rather than borrowing a plausible-looking number from a neighbouring
  population. A guess that looks like the honest feature is worse than nothing.
- **Every residual risk is flagged `is_lower_bound`**, because published
  detection rates belong to clinical panels that read more positions than
  DNAInsight does.
- **Joint reproductive risk** for autosomal recessive and X-linked recessive
  inheritance, returned as a range rather than a number, because the inputs are
  unverified population figures and false precision here would be worse than an
  honest interval.
- **ACMG secondary-findings coverage is reported honestly as near zero.** For an
  array the answer is zero in nearly every gene, and the report says zero in
  words rather than rendering an empty row a reader could mistake for a clean
  result. Carrying probes for three BRCA founder variants is not BRCA testing.

### Added: the grounded local assistant

- **`backend/assistant.py`** is refusal-first. Retrieval is limited to this
  profile's own findings, every claim must cite a finding id that was actually
  in the context, and a response citing anything else is rejected and replaced
  with the refusal.
- **Genotypes are stripped before anything leaves the process**, and the
  assembled prompt is re-scanned for genotype strings afterward. If any
  survived, nothing is sent.
- **It fails closed.** A rejected response returns the refusal text, never the
  model output alongside a warning, because that would put the hallucination on
  screen.
- Ollama is reached over loopback only, and like every other external tool it is
  unavailable until installed and licence-accepted.
- `GET /api/v3/assistant/contract` publishes the grounding contract, the refusal
  text, the redacted fields and the banned phrases. A safety rule nobody can
  read is a safety rule nobody can check.

### Added: the external tool registry and the v3 endpoint surface

- **`backend/external.py`** is the registry, licence gate and runner: 11 tools,
  5 permanently blocked entries, 2 reference panels, a three-state licence
  workflow and one `unavailable()` payload shape used by every adapter.
- **`available: False` plus `not_attempted: True` is the load-bearing pair.**
  "We looked and found nothing" and "we could not look at all" are different
  claims, and collapsing them is the exact failure this project already refuses
  for not-testable genosets, unverifiable strands and no-calls.
- **`backend/routes_v3.py`** adds 31 paths across 32 method bindings, registered
  as a third blueprint so a v3 subsystem that fails to import cannot take the v1
  or v2 application down with it. A missing subsystem returns 501, not 404: the
  path exists and the capability is real, it is this installation that cannot
  serve it.
- **`data/build_panel.py`** builds the 1000 Genomes plus public-tier SGDP
  reference panel, refusing the restricted SGDP tier by construction and placing
  HGDP behind a second opt-in that is separate from the licence gate.
- **`data/build_pgx_alleles.py`** reconciles `backend/diplotype.py` against the
  CPIC allele definition tables and writes the reconciliation report into its
  own output.
- **`data/tools_manifest.json`** is generated from `external.py` and checked
  against it by a test, so the published manifest cannot drift from the registry
  that actually gates execution.

### Changed

- Version 2.0.0 to 3.0.0.
- `backend/filters.py` gained the `/imputed`, `/typed`, `/provisional` flags and
  the `/dr2` comparison operator.
- `backend/interactive_report.py`'s `_KEEP` list carries the v3 fields into the
  offline file: `imputed`, `dr2`, `imputation_quality_band`,
  `imputation_capped`, `magnitude_ceiling`, `provenance`, `source_ids`,
  `conflicts`, `provisional`, `residual_risk`, `diplotype` and `phenotype`. An
  imputed call that loses its DR2 on the way out is exactly the opaque number
  this project refuses to ship, and the report is the artefact a clinician
  holds, so it is the artefact that has to carry its own evidence.
- `pipeline.available_subsystems()` reports the ten new modules plus every
  external tool capability under a `tool_` prefix.
- `app.py` registers the v3 blueprint and initialises the ledger and provenance
  schemas. Both are `CREATE TABLE IF NOT EXISTS` plus additive column migration,
  and both are wrapped defensively so a partial checkout degrades to a working
  v1 and v2 application rather than a server that will not boot.

### Fixed

- **External tool capability names collided with DNAInsight module names in
  `available_subsystems()`.** Beagle's capability is literally called
  `imputation` and Ollama's is `assistant`, which are also the names of the two
  DNAInsight modules that drive them. The first wiring attempt merged the tool
  capabilities into the subsystem map unprefixed, so the tool flag overwrote the
  module flag and **a user without Beagle installed was told DNAInsight's own
  imputation module was missing.** That is precisely the conflation the map
  exists to prevent: an absent optional tool was reported as a broken install.
  Every tool capability is now prefixed `tool_`, keeping the two namespaces
  apart. Caught by a test written for exactly that risk before the bug appeared,
  and the test now pins the prefix so the collision cannot return.

### Testing

- 3220 tests passing, up from 1929. New suites for ledger, provenance,
  sequencing, haplogroups, relatedness, imputation, ancestry, diplotype,
  carrier, assistant and the v3 builders.
- The builders suite compares `data/tools_manifest.json` against
  `backend/external.py` field by field, so the two cannot drift.
- Behavioural guarantees with their own tests rather than their own comments:
  an imputed call never reaches typed parity, no output path emits the bare
  phrase "not a carrier", no prescription-guard string contains imperative
  dosing language, and the `tool_` prefix keeps the two capability namespaces
  separate.

### Known gaps, listed rather than hidden

`docs/KNOWN_GAPS.md` is new and is the most important document in this release.
Everything in it is shipped and working, and every item rests on a figure that
was not machine-checked at source: all 49 Y markers, 15 of the 28 mtDNA nodes,
17 star alleles, 22 of 23 carrier variant mappings, all 25 carrier frequencies,
all 6 detection rates, the genetic map, the relationship bands, and every
external tool argument form, none of which has been executed against an
installed binary.

One defect is unresolved and named there rather than quietly carried:
`data/evidence_overlay.py` files **rs28371706 under CYP2C9** while the same
rsID is widely reported as the **CYP2D6\*17** defining variant c.1023C>T. Both
attributions cannot be correct. Separately, the CPIC builder run on 2026-08-04
found three direct base disagreements against `backend/diplotype.py`, both
sides stating the positive chromosomal strand, so they are genuine conflicts
rather than convention differences. Neither file was edited. Resolve at source
before relying on those calls.

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
  no longer exist. Its durable content, pipeline stage ordering, the Tier 2
  reference builder's parsing traps, and the design decisions that should not be
  quietly reversed, is now in `CONTRIBUTING.md` under Architecture Notes.
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
