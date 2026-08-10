---
created: 2026-08-04
modified: 2026-08-04
tags: [archivist, dnainsight, reference, decision]
aliases: [DNAInsight v3 API Contract, API_V3]
---

# DNAInsight v3.0 API Contract

Single source of truth for the v3 HTTP surface. The backend
(`backend/routes_v3.py`) and the single-page app (`frontend/index.html`) are
both built to this document. If they disagree, this document is right and the
code is wrong.

v3 is a third blueprint, additive to v1 and v2. Every v1 and v2 endpoint keeps
its exact behaviour. `docs/API_V2.md` remains the contract for those, and for the
finding object, which v3 extends rather than replaces.

Base URL when running locally: `http://127.0.0.1:5050`

31 paths, 32 method bindings.

## 0. Design constraints

The four from v1 and v2 still hold:

1. Offline first. **The running application makes zero network calls.** The
   optional MyVariant pass and the SNPedia harvest remain the only network paths
   and both are v2 endpoints. No v3 endpoint reaches the network. Builders and
   tool installers do, but they are separate programs the user runs themselves.
2. No genotype ever leaves the machine. The assistant is the new pressure point
   and is handled in section 4.11.
3. Nothing SNPedia-derived is ever written inside the repository.
4. A capability that cannot run says so, and says what would fix it.

v3 adds one:

5. **An external tool's absence is a normal state, not an error.** Six of the ten
   v3 capabilities need a third-party program the user installs themselves. None
   of them raises when it is missing. Every one returns the payload in section 1.

## 1. The degradation contract

This is the load-bearing shape of the whole v3 surface. Every adapter returns it
instead of raising when its tool is not ready.

```json
{
  "available": false,
  "capability": "ancestry_global",
  "tool": "fastmixture",
  "tool_id": "fastmixture",
  "state": "not_installed",
  "reason": "fastmixture is not installed. This analysis was not attempted, which is different from finding nothing.",
  "not_attempted": true,
  "results": [],
  "how_to_enable": {
    "tool": "fastmixture",
    "homepage": "https://github.com/Rosemeis/fastmixture",
    "install_to": "~/.dnainsight/tools/fastmixture",
    "expected_files": ["fastmixture"],
    "requires": [],
    "steps": ["..."],
    "licence": "GNU General Public License v3.0",
    "note": "",
    "alternative": "Or set DNAINSIGHT_TOOL_FASTMIXTURE to an existing install path."
  }
}
```

`not_attempted` is the field that matters. **"We looked and found nothing" and
"we could not look at all" are different claims.** Collapsing them is the exact
failure this project already refuses in three other places: a genoset over
positions the array never read is reported not testable rather than absent, a
strand that cannot be verified is badged rather than guessed, and a no-call
scores zero rather than counting as a negative finding. A UI that renders this
payload as an empty result set is broken.

`state` is one of:

| state | Meaning | Fix |
|---|---|---|
| `not_installed` | no binary found in any search location | install it |
| `runtime_missing` | binary found, its runtime is not, for example Java | install the runtime |
| `licence_not_accepted` | installed and runnable, consent not recorded | POST the licence endpoint |
| `blocked` | permanently excluded on licence grounds | use the named replacement |
| `ready` | all three conditions met | nothing |

An installed tool whose licence has not been accepted is deliberately reported
unavailable. The user is never in a state where DNAInsight silently ran
something they did not agree to.

Some adapters carry a `problem` key alongside, which distinguishes causes the
`state` cannot: `tool_missing`, `panel_missing`, `input_not_phased`, `no_input`,
`bad_mode` and `run_failed`. A missing panel and a missing tool have different
fixes and therefore different messages. `run_failed` carries the stderr tail,
because "imputation failed" with no detail is unactionable.

## 2. Additions to the finding object

Every v2 key is unchanged. v3 adds the following, all optional, all absent
rather than null on a v2-shaped finding.

### 2.1 Imputation quality

| Key | Type | Notes |
|---|---|---|
| `imputed` | boolean | true when the call was predicted, not measured. **An absent key means typed**, not unknown |
| `dr2` | float or null | Beagle's dosage r-squared, 0 to 1. Absent on a typed call |
| `imputation_quality_band` | string | `high`, `moderate`, `low`, `unusable`, `unknown` at cut points 0.9, 0.8, 0.3, 0.0 |
| `imputation_capped` | boolean | the magnitude ceiling actually bound |
| `magnitude_ceiling` | float | the highest magnitude this call may reach, always strictly below 10.0 |

An imputed finding always carries a named step in `magnitude_factors` beginning
`imputation cap`, **including when the cap did not bind**. A silent no-op would
leave the reader unable to tell the rule ran at all.

### 2.2 Provenance and conflicts

| Key | Type | Notes |
|---|---|---|
| `provenance` | object | per-field source attribution |
| `source_ids` | array of string | source ids from `provenance.SOURCES` |
| `conflicts` | array | disagreements between sources, each carrying both positions, `verdict` null and `resolved` false |

There is no code path that resolves a conflict. This mirrors pooled DNA files,
where two disagreeing arrays both keep their call.

### 2.3 Clinical extras

| Key | Type | Notes |
|---|---|---|
| `provisional` | boolean | the code itself marks this call as not yet trustworthy, for example a CYP2D6 diplotype or an unverified haplogroup marker |
| `diplotype` | string | for example `*1/*2` |
| `phenotype` | string | one of `Poor`, `Intermediate`, `Normal`, `Rapid`, `Ultrarapid`, `Indeterminate` |
| `residual_risk` | object or null | see 4.8 |

All of these travel into the interactive offline report. The report is the
artefact a clinician holds, so it is the artefact that has to carry its own
evidence. An imputed call that loses its DR2 on the way out is exactly the
opaque number this project refuses to ship.

## 3. Capabilities

`GET /api/v3/capabilities`

```json
{
  "version": "3.0.0",
  "subsystems": {"frequency": true, "genosets": true, "imputation": true,
                 "ancestry": true, "tool_imputation": false,
                 "tool_ancestry_global": false, "tool_assistant": false},
  "external": {"imputation": false, "ancestry_global": false,
               "haplogroup_y": false, "assistant": false,
               "panel_onekg_sgdp": false, "panel_hgdp": false},
  "panels": {"onekg_sgdp": false, "hgdp_optional": false},
  "offline": {"core": true, "note": "..."}
}
```

**Three axes, kept separate on purpose.**

- `subsystems` is code that ships with DNAInsight and works offline.
- `external` is third-party tools the user installed and accepted.
- `panels` is reference data the user built.

A control needs all three of its dependencies before it is worth showing.
Collapsing them into one boolean would make an unbuilt panel look like a broken
feature.

**Tool capabilities appear in `subsystems` under a `tool_` prefix and this is
not cosmetic.** Beagle's capability is literally named `imputation` and Ollama's
is `assistant`, which are also the names of the DNAInsight modules that drive
them. Merging them unprefixed made a user without Beagle installed see
DNAInsight's own imputation module reported missing. The prefix keeps the two
namespaces apart and a test pins it.

## 4. Endpoints

Every path below returns **501** when its subsystem module is absent from the
build, with `{available: false, not_attempted: true, error, capability}`. 501
rather than 404: the path exists and the capability is real, it is this
installation that cannot serve it. A 404 would tell the frontend the endpoint
was never part of the API, which is a different and wrong story.

### 4.1 External tools, panels and the licence gate

| Method | Path | Notes |
|---|---|---|
| GET | `/api/v3/tools` | every tool, every blocked entry, every panel, plus `tools_root`, `panel_root`, `policy` and `offline` |
| GET | `/api/v3/tools/<tool_id>` | full state plus `licence_notice` and `install` |
| POST | `/api/v3/tools/<tool_id>/licence` | `{accept_license: true}`. **403** with the full notice until it is true. **409** for a blocked tool |
| DELETE | `/api/v3/tools/<tool_id>/licence` | withdraw. The capability reports unavailable immediately |
| GET | `/api/v3/panels/<panel_id>` | build state, expected files, licence, exclusions, plus the ancestry `manifest` when available |
| GET | `/api/v3/licence-audit` | `{ok, violations, warnings, ...}` |

The 403 shape is deliberately identical to the SNPedia harvest gate in v2. One
rule to learn.

**409 for a blocked tool is not a missing consent and no body clears it.** The
five blocked entries are ADMIXTURE, RFMix v2, yhaplo, yallHap and DIYDodecad
with the community model files. Each response carries `reason` and
`replacement`.

`/api/v3/licence-audit` is the runtime half of `data/DATA_SOURCES.md`. A
non-empty `violations` list means something non-redistributable reached a
bundled artefact. Per-record licensing produces a **warning** rather than a
violation, because the PGS Catalog is legitimately bundled after filtering, and
the warning exists so nobody forgets the filtering is load bearing.

### 4.2 Sequencing ingest

| Method | Path | Body | Notes |
|---|---|---|---|
| POST | `/api/profiles/<pid>/sources/sequencing` | multipart `file`, optional `role`, `label`, `sample` | `.vcf`, `.gz`, `.bgz`, `.bam`, `.cram` |
| POST | `/api/v3/sequencing/inspect` | `{filename}` | format, build and liftover availability for an already-uploaded file, without ingesting it |

`role` is the same eight-value set as the v2 source API, so a sequencing file
can be a `self` file that pools, or a `mother` file that only compares.

Success:

```json
{"message": "Sequencing source added.", "source_id": 2, "format": "vcf",
 "build": "GRCh37", "build_confidence": "high", "build_evidence": [],
 "sample": "NA12878", "sample_count": 1, "snp_count": 612334,
 "skipped": {"indel": 0, "symbolic": 0, "no_id": 0,
             "multiallelic_complex": 0, "no_call": 0},
 "warnings": []}
```

`skipped` always carries all five reasons, zero-valued when nothing was skipped.
A user who uploads a four-million-record gVCF and gets 600,000 calls is owed the
arithmetic.

**Build detection ranks its evidence rather than pooling it.**

| Signal | Weight |
|---|---|
| contig `length=` | decides alone, confidence `high`. It is measured from the FASTA that produced the coordinates and is the only header field that cannot quietly disagree with the body |
| `assembly=` tag | used only when no length is present, confidence `medium`. It is a claim |
| `##reference` | same tier. Frequently a stale path |
| `chr` prefix | recorded and never voted on. It says the file came through UCSC tooling, which distributes both hg19 and hg38. It is evidence about the distributor |

Lengths from both builds return `build: null` with confidence `conflict`. Such a
file is concatenated or hand-edited and neither is safe to annotate.

**422 on a build mismatch**, with `detected_build`, `expected_build` and a hint.
Coordinates are never translated silently. **422** also when the file parsed but
yielded no usable genotypes.

Liftover is available only when the user has placed a UCSC chain file in
`~/.dnainsight/panels/chains/`. Without one, `liftover()` returns the section 1
payload. It does not pass coordinates through untranslated, which would produce
a file that is wrong at every position while looking converted.

BAM and CRAM need samtools, and extraction is targeted at the reference
positions DNAInsight already annotates. It is never full variant calling.

### 4.3 Ancestry

| Method | Path | Query | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/ancestry` | `panel`, `mode` | `mode` is `projection` (default), `supervised` or `unsupervised` |
| GET | `/api/profiles/<pid>/ancestry/painting` | `phased_vcf`, `panel` | needs phased input; without `phased_vcf` returns the FLARE unavailable payload with the fix |

```json
{"available": true,
 "proportions": [{"population": "CEU", "label": "...", "superpop": "EUR",
                  "proportion": 0.62,
                  "interval": {"low": 0.58, "high": 0.66,
                               "level": 0.95, "method": "wilson_marker_count"},
                  "informative_markers": 4820, "markers_read": 3910,
                  "coverage": 0.81}],
 "not_resolvable": [{"population": "...", "proportion": null,
                     "state": "not_resolvable", "reason": "..."}],
 "marker_coverage": {}, "panel_manifest": {}, "caveats": []}
```

**A population the array cannot resolve appears in `not_resolvable` with
`proportion: null`, never as zero percent.** Zero percent is a measurement: it
says we looked and found none. Not resolvable says we could not look. Reporting
one as the other is how every incumbent turns a model artefact into an apparent
fact about a person. The threshold is 20 percent marker coverage or 50 markers,
whichever binds.

**Every interval names the method that produced it.** With `bootstrap` at its
default of 0 the intervals are Wilson approximations from marker counts and are
labelled `wilson_marker_count`. Bootstrap replicates are opt-in because each one
is a full re-invocation of the external tool, and presenting an approximation as
a bootstrap would be a quiet overclaim.

`panel_manifest` publishes the model: populations, per-population sample counts,
source, licence, build, marker count and a content hash. Safe to call before
anything is built, in which case countable fields are `null` rather than `0`,
because a panel with zero populations and a panel that was never built are
different states.

`chromosome_painting` returns every chromosome, including ones with no segment.
An absent bar reads as "not analysed" and a present empty bar reads as
"analysed, nothing shared". Segments carry `x1` and `x2` as fractions so a
template can draw a rect without doing genomic arithmetic; a chromosome of
unknown length gets `null` rather than a fabricated scale, because a bar drawn
to an invented length is a lie in picture form.

### 4.4 Haplogroups

| Method | Path | Query | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/haplogroups` | `verified_only` | both systems, both ceilings, in one payload |

Returns the Y call, the mtDNA call, a `resolution_ceiling` for each, an
`unverified_markers` audit list and a `tree` stamp.

**The resolution ceiling is computed, not decorative.** `markers_available` and
`markers_in_tree` are counted from this genotype map and this tree, so the field
describes this array and this person rather than a marketing average. It also
carries one plain sentence saying where the data runs out.

Marker state is tri-state: `derived`, `ancestral`, `NOT_ON_ARRAY`, plus
`NO_CALL`, `UNUSABLE` and `DISCORDANT`. The third state is the point. An array
that never read a marker has not shown the marker is ancestral.

**`tree` is always present.** A haplogroup is meaningless without the tree
version that produced it. The bundled backbone stamps itself
`DNAInsight backbone 0.1`.

With Yleaf, HaploGrep 3 or Clade Finder installed the call is refined and the
source changes from `bundled_backbone`. **Where two tools disagree the
disagreement is surfaced, not reconciled.**

**`verified_only=true` refuses every unverified marker and yields an unresolved
call.** That is the honest state of this release: all 49 bundled Y markers and
15 of the 28 mtDNA nodes carry no verified defining position. See `docs/KNOWN_GAPS.md`. The
unverified count is on every response, because the difference between a caveat
somebody wrote once and a caveat the user actually sees is the whole point.

### 4.5 Household genomics

| Method | Path | Query | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/household` | | pairwise IBD across loaded kits |
| GET | `/api/profiles/<pid>/household/browser` | `a`, `b` | segment layout for one pair, SVG-ready. **404** for an unknown pair |
| GET | `/api/profiles/<pid>/phasing` | | which allele came from which parent |

Every response carries `scope`:

> Only the DNA files loaded into this profile are compared. DNAInsight has no
> matching database and cannot find relatives you have not already uploaded.
> That is a design decision, not a limitation to be fixed.

GEDmatch's profile count is a network effect, not a software moat, and a local
tool cannot copy it. What a local tool can do is the family case, which needs no
database and is immune to the opt-out circumvention documented at GEDmatch in
2023 and at MyHeritage in 2025.

**Relationships come back as ranges, never one answer.** Every band whose
published range contains the shared total is returned. A half sibling, a
grandparent and an aunt all share about a quarter of their DNA and no total
separates them. `longest_cm` and `segment_count` are carried through because
they are what discriminates within a band, but no rule is applied to them here.

**Every cM figure carries `cm_estimated: true`.** The genetic map is approximate
whole-chromosome averages, not a published map, so a number derived from it must
not read as a measurement. See `docs/KNOWN_GAPS.md`.

Role disagreement is reported plainly. A declared sibling sharing about 1,700 cM
is a half sibling, and the person deserves to be told that rather than shown a
number and left to work it out.

Phasing returns `{rsid: {"maternal": a, "paternal": b}}` for resolvable
positions and `{rsid: null}` for the rest, so a caller can tell "we worked it
out" from "we could not" without a second lookup. Where both parents are
heterozygous at a heterozygous child position, nothing is determined and the
position is reported ambiguous rather than guessed.

### 4.6 Cross-vendor concordance

| Method | Path | Query | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/concordance` | | how much your own kits agree, with every disagreement classified |

Every response carries `scope`:

> Only the DNA files you have loaded into this profile are compared, and only
> against each other. This is not an accuracy measurement against any truth set:
> DNAInsight has no reference genome for you, so it can say that two of your
> files disagree but never which of them is right.

`merge.py` has pooled kits and retained conflicts since v2.0 and deliberately
refuses to reconcile them. What it never said was WHICH two files disagreed or
out of how many shared positions, which is the number a user with two kits
actually wants.

**Every disagreement is classified before it is counted.** Most apparent vendor
disagreement is strand orientation, not error: one company reports a SNP on the
plus strand and another on the minus strand, so the same person reads AA in one
file and TT in the other. Publishing that as a vendor error rate would be a
false accusation about a named company in the most credible form one can take,
a statistic.

| Class | Meaning |
|---|---|
| `orientation_artifact` | one call is the exact complement of the other and neither is an irreducible heterozygote, so the complement explains it completely |
| `indeterminate` | at least one call is a palindromic A/T or C/G heterozygote, which reads the same on either strand, so a flip and a real difference cannot be told apart |
| `genuine` | the residue, and only once every other explanation is ruled out |
| `agreement` | the two calls match |
| `not_comparable` | at least one side is not a pair of ACGT base calls |

`genuine + orientation_artifact + indeterminate == conflicts`, always. **The
indeterminate bucket is never folded into either neighbour.** Folding it into
genuine overstates vendor disagreement; folding it into artifact hides real
disagreement. Same treatment "not testable on your array" gets in genosets and
"strand ambiguous" gets in scoring.

**No rate without its denominator.** Every rate travels with `shared`. A pair
that shares nothing reports `comparable: false` and `null` rates, never 0.0 and
never 100.0, because "they never agreed" and "we never compared them" are
different claims.

**`findings_covered` is `null` when no findings were supplied and `0` when they
were supplied and none were covered.** Null means nobody asked; zero means
someone asked and the answer was none.

**One kit returns `available: false`, `not_attempted: true` and `totals: null`.**
Totals of zero would claim a comparison ran and found nothing. One kit is an
absent comparison, not a failed one.

Same-provider pairs are compared and flagged with `same_provider`, never
dropped: two kits from one vendor years apart ran on different chips. Kits with
no declared provider each form their own coverage group, because pooling every
undeclared kit together would invent a vendor that agreed with itself.

Not to be confused with `/conflicts/sources` in 4.11, which is about ClinVar,
the GWAS Catalog and CPIC disagreeing over an interpretation. This endpoint is
about two of your own files disagreeing over a base call.

### 4.7 Imputation

| Method | Path | Body | Notes |
|---|---|---|---|
| POST | `/api/profiles/<pid>/imputation` | `{panel, dr2_threshold}` | starts a background job. **409** when one is already running |
| GET | `/api/profiles/<pid>/imputation/status` | | `{running, done, error, result}` |
| GET | `/api/profiles/<pid>/imputation/coverage` | | the three-way split |
| GET | `/api/profiles/<pid>/imputation/safety` | | `{ok, violations, checked}` |

`dr2_threshold` defaults to 0.8. Whole-genome imputation is minutes, not
milliseconds, so it cannot block a request thread.

Coverage reports typed, imputed with a usable DR2, and imputed without one, as
three separate counts. **A run that imputes 80 million variants of which half
are unusable is not an 80 million variant run**, and this report is what stops
the headline number from being quoted alone.

`/imputation/safety` is invariant 1 of this project, "do not alarm a
non-carrier", restated for imputation. A finding is a violation when it is
imputed AND pathogenic AND any of: `dr2` missing, quality band `low`,
`unusable` or `unknown`, or the cap never ran, detectable by a missing
`imputation_quality_band` or a missing cap step in `magnitude_factors`. An
uncapped imputed pathogenic call can sort above a typed one, which is the parity
failure this whole design prevents. The endpoint audits without raising.

Four mandatory caveats travel with every result and cannot be suppressed: an
imputed genotype was never measured; panel bias hits non-European ancestries
hardest and degrades the DR2 figures themselves for those people; accuracy falls
sharply below about 1 percent minor allele frequency, which is exactly where
clinical interest lives; and **no imputed call is confirmatory**.

### 4.8 Pharmacogenomics

| Method | Path | Body | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/pgx/diplotypes` | | all nine genes |
| POST | `/api/profiles/<pid>/pgx/prescription-guard` | `{medications: ["warfarin"]}` | only the pairs that apply |

Genes: CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, NUDT15, DPYD, UGT1A1, CYP2D6.

Each allele carries a tri-state: `present`, `absent`, `untestable`.
**Untestable is why this endpoint exists.** An array that does not carry
rs4244285 has not shown CYP2C19\*2 is absent, it has shown nothing, and
reporting that as \*1/\*1 is the most dangerous thing a consumer
pharmacogenomics tool can do.

Alleles are assigned most specific first, so an allele defined by three variants
is consumed before the one-variant allele that is a subset of it, and component
variants are consumed as they are used.

**`phenotype` defaults to `Indeterminate`, never `Normal`.** It is returned
whenever nothing could be called, whenever a chromosome was filled with the
reference allele while some allele of that gene was untestable, whenever an
activity value is unknown, and whenever the allele name is not in the table. The
cost of Indeterminate is a user who has to ask a pharmacist. The cost of a wrong
Normal is a user who does not.

`provisional_genes` is always `["CYP2D6"]`. **Arrays cannot see copy number or
hybrid alleles**, so \*2xN duplications and \*5 whole-gene deletions are
invisible and every CYP2D6 call is a partial view.

`unverified` lists the 17 star alleles this build could not corroborate,
including the three that directly conflict with CPIC's own tables. It is on
every response by design.

**The prescription guard never instructs.** No output string tells anyone to
begin, continue or change a medicine, `BANNED_PRESCRIPTIVE_PHRASES` names the
language, and `audit_language` exists so the claim is testable rather than
believed. Output is scoped to the medicines supplied, not a 500-row interaction
table the user has to search.

### 4.9 Carrier screening

| Method | Path | Body or query | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/carrier` | `population` | 11-gene panel with residual risk |
| POST | `/api/profiles/<pid>/carrier/joint` | `{gene, population}` | needs a loaded `mate` source. **409** without one |
| GET | `/api/profiles/<pid>/acmg` | | secondary-findings coverage |

Panel: CFTR, HEXA, SMN1, HBB, GJB2, PAH, ATP7B, GALT, ACADM, BTD, G6PD.

Three statuses, and the third is why the module exists:

| Status | Meaning |
|---|---|
| `carrier` | at least one tested variant detected |
| `not_carrier_for_tested_variants` | every variant that COULD be read was read, none detected |
| `untestable` | nothing in this gene could be read, so nothing was found and nothing was excluded |

**The bare phrase "not a carrier" is forbidden in code**, along with
"non-carrier", "no risk", "rules out" and five more. The negative is always
scoped to the number of variants actually read. If the panel lists five variants
and this file could read two, the answer is "not a carrier for the 2 variants
tested", not for the 5, and certainly not for the gene.

Residual risk after a negative result, where `f` is the population carrier
frequency and `DR` the panel detection rate:

```
residual risk = f (1 - DR) / (1 - f * DR)
```

- **When `f` or `DR` is unknown for this gene and population, `residual_risk` is
  `null` with a `reason` naming the missing input.** Substituting a
  plausible-looking number would produce a figure that looks like the honest
  feature while being a guess.
- **Every result carries `is_lower_bound: true`.** Published detection rates
  belong to clinical panels that read more positions than DNAInsight does.

Joint reproductive risk is a **range**, not a number, because both inputs rest
on unverified population figures. Autosomal recessive is
`risk_a * risk_b * 0.25`. X-linked recessive is `risk_a * 0.5 * 0.5` for an
affected son, and `risk_b` is explicitly reported as unused rather than silently
ignored, because silently ignoring an argument is how a caller ends up believing
the other parent was accounted for.

`/acmg` reports coverage across the encoded ACMG SF gene list and says **zero in
words** rather than rendering an empty row a reader could mistake for a clean
result. It also prints its own list discrepancy: 82 genes encoded against a
published count of 81, not reconciled item by item. Carrying probes for three
BRCA founder variants is not BRCA testing.

### 4.10 Reclassification ledger

| Method | Path | Query or body | Notes |
|---|---|---|---|
| GET | `/api/profiles/<pid>/snapshots` | | `{snapshots, latest}` |
| GET | `/api/profiles/<pid>/changes` | `since`, `limit` | every change across the whole history |
| POST | `/api/profiles/<pid>/addendum` | `{old_id, new_id, format}` | `format: "html"` also writes a report file |

A snapshot is written automatically at the end of every v2 scan.

**15 change kinds**, ranked by severity. `vus_resolved_pathogenic` is the
highest in the system by design: a variant of uncertain significance that became
pathogenic **and that this person carries** is the single highest-value event in
personal genomics, and no consumer product surfaces it as an event.

Two kinds exist beyond the obvious set. `sig_reclassified` covers a move that
cannot be placed on the benign to pathogenic axis at all, such as pathogenic to
conflicting: calling that an upgrade or a downgrade would be a lie and dropping
it would hide a real change. `carrier_status_changed` covers a reference update
that flips the recorded risk allele, so a person told they carry nothing is now
told they carry something.

**Direction is always explicit.** A change record never says only "changed": it
names the field, both values in display form, and `up`, `down` or `lateral`.

`since` compares against the LATER snapshot of each pair as an ISO 8601 string,
with a trailing `Z` stripped from both sides so `2026-08-01` and
`2026-08-01T00:00:00Z` behave the same. `limit` caps after sorting by severity,
so a caller asking for 10 gets the 10 that matter.

Diffs quote database versions recorded **at scan time**, so a change reads
"between ClinVar 2026-07-27 and ClinVar 2026-08-24" instead of the useless
"between two scans".

**The addendum is dated and additive.** `supersedes` is always `null`, snapshots
are insert-only, and no code path in the module writes to a prior snapshot or a
prior report. A clinician who acted on the January report has to be able to see
exactly what the January report said, because their decision is only defensible
against the evidence that existed at the time.

With one snapshot on record the addendum returns a well-formed baseline payload.
The first scan genuinely has nothing to compare against and that is not an error.

### 4.11 Provenance and manifests

| Method | Path | Body | Notes |
|---|---|---|---|
| POST | `/api/profiles/<pid>/manifest` | `{report_type, scan_parameters}` | returns `{manifest, text}` |
| POST | `/api/v3/verify-manifest` | `{manifest}` or the signed object | structured verdict. **400** on malformed input |
| GET | `/api/profiles/<pid>/conflicts/sources` | | `{count, findings, policy}` |

The manifest records the DNAInsight version, every database name with version,
retrieval date and licence, the sha256 of every input DNA file, the finding
count, the scan parameters and a UTC timestamp. A clinician holding a report can
ask which ClinVar release said this, and when. Without a manifest that question
has no answer.

An input file whose hash cannot be computed is recorded with `"sha256": null`
and `"present": false` rather than omitted. **A silently missing input is
indistinguishable from a report generated from nothing.**

**Read `signature.scope` before quoting the signature.** It is HMAC-SHA256 over
the canonical manifest, keyed by a secret generated on and stored only on this
machine. It proves the manifest was not altered after generation **on this
machine**. It is **not** a public-key attestation, it does not prove authorship
to a third party, and anyone holding the key file can sign anything.
Overclaiming that would be worse than not signing at all.

`verify-manifest` **never returns a bare boolean.** The verdict names the field
that failed, in check order: `signed`, `manifest`, `manifest.<name>`,
`signature`, `signature.algorithm`, `signature.signed_at`, `signature.value`.
"Verification failed" tells a clinician nothing they can act on;
"manifest.finding_count was altered" tells them exactly what to distrust. An
unsigned payload is a structured failure, not an exception, because being handed
a plain manifest instead of a signed one is a routine mistake.

`/conflicts/sources` answers `"policy": "Conflicts are surfaced, never
resolved."` Every conflict carries both positions with `verdict: null` and
`resolved: false`, and the ordering is stable so repeated calls produce
identical output.

### 4.12 Grounded local assistant

| Method | Path | Body | Notes |
|---|---|---|---|
| POST | `/api/profiles/<pid>/assistant` | `{question, model}` | answers from this profile's findings or refuses |
| GET | `/api/v3/assistant/contract` | | the grounding contract, refusal text, redacted fields, banned phrases |

The order of operations is the security design and is not negotiable:

1. `external.guard("ollama", "assistant")`. Ollama is an external tool like any
   other: unavailable until installed and licence-accepted, degrading to the
   section 1 payload.
2. Build context from the local evidence store only. No embedding model, no
   vector database, no external lookup.
3. **Refuse if the context is empty or the question is out of scope. Nothing is
   sent in that case.**
4. Redact, assemble the prompt, then **re-scan the finished prompt for genotype
   strings. If any survived, do not send.**
5. POST to the loopback Ollama endpoint, `http://127.0.0.1:11434/api/generate`
   by default. A non-loopback host is refused.
6. Validate the response post-hoc and **fail closed**.

Validation rejects an empty response, a response citing an ID that was not in
the context, a response containing banned medical-advice or dosing phrasing, and
a response that makes claims and cites nothing. **On rejection, `answer` is the
refusal text, never the model output.** Returning the model output alongside a
warning would put the hallucination on screen, which is what this exists to
prevent.

A response that is itself a refusal, saying the findings do not support an
answer and citing nothing, is **accepted**. That is the model doing exactly what
the contract asks.

`GET /api/v3/assistant/contract` publishes all of it. A safety rule nobody can
read is a safety rule nobody can check.

## 5. Filter grammar additions

Additive to the v2 grammar in `docs/API_V2.md` section 3.5.

| Form | Effect |
|---|---|
| `/imputed` | only imputed calls |
| `/typed` | only measured calls |
| `/provisional` | only calls the code marks as not yet trustworthy |
| `/dr2>=0.9` | imputation quality comparison, also `>`, `<=`, `<`, `=` |

Two rules that a naive implementation gets wrong:

- **`/typed` matches a finding carrying no `imputed` key at all.** The entire
  bundled reference predates imputation and every one of its findings is typed.
  Treating a missing key as "unknown" would hide all of it.
- **A finding with no DR2 never matches a DR2 comparison.** Treating a missing
  DR2 as zero would make `/dr2<0.5` return every typed call in the file.

## 6. Error shape

Unchanged from v2: `{"error": "human readable sentence"}` with a real status
code. v3 uses the following.

| Code | Meaning in v3 |
|---|---|
| 400 | bad input, missing required field, unreadable file |
| 403 | external tool licence not accepted |
| 404 | unknown profile, unknown pair of labels, file not in the upload directory |
| 409 | tool permanently blocked, job already running, or no `mate` source loaded |
| 413 | upload too large |
| 422 | genome build mismatch, or a file that parsed but yielded no usable genotypes |
| 500 | unexpected failure inside a subsystem |
| 501 | the subsystem module is absent from this build |

An absent external tool is **not** an error code. It is HTTP 200 carrying the
section 1 payload.

## 7. How the imputation cap works

| Constant | Value | Meaning |
|---|---|---|
| `TYPED_MAGNITUDE_CEILING` | 10.0 | the normal maximum |
| `IMPUTED_MAGNITUDE_CEILING` | 9.5 | the best an imputed call can ever reach |
| `PARITY_MARGIN` | 0.5 | the gap that guarantees separation |
| `IMPUTED_MAGNITUDE_CAP` | 3.0 | hard shelf below the DR2 threshold |
| `DEFAULT_DR2_THRESHOLD` | 0.8 | the working floor |

Two tiers. At or above the threshold an imputed call ceilings at 9.5, **strictly
below** the typed ceiling, so a prediction never ties a measurement no matter
how confident it is. Below the threshold, or with DR2 missing, it is capped at
3.0, because below the working floor the call is a lead and not a finding. The
return value is always strictly less than `TYPED_MAGNITUDE_CEILING` and the
tests assert that property.

The cap appends a step to `magnitude_factors` beginning `imputation cap`,
**including when the magnitude did not move**, and never touches a typed call.

## 8. Compatibility with v2

- Every v1 and v2 endpoint behaves exactly as `docs/API_V2.md` documents.
- The finding object gains keys only. No v2 key changed type or meaning, and no
  v2 consumer needs to know about the v3 keys.
- `GET /api/capabilities` is unchanged. `GET /api/v3/capabilities` is the
  three-axis superset.
- v3 is registered as a separate blueprint and app startup wraps it. **A build
  with no v3 modules serves a working v1 and v2 application**, printing a
  warning rather than failing to boot.
- The two v3 schemas are `CREATE TABLE IF NOT EXISTS` plus additive column
  migration. There is no DROP and no table rebuild on any path, so an existing
  database is never rewritten by an upgrade.
