# DNAInsight Data Sources and Licences

Every external dataset DNAInsight touches, what it is used for, and what may
legally be done with it. Referenced from `CHANGELOG.md`.

Last reviewed: 2026-08-09.

## Amendment 2026-08-09: attribution licences are now permitted

The rule below said CC0 and public domain only. It now also permits **CC-BY**,
on the owner's explicit instruction of 2026-08-09, with attribution recorded in
this document and in `NOTICE` at the repository root.

The reasoning is that CC-BY is not the hazard the original rule was written
against. That rule exists to stop a **share-alike or non-commercial** term
attaching itself to the repository and stripping the MIT grant from downstream
users. CC-BY does neither. It imposes an attribution obligation, which is a
documentation duty rather than a licence change, and MIT already carries a
notice-retention requirement of its own. Share-alike and non-commercial terms
remain refused, and nothing about SNPedia or PharmGKB changes.

**No CC-BY data is bundled today.** This amendment records an approved policy,
not an exercised one. The reason it was approved is that the Y phylogeny
supplements needed to verify `Y_BACKBONE` are almost all CC-BY, and the CC0-only
rule would have left that table permanently unverifiable.

`licence_audit()` in `backend/provenance.py` is deliberately **not** relaxed as
part of this amendment. It still refuses anything outside CC0 and public domain,
because loosening a guard before there is anything for it to guard costs a
protection and buys nothing. The first CC-BY artefact to actually land must
update `provenance.SOURCES` in the same change that adds the file, and that
change is where the audit gets extended.

## The project's position, in one paragraph

**The repository ships only CC0 and public domain data, so the whole thing stays
redistributable under its MIT licence.** Anything with a share-alike,
non-commercial, or no-commercial-sale term is never bundled. Where such a source
is genuinely useful it is made available as an opt-in local fetch that writes
outside the repository tree, so a user can enrich their own installation without
that licence ever attaching to this project or to anything built from it. This is
a deliberate constraint and it costs real coverage: SNPedia in particular carries
annotations nothing else has. The alternative is worse. A single CC-BY-NC-SA file
committed to `data/` would relicense the repository and silently strip the MIT
grant from every downstream user, including anyone who sells a service built on
it.

Two practical rules follow, and both are enforced in code rather than left to
good intentions:

1. **Licence contamination travels with the column, not with the file.** A CC0
   dump can contain columns that are not CC0. The CPIC database is the live
   example: it is CC0-1.0, but its `clinpgxlevel` and `pgxtesting` columns are
   PharmGKB-sourced and carry PharmGKB's terms. `data/build_full_reference.py`
   excludes those two columns in the HTTP request itself, so they are never even
   transferred.
2. **Non-free data is written outside the tree.** The SNPedia harvester writes to
   `~/.dnainsight/`, not to `data/`. `.gitignore` carries a second line of
   defence for the cache filenames anyway, because one careless `git add -A` is
   all it takes.

v3.0 adds a third rule, which came out of HGDP and is genuinely new:

3. **A consent objection is not answered by accepting a licence.** A dataset can
   be legally open and still carry an unresolved question about what its
   participants agreed to. Those sit behind a SECOND opt-in that is separate
   from the licence gate, because clicking "I accept the terms" does not address
   an objection that was never about terms. Section 13 is the live case.

Rule 1 is now enforced at runtime as well as at build time.
`backend/provenance.py` carries `SOURCES`, a machine-readable mirror of this
document, and `licence_audit()` checks every declared bundled artefact against
it. `GET /api/v3/licence-audit` returns the result. A document nobody re-reads
drifts; this one is checked by the application it governs.

## Summary

| Source | Licence | In repo? | Used for |
|---|---|---|---|
| CPIC | CC0-1.0 | Bundled | Pharmacogenomic actionability levels |
| CPIC allele definition tables | CC0-1.0 | Fetched locally | Star-allele definitions, reconciliation |
| ClinVar | US public domain | Bundled and fetched | Clinical significance, review stars |
| gnomAD | CC0-1.0 | Bundled | Population allele frequencies |
| 1000 Genomes via Ensembl | Open, no restriction | Bundled | Per-population frequencies |
| 1000 Genomes phase 3 genotypes | Open, no restriction | Fetched locally | Reference panel for ancestry, phasing, imputation |
| SGDP public tier, 279 samples | Unrestricted | Fetched locally | Reference panel breadth |
| SGDP restricted tier, 21 samples | Signed agreement, non-commercial | **Never fetched** | Nothing, excluded by construction |
| HGDP-CEPH | Open, but consent contested | **Second opt-in only** | Optional panel breadth |
| Allen Ancient DNA Resource | Unread at review time | **Not used** | Nothing |
| GWAS Catalog | EMBL-EBI terms | Fetched locally | Replicated trait associations |
| PGS Catalog | EMBL-EBI terms, per-score overrides | Fetched locally, filtered | Polygenic score definitions |
| MyVariant.info | Apache-2.0 code, mixed data | Live API, not stored | On-demand annotation |
| SNPedia | CC-BY-NC-SA-3.0-US | **Never bundled**, opt-in only | Optional local enrichment |
| ClinPGx / PharmGKB | CC-BY-SA-4.0 plus extra clause | **Not used** | Nothing |
| PhyloTree, ISOGG, YFull | Varies, not assessed | **Not used** | Nothing yet, the verification target |
| UCSC chain files | Open, UCSC terms | User-supplied | Optional liftover |
| Local language models via Ollama | Per-model, varies | User-installed | Optional grounded assistant |

"Bundled" means the data, or a derived table, is committed to this repository.
"Fetched locally" means the user's own machine downloads it into a gitignored
path at build time. "User-supplied" means DNAInsight never fetches it at all and
the user brings their own copy.

External programs are not data and are recorded separately, in
`docs/EXTERNAL_TOOLS.md`. The rule there is the same rule as here, applied to
executables: nothing whose licence forbids redistribution or commercial use is
bundled, and the subprocess boundary is the licence boundary.

## 1. CPIC

- **Name:** Clinical Pharmacogenetics Implementation Consortium database
- **URL:** https://cpicpgx.org/ and the API at https://api.cpicpgx.org/v1/
- **Used for:** the `cpic_level` field. CPIC's A to D actionability grading drives
  the pre-prescription silo and the `cpic_a` and `cpic_b` base score tiers in
  `backend/scoring.py`. In Tier 2 the gene and drug pair table supplies
  `cpiclevel` per gene.
- **Licence identifier:** `CC0-1.0`. CPIC places its database in the public
  domain. Verified against the CPIC site and the database release notes.
- **In repo:** Yes. Levels are baked into `data/evidence_overlay.py` for Tier 1
  and fetched into the gitignored `data/reference.db` for Tier 2.
- **Redistribution constraint:** None. CC0 imposes no conditions.
- **Caveat that matters:** the CPIC dump carries `clinpgxlevel` and `pgxtesting`
  columns that are **not** CC0. They originate from PharmGKB and inherit its
  no-commercial-sale clause. Only `genesymbol`, `drugname` and `cpiclevel` are
  ever read. See section 9.

## 2. ClinVar

- **Name:** ClinVar, NCBI
- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/ and
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
- **Used for:** `clinical_sig`, `clinvar_sig_code`, `review_status`,
  `review_stars` and `condition`. ClinVar is the backbone of Tier 2 and the
  source of the star ratings that gate the scoring tiers.
- **Licence identifier:** No SPDX identifier applies. ClinVar is a US Government
  work and is in the **public domain** in the United States (17 U.S.C. 105).
  NCBI's policy **requests** citation of the database and the relevant
  publication but does **not require** attribution as a licence condition.
- **In repo:** Both. A curated subset is committed in Tier 1. The full
  `variant_summary.txt.gz` derivative lives in the gitignored
  `data/reference.db`.
- **Redistribution constraint:** None binding. Two practical cautions rather than
  legal ones: individual submitter records may include data whose onward use is
  governed by that submitter's consent terms, and clinical assertions change
  between monthly releases, so a stale copy can be actively misleading. The
  `build_date` and `clinvar_version` keys in the Tier 2 `meta` table exist so a
  database found on disk states its own age.
- **Why it is not committed in full:** the derived database exceeds GitHub's
  100 MB file limit, and it is fully reproducible from the builder.

## 3. gnomAD

- **Name:** Genome Aggregation Database, Broad Institute
- **URL:** https://gnomad.broadinstitute.org/
- **Used for:** population allele frequencies behind the `freq` field and the
  frequency facet, via `backend/frequency.py` and `data/frequencies.json`.
- **Licence identifier:** `CC0-1.0`. gnomAD releases its aggregate frequency
  data into the public domain and explicitly permits use without restriction,
  including in commercial products.
- **In repo:** Yes, as the derived frequency table.
- **Redistribution constraint:** None. The Broad requests citation as a courtesy.
  Note that only **aggregate** frequencies are CC0 and only aggregates are used;
  individual-level gnomAD data is not public and is not touched.

## 4. 1000 Genomes via Ensembl

- **Name:** 1000 Genomes Project phase 3, accessed through Ensembl
- **URL:** https://www.internationalgenome.org/ and https://rest.ensembl.org/
- **Used for:** per-population frequency breakdowns (AFR, AMR, EAS, EUR, SAS)
  that the population selector offers, complementing gnomAD.
- **Licence identifier:** No formal SPDX identifier. The 1000 Genomes Project
  data is **open with no restriction on use**, per the project's data reuse
  statement, and Ensembl distributes its own annotation with no usage
  restrictions and asks only for citation.
- **In repo:** Yes, as derived per-population frequencies.
- **Redistribution constraint:** None. Citation requested.

## 5. GWAS Catalog

- **Name:** NHGRI-EBI GWAS Catalog
- **URL:** https://www.ebi.ac.uk/gwas/ . The documented bulk endpoint
  `https://www.ebi.ac.uk/gwas/api/search/downloads/alternative` returned 404 on
  2026-07-27; the live equivalent is
  `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip`,
  whose member is `gwas-catalog-download-associations-alt-full.tsv`. The builder
  tries the documented URL first and falls back.
- **Used for:** `gwas_traits` and `gwas_studies` in Tier 2, and the
  `gwas_replicated` scoring tier. Only associations at `PVALUE_MLOG >= 7.3` with
  two or more independent studies are kept.
- **Licence identifier:** No SPDX identifier. Governed by the **EMBL-EBI terms
  of use**: the catalogue is freely available, may be used for research and
  commercial purposes, and may be redistributed. Citation of the catalogue
  publication is requested.
- **In repo:** No. Fetched locally into the gitignored `data/reference.db`. This
  is a size decision, not a licence decision.
- **Redistribution constraint:** None binding. EMBL-EBI disclaims warranty and
  asks that the catalogue be cited and not misrepresented as EBI-endorsed
  analysis.

## 6. PGS Catalog

- **Name:** Polygenic Score Catalog
- **URL:** https://www.pgscatalog.org/
- **Used for:** polygenic score definitions in `backend/prs.py` and
  `data/prs_models.json`: the variant list, effect weights and per-score
  metadata.
- **Licence identifier:** No single identifier. The catalogue as a whole is
  distributed under the **EMBL-EBI terms of use**, but individual scores carry
  **per-score licence overrides** set by the submitting authors, and those
  overrides win.
- **In repo:** Partly, and only after filtering.
- **Redistribution constraint, and this one is load-bearing:** roughly **31
  scores are CC BY-NC-ND**, which forbids both commercial use and derivative
  works. A no-derivatives term is incompatible with reweighting or subsetting a
  score, which is exactly what a scoring engine does. Those scores **must be
  filtered out** before any PGS data is committed or shipped. Do not treat the
  catalogue's default terms as applying to every score in it. Check each score's
  `license` field, and when it is anything other than the permissive default,
  exclude the score rather than trying to reason about whether the particular use
  qualifies.

## 7. MyVariant.info

- **Name:** MyVariant.info, Su and Wu labs, Scripps Research
- **URL:** https://myvariant.info/
- **Used for:** live on-demand annotation lookups in `backend/scanner.py` when a
  variant is not covered by Tier 1 or Tier 2 and the user has enabled network
  calls.
- **Licence identifier:** The **service code is `Apache-2.0`**. The **data is
  not** under a single licence: MyVariant.info is an aggregator, and each field
  keeps the licence of the upstream source it came from. A single response can
  blend US public domain ClinVar fields, CC0 gnomAD fields, and fields from
  sources with more restrictive terms.
- **In repo:** No. Responses are used transiently for display and are not written
  into any committed file.
- **Redistribution constraint:** Because the licence is per field, **no
  MyVariant.info response should ever be persisted into a bundled artefact**. Use
  it live, show it to the user, do not bake it in. If a field from it is wanted
  permanently, go to that field's upstream source directly and take the licence
  from there.

## 8. SNPedia

- **Name:** SNPedia
- **URL:** https://www.snpedia.com/
- **Used for:** optional enrichment only: plain-language variant summaries,
  genoset definitions and magnitude values, through `backend/snpedia.py`.
- **Licence identifier:** `CC-BY-NC-SA-3.0-US`. All three terms bite:
  attribution, non-commercial, and share-alike.
- **In repo: NEVER.** Not in `data/`, not in a fixture, not in a test file, not
  in a docstring example. A harvested copy inside this repository would trigger
  the share-alike term and relicense the project away from MIT, and the
  non-commercial term would strip the right to sell anything built on it.
- **How it is offered instead:** an **opt-in local harvest**. The
  `/api/admin/snpedia/harvest` endpoint refuses with HTTP 403 and a licence
  notice unless the caller passes `accept_license: true`, and the harvester
  writes its cache to `~/.dnainsight/`, deliberately outside the repository tree.
  `.gitignore` also blocks `snpedia_cache.db` and `*snpedia_cache*` as a second
  line of defence.
- **Redistribution constraint:** Full CC-BY-NC-SA obligations fall on the user
  who harvests, for their own copy. They do not touch this repository because
  this repository never holds the data.

## 9. ClinPGx and PharmGKB

- **Name:** PharmGKB, now consolidated under ClinPGx
- **URL:** https://www.pharmgkb.org/ and https://www.clinpgx.org/
- **Used for:** **nothing.** Listed here so the exclusion is a recorded decision
  rather than an oversight, and so nobody adds it later assuming it is fine.
- **Licence identifier:** Claims `CC-BY-SA-4.0`, but the accompanying data use
  agreement **adds a term prohibiting sale of the data or of products containing
  it**. That extra restriction is not part of CC-BY-SA-4.0 and cannot be added to
  it, so the effective licence is a bespoke non-commercial one and is **not** the
  standard identifier it appears to be.
- **In repo: no, and it must stay that way.** Bundling it would forbid downstream
  commercial use, which conflicts with the MIT grant this project makes.
- **Where it leaks in from:** the CPIC dump. `clinpgxlevel` and `pgxtesting`
  arrive inside an otherwise CC0 CPIC download and look safe by association. They
  are not. `data/build_full_reference.py` excludes them at the wire level by
  naming its columns in the request, and `data/evidence_overlay.py` states in its
  own header that nothing in it derives from PharmGKB bulk downloads.
- **What replaces it:** CPIC for actionability levels, ClinVar for clinical
  significance, and FDA drug label tiers for pharmacogenomic labelling. Between
  them the coverage loss is small and the licence position is clean.

## 10. CPIC allele definition tables

- **Name:** CPIC allele definition tables, the star-allele half of the CPIC
  database
- **URL:** https://api.cpicpgx.org/v1/ , the `allele_definition` and
  `allele_location_value` resources
- **Used for:** `data/build_pgx_alleles.py` builds `data/pgx_alleles.json` and
  reconciles it against the hand-written allele table in
  `backend/diplotype.py`. The reconciliation report is embedded in the output.
- **Licence identifier:** `CC0-1.0`, same as section 1, and the same column
  caveat applies. `clinpgxlevel` and `pgxtesting` are never requested.
- **In repo:** No. `data/pgx_alleles.json` is a local build artefact.
- **Three parsing traps, verified on 2026-08-04:**
  - **CPIC positions are GRCh38. `backend/diplotype.py` records GRCh37.** The
    builder therefore compares **bases only and never positions**. Comparing
    positions across builds would produce a wall of false conflicts and bury the
    real ones.
  - **CPIC uses IUPAC ambiguity codes in `variantallele`.** An ambiguity code is
    never treated as a match. "R matches A" is how a reconciliation tool talks
    itself into agreement it does not have.
  - Both sides state the **positive chromosomal strand**, so a base
    disagreement is a genuine conflict rather than a convention difference. The
    2026-08-04 run found three. They are recorded in `docs/KNOWN_GAPS.md` and
    neither file was edited to hide them.
- **Redistribution constraint:** None. CC0.

## 11. 1000 Genomes phase 3 genotypes, for the reference panel

- **Name:** 1000 Genomes Project phase 3, full genotype callset
- **URL:** IGSR at https://ftp.1000genomes.ebi.ac.uk/ , and the Beagle-hosted
  phase 3 v5a distribution at https://bochet.gcc.biostat.washington.edu/beagle/
- **Used for:** the reference panel behind ancestry, phasing and imputation,
  built by `data/build_panel.py` into `~/.dnainsight/panels/`.
- **Licence identifier:** No formal SPDX identifier. **Open with no restriction
  on use**, per the project's data reuse statement. Commercial use permitted,
  citation requested.
- **In repo:** No, and it never will be. The panel is tens of GB and is built on
  the user's own machine.
- **The trap that decided the default, verified on 2026-08-04:** the IGSR
  20130502 release genotype VCFs publish `ID` as `.` for **every record**,
  confirmed across all 1,103,547 chr22 records. An rsID-keyed feature gets
  nothing from them, which means `--array-file` matches nothing and
  `informative_markers.tsv` cannot be written. `build_panel.py` therefore
  defaults to `--onekg-source beagle`, the same phase 3 v5a callset with IDs
  populated. `igsr` remains selectable and prints the caveat.

## 12. Simons Genome Diversity Project

- **Name:** SGDP, Simons Genome Diversity Project
- **URL:** https://reichdata.hms.harvard.edu/pub/datasets/sgdp/
- **Used for:** population breadth in the reference panel, alongside 1000
  Genomes.
- **Licence identifier:** No SPDX identifier. **Two tiers with different terms,
  and the split is the whole point of this entry.**
  - The **279-sample public tier is unrestricted** and is what DNAInsight uses.
  - The **21-sample restricted tier** sits behind a signed agreement whose terms
    include "I will not use the data for any commercial purposes".
- **In repo:** No. Fetched into `~/.dnainsight/panels/` at build time.
- **Redistribution constraint:** **The restricted tier is excluded by
  construction.** The builder does not fetch it, does not offer a flag to fetch
  it, and there is no consent path that enables it. A non-commercial term inside
  a bundled panel is the same failure as a non-commercial term inside `data/`.
- **Known unverified detail:** the per-sample VCF filename template is
  **assumed, not verified**, because the Harvard host's certificate chain did not
  validate at review time. `--sgdp-vcf-template`, `--sgdp-dir` and
  `--sgdp-metadata` all override it. Recorded in `docs/KNOWN_GAPS.md`.

## 13. HGDP-CEPH

- **Name:** Human Genome Diversity Project, CEPH panel, accessed via IGSR
- **URL:** https://www.internationalgenome.org/data-portal/data-collection/hgdp
- **Used for:** optional additional population breadth in the reference panel.
  Off by default.
- **Licence identifier:** No SPDX identifier. **Open access under the Fort
  Lauderdale Principles. The licence is not the problem.**
- **Why it is gated anyway:** two published events, both recorded with dates so
  the decision can be re-examined rather than inherited.
  - **Nature Genetics, 24 November 2025** concluded that broad reuse of HGDP may
    diverge from what participants consented to.
  - **The PRIMED Consortium voted on 21 August 2024** to keep permitting its
    use, while acknowledging "failure to obtain informed consent consistent with
    current standards from many participants".
  Those two facts are recorded together because they point in opposite
  directions and a reader is entitled to both.
- **In repo:** No, and not in the default panel either.
- **How it is offered instead:** a **second opt-in**, separate from the licence
  gate. `data/build_panel.py --include-hgdp` is refused unless
  `--accept-consent-caveat` is also passed, and `--accept-terms` does not imply
  it. This is rule 3 in practice: accepting a licence does not answer a consent
  objection, so the two acknowledgements are not allowed to be the same click.

## 14. Allen Ancient DNA Resource

- **Name:** AADR, Allen Ancient DNA Resource
- **Used for:** **nothing.** Listed so the exclusion is a recorded decision.
- **Licence identifier:** **CC0 1.0**, read 2026-08-10. The Harvard Dataverse
  record for the AADR states `License/Data Use Agreement: CC0 1.0`. CC0 already
  sat inside the original bundling rule, so no part of the CC-BY amendment is
  needed to accept it.
- **Citation:** Mallick S, Micco A, Mah M, Ringbauer H, Lazaridis I, Olalde I,
  Patterson N, Reich D (2024). The Allen Ancient DNA Resource (AADR): a curated
  compendium of ancient human genomes. *Scientific Data* 11, 182.
- **In repo:** No, and it is not fetchable through any flag.
- **Why the refusal STANDS anyway, decided 2026-08-10:** the exclusion rested on
  two grounds and reading the licence answered only the first.

  1. *The terms were not readable.* **Answered.** They are CC0 1.0.
  2. *A compendium cannot grant rights its components did not grant.*
     **Unanswered.** The AADR aggregates many published datasets, each with its
     own upstream terms, and the Dataverse record speaks for the compendium
     rather than for the constituent studies. Nothing read so far addresses
     whether every contributing study permitted redistribution under CC0.

  A licence declaration by an aggregator is evidence about the aggregator's
  intent, not proof that the underlying rights existed to be granted. Ground 2
  is the load-bearing one and it is still open, so `REFUSAL_AADR` remains in
  force. This is recorded as a decision with a reason rather than an absence, so
  that the next reader does not re-litigate ground 1 and mistake it for the
  whole question.
- **What would settle it:** an audit of the constituent studies' own terms, or a
  statement from the compendium's maintainers that redistribution rights were
  obtained for every included dataset.

## 15. PhyloTree, ISOGG and YFull

- **Name:** PhyloTree mtDNA build, the ISOGG Y-DNA haplogroup tree, YFull YTree
- **Used for:** **nothing yet, and that is a defect rather than a policy.**
  `backend/haplogroups.py` ships a backbone tree of 49 Y markers and 28 mtDNA
  nodes that was written from recall and **not machine-checked against any of
  these three**. Every Y marker and 11 of the 28 mtDNA nodes are flagged
  `verified: false`, `unverified_markers()` returns the audit list, and
  `verified_only=true` refuses them all.
- **Licence identifier:** Not assessed. Each has its own terms and none has been
  read against the five questions below.
- **In repo:** No data from any of them.
- **What has to happen before they can be used:** answer the questions in
  "Adding a new source" for each tree, in this file, then verify the backbone
  row by row and flip each `verified` flag with the source recorded. Until then
  the haplogroup call ships flagged provisional rather than shipping wrong and
  silent. See `docs/KNOWN_GAPS.md`.

## 16. UCSC chain files

- **Name:** UCSC liftOver chain files, for example hg38ToHg19.over.chain.gz
- **URL:** https://hgdownload.soe.ucsc.edu/goldenPath/
- **Used for:** optional coordinate liftover in `backend/sequencing.py`. The
  chain parser is DNAInsight's own and is MIT; only the chain data is external.
- **Licence identifier:** UCSC genome data is free for all uses; the Genome
  Browser software itself carries a separate licence which DNAInsight does not
  use or ship.
- **In repo:** No. **User-supplied.** DNAInsight never downloads a chain file.
  Drop one into `~/.dnainsight/panels/chains/`.
- **Behaviour without one:** `liftover()` returns the standard unavailable
  payload. It does **not** pass coordinates through untranslated, which would
  produce a file that is wrong at every position while looking converted.

## 17. Local language models via Ollama

- **Name:** whatever model the user pulls, for example `llama3.1:8b`
- **Used for:** the optional grounded assistant in `backend/assistant.py`.
- **Licence identifier:** **Per model, and DNAInsight asserts nothing about
  any of them.** Ollama itself is MIT. Model weights carry their own terms,
  several of which are neither open source nor commercial-use-clean.
- **In repo:** No weights, no model files, no prompt caches, nothing.
- **Redistribution constraint:** falls entirely on the user who pulls the model.
  DNAInsight sends it finding text and citations over loopback and stores none
  of its output as a bundled artefact.

## 18. NCBI dbSNP

- **Name:** NCBI dbSNP, queried through E-utilities `esummary`, `db=snp`
- **Used for:** auditing the Y backbone in `backend/haplogroups.py`. Variant
  class, chromosome, GRCh38 and GRCh37 positions, and the reference/alternate
  allele set. `tools/audit_y_dbsnp.py` runs it and `docs/Y_BACKBONE_AUDIT.md`
  records the result.
- **Licence identifier:** US Government work, **public domain**. NCBI databases
  carry no copyright restriction on their content, which puts dbSNP inside the
  bundling rule without any amendment.
- **In repo:** derived values only, and only where they are facts dbSNP can
  actually settle: `assembly`, `ref_carries`, `variant_type`, `ancestral_seq`
  and `derived_seq` on audited backbone rows, plus the positions quoted in
  `docs/Y_BACKBONE_AUDIT.md`. No bulk extract of dbSNP is stored.
- **What it cannot settle, recorded because the distinction is load-bearing:**
  dbSNP reports **reference over alternate**. The backbone records **ancestral
  over derived**. On the Y these routinely disagree, because the GRCh38
  reference Y descends from a lineage carrying the derived allele at many
  backbone nodes. The audit measured it: the reference carries the derived
  allele at **10 of the 17** nodes where the state is determinable. So no dbSNP
  value is written into an `ancestral` or `derived` field, and no row was marked
  `verified` on the strength of the audit.
- **What it cannot reach:** dbSNP does not index Y marker NAMES. `esearch` for
  M91, M175 and M267 each returns zero hits, so the 31 backbone markers that
  carry no rsID could not be audited at all and are reported as unaudited rather
  than assumed correct.

## Adding a new source

Before committing any new dataset, answer all six in writing in this file:

1. What is the **exact** licence identifier, and does the accompanying data use
   agreement add terms the identifier does not carry? Section 9 exists because
   of that second half.
2. Is it **per file or per record**? Section 6 exists because of that.
3. Are any **columns** licensed differently from the file? Section 1 exists
   because of that.
4. Does it permit **commercial use and derivative works**? If either answer is
   no, it is a local opt-in fetch, not a bundled file.
5. Is the derived artefact under 100 MB? If not it is gitignored and rebuilt by a
   script, regardless of licence.
6. Is there an **objection to reuse that the licence does not settle**, such as
   a documented consent problem? Sections 12, 13 and 14 exist because of that.
   A tiered dataset gets only its clean tier, a contested one gets a second
   opt-in of its own, and one whose terms nobody has read gets excluded until
   somebody reads them.

Then add it to `SOURCES` in `backend/provenance.py` with the same answers, and
list any bundled artefact it feeds in `BUNDLED_ARTEFACTS`. An artefact fed by a
source id that is not in `SOURCES` is a licence audit **violation**, not a
warning, because an unassessed source is exactly how contamination enters.
