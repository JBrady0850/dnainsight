# DNAInsight Data Sources and Licences

Every external dataset DNAInsight touches, what it is used for, and what may
legally be done with it. Referenced from `CHANGELOG.md`.

Last reviewed: 2026-07-27.

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

## Summary

| Source | Licence | In repo? | Used for |
|---|---|---|---|
| CPIC | CC0-1.0 | Bundled | Pharmacogenomic actionability levels |
| ClinVar | US public domain | Bundled and fetched | Clinical significance, review stars |
| gnomAD | CC0-1.0 | Bundled | Population allele frequencies |
| 1000 Genomes via Ensembl | Open, no restriction | Bundled | Per-population frequencies |
| GWAS Catalog | EMBL-EBI terms | Fetched locally | Replicated trait associations |
| PGS Catalog | EMBL-EBI terms, per-score overrides | Fetched locally, filtered | Polygenic score definitions |
| MyVariant.info | Apache-2.0 code, mixed data | Live API, not stored | On-demand annotation |
| SNPedia | CC-BY-NC-SA-3.0-US | **Never bundled**, opt-in only | Optional local enrichment |
| ClinPGx / PharmGKB | CC-BY-SA-4.0 plus extra clause | **Not used** | Nothing |

"Bundled" means the data, or a derived table, is committed to this repository.
"Fetched locally" means the user's own machine downloads it into a gitignored
path at build time.

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

## Adding a new source

Before committing any new dataset, answer all five in writing in this file:

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
