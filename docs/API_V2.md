---
created: 2026-07-26
modified: 2026-07-26
tags: [archivist, dnainsight, reference, decision]
aliases: [DNAInsight v2 API Contract, API_V2]
---

# DNAInsight v2.0 API Contract

Single source of truth for the v2 HTTP surface and the finding object. The
backend (`backend/routes.py`) and the single-page app (`frontend/index.html`)
are both built to this document. If they disagree, this document is right and
the code is wrong.

Base URL when running locally: `http://127.0.0.1:5050`

## 0. Design constraints carried over from v1

1. Offline first. Every endpoint below works with no internet except the ones
   explicitly marked NETWORK.
2. No genotype ever leaves the machine. The MyVariant path sends rsIDs only.
   The SNPedia harvester sends page titles only.
3. Nothing SNPedia-derived is ever written inside the repository. The harvest
   cache lives at `~/.dnainsight/snpedia_cache.db`.

## 1. Entity types

Every finding carries `entity_type`, which tells the UI how to render it and
which filters legitimately apply.

| entity_type | rsid field holds | Frequency applies | Publications apply | Repute applies |
|---|---|---|---|---|
| `snp`     | an rsID, e.g. `rs1801133` | yes | yes | yes |
| `genoset` | a genoset name, e.g. `dgs004` | no  | no  | yes |
| `trait`   | a trait key, e.g. `lactase_persistence` | no | no | no, always empty |
| `prs`     | a model id, e.g. `t2d` | no | no | no, always empty |

Genosets, traits and polygenic scores are exempt from the frequency and
publication filters by design. A UI that filters them out with a frequency
slider is broken, because they have no single position to have a frequency at.

## 2. The finding object

Returned by `/findings`, and embedded in reports. Every key is always present.
Absent data is `null`, `""`, `[]` or `false`, never a missing key, so template
code cannot raise.

### 2.1 Identity and genotype

| Key | Type | Notes |
|---|---|---|
| `rsid` | string | lowercase; or the genoset / trait / model id |
| `entity_type` | string | one of the four above |
| `gene` | string | HUGO symbol, `""` for genosets |
| `chromosome` | string | `""` for genosets |
| `position` | integer | GRCh37, `0` when not applicable |
| `allele1`, `allele2` | string | oriented to the reference, `"N"` for a no-call |
| `genotype` | string | two characters, e.g. `"CT"` |
| `token` | string | SNPedia-style, e.g. `"(C;T)"` |
| `zygosity` | string | `homozygous`, `heterozygous`, `hemizygous`, `no_call` |

### 2.2 Interest and direction (the Magnitude and Repute equivalents)

Computed locally from CC0 and public-domain evidence. See section 6 for the
formula. When a local SNPedia cache exists and the user has opted in, the
cached values override the computed ones and `magnitude_source` says so.

| Key | Type | Notes |
|---|---|---|
| `magnitude` | float or null | 0 to 10. `null` means unscored, which the UI must treat as 1 for sorting, matching the documented convention that a blank magnitude sorts as 1 |
| `magnitude_source` | string | `computed`, `snpedia`, or `""` |
| `repute` | string | `Good`, `Bad` or `""`. Never set for `trait` or `prs` |
| `summary` | string | one line |
| `interpretation` | string | the longer body text |
| `confidence` | string | `high`, `moderate`, `low`, `none` |

### 2.3 Evidence

| Key | Type | Notes |
|---|---|---|
| `clinical_sig` | string | ClinVar germline classification, lowercase |
| `clinvar_sig_code` | integer or null | 5 pathogenic, 4 likely pathogenic, 3 likely benign, 2 benign, 1 uncertain, 6 drug response, 7 histocompatibility, 255 other |
| `review_status` | string | verbatim ClinVar review status string |
| `review_stars` | integer | 0 to 4 |
| `cpic_level` | string | `A`, `A/B`, `B`, `B/C`, `C`, `C/D`, `D`, `Retired` or `""` |
| `pgx_level` | string | FDA label tier or `""` |
| `evidence` | string | short human label, e.g. `CPIC Level A` |
| `publications` | integer | citation count, `0` when unknown |
| `conditions` | string | joined condition names |
| `conditions_list` | array of string | |
| `sources` | array of string | e.g. `["bundled_reference", "clinvar"]` |

### 2.4 Population frequency

Set by `backend/frequency.py`. `freq` is a percentage 0 to 100 for the user's
exact genotype in the selected population.

| Key | Type | Notes |
|---|---|---|
| `freq` | float or null | `0.0` means not observed in that panel. `null` means unknown. These are different facts and the UI must render them differently |
| `freq_population` | string | population code the value came from |
| `freq_band` | string | `very_rare`, `rare`, `uncommon`, `common`, `majority`, `unknown` |
| `freq_color` | string | 7-character hex for the rarity heat |
| `freq_derived` | boolean | true when Hardy-Weinberg was used |
| `freq_method` | string | `observed`, `hardy_weinberg`, `unavailable` |
| `freq_flipped` | boolean | the alleles were complemented to match the table |
| `freq_ambiguous` | boolean | palindromic site, the frequency may belong to the other reading |
| `gmaf` | float or null | global minor allele frequency, fraction 0 to 1 |
| `minor_allele` | string | |
| `population_series` | array | `[{code, label, brief, frequency, yours}]` |

### 2.5 Strand and quality

The single largest correctness risk in this class of tool. Surface it, do not
hide it.

| Key | Type | Notes |
|---|---|---|
| `orientation` | string | `plus`, `minus` or `""` |
| `stabilized_orientation` | string | governs genotype matching |
| `flipped` | boolean | a complement was applied during annotation |
| `ambiguous` | boolean | A/T or C/G heterozygote, strand cannot be settled |
| `dubious` | boolean | the call is suspect for any reason |
| `variant_allele` | string | authoritative alternate allele when known |
| `variant_copies` | integer or null | copies the person carries, 0, 1 or 2 |
| `carrier` | boolean or null | `false` means the classification does not apply to them |

### 2.6 Multi-file, family and grouping

| Key | Type | Notes |
|---|---|---|
| `count` | integer | how many pooled source files produced a real call here |
| `labels` | array of string | which files contributed |
| `conflict` | boolean | pooled files disagree at this position |
| `calls` | array | `[{label, allele1, allele2, genotype}]`, all retained, never reconciled |
| `comparison` | array | `[{label, role, genotype, shared}]` for non-self sources |
| `probability` | float or null | offspring probability when both parents are loaded |
| `mendelian_ok` | boolean or null | |
| `topics` | array of string | |
| `medicines` | array of string | |

### 2.7 Genoset, trait and PRS extras

| Key | Type | Applies to |
|---|---|---|
| `criteria` | string | genoset, the rule text |
| `matched_rsids` | array of string | genoset |
| `coverage` | float | genoset and prs, 0 to 1 |
| `percentile` | float or null | prs |
| `band` | string | prs, `low` through `high` |
| `reliable` | boolean | prs, false below 0.90 coverage |
| `caveats` | array of string | prs and trait |

## 3. Endpoints

### 3.1 System

| Method | Path | Notes |
|---|---|---|
| GET | `/api/status` | `{status, version}` |
| GET | `/api/version` | `{version, snp_count, genoset_count, prs_count, trait_count, frequency_rsids}` |
| GET | `/api/populations` | `[{code, label, brief, superpop, available}]` plus `default` and the `MAX`, `AVG`, `MIN`, `GLOBAL` aggregate modes |
| GET | `/api/capabilities` | which subsystems have data: `{frequency, genosets, prs, traits, snpedia, api}`. The UI hides controls whose data is absent rather than showing dead sliders |

### 3.2 Profiles and sources

| Method | Path | Body | Notes |
|---|---|---|---|
| GET | `/api/profiles` | | includes `findings_summary` and `source_count` |
| POST | `/api/profiles` | multipart `name`, `dob`, `sex`, `file`, optional `label` | creates the profile with one `self` source |
| GET | `/api/profiles/<pid>` | | |
| DELETE | `/api/profiles/<pid>` | | |
| GET | `/api/profiles/<pid>/sources` | | `[{id, label, role, provider, snp_count, contributed, overlapped, conflicting, uploaded_at}]` |
| POST | `/api/profiles/<pid>/sources` | multipart `file`, `role`, optional `label` | `role` is one of `self`, `mother`, `father`, `mate`, `child`, `sibling`, `other`, `ignore`. Adding a `self` source POOLS it. Any other role adds comparison rows only |
| PATCH | `/api/profiles/<pid>/sources/<sid>` | `{role}` or `{label}` | re-role without re-uploading |
| DELETE | `/api/profiles/<pid>/sources/<sid>` | | |

Pooling contract, inherited from how the reference product actually behaves:
conflicting calls between two `self` files are BOTH retained and surfaced. There
is no voting, no confidence weighting and no automatic winner.

### 3.3 Scan

`POST /api/profiles/<pid>/scan`

```json
{
  "use_api": true,
  "population": "CEU",
  "include_genosets": true,
  "include_traits": true,
  "include_prs": true,
  "use_snpedia": false
}
```

`use_api` is NETWORK. `use_snpedia` reads the local cache only and is ignored
when the cache is absent.

`GET /api/profiles/<pid>/scan/status` returns
`{running, done, phase, processed, total, findings, error}` where `phase` is one
of `bundled`, `orientation`, `frequency`, `genosets`, `traits`, `prs`, `api`,
`writing`, `complete`.

### 3.4 Findings and filtering

`GET /api/profiles/<pid>/findings`

Every filter is a query parameter and every one is optional. Omitting a
parameter means "do not filter on this".

| Parameter | Type | Notes |
|---|---|---|
| `silo` | string | `pre_prescription`, `actionable`, `informational` |
| `entity_type` | csv | e.g. `snp,genoset` |
| `min_magnitude`, `max_magnitude` | float | a null magnitude is treated as 1 |
| `repute` | csv | any of `Good`, `Bad`, `unset`. Default all three |
| `min_publications`, `max_publications` | int | NOT applied to genoset, trait or prs |
| `min_freq`, `max_freq` | float | percentage. NOT applied to genoset, trait or prs |
| `require_frequency` | bool | drop findings with no frequency. Default false |
| `clinvar_sig` | csv of codes | e.g. `5,4`. Default when `clinvar_only=true` is `5,4` |
| `clinvar_only` | bool | restrict to findings with any ClinVar record |
| `min_stars` | int | 0 to 4 |
| `gene` | csv | |
| `topic`, `medicine`, `condition` | csv | |
| `zygosity` | csv | `homozygous`, `heterozygous`, `hemizygous`, `no_call` |
| `carrier_only` | bool | drop non-carrier reference genotypes |
| `conflicts_only` | bool | pooled disagreements only |
| `ambiguous_only` | bool | palindromic or flipped calls only, for QC |
| `q` | string | free text, see 3.5 |
| `population` | string | recomputes frequency in this population |
| `sort` | string | see below |
| `order` | string | `asc` or `desc` |
| `limit`, `offset` | int | `limit` default 200, `0` means all |

`sort` accepts `magnitude`, `frequency`, `publications`, `location`, `modified`,
`gmaf`, `stars`, `gene`, `rsid`. Each pairs with `order` to give both directions,
which is how the reference product exposed 14 separate sort options.

Response:

```json
{
  "findings": [],
  "total": 412,
  "returned": 200,
  "offset": 0,
  "summary": {"pre_prescription": 15, "actionable": 23, "informational": 14},
  "filtered_summary": {},
  "ranges": {"magnitude": [0, 8.5], "publications": [0, 412], "frequency": [0, 99.4]},
  "population": "CEU",
  "qc": {"total": 412, "flipped": 37, "ambiguous": 13, "unknown_orientation": 4}
}
```

`ranges` exists so the UI can set slider bounds from the data rather than
hard-coding them, and so Reset restores data-derived bounds.

### 3.5 Free-text query grammar for `q`

Operators are stripped from the text before the remainder is matched as a
case-insensitive regular expression across rsid, gene, summary, interpretation,
conditions, topics, medicines and the genotype token.

| Form | Effect |
|---|---|
| `chr7` | everything on chromosome 7 |
| `chr7:1234` | that exact position |
| `chr7:1000-2000` | inclusive range |
| `chr7:1000+500` | position plus offset |
| `/CLNSIG=5,4` | ClinVar significance whitelist |
| `/STARS>=2` | minimum review stars |
| `/MAG>=3` | minimum magnitude |
| `/COUNT>=2` | minimum pooled call count |
| `/dubious` | only suspect calls |
| `/flipped` | only strand-flipped calls |
| `/ambiguous` | only palindromic calls |

### 3.6 Facets

`GET /api/profiles/<pid>/facets` returns every filter value with its count, so
each dropdown can render `Name (n)`:

```json
{"genes": [{"value": "MTHFR", "count": 2}],
 "topics": [], "medicines": [], "conditions": [],
 "clinvar_diseases": [], "silos": [], "entity_types": [],
 "zygosity": [], "categories": []}
```

### 3.7 Subsystem views

| Method | Path | Notes |
|---|---|---|
| GET | `/api/profiles/<pid>/genosets` | `{matched, unmatched, incomplete}`. `incomplete` means required rsIDs were not on the array, which is "not testable" and must not be shown as "absent" |
| GET | `/api/profiles/<pid>/traits` | `{traits, blood_type}` |
| GET | `/api/profiles/<pid>/prs` | `{results, disclaimer}` |
| GET | `/api/profiles/<pid>/pgx` | pharmacogenomic view grouped by drug |
| GET | `/api/profiles/<pid>/conflicts` | pooled disagreements |
| GET | `/api/profiles/<pid>/trio` | `{trio_available, compared, violations, violation_rsids, note}` |
| GET | `/api/profiles/<pid>/qc` | strand, flip, ambiguity and no-call report |
| GET | `/api/profiles/<pid>/lookup/<rsid>` | single-SNP lookup, no scan required |

### 3.8 Export and reports

| Method | Path | Notes |
|---|---|---|
| GET | `/api/profiles/<pid>/export/json` | honours every `/findings` filter |
| GET | `/api/profiles/<pid>/export/csv` | as above |
| GET | `/api/profiles/<pid>/export/tsv` | as above |
| POST | `/api/profiles/<pid>/reports` | `{type}` one of `genetic`, `doctor`, `interactive` |
| GET | `/api/reports/<rid>/view` | serves the HTML |
| GET | `/api/reports/<rid>/download` | as an attachment |

`interactive` produces a single self-contained HTML file with the filter engine
and the data embedded, openable with no server and no internet.

### 3.9 Admin

| Method | Path | Notes |
|---|---|---|
| GET | `/api/admin/db-status` | bundled reference metadata and staleness |
| POST | `/api/admin/update-databases` | NETWORK, refresh ClinVar significance and review stars |
| GET | `/api/admin/update-databases/status` | |
| GET | `/api/admin/snpedia/status` | `{available, path, snps, genotypes, genosets, last_harvest, license, notice}` |
| POST | `/api/admin/snpedia/harvest` | NETWORK, `{accept_license: true, scope}`. Returns 403 with the licence notice when `accept_license` is not true |
| GET | `/api/admin/snpedia/harvest/status` | |
| DELETE | `/api/admin/snpedia/cache` | purge the local cache |

## 4. Error shape

Every error is `{"error": "human readable sentence"}` with a real status code.
`400` bad input, `403` licence not accepted, `404` unknown id, `409` an
operation is already running, `413` upload too large, `500` unexpected.

## 5. Frontend responsibilities

The findings view is ONE virtualised, filtered, sorted list of cards, not a
chaptered document. Specifically required:

1. Card border colour encodes repute: green `#60B060` Good, red `#FF9090` Bad,
   grey `#C0C0C0` unset. A colourblind mode swaps to `#998EC3` and `#F1A340`.
2. Three dual-handle range sliders: magnitude, publications, frequency. Each
   with plus and minus nudge buttons and data-derived reset bounds.
3. Repute tri-state checkboxes. All three on by default.
4. ClinVar mode: All, ClinVar only, plus a dropdown of the nine significance
   codes. Filtering defaults to pathogenic and likely pathogenic only.
5. Searchable multi-select facets for gene, topic, medicine and condition, each
   showing its own count.
6. Sort dropdown with both directions for every key in 3.4.
7. Reference population selector that recomputes frequency live, plus the
   aggregate modes.
8. Show and Require checkbox rows: genosets, homozygous, heterozygous,
   conflicts, require frequency.
9. Visible counters: allowed, visible, offscreen, with a "show twice as many"
   button.
10. A repute distribution chart and a per-category stacked chart.
11. Table view with CSV and TSV export honouring active filters.
12. A QC banner when flipped or ambiguous counts are non-zero, because a user
    reading a palindromic site as settled fact is the main way this class of
    tool misleads people.
13. Reset clears every filter. Escape also resets. Ctrl+H opens help.

## 6. How magnitude and repute are computed without SNPedia

SNPedia's Magnitude and Repute are hand-curated and licensed CC-BY-NC-SA-3.0-US,
so they cannot ship in this repository. DNAInsight computes its own from CC0 and
public-domain inputs. The scale is deliberately the same 0 to 10 shape so the
numbers are legible to anyone who has read a Promethease report, but they are
NOT the same numbers, and the UI labels them as DNAInsight magnitude.

Base, by the strongest available evidence:

| Evidence | Base |
|---|---|
| CPIC Level A, or ClinVar pathogenic at 3 stars or better | 6.0 |
| CPIC Level B, or ClinVar pathogenic at 2 stars | 4.5 |
| FDA label tier Testing Required or Testing Recommended | 4.0 |
| ClinVar likely pathogenic at 2 stars or better | 3.5 |
| Replicated GWAS association | 2.5 |
| ClinVar single submitter, or uncertain significance | 1.5 |
| Everything else | 1.0 |

Adjustments, applied in order and clamped to 0 through 10:

1. Carrier status. Not a carrier of the reported variant: multiply by 0.25.
   Homozygous for it: multiply by 1.3. This is the single biggest honesty
   improvement over a naive report, which will happily alarm someone about a
   variant they do not carry.
2. Rarity. A `freq_band` of `very_rare` adds 0.5, `rare` adds 0.25,
   `majority` subtracts 0.5.
3. Publications. Add `min(1.0, log10(1 + publications) / 2)`.
4. No-call: force to 0.0.
5. Strand ambiguous: cap at 2.0 and set `dubious` true, because an unverifiable
   call must not outrank a verifiable one.

Repute:

- `Bad` when the direction of effect is harmful and the person carries the
  allele: pathogenic, likely pathogenic, risk factor, increased toxicity, poor
  metabolizer, contraindicated, loss of function.
- `Good` when protective, normal function, favourable response, or a documented
  reduced-risk genotype.
- `""` for anything neutral, ancestral, trait-like or conflicting, and ALWAYS
  for `trait` and `prs`. A trait is not good or bad.

`confidence` is `high` at 3 stars or better or CPIC A, `moderate` at 2 stars or
CPIC B, `low` at 1 star or a single submitter, and `none` for a no-call or an
unresolvable strand.

## 7. Compatibility with v1.2

`GET /api/profiles/<pid>/findings` still returns a `findings` array whose objects
contain every v1.2 key, so an existing consumer keeps working. The v2 keys are
additive. `POST /api/profiles/<pid>/scan` still accepts a body of only
`{"use_api": bool}`.
