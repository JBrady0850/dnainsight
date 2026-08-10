# Y_BACKBONE audit against NCBI dbSNP

Run 2026-08-10 against `backend/haplogroups.py` at `b9ef689` (v3.2.0).
Source: NCBI dbSNP via E-utilities `esummary`, `db=snp`, retrieved live.
dbSNP is a US Government work and is public domain, so its values may be
bundled. See `data/DATA_SOURCES.md`.

**No row of `Y_BACKBONE` was modified. Every entry remains `verified: false`.**

## What this audit can and cannot settle

dbSNP settles **variant class**, **chromosome**, **GRCh38 and GRCh37 position**
and the **reference and alternate allele set**.

dbSNP cannot settle **ancestral against derived**. It reports reference over
alternate; the backbone records ancestral over derived. On the Y these routinely
disagree, because the reference Y descends from a lineage carrying the derived
allele at many backbone nodes. That is precisely what the `ref_carries` field
shipped in v3.1.1 exists to guard, and this audit therefore reports which state
the reference carries rather than deciding the assignment.

## Coverage: 18 of 49 markers

`Y_BACKBONE` holds 50 entries, one of which is the marker-free root. Of the 49
real markers, **18 carry an rsID and 31 do not**. Only the 18 could be audited.

dbSNP cannot resolve a marker name to an rsID. `esearch` on `db=snp` for `M91`,
`M175` and `M267` each returned **0 hits**. The name to rsID mapping lives in
sources this project cannot bundle or does not yet hold, so the remaining 31 are
recorded as unaudited rather than guessed at.

## Result summary, 18 markers

| Verdict | Count | Markers |
|---|---|---|
| consistent | 12 | M168, M96, M89, M201, M170, M253, M304, M9, M242, M207, M173, U106 |
| consistent, multi-allelic in dbSNP | 4 | M45, M343, M269, P312 |
| **CLASS ERROR** | **1** | **M17** |
| **CONFLICT** | **1** | **M20** |

## Finding 1: M17 is an indel, not a base substitution. Confirmed.

`Y_BACKBONE["R1a1a"]` records marker M17 as ancestral **G**, derived **A**.
dbSNP records rs3908 as:

```
snp_class : delins
SPDI      : NC_000024.10:19571278:GGGG:GGG
HGVS      : NC_000024.10:g.19571282del, NC_000024.9:g.21733168del
SEQ       : [G/-]        LEN: 4        GENE: TXLNGY
GRCh38    : Y:19571279          GRCh37: Y:21733165
```

M17 is a **single-base deletion inside a four-base G homopolymer**, GGGG to GGG.
It is not a G to A substitution and there is no A allele at that site. The entry
is wrong in kind, not merely in value.

This matters more than one row. **M17 defines R1a1a**, one of the most common Y
haplogroups in Europe and South Asia. A genotyping rule expecting a G or A base
call at that position cannot fire correctly against a deletion, so the node is
currently unreachable by the logic that is supposed to reach it.

This is the first of the four indel leads to be corroborated against a citable
accession rather than recall. The Karafet 2008 article text already corroborated
M91 independently.

## Finding 2: M20 disagrees with dbSNP on the allele pair.

`Y_BACKBONE["L"]` records marker M20 as ancestral **A**, derived **C**.
dbSNP records rs3911 as:

```
snp_class : snv
SPDI      : NC_000024.10:19571567:A:G     (single, not multi-allelic)
HGVS      : NC_000024.10:g.19571568A>G
GRCh38    : Y:19571568          GRCh37: Y:21733454
```

dbSNP carries **A/G**. The table carries **A/C**. There is no C allele in the
dbSNP record, and complementing does not reconcile the two either. One of the
two is wrong and this audit cannot say which.

**Inference, based on available evidence:** the rsID assignment itself may be
the error rather than the alleles. rs3911 sits 289 bases from rs3908 in the same
gene, and both were recorded from literature recall in the same pass. This may
vary and needs a primary source to settle.

## Finding 3: the reference Y carries the derived allele at 10 of 17 nodes.

Of the 17 markers where the reference state could be determined, M17 excluded as
an indel:

| Reference carries | Count | Markers |
|---|---|---|
| **derived** | **10** | M168, M96, M89, M253, M9, M45, M207, M173, M343, M269 |
| ancestral | 7 | M201, M170, M304, M20, M242, U106, P312 |

The v3.1.1 CHANGELOG claimed a builder mapping `ref` onto `ancestral` "would
have inverted roughly half the tree with the entire suite still green". That was
an argument at the time. It is now a measurement: **59 percent of the auditable
backbone, 10 of 17 nodes**. The `ref_carries` guard is doing exactly the work it
was added to do, and it must stay in force before any dbSNP-derived value is
written into an `ancestral` or `derived` field.

## Finding 4: four markers are multi-allelic in dbSNP.

M45, M343, M269 and P312 each carry more than two alleles at the site:

| Marker | rsID | dbSNP alleles | Table pair |
|---|---|---|---|
| M45 | rs2032631 | A / G / T | G > A |
| M343 | rs9786184 | A / C / G | C > A |
| M269 | rs9786153 | C / G / T | T > C |
| P312 | rs34276300 | A / C / G | A > G |

In every case the table's pair is a subset of the observed set, so none is a
conflict. This is recorded because a first pass of this audit parsed only the
first SPDI, blanked the allele pair for all four, and reported them as
conflicts. That draft was wrong and was discarded before it was acted on. Any
future builder reading the `spdi` field must parse the whole comma-separated
list.

## Full table

| Node | Marker | rsID | Class | GRCh38 | GRCh37 | dbSNP ref/alt | Table anc>der | Verdict | Ref carries |
|---|---|---|---|---|---|---|---|---|---|
| CT | M168 | rs2032595 | snv | Y:12702062 | Y:14813991 | T/C | C>T | consistent | derived |
| E | M96 | rs9306841 | snv | Y:19617112 | Y:21778998 | C/G | G>C | consistent | derived |
| F | M89 | rs2032652 | snv | Y:19755427 | Y:21917313 | T/C | C>T | consistent | derived |
| G | M201 | rs2032636 | snv | Y:12915617 | Y:15027529 | G/T | G>T | consistent | ancestral |
| I | M170 | rs2032597 | snv | Y:12735858 | Y:14847792 | A/C | A>C | consistent | ancestral |
| I1 | M253 | rs17307677 | snv | Y:16050954 | Y:18162834 | T/C | C>T | consistent | derived |
| J | M304 | rs13447352 | snv | Y:20587967 | Y:22749853 | A/C | A>C | consistent | ancestral |
| K | M9 | rs3900 | snv | Y:19568371 | Y:21730257 | G/C | C>G | consistent | derived |
| L | M20 | rs3911 | snv | Y:19571568 | Y:21733454 | A/G | A>C | **CONFLICT** | ancestral |
| P | M45 | rs2032631 | snv | Y:19705901 | Y:21867787 | A/G/T | G>A | consistent, multi | derived |
| Q | M242 | rs8179021 | snv | Y:12906671 | Y:15018582 | C/T | C>T | consistent | ancestral |
| R | M207 | rs2032658 | snv | Y:13470103 | Y:15581983 | G/A | A>G | consistent | derived |
| R1 | M173 | rs2032624 | snv | Y:12914512 | Y:15026424 | C/A | A>C | consistent | derived |
| R1a1a | M17 | rs3908 | **delins** | Y:19571279 | Y:21733165 | GGGG/GGG | G>A | **CLASS** | undetermined |
| R1b | M343 | rs9786184 | snv | Y:3019783 | Y:2887824 | A/C/G | C>A | consistent, multi | derived |
| R-M269 | M269 | rs9786153 | snv | Y:20577481 | Y:22739367 | C/G/T | T>C | consistent, multi | derived |
| R-U106 | U106 | rs16981293 | snv | Y:8928037 | Y:8796078 | C/T | C>T | consistent | ancestral |
| R-P312 | P312 | rs34276300 | snv | Y:19995425 | Y:22157311 | A/C/G | A>G | consistent, multi | ancestral |

## The 31 markers that could not be audited

M31, M91, M60, M145, M174, M2, M35, P143, M130, M217, F1329, P15, F929, M69,
L15, M429, M438, P37.2, M267, M172, M410, L298, M184, M214, M231, M178, M175,
M122, P331, M3, M420.

Three of the four indel leads sit in this set: **M91, M60 and M175**. So does
**M267**, the marker the retracted audit named as reversed. None of them can be
checked without an rsID, and none of those claims is supported by anything
gathered so far except the Karafet article text on M91.

## What is needed next, in priority order

1. **A name to rsID source for the 31.** Karafet 2008 **Supplemental Table 1**
   remains the target artifact. The article body supplied so far does not contain
   it, and neither ISOGG page supplied so far contains rsIDs: the 2018 index page
   carries 2 table rows and 0 rsIDs, and the 2014 Source page carries 4,129 rows
   with columns `SNP | Haplogrp | Other Names | Sources` and no rsID, position,
   allele or mutation-type column.
2. **A primary source for M20**, to decide whether the rsID or the allele pair is
   the error.
3. **A representation decision for indels.** `Y_BACKBONE` has no way to express
   "9T to 8T" or "GGGG to GGG". Until it does, M17 cannot be recorded correctly
   even though its true state is now known, and the same will apply to M91.

## Reproducing this audit

```
python tools/audit_y_dbsnp.py
python tools/audit_y_dbsnp.py --json
```

One batched E-utilities request. It needs the network and nothing else. It reads
`Y_BACKBONE` and writes nothing back to it.

`tests/test_y_dbsnp_audit.py` covers the classification logic offline against
real dbSNP records captured on the run date, including the multi-allelic SPDI
parsing whose absence produced the discarded draft described above.
