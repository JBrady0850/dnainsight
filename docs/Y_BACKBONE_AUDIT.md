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


---

# Re-run 2026-08-10, v3.4.0

Run after the Karafet 2008 supplement was folded in. **35 of 49 markers are now
reachable by dbSNP, up from 17.** No row of `Y_BACKBONE` was modified by this
run and every entry remains `verified: false`.

## What this run independently confirms

The audit is a report, so its value is what it corroborates from a different
direction than the publication did.

**M20 now returns `consistent`.** dbSNP gives rs3911 as A/G. The previous run
could not comment, because the table recorded a derived C that dbSNP has no
allele for. The supplement gives A->G, the table was corrected to A->G, and
dbSNP agrees. That conflict is closed by two independent sources rather than by
a judgement call.

**M267 now returns `consistent`.** dbSNP gives rs9341313 as T/G against the
corrected T->G. The pair had been stored C->A, the source pair complemented and
reversed.

**M145 and P143 return `strand`, which is the expected answer.** Both are the
same call read on the opposite strand, and the table records the strand rather
than rewriting the alleles.

**Five markers return `CLASS`, and four of those are the audit working.** M60
comes back `ins`, M175 `delins CTTCTCTTCTC/CTTCTC`, M17 `delins GGGG/GGG`.
dbSNP is reporting that these sites are length polymorphisms, which is exactly
what they are now recorded as. The verdict fires because the entry no longer
carries a base pair to compare, so the comparison is skipped rather than faked.

**M2 (rs3893) and M91 (rs2032651) come back empty**: no class, no positions, no
alleles. Both are 1994-2001 era accessions and the likeliest explanation is that
they have been merged into newer ones, which `esummary` does not follow. The
audit reports them as `CLASS` because it has nothing to compare, not because a
class conflict was found. Neither carries `dbsnp_checked`. **Following those
merges is the next concrete change to this tool**, and it is a tool change
rather than a data change.

## What still cannot be settled here

14 markers carry no rsID and are unreachable: M31, M35, P15, F1329, F929, P37.2,
M429, L15, M410, P331, L298, M178, M122, M420. Six post-date the 2008 paper and
are recorded as permanently unresolvable from it. Five were genotyped by that
survey and assigned no RefSNP ID. Two, M31 and M429, are held because the
supplement and the stored pair conflict.

The reference-orientation counts below remain the material result. **10 of 30
determinable markers have a reference carrying the DERIVED allele**, plus one
more on the opposite strand. A builder that mapped `ref` onto `ancestral` would
invert every one of those and nothing in the suite would catch it.

## Full output

```
node    marker rsid        class   GRCh38        GRCh37        ref/alt       anc>der   verdict                     ref_carries
------------------------------------------------------------------------------------------------------------------------------
B       M60    rs2032623   ins     Y:19716186    Y:21878073    /T            >        CLASS                       undetermined
BT      M91    rs2032651   None    None          None          /             >        CLASS                       undetermined
C       M130   rs35284970  snv     Y:2866813     Y:2734854     C/T           C>T       consistent                  ancestral
C-M217  M217   rs2032668   snv     Y:13325453    Y:15437333    A/C           A>C       consistent                  ancestral
CF      P143   rs4141886   snv     Y:12077161    Y:14197867    A/G           C>T       strand                      derived (opposite strand)
CT      M168   rs2032595   snv     Y:12702062    Y:14813991    T/C           C>T       consistent                  derived
D       M174   rs2032602   snv     Y:12842354    Y:14954280    T/A/C         T>C       consistent (multi-allelic)  ancestral
DE      M145   rs3848982   snv     Y:19555322    Y:21717208    C/T           G>A       strand                      ancestral (opposite strand)
E       M96    rs9306841   snv     Y:19617112    Y:21778998    C/G           G>C       consistent                  derived
E-M2    M2     rs3893      None    None          None          /             A>G       CLASS                       undetermined
F       M89    rs2032652   snv     Y:19755427    Y:21917313    T/C           C>T       consistent                  derived
G       M201   rs2032636   snv     Y:12915617    Y:15027529    G/T           G>T       consistent                  ancestral
H       M69    rs2032673   snv     Y:19732172    Y:21894058    T/C           T>C       consistent                  ancestral
I       M170   rs2032597   snv     Y:12735858    Y:14847792    A/C           A>C       consistent                  ancestral
I1      M253   rs17307677  snv     Y:16050954    Y:18162834    T/C           C>T       consistent                  derived
I2      M438   rs17307294  snv     Y:14526924    Y:16638804    A/G           A>G       consistent                  ancestral
J       M304   rs13447352  snv     Y:20587967    Y:22749853    A/C           A>C       consistent                  ancestral
J1      M267   rs9341313   snv     Y:20579932    Y:22741818    T/G           T>G       consistent                  ancestral
J2      M172   rs2032604   snv     Y:12857709    Y:14969634    T/G           T>G       consistent                  ancestral
K       M9     rs3900      snv     Y:19568371    Y:21730257    G/C           C>G       consistent                  derived
L       M20    rs3911      snv     Y:19571568    Y:21733454    A/G           A>G       consistent                  ancestral
N       M231   rs9341278   snv     Y:13357844    Y:15469724    G/A/T         G>A       consistent (multi-allelic)  ancestral
NO      M214   rs2032674   snv     Y:13360045    Y:15471925    T/C           T>C       consistent                  ancestral
O       M175   rs2032678   delins  Y:13396820    Y:15508700    CTTCTCTTCTC/CTTCTC>        CLASS                       undetermined
P       M45    rs2032631   snv     Y:19705901    Y:21867787    A/G/T         G>A       consistent (multi-allelic)  derived
Q       M242   rs8179021   snv     Y:12906671    Y:15018582    C/T           C>T       consistent                  ancestral
Q-M3    M3     rs3894      snv     Y:16984483    Y:19096363    G/A/T         C>T       strand                      ancestral (opposite strand)
R       M207   rs2032658   snv     Y:13470103    Y:15581983    G/A           A>G       consistent                  derived
R-M269  M269   rs9786153   snv     Y:20577481    Y:22739367    C/G/T         T>C       consistent (multi-allelic)  derived
R-P312  P312   rs34276300  snv     Y:19995425    Y:22157311    A/C/G         A>G       consistent (multi-allelic)  ancestral
R-U106  U106   rs16981293  snv     Y:8928037     Y:8796078     C/T           C>T       consistent                  ancestral
R1      M173   rs2032624   snv     Y:12914512    Y:15026424    C/A           A>C       consistent                  derived
R1a1a   M17    rs3908      delins  Y:19571279    Y:21733165    GGGG/GGG      >        CLASS                       undetermined
R1b     M343   rs9786184   snv     Y:3019783     Y:2887824     A/C/G         C>A       consistent (multi-allelic)  derived
T       M184   rs20320     snv     Y:12786229    Y:14898163    G/A           G>A       consistent                  ancestral

consistent                   21   M130, M217, M168, M96, M89, M201, M69, M170, M253, M438, M304, M267, M172, M9, M20, M214, M242, M207, U106, M173, M184
consistent (multi-allelic)    6   M174, M231, M45, M269, P312, M343
strand                        3   P143, M145, M3
CLASS                         5   M60, M91, M2, M175, M17

what the reference carries:
  ancestral                      17   M130, M217, M174, M201, M69, M170, M438, M304, M267, M172, M20, M231, M214, M242, P312, U106, M184
  ancestral (opposite strand)     2   M145, M3
  derived                        10   M168, M96, M89, M253, M9, M45, M207, M269, M173, M343
  derived (opposite strand)       1   P143
  undetermined                    5   M60, M91, M2, M175, M17

audited 35 of 49 markers.
14 carry no rsID and cannot be reached by dbSNP:
  M31, M35, P15, F1329, F929, P37.2, M429, L15, M410, P331, L298, M178, M122, M420

dbSNP cannot settle ancestral against derived. No value above may be
copied into an ancestral or derived field on the strength of this run.
```
