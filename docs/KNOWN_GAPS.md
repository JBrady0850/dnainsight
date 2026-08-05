---
created: 2026-08-04
modified: 2026-08-04
tags: [archivist, dnainsight, reference, decision]
aliases: [DNAInsight Known Gaps, KNOWN_GAPS]
---

# DNAInsight v3.0 Known Gaps

## Why this file exists

A figure nobody re-verifies drifts. It gets copied into a second file, then into
a report, then into somebody's decision, and by then nobody remembers that the
first copy was written from memory. The only defence is to write down which
figures were checked at source and which were not, and to keep that list where
the people using the software can read it.

**Everything listed below is shipped and working.** None of it is broken, none
of it is a stub, and none of it will crash. Every item rests on a figure that was
**not machine-checked at source**, which is a different and quieter problem. A
tool that hides its unverified figures is worse than one that lists them, because
the first kind lets a reader assume that everything unmentioned was verified.

This project already refuses to conflate "we looked and found nothing" with "we
could not look at all". This document applies the same rule one level up, to the
software itself: here is what we checked, and here is what we did not.

Every item names the file and the field, so closing a gap is a bounded task
rather than an investigation.

## Summary

| Area | Unverified | Where |
|---|---|---|
| **Open defect** | rs28371706 filed under two different genes | `data/evidence_overlay.py:79`, `backend/diplotype.py` |
| CPIC base conflicts | 3 direct disagreements, both sides plus strand | `backend/diplotype.py` |
| Y backbone | **49 of 49 markers**, 31 with no rsID at all | `backend/haplogroups.py` |
| mtDNA backbone | 33 of 50 defining positions, 15 of 28 nodes fully | `backend/haplogroups.py` |
| Star alleles | 17 alleles across 7 genes | `backend/diplotype.py` |
| Carrier variant mappings | 22 of 23 | `backend/carrier.py` |
| Carrier frequencies | **25 of 25** | `backend/carrier.py` |
| Detection rates | 6 of 8 | `backend/carrier.py` |
| Total pathogenic counts | 11 of 11, order of magnitude only | `backend/carrier.py` |
| ACMG SF gene list | 82 encoded against 81 published | `backend/carrier.py` |
| Genetic map | approximate, 3,639 cM against about 3,550 | `backend/relatedness.py` |
| Relationship bands | DNAInsight's own approximation | `backend/relatedness.py` |
| External tool arguments | **every tool, none executed against a binary** | 6 modules |
| SGDP filename template | assumed, host certificate did not validate | `data/build_panel.py` |

---

## 1. Open defect: rs28371706 is filed under two different genes

**Status: UNRESOLVED. Highest priority in this release.**

`data/evidence_overlay.py` line 79 files **rs28371706 under CYP2C9**, with
medicines warfarin and phenytoin:

```python
"rs28371706":  ("T", "A", 3, 90,  ["CYP2C9"], ["warfarin", "phenytoin"]),
```

**rs28371706 is widely reported as the CYP2D6\*17 defining variant, c.1023C>T.**
`backend/diplotype.py` uses it that way.

**Both attributions cannot be correct.** One of the two files is wrong, and this
build did not establish which. Neither file was edited, because guessing which
one to change would replace a visible conflict with an invisible one.

Consequences, so nobody has to work them out under time pressure:

- If the overlay is wrong, a CYP2C9 finding is being surfaced for a variant that
  has nothing to do with warfarin dosing, in the pre-prescription silo.
- If `diplotype.py` is wrong, every CYP2D6\*17 call is wrong, and \*17 is a
  reduced-function allele that is common in African ancestry, so the error would
  fall unevenly.

**Resolve at source**, against dbSNP and the CPIC allele definition tables, then
fix whichever file is wrong and record the source in the commit.

### The three CPIC base conflicts

The CPIC allele definition builder was run on **2026-08-04**
(`python data/build_pgx_alleles.py --accept-terms`) and reconciled
`backend/diplotype.py` against CPIC's own tables. Result: **12 confirmed, 3
conflicting, 2 not present in CPIC's current table.**

**Both sides state the positive chromosomal strand.** These are therefore
genuine disagreements about the base, not a strand convention difference, and
they cannot be explained away.

| Gene | Allele | rsID | DNAInsight | CPIC |
|---|---|---|---|---|
| CYP2D6 | \*2 | rs1135840 | C | **G** |
| CYP2D6 | \*17 | rs28371706 | T | **A** |
| TPMT | \*2 | rs1800462 | C | **G** |

Not present in CPIC's current table at all: **SLCO1B1 \*1B** and **SLCO1B1
\*17**. That is not necessarily an error in either place. CPIC's SLCO1B1
nomenclature has moved, and an allele can be retired rather than wrong. It does
mean those two calls have no current upstream definition to check against.

**Nothing was edited.** The builder reports the conflict, `diplotype.py` still
carries its original bases, and both are recorded here. Resolve at source before
relying on these calls.

---

## 2. Haplogroup backbone

`backend/haplogroups.py` ships a bundled backbone so that a haplogroup call is
possible with no external tool installed. The tree stamps itself
**`DNAInsight backbone 0.1`**, and the version number is honest.

`unverified_markers()` returns the full audit list and it is attached to every
`/haplogroups` response. `verified_only=true` refuses every unverified marker
and yields an **unresolved call**, which is the honest state of this release.

### 2.1 Y chromosome: 49 of 49 markers unverified

**Every single Y marker is `verified: false`.** Every rsID and every
ancestral/derived allele pair is a **literature-recall candidate**, not
machine-checked against dbSNP, ISOGG or YFull.

**31 of the 49 have no rsID at all** and fall back to the marker name, which
means they read as **not testable on an rsID-keyed array**. That is not a bug in
the caller; it is a hole in the data. A marker with no rsID cannot be looked up
in a 23andMe or AncestryDNA export, so those 31 nodes contribute nothing to a
real call today.

Candidate rsIDs are recorded, and are candidates only, for these 18 nodes:

```
CT  E  F  G  I  I1  J  K  L  P  Q  R  R1  R1a1a  R1b  R-M269  R-U106  R-P312
```

One of them carries an extra warning in the data: **CT / M168 / rs2032595** was
supplied as an example pairing in the Wave 3 brief and is still unconfirmed.

Verification target: ISOGG or YFull, marker by marker, checking both the rsID
and the derived and ancestral alleles on the plus strand. Flip `verified` per row
and record the source.

### 2.2 mtDNA: 17 of 50 defining positions verified

**Verified, 17 positions across 13 nodes:**

| Node | Positions |
|---|---|
| M | 10400 |
| C | 13263 |
| D | 5178 |
| A | 663 |
| I | 10034 |
| W | 8994 |
| X | 13966 |
| U | 11467, 12308, 12372 |
| K | 9055 |
| J | 13708, 16069 |
| T | 4917 |
| H | 2706, 7028 |
| V | 4580 |

**Not verified: 33 positions across 22 nodes.** Fifteen nodes carry **no
verified defining position at all**, so those branches rest entirely on recall:

```
L0   L1   L1'2'3'4'5'6   L2   L2'3'4'5'6   L3   L4   L5   L6
N    R    R0    HV    JT    B
```

**The entire L lineage is unverified.** That falls unevenly: L0 through L6 are
the African macro-haplogroups, so the branch of the tree with no verified
positions is the branch that matters most for African maternal ancestry.

**B deserves a specific warning.** B's real defining feature is the 9-base-pair
deletion at 8281 to 8289, which **no consumer array calls**. The bundled tree
uses a 16189 proxy instead, which is weak, and the node is flagged. A B call from
this build should be treated as a hypothesis.

Two other proxies are flagged in the data and worth knowing about: **HV / 14766**
is a back mutation to the rCRS state and only discriminates below R0, and **H**
is defined by carrying the rCRS base at 2706 and 7028, because the rCRS is itself
an H sequence.

Verification target: PhyloTree, node by node.

---

## 3. Star alleles: 17 unverified across 7 genes

`backend/diplotype.py`. `unverified_entries()` returns this list and it is on
every `/pgx/diplotypes` response.

| Gene | Alleles |
|---|---|
| CYP2C19 | \*6, \*9, \*35 |
| CYP2C9 | \*5, \*8, \*11 |
| CYP2D6 | \*2, \*17, \*41 |
| DPYD | HapB3 |
| NUDT15 | \*5 |
| SLCO1B1 | \*1B, \*15, \*17 |
| TPMT | \*2, \*3A, \*3B |

Reconciled by the CPIC builder on 2026-08-04: **12 confirmed, 3 conflicting**
(section 1), **2 absent from CPIC's current table** (SLCO1B1 \*1B and \*17).

Three further facts about this table that are not gaps but are frequently
misread, so they are recorded here rather than in a docstring nobody opens.

**UGT1A1\*28 is called through the rs887829 tag, not the (TA)7 repeat.** The
real \*28 allele is a dinucleotide repeat length polymorphism and an array cannot
count repeats. rs887829 is in strong linkage with it, but **linkage is
ancestry-dependent**, so the tag can be right about the population and wrong
about the person. The call is a proxy and is labelled one.

**Six genes are on the minus strand**: TPMT, NUDT15, DPYD, VKORC1, CYP2D6 and
HBB. The base recorded in the allele table therefore differs from the c. notation
in the literature. That is expected and correct, and it is also the single
easiest way to introduce a wrong "fix" while verifying. Check the strand before
changing a base.

**CYP2D6 is permanently provisional**, independent of any verification work.
Arrays cannot see copy number or hybrid alleles, so \*2xN duplications, \*5
whole-gene deletions and CYP2D6-CYP2D7 hybrids are invisible. Verifying the
bases will not fix that, and nothing short of a different assay will.

---

## 4. Carrier data

`backend/carrier.py`. `unverified_figures()` returns all of it and it is
attached to every `/carrier` response.

### 4.1 Variant mappings: 22 of 23 unverified

**The only corroborated mapping in the panel is G6PD rs1050828 to T.**

Every other rsID-to-variant mapping across CFTR, HEXA, SMN1, HBB, GJB2, PAH,
ATP7B, GALT, ACADM and BTD was recalled, not corroborated in-tree.

Two entries carry warnings beyond the general one:

- **CFTR F508del (rs113993960)** is a three-base deletion. Arrays report
  substitutions. Vendors that carry this position use proprietary i-numbered
  probes with inconsistent D/I tokens, so the module treats it as **not
  readable** rather than risk the exact failure mode `backend/traits.py`
  documents for rs8176719, where a failed probe read as a deletion produced a
  confident wrong answer. F508del is roughly 70 percent of CF alleles in European
  ancestry, so this is the most consequential position the panel cannot read.
- **CFTR R117H (rs78655421)** has a consequence that depends on the poly-T tract
  in cis, which an array cannot read. **Even a positive result here is not
  interpretable** without clinical testing.

### 4.2 Carrier frequencies: 25 of 25 unverified

Every figure below was recalled from the published literature and **none was
re-verified at source in this build.**

| Gene | Population | Frequency |
|---|---|---|
| CFTR | European | 1 in 25 |
| CFTR | Ashkenazi | 1 in 24 |
| CFTR | Hispanic | 1 in 58 |
| CFTR | African | 1 in 61 |
| CFTR | East Asian | 1 in 94 |
| HEXA | Ashkenazi | 1 in 27 |
| HEXA | European | 1 in 300 |
| HEXA | General | 1 in 300 |
| SMN1 | European | 1 in 47 |
| SMN1 | Ashkenazi | 1 in 67 |
| SMN1 | African | 1 in 72 |
| SMN1 | East Asian | 1 in 59 |
| SMN1 | Hispanic | 1 in 68 |
| HBB | African | 1 in 13 |
| HBB | Middle Eastern | 1 in 30 |
| GJB2 | European | 1 in 33 |
| GJB2 | Ashkenazi | 1 in 25 |
| GJB2 | East Asian | 1 in 50 |
| PAH | European | 1 in 50 |
| ATP7B | General | 1 in 90 |
| GALT | European | 1 in 107 |
| ACADM | European | 1 in 65 |
| BTD | General | 1 in 120 |
| G6PD | African | 0.11 |
| G6PD | Middle Eastern | 0.05 |

**The two G6PD figures are allele frequencies, not carrier frequencies**,
because G6PD is X-linked. The proportion of heterozygous females is higher than
the proportion of hemizygous males, and the two are not interchangeable. The
module records the distinction; a reader skimming the table might not.

### 4.3 Detection rates: 6 unverified, 2 verified

Unverified, all from the 23-variant ACMG/ACOG CFTR panel except the last:

| Gene | Population | Rate |
|---|---|---|
| CFTR | European | 0.88 |
| CFTR | Ashkenazi | 0.94 |
| CFTR | Hispanic | 0.72 |
| CFTR | African | 0.65 |
| CFTR | East Asian | 0.49 |
| HEXA | Ashkenazi | 0.94 |

**Two detection rates are `verified: true`, and the reason matters.** They are
facts about the assay rather than about a population, so they can be asserted
without a citation to a study:

- **SMN1 = 0.0.** SMN1 carrier status is a copy number state. A SNP array
  measures bases, not copy number, so SMN1 carrier status is undetectable from
  array data **by construction**.
- **GJB2 European = 0.0.** The dominant European allele c.35delG is a deletion
  the array cannot read, so the array tests none of the common alleles in that
  population.

**Detection rate is deliberately ABSENT for HBB, PAH, ATP7B, GALT, ACADM, BTD
and G6PD, and for every population not listed above.** In those cases
`residual_risk` returns `None` with a reason naming the missing input, rather
than borrowing a number from a neighbouring population. A guess that looks like
the honest feature is worse than nothing at all.

### 4.4 Total pathogenic variant counts: 11 of 11 order-of-magnitude only

CFTR 2100, HEXA 130, HBB 400, GJB2 100, PAH 1000, ATP7B 800, GALT 330, ACADM
100, BTD 250, G6PD 200, plus SMN1 recorded as not applicable.

These exist to make the sentence "not a carrier for the N variants tested" mean
something by comparison. They are **orders of magnitude, not counts**, and should
never be quoted as figures.

### 4.5 Every residual risk is a lower bound

`is_lower_bound` is `true` on every result and this is a permanent property of
the design, not a gap that verification closes. **Published detection rates
belong to clinical panels that read more positions than DNAInsight does.** The
true residual risk after a negative DNAInsight result is therefore higher than
the number shown, and the payload says so.

### 4.6 ACMG secondary findings list

`ACMG_SF_GENES` is **hand-encoded as v3.2 (2023)** and holds **82 genes against a
published count of 81**. It was **not reconciled item by item.**
`ACMG_SF_LIST_VERIFIED` is `False`, and `acmg_coverage_report` prints the
discrepancy in its own output rather than leaving it in a comment:

> The published ACMG SF v3.2 table lists 81 genes. This hand-encoded copy holds
> 82 and was not reconciled item by item, so a gene missing from this list is not
> evidence that ACMG does not list it.

The four recorded array probe entries, including the three Ashkenazi BRCA1 and
BRCA2 founder variants, are also unverified rsIDs. **Carrying probes for three
BRCA variants is not BRCA testing.** Thousands of pathogenic BRCA variants exist
and a negative array result excludes three of them.

---

## 5. Tool integration: assumptions never executed against a binary

**No external tool's command-line arguments or output format has been run
against an installed binary.** Everything below is read from documentation. Each
one is a plausible assumption and each one is a place where a first real run will
probably fail with an obvious error, which is the good case. The bad case is a
tool that accepts the arguments and produces something DNAInsight parses wrongly.

### Beagle

- Argument form `gt= ref= map= out= chrom= nthreads=`.
- Output at `<out>.vcf.gz`.
- `DR2=`, `AF=` and `IMP` present in the INFO field.

### fastmixture

- Flags `--bfile --out --seed --projection --reference`.
- Output as `ancestry*.Q`.
- **Without `q_columns.tsv` the columns are reported as `component_N` with
  `components_labelled: false`**, because a numbered component is not a
  population. This is correct behaviour rather than a gap, and it is recorded
  here so nobody "fixes" it by inventing labels.

### FLARE

**This is the weakest assumption in the release.** Per-haplotype ancestry is
assumed to arrive as **`AN1` and `AN2` FORMAT fields holding integer indices into
the reference panel population order**. If the field names differ, painting
produces nothing. If the field names match but the **indexing base or the
population order differs**, painting produces a plausible-looking picture with
the ancestries assigned to the wrong populations, which is far worse. Verify the
ordering explicitly against a known sample before trusting any painted
chromosome.

### Yleaf, HaploGrep 3, Clade Finder, IBIS

CLI flags are deliberately isolated in named constants so correcting one is a
one-line change in a known place:

| Tool | Constant | File |
|---|---|---|
| Yleaf | `YLEAF_ARGS` | `backend/haplogroups.py` |
| HaploGrep 3 | `HAPLOGREP_ARGS` | `backend/haplogroups.py` |
| Clade Finder | `CLADEFINDER_ARGS` | `backend/haplogroups.py` |
| IBIS | `IBIS_ARGS` | `backend/relatedness.py` |

**IBIS has a known input mismatch, not just an assumption.** The module writes
`.ped` and `.map` files, but `IBIS_ARGS` references bed/bim/fam. **A PLINK
conversion step, or a switch to IBIS text-input mode, is still needed.** This one
will fail on first run.

---

## 6. Genetic map and relationship bands

`backend/relatedness.py`.

### The genetic map is approximate

`AVERAGE_CM_PER_MB` holds **approximate whole-chromosome averages**, not a
published genetic map. They total **3,639 cM autosomal against a real figure near
3,550**, so the map is about 2.5 percent long overall, and per-chromosome error
will be larger than that because recombination is not uniform along a chromosome.

**Every cM this produces is flagged `cm_estimated: true`.** That flag is the
whole reason the approximation is acceptable to ship: a number that announces
itself as an estimate is usable, and a number that does not is not.

Verification target: a public-domain or permissively licensed genetic map, at
sufficient resolution to replace whole-chromosome averages with per-interval
rates.

### The relationship bands are DNAInsight's own

`RELATIONSHIP_BANDS` is **DNAInsight's own approximate table**, written to the
shape of the published literature and **deliberately not a copy of the Shared cM
Project dataset**, for licence reasons. It is 15 bands from identical twin to
fourth cousin.

Both the map and the bands need confirming against a public-domain source. Note
the constraint before starting: the obvious reference dataset is the one that
cannot be copied, so this is a sourcing problem before it is a data-entry
problem.

Two properties that are correct and should survive any revision: **bands overlap
heavily by design**, and `predict_relationship` returns **every** band whose
range contains the total. Narrowing that to one answer would be inventing
certainty that no shared total contains.

---

## 7. Builder endpoint notes

`data/build_panel.py` and `data/build_pgx_alleles.py`.

### The IGSR ID column is empty, verified

**The IGSR 20130502 release genotype VCFs publish `ID` as `.` for every record.**
Verified across **all 1,103,547 chr22 records**. An rsID-keyed feature gets
nothing from them: `--array-file` matches nothing, and `informative_markers.tsv`
cannot be written.

`build_panel.py` therefore defaults to **`--onekg-source beagle`**, the same
phase 3 v5a callset with IDs populated. `igsr` remains selectable and prints the
caveat. This one is verified rather than assumed, and is recorded here because a
future maintainer will otherwise "restore" the canonical source and silently
break the join.

### The SGDP filename template is assumed

**The SGDP per-sample VCF filename template is assumed, not verified**, because
the Harvard host's certificate chain **did not validate at review time**. Three
flags override it, so a working template can be supplied without touching code:

```
--sgdp-vcf-template   the per-sample URL template, must contain {sample}
--sgdp-dir            read per-sample VCFs from a local directory instead
--sgdp-metadata       read the public-tier metadata CSV from a local path or another URL
```

`--sgdp-metadata` exists specifically because the host was unreachable, and is
needed alongside `--sgdp-dir`.

### CPIC build and strand handling

Three properties of `build_pgx_alleles.py` that look like bugs and are not:

- **CPIC positions are GRCh38 while `diplotype.py` records GRCh37.** The builder
  therefore compares **bases only and never positions**. Comparing positions
  across builds would produce a wall of false conflicts and bury the three real
  ones.
- **CPIC uses IUPAC ambiguity codes in `variantallele`, and these are never
  treated as a match.** "R matches A" is how a reconciliation tool talks itself
  into agreement it does not have.
- Both sides state the **positive chromosomal strand**, which is what makes the
  three conflicts in section 1 real rather than a convention artefact.

---

## 8. Closing a gap

**Verification is the most valuable contribution to v3.0.** It is worth more than
a new feature, and it is bounded: every item above names a file and a field.

1. Check the figure at source. ISOGG or YFull for Y markers, PhyloTree for mtDNA,
   the CPIC allele definition tables for star alleles, dbSNP for rsID-to-variant
   mappings, a citable publication for a carrier frequency or detection rate.
2. Fix the value if it is wrong. **Record what changed and why.**
3. Flip the `verified` flag on that row, and only that row. The flags are
   per-row precisely so partial progress is possible and visible.
4. Put the source in the commit message and the pull request.
5. Update this file. An item that has been verified should leave this list, and
   the summary table at the top should stop counting it.

If you verify a figure and it turns out to be **right**, that is a full
contribution and it should still be committed. "Checked on this date against this
source, unchanged" is exactly the record whose absence created this document.
