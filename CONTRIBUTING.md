# Contributing to DNAInsight

Thank you for your interest in contributing. DNAInsight is a privacy-first personal genomics tool and we welcome evidence-based improvements.

---

## Ways to Contribute

- Add or improve SNP entries in the bundled reference
- Fix bugs in parsing, scanning, or report generation
- Improve the frontend UI or report templates
- Add support for new DNA file formats
- Write or improve tests
- Update documentation

---

## Adding SNPs to the Reference

The bundled reference lives in `data/build_reference.py`. Each entry follows this format:

```python
("rsID", "GENE", "CATEGORY", "clinical_sig", "Plain-English interpretation text."),
```

**Categories:** `PHARM`, `METAB`, `INFLAM`, `NEURO`, `DETOX`, `CARDIO`

**clinical_sig values:** `drug response`, `risk factor`, `informational`, `pathogenic`

### Inclusion Criteria

All submitted SNPs must meet ALL of the following:

1. **Array coverage** — Confirmed present on at least one major consumer array (23andMe v4/v5, AncestryDNA v2/v3). Check the UCSC Genome Browser or SNPedia.
2. **Evidence level** — CPIC Level A or B, PharmGKB Level 1A/1B/2A, high ClinVar significance, or a replicated GWAS association (p < 5e-8, multiple cohorts).
3. **Actionability** — The finding must have a clear lifestyle, supplement, monitoring, or physician-discussion implication. Variants of uncertain significance (VUS) are not eligible.
4. **Plain-English interpretation** — The interpretation text must be written for a non-clinician. Avoid unexplained acronyms. Include the specific action or discussion point.
5. **No duplication** — Search `build_reference.py` for the rsID before submitting.

### Submission Process

1. Fork the repository and create a branch: `git checkout -b snp/add-GENE-rsXXXXXX`
2. Edit `data/build_reference.py` — add your tuple in the correct category section.
3. Rebuild the reference: `python data/build_reference.py`
4. Test by running DNAInsight locally and scanning a sample DNA file.
5. Open a Pull Request with:
   - The rsID, gene, and category
   - A link to the supporting evidence (CPIC guideline URL, PharmGKB entry, ClinVar page, or GWAS catalog entry)
   - The plain-English interpretation text you used and why it was worded that way

### What Gets Rejected

- Variants with no consumer array coverage
- VUS or conflicting evidence
- SNPs already in the reference
- Interpretations written in clinical jargon without plain-English explanation
- Population-specific variants without noting the population

---

## Code Contributions

### Setup

```bash
git clone https://github.com/yourusername/dnainsight.git
cd dnainsight
pip install -r requirements.txt
python app.py
```

### Project Layout

```
backend/        Flask API, parsers, scanner, report generators
data/           SNP reference build script and output JSON
frontend/       Single-page HTML/JS application
tests/          Unit tests (pytest)
grok/           AI analysis prompt
```

### Code Standards

- Python 3.10+, no external dependencies beyond `requirements.txt`
- Follow existing code style (PEP 8, 4-space indent)
- Add or update tests in `tests/` for any logic change
- Run `python -m pytest tests/` before submitting
- Do not commit `dnainsight.db`, uploaded DNA files, or generated reports

### Pull Request Checklist

- [ ] Tests pass locally (`python -m pytest tests/`)
- [ ] No DNA data files committed
- [ ] `requirements.txt` updated if new dependencies added (with pinned versions)
- [ ] README updated if user-facing behavior changed

---

## Architecture Notes

These constraints exist because a fix once made the mistake explicit. Read them before touching the affected code.

### Pipeline stage order

`backend/pipeline.py` runs merge, then bundled and API annotation, then orientation, frequency, snpedia, scoring, then genosets, traits and polygenic scores. Two adjacencies are load-bearing, not stylistic:

- Strand resolution must precede frequency lookup, or minus-strand frequencies silently return `unavailable`.
- Frequency must precede scoring, or the rarity adjustment has nothing to read and every magnitude flattens toward the default of 1.0.

### Tier 2 reference builder parsing traps

`data/build_full_reference.py` produces the gitignored, locally-built `data/reference.db` from ClinVar, the GWAS Catalog and CPIC. Known traps, verified against the live files:

- The ClinVar rsID column is literally `RS# (dbSNP)`, with a space and a hash.
- The gene column differs between ClinVar tables: `variant_summary.txt.gz` calls it `GeneSymbol`, `gene_specific_summary.txt` calls it `#Symbol`. The builder resolves every column by name, never by position. `PhenotypeIDS` ends in a capital S.
- Use `PositionVCF`, not `Start`/`Stop`, which are right-shifted.
- Filter to `Assembly == "GRCh37"`. `na` and NCBI36 rows exist. A `-1` rsID means no dbSNP mapping.
- Review status: the docs describe "no classification for the individual variant" but the data emits `single`. Code against the data string; the full mapping is `REVIEW_STATUS_STARS` in `backend/scoring.py`.
- GWAS `P-VALUE` underflows to 0 for thousands of rows; use `PVALUE_MLOG`. The sign of beta is in the `95% CI (TEXT)` column, not the numeric one.
- The GWAS Catalog's documented download endpoint has intermittently returned 404 in favour of `gwas-catalog-associations_ontology-annotated-full.zip` on its FTP mirror. The builder tries the documented URL first, then falls through to the FTP mirror, so a restored primary endpoint is picked up automatically.
- CPIC is CC0 and safe to bundle, but never copy the `clinpgxlevel` or `pgxtesting` columns out of the CPIC dump. They arrive inside a CC0 file, which makes them look safe, but they are PharmGKB-sourced and carry a no-commercial-sale clause. They are excluded at the wire level by `CPIC_SELECT`, not filtered afterward.

### Design decisions that should not be quietly reversed

1. **No SNPedia-derived data or PharmGKB/ClinPGx bulk data ships in this repository.** See the Licensing section of `CHANGELOG.md` for why.
2. **No `localStorage`, `sessionStorage` or any browser storage API** in the frontend or the interactive report. UI state lives in a plain object.
3. **No personal genetic data in the repository.** Never commit `uploads/` content or any `.db` file.
4. **Conflicting pooled calls from multiple DNA files are both retained.** No voting, no winner.
5. **A missing SNP evaluates false in a genoset**, and the genoset is reported as not testable rather than absent.
6. **Traits and polygenic scores never get a repute.** A trait is not good or bad.
7. **Non-carriers are down-weighted and their repute cleared.** A ClinVar classification describes an allele, not a position.
8. **Palindromic A/T and C/G sites are flagged, not guessed**, and capped at magnitude 2.
9. **`0.0` frequency and `null` frequency are different facts** and must render differently.
10. **A null magnitude sorts and filters as 1**, not 0.

---

## Reporting Bugs

Open a GitHub Issue with:
- DNAInsight version (shown in the UI footer)
- DNA provider and array version (e.g., 23andMe v5)
- Steps to reproduce
- Expected vs. actual behavior
- Any error output from the terminal

---

## Legal

By contributing, you agree that your contributions will be licensed under the MIT License. Do not submit real DNA data in issues, PRs, or example files.
