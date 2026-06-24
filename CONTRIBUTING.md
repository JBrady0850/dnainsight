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
