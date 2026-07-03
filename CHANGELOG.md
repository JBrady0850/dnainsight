# Changelog

All notable changes to DNAInsight are documented here. This project follows
[Semantic Versioning](https://semver.org/).

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
