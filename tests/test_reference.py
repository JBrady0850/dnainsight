"""Unit tests for the reference builders.

Tier 1 is data/build_reference.py, the curated in-repo table. Tier 2 is
data/build_full_reference.py, the locally built full-array database.

Every Tier 2 test here runs with no network and no built database. The schema,
the CLI, the column-name constants and the read-side helpers are all verifiable
offline, and they are the parts that break silently.
"""

import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from data.build_reference import REFERENCE, build_reference, find_duplicates

from backend import scoring
from data import build_full_reference as bfr


class TestReferenceIntegrity:
    def test_no_duplicate_rsids(self):
        """Every rsID in the curated table must be unique; duplicates silently
        drop annotations at build time."""
        dupes = find_duplicates()
        assert dupes == [], f"Duplicate rsIDs in REFERENCE: {sorted(set(dupes))}"

    def test_build_count_matches_table(self):
        """With no duplicates, the built dict size equals the table length."""
        ref = build_reference()
        assert len(ref) == len(REFERENCE)

    def test_all_rows_have_five_fields(self):
        for row in REFERENCE:
            assert len(row) == 5, f"Malformed row: {row!r}"

    def test_all_rsids_well_formed(self):
        for rsid, *_ in REFERENCE:
            assert rsid.startswith("rs"), f"Bad rsID: {rsid}"

    def test_category_values_valid(self):
        valid = {"PHARM", "METAB", "INFLAM", "NEURO", "DETOX", "CARDIO"}
        for rsid, gene, category, *_ in REFERENCE:
            assert category in valid, f"{rsid}: unknown category {category}"


# ---------------------------------------------------------------------------
# Tier 2: the star mapping is borrowed, never redefined
# ---------------------------------------------------------------------------

# Every review-status string the builder can legitimately meet, with the star
# count it must produce. This covers the published documentation vocabulary, the
# different strings the DATA actually emits, the somatic aggregate vocabulary and
# the submitted-record vocabulary.
DOCUMENTED_STATUSES: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no assertion provided": 0,
    "no classification for the single variant": 0,
    "no classification for the individual variant": 0,
    "no classifications from unflagged records": 0,
    "flagged submission": 0,
    "-": 0,
    "": 0,
}


class TestTier2StarMappingIsShared:
    def test_mapping_is_the_same_object_as_scoring(self):
        """The builder must import REVIEW_STATUS_STARS, not carry its own copy.

        Identity is the assertion that matters. An equal-but-separate dict would
        pass a value comparison and then drift the moment either copy is edited.
        """
        assert bfr.REVIEW_STATUS_STARS is scoring.REVIEW_STATUS_STARS

    def test_mapping_is_not_redefined_in_the_builder_source(self):
        """No assignment to REVIEW_STATUS_STARS may appear in the builder.

        Checked against the source text because a re-assignment placed after the
        import would be invisible to an identity check performed at any other
        moment.
        """
        source = (Path(bfr.__file__)).read_text(encoding="utf-8")
        for forbidden in ("REVIEW_STATUS_STARS = {",
                          "REVIEW_STATUS_STARS: dict",
                          "REVIEW_STATUS_STARS={"):
            assert forbidden not in source, f"builder redefines: {forbidden!r}"

    def test_every_documented_status_is_covered(self):
        for status in DOCUMENTED_STATUSES:
            assert status in bfr.REVIEW_STATUS_STARS, f"unmapped status: {status!r}"

    def test_every_documented_status_maps_to_its_star_count(self):
        for status, stars in DOCUMENTED_STATUSES.items():
            assert bfr.REVIEW_STATUS_STARS[status] == stars, status
            assert scoring.review_stars(status) == stars, status

    def test_docs_spelling_and_data_spelling_agree(self):
        """The docs say "individual variant", the data says "single variant".

        Both must resolve, and to the same thing, because the builder feeds it
        whatever the file contained.
        """
        docs = "no classification for the individual variant"
        data = "no classification for the single variant"
        assert scoring.review_stars(docs) == scoring.review_stars(data) == 0

    def test_unknown_status_degrades_to_zero_without_raising(self):
        """ClinVar adds status strings between releases."""
        assert scoring.review_stars("reviewed by a committee of one") == 0
        assert scoring.review_stars(None) == 0

    def test_star_range_is_zero_to_four(self):
        for status, stars in bfr.REVIEW_STATUS_STARS.items():
            assert 0 <= stars <= 4, f"{status!r} maps outside 0 to 4"


# ---------------------------------------------------------------------------
# Tier 2: ClinVar column-name constants
# ---------------------------------------------------------------------------

class TestClinVarColumnNames:
    """The exact literals. Every one of these has been mistyped at least once."""

    def test_rsid_column_has_a_space_and_a_hash(self):
        assert bfr.COL_RSID == "RS# (dbSNP)"
        assert "#" in bfr.COL_RSID
        assert " " in bfr.COL_RSID

    def test_gene_column_carries_the_leading_hash(self):
        assert bfr.COL_GENE == "#Symbol"
        assert bfr.COL_GENE != "#GeneSymbol"
        assert bfr.COL_GENE.startswith("#")

    def test_phenotype_ids_column_ends_in_capital_s(self):
        assert bfr.COL_PHENOTYPE_IDS == "PhenotypeIDS"
        assert bfr.COL_PHENOTYPE_IDS.endswith("S")
        assert bfr.COL_PHENOTYPE_IDS != "PhenotypeIDs"

    def test_remaining_column_literals(self):
        assert bfr.COL_POSITION_VCF == "PositionVCF"
        assert bfr.COL_ASSEMBLY == "Assembly"
        assert bfr.COL_CHROMOSOME == "Chromosome"
        assert bfr.COL_CLINICAL_SIG == "ClinicalSignificance"
        assert bfr.COL_REVIEW_STATUS == "ReviewStatus"
        assert bfr.COL_PHENOTYPE_LIST == "PhenotypeList"

    def test_coordinates_come_from_position_vcf_not_start_or_stop(self):
        """Start and Stop are shifted relative to PositionVCF."""
        assert bfr.CLINVAR_COLUMNS["position"] == (bfr.COL_POSITION_VCF,)
        for candidates in bfr.CLINVAR_COLUMNS.values():
            assert "Start" not in candidates
            assert "Stop" not in candidates

    def test_assembly_target_is_grch37(self):
        assert bfr.TARGET_ASSEMBLY == "GRCh37"

    def test_absent_dbsnp_mapping_sentinel(self):
        assert bfr.NO_DBSNP_MAPPING == "-1"

    def test_documented_literals_are_offered_first(self):
        """The documented spelling leads each candidate list."""
        assert bfr.CLINVAR_COLUMNS["rsid"][0] == bfr.COL_RSID
        assert bfr.CLINVAR_COLUMNS["gene"][0] == bfr.COL_GENE


class TestColumnResolution:
    """Positions are resolved by name, never hard-coded."""

    HEADER = [
        "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID",
        "ClinicalSignificance", "ClinSigSimple", "LastEvaluated", "RS# (dbSNP)",
        "nsv/esv (dbVar)", "RCVaccession", "PhenotypeIDS", "PhenotypeList",
        "Origin", "OriginSimple", "Assembly", "ChromosomeAccession",
        "Chromosome", "Start", "Stop", "ReferenceAllele", "AlternateAllele",
        "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Guidelines",
        "TestedInGTR", "OtherIDs", "SubmitterCategories", "VariationID",
        "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
    ]

    def test_resolves_every_field_against_a_real_header(self):
        cols = bfr.resolve_columns(self.HEADER, bfr.CLINVAR_COLUMNS,
                                   bfr.CLINVAR_MANDATORY)
        assert len(cols) == len(bfr.CLINVAR_COLUMNS)
        assert cols["rsid"] == self.HEADER.index("RS# (dbSNP)")
        assert cols["position"] == self.HEADER.index("PositionVCF")
        assert cols["phenotype_ids"] == self.HEADER.index("PhenotypeIDS")

    def test_gene_falls_back_to_the_spelling_the_live_file_uses(self):
        """variant_summary.txt.gz publishes GeneSymbol, not #Symbol."""
        cols = bfr.resolve_columns(self.HEADER, bfr.CLINVAR_COLUMNS,
                                   bfr.CLINVAR_MANDATORY)
        assert cols["gene"] == self.HEADER.index("GeneSymbol")

    def test_gene_resolves_when_the_documented_spelling_is_present(self):
        header = ["#Symbol", "Assembly", "RS# (dbSNP)", "ClinicalSignificance",
                  "ReviewStatus"]
        cols = bfr.resolve_columns(header, bfr.CLINVAR_COLUMNS,
                                   bfr.CLINVAR_MANDATORY)
        assert cols["gene"] == 0

    def test_positions_are_not_assumed(self):
        """A reordered header must still resolve correctly."""
        header = list(reversed(self.HEADER))
        cols = bfr.resolve_columns(header, bfr.CLINVAR_COLUMNS,
                                   bfr.CLINVAR_MANDATORY)
        assert cols["rsid"] == header.index("RS# (dbSNP)")
        assert cols["assembly"] == header.index("Assembly")

    def test_missing_mandatory_column_raises_and_names_the_field(self):
        with pytest.raises(KeyError) as excinfo:
            bfr.resolve_columns(["Assembly", "Chromosome"], bfr.CLINVAR_COLUMNS,
                                bfr.CLINVAR_MANDATORY)
        assert "rsid" in str(excinfo.value)

    def test_missing_optional_column_is_simply_absent(self):
        header = ["RS# (dbSNP)", "Assembly", "ClinicalSignificance", "ReviewStatus"]
        cols = bfr.resolve_columns(header, bfr.CLINVAR_COLUMNS,
                                   bfr.CLINVAR_MANDATORY)
        assert "gene" not in cols
        assert cols["rsid"] == 0

    def test_gwas_header_resolves_including_the_awkward_names(self):
        header = ["DATE ADDED TO CATALOG", "PUBMEDID", "DISEASE/TRAIT", "SNPS",
                  "STRONGEST SNP-RISK ALLELE", "P-VALUE", "PVALUE_MLOG",
                  "OR or BETA", "95% CI (TEXT)", "MAPPED_TRAIT",
                  "MAPPED_TRAIT_URI", "STUDY ACCESSION"]
        cols = bfr.resolve_columns(header, bfr.GWAS_COLUMNS, bfr.GWAS_MANDATORY)
        assert len(cols) == len(bfr.GWAS_COLUMNS)
        assert cols["ci_text"] == header.index("95% CI (TEXT)")
        assert cols["study"] == header.index("STUDY ACCESSION")

    def test_gwas_significance_uses_mlog_and_never_p_value(self):
        """P-VALUE underflows to 0 for thousands of rows."""
        assert bfr.GWAS_COLUMNS["mlog"] == ("PVALUE_MLOG",)
        for candidates in bfr.GWAS_COLUMNS.values():
            assert "P-VALUE" not in candidates
        assert bfr.GWAS_MLOG_THRESHOLD == pytest.approx(7.3)
        assert bfr.GWAS_MIN_STUDIES == 2


# ---------------------------------------------------------------------------
# Tier 2: schema
# ---------------------------------------------------------------------------

class TestSchemaDDL:
    """The DDL is executed for real, in memory, so a typo cannot survive."""

    @staticmethod
    def _built() -> sqlite3.Connection:
        conn = sqlite3.connect(":memory:")
        bfr.create_schema(conn)
        return conn

    def test_ddl_executes_without_error(self):
        conn = self._built()
        conn.close()

    def test_variants_table_has_exactly_the_declared_columns(self):
        conn = self._built()
        try:
            info = conn.execute("PRAGMA table_info(variants)").fetchall()
        finally:
            conn.close()
        assert tuple(row[1] for row in info) == bfr.VARIANT_COLUMNS

    def test_variants_is_keyed_on_rsid(self):
        conn = self._built()
        try:
            info = conn.execute("PRAGMA table_info(variants)").fetchall()
        finally:
            conn.close()
        primary = [row[1] for row in info if row[5]]
        assert primary == ["rsid"]

    def test_required_columns_are_present(self):
        conn = self._built()
        try:
            info = conn.execute("PRAGMA table_info(variants)").fetchall()
        finally:
            conn.close()
        names = {row[1] for row in info}
        for column in ("rsid", "gene", "chromosome", "position", "clinical_sig",
                       "clinvar_sig_code", "review_status", "review_stars",
                       "condition", "cpic_level", "gwas_traits", "gwas_studies",
                       "publications", "risk_allele", "source"):
            assert column in names, f"missing column: {column}"

    def test_indices_on_gene_and_review_stars_exist(self):
        conn = self._built()
        try:
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='index' "
                "AND tbl_name='variants'").fetchall()
        finally:
            conn.close()
        names = {row[0] for row in rows}
        assert "idx_variants_gene" in names
        assert "idx_variants_review_stars" in names

    def test_indices_cover_the_intended_columns(self):
        conn = self._built()
        try:
            gene = conn.execute("PRAGMA index_info(idx_variants_gene)").fetchall()
            stars = conn.execute(
                "PRAGMA index_info(idx_variants_review_stars)").fetchall()
        finally:
            conn.close()
        assert [row[2] for row in gene] == ["gene"]
        assert [row[2] for row in stars] == ["review_stars"]

    def test_meta_table_exists_and_is_key_value(self):
        conn = self._built()
        try:
            info = conn.execute("PRAGMA table_info(meta)").fetchall()
        finally:
            conn.close()
        assert [row[1] for row in info] == ["key", "value"]

    def test_create_schema_is_idempotent(self):
        """A --resume build runs the DDL against an existing database."""
        conn = sqlite3.connect(":memory:")
        try:
            bfr.create_schema(conn)
            bfr.create_schema(conn)
            assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 0
        finally:
            conn.close()

    def test_clinvar_upsert_statement_is_valid_sql(self):
        """Built by string assembly, so it is compiled against the real schema."""
        conn = self._built()
        try:
            conn.execute("EXPLAIN " + bfr._clinvar_upsert_sql(),
                         tuple([None] * len(bfr.VARIANT_COLUMNS)))
        finally:
            conn.close()

    def test_upsert_keeps_the_better_attested_row(self):
        """The same rsID arrives repeatedly; PAR variants arrive four times."""
        conn = self._built()
        sql = bfr._clinvar_upsert_sql()
        weak = ("rs1", "AAA", "1", 100, "uncertain significance", 1,
                "criteria provided, single submitter", 1, "cond", None, None,
                0, 0, None, "clinvar")
        strong = ("rs1", "BBB", "1", 100, "pathogenic", 5,
                  "reviewed by expert panel", 3, "cond", None, None,
                  0, 0, None, "clinvar")
        try:
            conn.execute(sql, weak)
            conn.execute(sql, strong)
            conn.execute(sql, weak)
            rows = conn.execute(
                "SELECT rsid, gene, review_stars FROM variants").fetchall()
        finally:
            conn.close()
        assert rows == [("rs1", "BBB", 3)]

    def test_other_significance_never_outranks_pathogenic(self):
        """Code 255 means "other" and is not a higher grade than code 5."""
        conn = self._built()
        sql = bfr._clinvar_upsert_sql()
        pathogenic = ("rs2", "GENE1", "2", 200, "pathogenic", 5,
                      "no assertion criteria provided", 0, "c", None, None,
                      0, 0, None, "clinvar")
        other = ("rs2", "GENE2", "2", 200, "other", 255,
                 "no assertion criteria provided", 0, "c", None, None,
                 0, 0, None, "clinvar")
        try:
            conn.execute(sql, pathogenic)
            conn.execute(sql, other)
            gene = conn.execute("SELECT gene FROM variants").fetchone()[0]
        finally:
            conn.close()
        assert gene == "GENE1"


# ---------------------------------------------------------------------------
# Tier 2: CLI
# ---------------------------------------------------------------------------

class TestCommandLine:
    def test_parser_builds(self):
        assert bfr.build_parser() is not None

    def test_every_documented_flag_is_accepted(self):
        args = bfr.build_parser().parse_args([
            "--limit", "20000",
            "--skip-clinvar",
            "--skip-gwas",
            "--skip-cpic",
            "--array-file", "uploads/mine.txt",
            "--out", "data/other.db",
            "--stats",
            "--resume",
        ])
        assert args.limit == 20000
        assert args.skip_clinvar is True
        assert args.skip_gwas is True
        assert args.skip_cpic is True
        assert args.array_file == "uploads/mine.txt"
        assert args.out == "data/other.db"
        assert args.stats is True
        assert args.resume is True

    def test_defaults_are_a_full_build(self):
        args = bfr.build_parser().parse_args([])
        assert args.limit == 0
        assert args.skip_clinvar is False
        assert args.skip_gwas is False
        assert args.skip_cpic is False
        assert args.array_file is None
        assert args.out is None
        assert args.stats is False
        assert args.resume is False

    def test_each_flag_parses_on_its_own(self):
        for flag, attribute in (("--skip-clinvar", "skip_clinvar"),
                                ("--skip-gwas", "skip_gwas"),
                                ("--skip-cpic", "skip_cpic"),
                                ("--stats", "stats"),
                                ("--resume", "resume")):
            args = bfr.build_parser().parse_args([flag])
            assert getattr(args, attribute) is True, flag

    def test_limit_must_be_an_integer(self):
        with pytest.raises(SystemExit):
            bfr.build_parser().parse_args(["--limit", "many"])

    def test_unknown_flag_is_rejected(self):
        with pytest.raises(SystemExit):
            bfr.build_parser().parse_args(["--skip-snpedia"])


# ---------------------------------------------------------------------------
# Tier 2: read-side helpers degrade to None, never raise
# ---------------------------------------------------------------------------

class TestReadSideWithoutADatabase:
    def test_open_returns_none_when_the_file_is_absent(self, tmp_path):
        assert bfr.open_reference_db(tmp_path / "reference.db") is None

    def test_open_returns_none_for_the_default_path_when_absent(self, tmp_path,
                                                                monkeypatch):
        """A fresh clone has no data/reference.db and must still scan."""
        monkeypatch.setattr(bfr, "DEFAULT_DB_PATH", tmp_path / "reference.db")
        assert bfr.open_reference_db() is None

    def test_lookup_returns_none_rather_than_raising(self, tmp_path, monkeypatch):
        monkeypatch.setattr(bfr, "DEFAULT_DB_PATH", tmp_path / "reference.db")
        assert bfr.lookup("rs1801133") is None
        assert bfr.lookup("rs9999999999") is None
        assert bfr.lookup("") is None

    def test_coverage_stats_is_empty_rather_than_raising(self, tmp_path,
                                                         monkeypatch):
        monkeypatch.setattr(bfr, "DEFAULT_DB_PATH", tmp_path / "reference.db")
        assert bfr.coverage_stats() == {}

    def test_lookup_returns_none_for_a_missing_row_in_a_real_database(self,
                                                                     tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        bfr.create_schema(conn)
        conn.close()
        opened = bfr.open_reference_db(target)
        assert opened is not None
        try:
            assert bfr.lookup("rs404", opened) is None
        finally:
            opened.close()

    def test_lookup_finds_a_row_and_normalises_the_rsid(self, tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        bfr.create_schema(conn)
        conn.execute(
            "INSERT INTO variants (rsid, gene, review_stars) VALUES (?, ?, ?)",
            ("rs1801133", "MTHFR", 3))
        conn.commit()
        conn.close()
        opened = bfr.open_reference_db(target)
        try:
            row = bfr.lookup("  RS1801133  ", opened)
        finally:
            opened.close()
        assert row is not None
        assert row["gene"] == "MTHFR"
        assert row["review_stars"] == 3

    def test_open_reference_db_is_read_only(self, tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        bfr.create_schema(conn)
        conn.close()
        opened = bfr.open_reference_db(target)
        try:
            with pytest.raises(sqlite3.OperationalError):
                opened.execute("INSERT INTO variants (rsid) VALUES ('rs1')")
        finally:
            opened.close()


# ---------------------------------------------------------------------------
# Tier 2: free-text parsing traps
# ---------------------------------------------------------------------------

class TestGwasFreeTextParsing:
    def test_beta_direction_reads_the_ci_text_column(self):
        """The sign of an effect size is only ever stated as free text."""
        assert bfr._beta_direction("[0.09-0.15] unit decrease") == "decrease"
        assert bfr._beta_direction("[1.01-1.05] unit increase") == "increase"
        assert bfr._beta_direction("[0.2-0.4] cm lower") == "decrease"
        assert bfr._beta_direction("[0.2-0.4] cm higher") == "increase"

    def test_beta_direction_is_blank_when_unstated(self):
        assert bfr._beta_direction("") == ""
        assert bfr._beta_direction("[NR]") == ""
        assert bfr._beta_direction("[1.01-1.05]") == ""

    def test_risk_allele_is_extracted_from_the_strongest_snp_cell(self):
        assert bfr._risk_allele("rs1234567-A") == "A"
        assert bfr._risk_allele("rs76418789-G") == "G"

    def test_undetermined_risk_allele_is_blank_never_guessed(self):
        """A wrong risk allele inverts carrier status."""
        assert bfr._risk_allele("rs1234567-?") == ""
        assert bfr._risk_allele("rs1234567") == ""
        assert bfr._risk_allele("") == ""

    def test_snps_cell_splits_on_semicolons_and_interactions(self):
        assert bfr._split_gwas_rsids("rs123") == ["rs123"]
        assert bfr._split_gwas_rsids("rs1; rs2") == ["rs1", "rs2"]
        assert bfr._split_gwas_rsids("rs1 x rs2") == ["rs1", "rs2"]
        assert bfr._split_gwas_rsids("rs3 x rs4; rs5") == ["rs3", "rs4", "rs5"]

    def test_non_rsid_loci_are_dropped(self):
        """A chr:pos locus cannot be joined against an array export."""
        assert bfr._split_gwas_rsids("chr7:1234") == []
        assert bfr._split_gwas_rsids("HLA-DRB1*15:01") == []
        assert bfr._split_gwas_rsids("") == []


class TestCellHelpers:
    HEADER = ["RS# (dbSNP)", "Assembly", "PositionVCF"]

    def test_absent_column_yields_empty_string(self):
        cols = {"rsid": 0}
        assert bfr._cell(["rs1"], cols, "gene") == ""

    def test_short_row_yields_empty_string(self):
        cols = {"position": 9}
        assert bfr._cell(["rs1"], cols, "position") == ""

    def test_placeholder_values_are_treated_as_blank(self):
        """ClinVar writes na and - for absent values."""
        cols = {"gene": 0}
        for placeholder in ("na", "NA", "-", "NULL"):
            assert bfr._cell([placeholder], cols, "gene") == ""

    def test_numeric_coercion_rejects_junk(self):
        assert bfr._int_or_none("11856378") == 11856378
        assert bfr._int_or_none("") is None
        assert bfr._int_or_none("na") is None
        assert bfr._float_or_none("7.3") == pytest.approx(7.3)
        assert bfr._float_or_none("NR") is None


# ---------------------------------------------------------------------------
# Tier 2: array coverage set
# ---------------------------------------------------------------------------

class TestArrayCoverage:
    def test_reads_the_rsid_column_of_a_23andme_export(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_text(
            "# This data file generated by 23andMe\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs1801133\t1\t11856378\tAA\n"
            "rs4680\t22\t19951271\tAG\n",
            encoding="utf-8")
        assert bfr.load_array_rsids(str(raw)) == {"rs1801133", "rs4680"}

    def test_vendor_internal_identifiers_are_not_rsids(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_text("rs1\t1\t1\tAA\ni3000001\t1\t2\tGG\n", encoding="utf-8")
        assert bfr.load_array_rsids(str(raw)) == {"rs1"}

    def test_comma_separated_exports_are_accepted(self, tmp_path):
        raw = tmp_path / "genome.txt"
        raw.write_text('rsid,chromosome,position,allele1,allele2\n'
                       '"rs1801133",1,11856378,A,A\n', encoding="utf-8")
        assert bfr.load_array_rsids(str(raw)) == {"rs1801133"}

    def test_a_directory_is_unioned(self, tmp_path):
        (tmp_path / "a.txt").write_text("rs1\t1\t1\tAA\n", encoding="utf-8")
        (tmp_path / "b.txt").write_text("rs2\t1\t2\tGG\n", encoding="utf-8")
        assert bfr.load_array_rsids(str(tmp_path)) == {"rs1", "rs2"}

    def test_missing_path_yields_an_empty_set_not_an_error(self, tmp_path):
        assert bfr.load_array_rsids(str(tmp_path / "nope.txt")) == set()

    def test_empty_directory_yields_an_empty_set(self, tmp_path):
        """An empty coverage set means build without an array filter."""
        assert bfr.load_array_rsids(str(tmp_path)) == set()

    def test_no_vendor_manifest_url_is_referenced(self):
        """Array manifests are not redistributable and must not be fetched.

        Hosts and filenames only. The builder names Illumina and Affymetrix in
        prose to explain why it will not download from them, and that sentence
        is the point rather than a violation.
        """
        source = Path(bfr.__file__).read_text(encoding="utf-8").lower()
        for vendor in ("illumina.com", "thermofisher.com", "affymetrix.com",
                       "23andme.com/download", "manifest.csv", ".bpm", ".egt"):
            assert vendor not in source, f"builder fetches from {vendor}"


# ---------------------------------------------------------------------------
# Tier 2: licence hygiene and metadata
# ---------------------------------------------------------------------------

class TestLicenceHygiene:
    def test_cpic_request_asks_only_for_cc0_columns(self):
        assert bfr.CPIC_SELECT == "genesymbol,drugname,cpiclevel"

    def test_pharmgkb_sourced_cpic_columns_are_never_requested(self):
        """They arrive inside a CC0 dump but carry a no-commercial-sale clause."""
        for column in ("clinpgxlevel", "pgxtesting"):
            assert column not in bfr.CPIC_SELECT

    def test_no_pharmgkb_column_reaches_the_schema(self):
        for column in ("clinpgxlevel", "pgxtesting", "pharmgkb"):
            assert column not in bfr.VARIANT_COLUMNS

    def test_no_snpedia_endpoint_is_referenced(self):
        """SNPedia is CC-BY-NC-SA and would relicense the repository."""
        source = Path(bfr.__file__).read_text(encoding="utf-8").lower()
        assert "snpedia.com" not in source
        assert "bots.snpedia" not in source

    def test_no_pharmgkb_bulk_download_is_referenced(self):
        source = Path(bfr.__file__).read_text(encoding="utf-8").lower()
        for host in ("api.pharmgkb.org", "s3.pgkb.org", "clinpgx.org/download"):
            assert host not in source, f"builder references {host}"

    def test_every_source_has_a_recorded_licence(self):
        for source in ("clinvar", "gwas", "cpic", "array"):
            assert source in bfr.SOURCE_LICENCES
            assert bfr.SOURCE_LICENCES[source].strip()

    def test_cpic_licence_statement_names_cc0(self):
        assert "CC0-1.0" in bfr.SOURCE_LICENCES["cpic"]

    def test_clinvar_licence_statement_names_public_domain(self):
        statement = bfr.SOURCE_LICENCES["clinvar"].lower()
        assert "public domain" in statement

    def test_data_sources_document_exists_and_covers_every_source(self):
        """The changelog references this file, so it has to be there."""
        doc = Path(bfr.__file__).parent / "DATA_SOURCES.md"
        assert doc.exists(), "data/DATA_SOURCES.md is missing"
        text = doc.read_text(encoding="utf-8")
        for name in ("CPIC", "ClinVar", "gnomAD", "1000 Genomes",
                     "GWAS Catalog", "PGS Catalog", "SNPedia",
                     "MyVariant.info", "PharmGKB"):
            assert name in text, f"DATA_SOURCES.md does not cover {name}"
        for identifier in ("CC0-1.0", "CC-BY-NC-SA-3.0-US", "Apache-2.0"):
            assert identifier in text, f"DATA_SOURCES.md omits {identifier}"


class TestMetaAndStats:
    def test_meta_records_the_build_date_counts_and_licences(self, tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        try:
            bfr.create_schema(conn)
            bfr.write_meta(conn, {"clinvar": 12, "gwas": 34, "cpic": 56},
                           coverage_size=640000, skipped=["gwas (--skip-gwas)"])
            meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
        finally:
            conn.close()
        assert meta["clinvar_rows"] == "12"
        assert meta["gwas_rows"] == "34"
        assert meta["cpic_rows"] == "56"
        assert meta["array_coverage_rsids"] == "640000"
        assert meta["assembly"] == "GRCh37"
        assert meta["build_date"].endswith("Z")
        assert "clinvar_licence" in meta
        assert meta["skipped_sources"] == "gwas (--skip-gwas)"

    def test_build_date_is_timezone_aware_utc_not_utcnow(self):
        """utcnow returns a naive datetime and is deprecated."""
        source = Path(bfr.__file__).read_text(encoding="utf-8")
        assert "utcnow" not in source
        assert "datetime.now(timezone.utc)" in source

    def test_write_meta_is_idempotent(self, tmp_path):
        """A --resume build rewrites meta over an existing row set."""
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        try:
            bfr.create_schema(conn)
            bfr.write_meta(conn, {"clinvar": 1}, coverage_size=0)
            bfr.write_meta(conn, {"clinvar": 2}, coverage_size=0)
            rows = conn.execute(
                "SELECT value FROM meta WHERE key='clinvar_rows'").fetchall()
        finally:
            conn.close()
        assert rows == [("2",)]

    def test_coverage_stats_returns_the_row_counts(self, tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        try:
            bfr.create_schema(conn)
            bfr.write_meta(conn, {"clinvar": 7}, coverage_size=3)
        finally:
            conn.close()
        opened = bfr.open_reference_db(target)
        try:
            stats = bfr.coverage_stats(opened)
        finally:
            opened.close()
        assert stats["clinvar_rows"] == 7
        assert stats["variants_rows"] == 0
        assert stats["array_coverage_rsids"] == 3

    def test_array_filter_flag_records_whether_one_was_applied(self, tmp_path):
        target = tmp_path / "reference.db"
        conn = sqlite3.connect(str(target))
        try:
            bfr.create_schema(conn)
            bfr.write_meta(conn, {}, coverage_size=0)
            unfiltered = conn.execute(
                "SELECT value FROM meta WHERE key='array_filter_applied'"
            ).fetchone()[0]
            bfr.write_meta(conn, {}, coverage_size=500)
            filtered = conn.execute(
                "SELECT value FROM meta WHERE key='array_filter_applied'"
            ).fetchone()[0]
        finally:
            conn.close()
        assert unfiltered == "no"
        assert filtered == "yes"


class TestGitignoreCoversTheDatabase:
    def test_reference_db_and_its_sidecars_are_ignored(self):
        """GitHub rejects files over 100 MB, so this must never be committed."""
        gitignore = Path(bfr.__file__).parent.parent / ".gitignore"
        text = gitignore.read_text(encoding="utf-8")
        for pattern in ("data/reference.db", "data/reference.db-wal",
                        "data/reference.db-shm", "data/*.txt.gz",
                        "data/_cache/"):
            assert pattern in text, f".gitignore is missing {pattern}"
