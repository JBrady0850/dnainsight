"""Tests for backend.sequencing: VCF and gVCF ingest, build detection, liftover.

Every test here runs with no external tools installed and no network. samtools is
never invoked; the alignment paths are exercised only through their degraded
payloads, which is exactly the state a normal user is in.
"""

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import external
from backend.parsers import ParseError
from backend.sequencing import (
    CONTIG_LENGTHS,
    MIN_PILEUP_DEPTH,
    REFERENCE_BUILD,
    SAMTOOLS_TOOL_ID,
    SAMTOOLS_CAPABILITY,
    SKIP_REASONS,
    BuildMismatch,
    ChainFile,
    alleles_from_gt,
    assert_build_compatible,
    call_genotype_from_counts,
    chain_dir,
    detect_build,
    detect_format,
    extract_positions,
    find_chain_file,
    is_simple_snv,
    iter_vcf_records,
    liftover,
    liftover_available,
    normalize_contig,
    parse_header_lines,
    parse_mpileup_line,
    parse_record,
    parse_sequencing_file,
    read_header,
    read_vcf,
    select_sample,
    split_gt,
    write_positions_bed,
)


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------

GRCH37_CONTIGS = [
    "##contig=<ID=1,length=249250621>",
    "##contig=<ID=2,length=243199373>",
    "##contig=<ID=3,length=198022430>",
    "##contig=<ID=X,length=155270560>",
]

GRCH38_CONTIGS = [
    "##contig=<ID=1,length=248956422>",
    "##contig=<ID=2,length=242193529>",
    "##contig=<ID=3,length=198295559>",
    "##contig=<ID=X,length=156040895>",
]

COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def build_vcf(records, contigs=None, samples=("NA00001",), meta=()):
    """Assemble a small VCF as a string."""
    lines = ["##fileformat=VCFv4.2"]
    lines.extend(meta)
    lines.extend(GRCH37_CONTIGS if contigs is None else contigs)
    lines.append("\t".join(list(COLUMNS) + list(samples)))
    lines.extend(records)
    return "\n".join(lines) + "\n"


def record(chrom="1", pos=100, rsid="rs1", ref="A", alt="G",
           qual="50", filt="PASS", info=".", fmt="GT", *calls):
    fields = [chrom, str(pos), rsid, ref, alt, qual, filt, info, fmt]
    fields.extend(calls or ("0/1",))
    return "\t".join(fields)


def write_vcf(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def write_gzipped_vcf(tmp_path, name, text):
    path = tmp_path / name
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)
    return path


@pytest.fixture
def isolated_home(tmp_path, monkeypatch):
    """Point ~/.dnainsight at a scratch directory so no real install is read."""
    home = tmp_path / "dnainsight_home"
    home.mkdir()
    monkeypatch.setenv("DNAINSIGHT_HOME", str(home))
    external.reset_cache()
    yield home
    external.reset_cache()


SIMPLE_CHAIN = """chain 4900 chr1 249250621 + 100 210 chr1 248956422 + 500 620 1
50\t10\t20
50

chain 3000 chr2 243199373 + 1000 1010 chr2 242193529 - 2000 2010 2
10

"""


# ---------------------------------------------------------------------------
# Contig naming
# ---------------------------------------------------------------------------

class TestNormalizeContig:
    def test_bare_number_is_unchanged(self):
        assert normalize_contig("1") == "1"

    def test_chr_prefix_is_stripped(self):
        assert normalize_contig("chr1") == "1"

    def test_chr_prefix_is_case_insensitive(self):
        assert normalize_contig("CHR17") == "17"

    def test_chrm_becomes_mt(self):
        assert normalize_contig("chrM") == "MT"

    def test_m_becomes_mt(self):
        assert normalize_contig("M") == "MT"

    def test_x_is_uppercased(self):
        assert normalize_contig("chrx") == "X"

    def test_numeric_23_is_left_alone(self):
        # Translating 23 to X would be a guess about the file's provenance.
        assert normalize_contig("23") == "23"

    def test_empty_name_is_empty(self):
        assert normalize_contig(None) == ""


# ---------------------------------------------------------------------------
# Header parsing
# ---------------------------------------------------------------------------

class TestHeaderParsing:
    def test_fileformat_is_read(self):
        header = parse_header_lines(["##fileformat=VCFv4.2", "\t".join(COLUMNS + ["S1"])])
        assert header.fileformat == "VCFv4.2"

    def test_reference_line_is_read(self):
        header = parse_header_lines([
            "##reference=file:///refs/human_g1k_v37.fasta",
            "\t".join(COLUMNS + ["S1"]),
        ])
        assert header.reference == "file:///refs/human_g1k_v37.fasta"

    def test_contig_id_and_length_are_read(self):
        header = parse_header_lines(GRCH37_CONTIGS + ["\t".join(COLUMNS + ["S1"])])
        assert header.contig_length("1") == 249250621

    def test_contig_assembly_tag_is_read(self):
        header = parse_header_lines([
            "##contig=<ID=1,length=249250621,assembly=b37>",
            "\t".join(COLUMNS + ["S1"]),
        ])
        assert header.contigs[0]["assembly"] == "b37"

    def test_quoted_commas_do_not_split_attributes(self):
        header = parse_header_lines([
            '##contig=<ID=1,length=249250621,description="one, two, three">',
            "\t".join(COLUMNS + ["S1"]),
        ])
        assert header.contigs[0]["length"] == 249250621

    def test_sample_names_are_read(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM", "DAD", "KID"])])
        assert header.samples == ["MUM", "DAD", "KID"]
        assert header.sample_count == 3

    def test_sites_only_header_has_no_samples(self):
        header = parse_header_lines(["\t".join(COLUMNS[:8])])
        assert header.samples == []

    def test_chr_prefix_is_recorded(self):
        header = parse_header_lines([
            "##contig=<ID=chr1,length=248956422>",
            "\t".join(COLUMNS + ["S1"]),
        ])
        assert header.chr_prefixed is True

    def test_non_ref_alt_line_sets_gvcf_hint(self):
        header = parse_header_lines([
            "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele\">",
            "\t".join(COLUMNS + ["S1"]),
        ])
        assert header.gvcf_hint is True

    def test_read_header_from_disk(self, tmp_path):
        path = write_vcf(tmp_path, "a.vcf", build_vcf([record()]))
        header = read_header(path)
        assert header.samples == ["NA00001"]

    def test_missing_chrom_line_raises(self, tmp_path):
        path = write_vcf(tmp_path, "b.vcf", "##fileformat=VCFv4.2\n")
        with pytest.raises(ParseError, match="#CHROM"):
            read_header(path)

    def test_data_before_chrom_line_raises(self, tmp_path):
        path = write_vcf(tmp_path, "c.vcf", "##fileformat=VCFv4.2\n1\t100\trs1\tA\tG\t.\t.\t.\n")
        with pytest.raises(ParseError, match="data line"):
            read_header(path)


# ---------------------------------------------------------------------------
# Build detection
# ---------------------------------------------------------------------------

class TestBuildDetection:
    def test_grch37_detected_from_chromosome_1_length(self):
        header = parse_header_lines(["##contig=<ID=1,length=249250621>"])
        result = detect_build(header)
        assert result["build"] == "GRCh37"
        assert result["confidence"] == "high"

    def test_grch38_detected_from_chromosome_1_length(self):
        header = parse_header_lines(["##contig=<ID=1,length=248956422>"])
        result = detect_build(header)
        assert result["build"] == "GRCh38"
        assert result["confidence"] == "high"

    def test_grch37_detected_from_chromosomes_2_and_3(self):
        header = parse_header_lines([
            "##contig=<ID=2,length=243199373>",
            "##contig=<ID=3,length=198022430>",
        ])
        assert detect_build(header)["build"] == "GRCh37"

    def test_grch38_detected_from_x_chromosome_length(self):
        header = parse_header_lines(["##contig=<ID=X,length=156040895>"])
        assert detect_build(header)["build"] == "GRCh38"

    def test_grch37_detected_through_chr_prefixed_contigs(self):
        header = parse_header_lines(["##contig=<ID=chr1,length=249250621>"])
        assert detect_build(header)["build"] == "GRCh37"

    def test_unknown_build_returns_none_when_header_is_bare(self):
        header = parse_header_lines(["##fileformat=VCFv4.2"])
        result = detect_build(header)
        assert result["build"] is None
        assert result["confidence"] == "none"

    def test_unrecognised_contig_length_returns_none(self):
        header = parse_header_lines(["##contig=<ID=1,length=123456789>"])
        result = detect_build(header)
        assert result["build"] is None
        assert result["confidence"] == "none"

    def test_contig_without_length_is_not_evidence(self):
        header = parse_header_lines(["##contig=<ID=1>"])
        assert detect_build(header)["build"] is None

    def test_mixed_build_lengths_report_conflict(self):
        header = parse_header_lines([
            "##contig=<ID=1,length=249250621>",
            "##contig=<ID=2,length=242193529>",
        ])
        result = detect_build(header)
        assert result["build"] is None
        assert result["confidence"] == "conflict"

    def test_mitochondrial_length_alone_decides_nothing(self):
        # 16569 bases of rCRS is identical in both builds.
        header = parse_header_lines(["##contig=<ID=MT,length=16569>"])
        result = detect_build(header)
        assert result["build"] is None
        assert CONTIG_LENGTHS["GRCh37"]["MT"] == CONTIG_LENGTHS["GRCh38"]["MT"]

    def test_reference_header_names_build_when_no_lengths(self):
        header = parse_header_lines(["##reference=file:///refs/human_g1k_v37.fasta"])
        result = detect_build(header)
        assert result["build"] == "GRCh37"
        assert result["confidence"] == "medium"

    def test_reference_header_names_grch38(self):
        header = parse_header_lines(["##reference=/data/GRCh38/genome.fa"])
        assert detect_build(header)["build"] == "GRCh38"

    def test_assembly_tag_names_build_when_no_lengths(self):
        header = parse_header_lines(["##contig=<ID=1,assembly=GRCh38>"])
        result = detect_build(header)
        assert result["build"] == "GRCh38"
        assert result["confidence"] == "medium"

    def test_reference_naming_both_builds_names_neither(self):
        header = parse_header_lines(["##reference=/refs/GRCh37/lifted_from_hg38.fa"])
        assert detect_build(header)["build"] is None

    def test_contig_lengths_outrank_contradictory_reference_text(self):
        header = parse_header_lines([
            "##reference=/refs/hg38.fa",
            "##contig=<ID=1,length=249250621>",
        ])
        result = detect_build(header)
        assert result["build"] == "GRCh37"
        assert result["confidence"] == "medium"

    def test_chr_prefix_never_votes_for_a_build(self):
        header = parse_header_lines(["##contig=<ID=chr1,length=249250621>"])
        naming = [e for e in detect_build(header)["evidence"] if e["signal"] == "contig_naming"]
        assert naming and all(e["supports"] is None for e in naming)

    def test_evidence_entries_carry_signal_detail_and_supports(self):
        header = parse_header_lines(["##contig=<ID=1,length=249250621>"])
        for entry in detect_build(header)["evidence"]:
            assert set(entry) == {"signal", "detail", "supports"}


# ---------------------------------------------------------------------------
# Build refusal
# ---------------------------------------------------------------------------

class TestAssertBuildCompatible:
    def test_matching_build_passes(self):
        assert assert_build_compatible("GRCh37", "GRCh37") is None

    def test_mismatched_build_raises(self):
        with pytest.raises(BuildMismatch):
            assert_build_compatible("GRCh38", "GRCh37")

    def test_build_mismatch_is_a_parse_error(self):
        # Callers keep one except clause.
        with pytest.raises(ParseError):
            assert_build_compatible("GRCh38", "GRCh37")

    def test_mismatch_message_names_both_builds(self):
        with pytest.raises(BuildMismatch) as exc:
            assert_build_compatible("GRCh38", "GRCh37")
        assert "GRCh38" in str(exc.value) and "GRCh37" in str(exc.value)

    def test_mismatch_carries_detected_and_expected(self):
        with pytest.raises(BuildMismatch) as exc:
            assert_build_compatible("GRCh38", "GRCh37")
        assert exc.value.detected == "GRCh38"
        assert exc.value.expected == "GRCh37"

    def test_unknown_build_raises_by_default(self):
        with pytest.raises(BuildMismatch):
            assert_build_compatible(None, "GRCh37")

    def test_unknown_build_is_allowed_when_flagged(self):
        assert assert_build_compatible(None, "GRCh37", allow_unknown=True) is None

    def test_wrong_build_still_raises_even_when_unknown_is_allowed(self):
        with pytest.raises(BuildMismatch):
            assert_build_compatible("GRCh38", "GRCh37", allow_unknown=True)

    def test_expected_none_skips_the_check(self):
        assert assert_build_compatible("GRCh38", None) is None

    def test_reference_build_constant_is_grch37(self):
        assert REFERENCE_BUILD == "GRCh37"


# ---------------------------------------------------------------------------
# Genotype decoding
# ---------------------------------------------------------------------------

class TestSplitGt:
    def test_unphased_separator(self):
        assert split_gt("0/1") == (["0", "1"], False)

    def test_phased_separator(self):
        assert split_gt("0|1") == (["0", "1"], True)

    def test_haploid_has_no_separator(self):
        assert split_gt("1") == (["1"], False)

    def test_empty_gt_yields_no_indexes(self):
        assert split_gt("") == ([], False)

    def test_mixed_separators_are_not_phased(self):
        indexes, phased = split_gt("0|1/2")
        assert indexes == ["0", "1", "2"]
        assert phased is False


class TestAllelesFromGt:
    def test_unphased_heterozygote(self):
        out = alleles_from_gt("0/1", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("A", "G")
        assert out["phased"] is False

    def test_phased_heterozygote_is_flagged(self):
        out = alleles_from_gt("0|1", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("A", "G")
        assert out["phased"] is True

    def test_phased_order_is_preserved(self):
        out = alleles_from_gt("1|0", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("G", "A")

    def test_homozygous_reference(self):
        out = alleles_from_gt("0/0", "C", ["T"])
        assert (out["allele1"], out["allele2"]) == ("C", "C")

    def test_homozygous_alternate(self):
        out = alleles_from_gt("1/1", "C", ["T"])
        assert (out["allele1"], out["allele2"]) == ("T", "T")

    def test_haploid_alt_duplicates_the_allele(self):
        out = alleles_from_gt("1", "A", ["T"])
        assert (out["allele1"], out["allele2"]) == ("T", "T")
        assert out["ploidy"] == 1

    def test_haploid_reference_duplicates_the_reference(self):
        out = alleles_from_gt("0", "A", ["T"])
        assert (out["allele1"], out["allele2"]) == ("A", "A")

    def test_haploid_missing_is_a_no_call(self):
        out = alleles_from_gt(".", "A", ["T"])
        assert (out["allele1"], out["allele2"]) == ("N", "N")
        assert out["no_call"] is True

    def test_missing_diploid_gt_becomes_n_over_n(self):
        out = alleles_from_gt("./.", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("N", "N")
        assert out["no_call"] is True

    def test_phased_missing_gt_becomes_n_over_n(self):
        out = alleles_from_gt(".|.", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("N", "N")

    def test_empty_gt_becomes_n_over_n(self):
        out = alleles_from_gt("", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("N", "N")
        assert out["ploidy"] == 0

    def test_partial_no_call_keeps_the_called_allele(self):
        out = alleles_from_gt("./1", "A", ["G"])
        assert (out["allele1"], out["allele2"]) == ("N", "G")
        assert out["partial_no_call"] is True
        assert out["no_call"] is False

    def test_multiallelic_selects_the_second_alt(self):
        out = alleles_from_gt("1/2", "A", ["G", "T"])
        assert (out["allele1"], out["allele2"]) == ("G", "T")

    def test_multiallelic_selects_reference_and_second_alt(self):
        out = alleles_from_gt("0/2", "A", ["G", "T"])
        assert (out["allele1"], out["allele2"]) == ("A", "T")

    def test_gt_index_beyond_alt_list_raises(self):
        with pytest.raises(ParseError, match="ALT allele"):
            alleles_from_gt("0/3", "A", ["G"])

    def test_non_numeric_gt_index_raises(self):
        with pytest.raises(ParseError, match="GT allele index"):
            alleles_from_gt("0/x", "A", ["G"])

    def test_alleles_are_uppercased(self):
        out = alleles_from_gt("0/1", "a", ["g"])
        assert (out["allele1"], out["allele2"]) == ("A", "G")


class TestIsSimpleSnv:
    def test_single_base_pair_is_simple(self):
        assert is_simple_snv("A", "G") is True

    def test_insertion_is_not_simple(self):
        assert is_simple_snv("A", "AGG") is False

    def test_deletion_is_not_simple(self):
        assert is_simple_snv("AGG", "A") is False

    def test_ambiguous_base_is_not_simple(self):
        assert is_simple_snv("N", "A") is False


# ---------------------------------------------------------------------------
# Record parsing and skip classification
# ---------------------------------------------------------------------------

class TestParseRecord:
    def test_snv_produces_the_parsers_dict_shape(self):
        parsed = parse_record(record(rsid="rs4988235", ref="G", alt="A"))
        assert parsed["skip"] is None
        assert set(parsed["call"]) == {"rsid", "chromosome", "position", "allele1", "allele2"}

    def test_snv_call_values(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["call"] == {
            "rsid": "rs1", "chromosome": "1", "position": 100,
            "allele1": "A", "allele2": "G",
        }

    def test_chr_prefix_is_normalised_on_records(self):
        parsed = parse_record("chr7\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["call"]["chromosome"] == "7"

    def test_insertion_is_skipped_as_indel(self):
        parsed = parse_record("1\t100\trs1\tA\tAGG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "indel"

    def test_deletion_is_skipped_as_indel(self):
        parsed = parse_record("1\t100\trs1\tAGG\tA\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "indel"

    def test_indel_without_rsid_is_still_counted_as_indel(self):
        parsed = parse_record("1\t100\t.\tA\tAGG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "indel"

    def test_symbolic_deletion_is_skipped_as_symbolic(self):
        parsed = parse_record("1\t100\trs1\tA\t<DEL>\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "symbolic"

    def test_breakend_is_skipped_as_symbolic(self):
        parsed = parse_record("1\t100\trs1\tA\tA[2:200[\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "symbolic"

    def test_spanning_deletion_star_is_skipped_as_symbolic(self):
        parsed = parse_record("1\t100\trs1\tA\tG,*\t50\tPASS\t.\tGT\t1/2")
        assert parsed["skip"] == "symbolic"

    def test_gvcf_reference_block_is_flagged_and_skipped(self):
        line = "1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=2000\tGT:DP\t0/0:30"
        parsed = parse_record(line)
        assert parsed["is_ref_block"] is True
        assert parsed["skip"] == "symbolic"

    def test_gvcf_variant_record_with_non_ref_padding_still_parses(self):
        line = "1\t100\trs1\tA\tG,<NON_REF>\t50\tPASS\t.\tGT\t0/1"
        parsed = parse_record(line)
        assert parsed["skip"] is None
        assert (parsed["allele1"], parsed["allele2"]) == ("A", "G")

    def test_gt_landing_on_non_ref_is_symbolic(self):
        line = "1\t100\trs1\tA\tG,<NON_REF>\t50\tPASS\t.\tGT\t0/2"
        assert parse_record(line)["skip"] == "symbolic"

    def test_missing_id_is_skipped_as_no_id(self):
        parsed = parse_record("1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "no_id"
        assert parsed["rsid"] is None

    def test_non_rs_identifier_is_treated_as_no_id(self):
        parsed = parse_record("1\t100\tCOSM123\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] == "no_id"

    def test_multiple_ids_pick_the_rs_one(self):
        parsed = parse_record("1\t100\tCOSM123;rs99\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["call"]["rsid"] == "rs99"

    def test_uppercase_rs_prefix_is_normalised(self):
        parsed = parse_record("1\t100\tRS99\tA\tG\t50\tPASS\t.\tGT\t0/1")
        assert parsed["call"]["rsid"] == "rs99"

    def test_no_call_is_skipped_as_no_call(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t./.")
        assert parsed["skip"] == "no_call"
        assert (parsed["allele1"], parsed["allele2"]) == ("N", "N")

    def test_partial_no_call_is_skipped_as_no_call(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t./1")
        assert parsed["skip"] == "no_call"

    def test_multiallelic_snv_selection_is_kept(self):
        parsed = parse_record("1\t100\trs1\tA\tG,T\t50\tPASS\t.\tGT\t1/2")
        assert parsed["skip"] is None
        assert (parsed["allele1"], parsed["allele2"]) == ("G", "T")

    def test_multiallelic_with_indel_selected_is_multiallelic_complex(self):
        parsed = parse_record("1\t100\trs1\tA\tG,AT\t50\tPASS\t.\tGT\t0/2")
        assert parsed["skip"] == "multiallelic_complex"

    def test_multiallelic_with_indel_not_selected_is_kept(self):
        parsed = parse_record("1\t100\trs1\tA\tG,AT\t50\tPASS\t.\tGT\t0/1")
        assert parsed["skip"] is None

    def test_polyploid_gt_is_multiallelic_complex(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1/1")
        assert parsed["skip"] == "multiallelic_complex"

    def test_haploid_chry_record_is_kept(self):
        parsed = parse_record("Y\t2655180\trs11575897\tG\tA\t50\tPASS\t.\tGT\t1")
        assert parsed["call"]["allele1"] == "A"
        assert parsed["call"]["allele2"] == "A"

    def test_gt_is_located_by_format_key_not_position(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tDP:GT\t30:0/1")
        assert (parsed["allele1"], parsed["allele2"]) == ("A", "G")

    def test_format_without_gt_is_a_no_call(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.\tDP\t30")
        assert parsed["skip"] == "no_call"

    def test_second_sample_is_read_when_index_given(self):
        line = "1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1"
        parsed = parse_record(line, sample_index=1)
        assert (parsed["allele1"], parsed["allele2"]) == ("G", "G")

    def test_sample_index_beyond_columns_raises(self):
        line = "1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0"
        with pytest.raises(ParseError, match="sample index"):
            parse_record(line, sample_index=4)

    def test_short_line_raises(self):
        with pytest.raises(ParseError, match="8 tab-separated"):
            parse_record("1\t100\trs1\tA")

    def test_unreadable_position_raises(self):
        with pytest.raises(ParseError, match="POS"):
            parse_record("1\tNOTANUMBER\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1")

    def test_sites_only_record_raises(self):
        with pytest.raises(ParseError, match="no FORMAT or sample"):
            parse_record("1\t100\trs1\tA\tG\t50\tPASS\t.")

    def test_filter_value_is_carried(self):
        parsed = parse_record("1\t100\trs1\tA\tG\t50\tLowQual\t.\tGT\t0/1")
        assert parsed["filter"] == "LowQual"


# ---------------------------------------------------------------------------
# Sample selection
# ---------------------------------------------------------------------------

class TestSelectSample:
    def test_default_is_the_first_sample(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM", "DAD", "KID"])])
        assert select_sample(header, None) == (0, "MUM")

    def test_selection_by_name(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM", "DAD", "KID"])])
        assert select_sample(header, "KID") == (2, "KID")

    def test_selection_by_index(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM", "DAD", "KID"])])
        assert select_sample(header, 1) == (1, "DAD")

    def test_unknown_name_raises_and_lists_the_samples(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM", "DAD"])])
        with pytest.raises(ParseError, match="MUM, DAD"):
            select_sample(header, "GRAN")

    def test_out_of_range_index_raises(self):
        header = parse_header_lines(["\t".join(COLUMNS + ["MUM"])])
        with pytest.raises(ParseError, match="out of range"):
            select_sample(header, 7)

    def test_sites_only_header_raises(self):
        header = parse_header_lines(["\t".join(COLUMNS[:8])])
        with pytest.raises(ParseError, match="no sample columns"):
            select_sample(header, None)


# ---------------------------------------------------------------------------
# Streaming readers
# ---------------------------------------------------------------------------

class TestReadVcf:
    def test_snps_are_collected(self, tmp_path):
        text = build_vcf([
            record(rsid="rs1", pos=100),
            record(rsid="rs2", pos=200, ref="C", alt="T"),
        ])
        result = read_vcf(write_vcf(tmp_path, "a.vcf", text))
        assert [s["rsid"] for s in result["snps"]] == ["rs1", "rs2"]

    def test_skipped_dict_always_carries_every_reason(self, tmp_path):
        result = read_vcf(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert set(result["skipped"]) == set(SKIP_REASONS)

    def test_indel_is_counted_in_skipped(self, tmp_path):
        text = build_vcf([
            record(rsid="rs1"),
            "1\t200\trs2\tA\tATT\t50\tPASS\t.\tGT\t0/1",
        ])
        result = read_vcf(write_vcf(tmp_path, "a.vcf", text))
        assert result["skipped"]["indel"] == 1
        assert result["snps"][0]["rsid"] == "rs1"

    def test_reference_blocks_are_counted_separately(self, tmp_path):
        text = build_vcf([
            "1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=5000\tGT:DP\t0/0:30",
            record(rsid="rs1", pos=6000),
        ])
        result = read_vcf(write_vcf(tmp_path, "g.vcf", text))
        assert result["ref_blocks"] == 1
        assert result["gvcf"] is True

    def test_non_pass_records_are_kept_and_counted(self, tmp_path):
        text = build_vcf(["1\t100\trs1\tA\tG\t50\tLowQual\t.\tGT\t0/1"])
        result = read_vcf(write_vcf(tmp_path, "a.vcf", text))
        assert result["non_pass"] == 1
        assert len(result["snps"]) == 1
        assert any("FILTER other than PASS" in w for w in result["warnings"])

    def test_blank_lines_in_the_body_are_ignored(self, tmp_path):
        text = build_vcf([record(rsid="rs1"), "", record(rsid="rs2", pos=300)])
        result = read_vcf(write_vcf(tmp_path, "a.vcf", text))
        assert len(result["snps"]) == 2

    def test_multi_sample_defaults_to_first_and_warns(self, tmp_path):
        text = build_vcf(
            ["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1"],
            samples=("MUM", "DAD"),
        )
        result = read_vcf(write_vcf(tmp_path, "a.vcf", text))
        assert result["sample"] == "MUM"
        assert result["sample_count"] == 2
        assert any("samples are present" in w for w in result["warnings"])

    def test_named_sample_selection_changes_the_genotype(self, tmp_path):
        text = build_vcf(
            ["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1"],
            samples=("MUM", "DAD"),
        )
        path = write_vcf(tmp_path, "a.vcf", text)
        mum = read_vcf(path, sample="MUM")["snps"][0]
        dad = read_vcf(path, sample="DAD")["snps"][0]
        assert (mum["allele1"], mum["allele2"]) == ("A", "A")
        assert (dad["allele1"], dad["allele2"]) == ("G", "G")

    def test_single_sample_file_produces_no_ambiguity_warning(self, tmp_path):
        result = read_vcf(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert not any("samples are present" in w for w in result["warnings"])

    def test_iter_vcf_records_streams_every_data_line(self, tmp_path):
        text = build_vcf([
            record(rsid="rs1"),
            "1\t200\trs2\tA\tATT\t50\tPASS\t.\tGT\t0/1",
        ])
        path = write_vcf(tmp_path, "a.vcf", text)
        records = list(iter_vcf_records(path))
        assert len(records) == 2
        assert [r["skip"] for r in records] == [None, "indel"]

    def test_iter_vcf_records_on_missing_file_raises(self, tmp_path):
        with pytest.raises(ParseError, match="File not found"):
            list(iter_vcf_records(tmp_path / "nope.vcf"))


# ---------------------------------------------------------------------------
# Container detection and gzip
# ---------------------------------------------------------------------------

class TestDetectFormat:
    def test_plain_vcf(self, tmp_path):
        path = write_vcf(tmp_path, "a.vcf", build_vcf([record()]))
        assert detect_format(path) == "vcf"

    def test_gzipped_vcf(self, tmp_path):
        path = write_gzipped_vcf(tmp_path, "a.vcf.gz", build_vcf([record()]))
        assert detect_format(path) == "vcf.gz"

    def test_gzip_is_detected_by_magic_not_by_suffix(self, tmp_path):
        # A gzipped file named .vcf is still gzip, and guessing from the suffix
        # produces a decode error that tells the user nothing.
        path = write_gzipped_vcf(tmp_path, "misnamed.vcf", build_vcf([record()]))
        assert detect_format(path) == "vcf.gz"

    def test_uncompressed_bam_magic(self, tmp_path):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        assert detect_format(path) == "bam"

    def test_cram_magic(self, tmp_path):
        path = tmp_path / "a.cram"
        path.write_bytes(b"CRAM\x03\x00rubbish")
        assert detect_format(path) == "cram"

    def test_unrelated_text_is_unknown(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_text("rsid\tchromosome\tposition\tgenotype\n", encoding="utf-8")
        assert detect_format(path) == "unknown"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ParseError, match="File not found"):
            detect_format(tmp_path / "nope.vcf")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

class TestParseSequencingFile:
    def test_returns_the_documented_key_set(self, tmp_path):
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert set(result) == {
            "provider", "format", "build", "build_confidence", "build_evidence",
            "sample", "sample_count", "snp_count", "snps", "skipped", "warnings",
        }

    def test_provider_is_vcf(self, tmp_path):
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert result["provider"] == "vcf"
        assert result["format"] == "vcf"

    def test_build_is_detected_from_contig_lengths(self, tmp_path):
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert result["build"] == "GRCh37"
        assert result["build_confidence"] == "high"

    def test_snp_count_matches_snps_length(self, tmp_path):
        text = build_vcf([record(rsid="rs1"), record(rsid="rs2", pos=200)])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["snp_count"] == len(result["snps"]) == 2

    def test_snps_match_the_parsers_dict_shape(self, tmp_path):
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", build_vcf([record()])))
        assert set(result["snps"][0]) == {
            "rsid", "chromosome", "position", "allele1", "allele2"
        }

    def test_skipped_reasons_are_all_reported(self, tmp_path):
        text = build_vcf([
            record(rsid="rs1"),
            "1\t200\trs2\tA\tATT\t50\tPASS\t.\tGT\t0/1",
            "1\t300\trs3\tA\t<DEL>\t50\tPASS\t.\tGT\t0/1",
            "1\t400\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
            "1\t500\trs5\tA\tG,AT\t50\tPASS\t.\tGT\t0/2",
            "1\t600\trs6\tA\tG\t50\tPASS\t.\tGT\t./.",
        ])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["skipped"] == {
            "indel": 1, "symbolic": 1, "no_id": 1,
            "multiallelic_complex": 1, "no_call": 1,
        }
        assert result["snp_count"] == 1

    def test_gzip_round_trip_matches_plain_text(self, tmp_path):
        text = build_vcf([record(rsid="rs1"), record(rsid="rs2", pos=200, ref="C", alt="T")])
        plain = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        gzipped = parse_sequencing_file(write_gzipped_vcf(tmp_path, "a.vcf.gz", text))
        assert plain["snps"] == gzipped["snps"]
        assert gzipped["format"] == "vcf.gz"

    def test_gvcf_provider_is_reported(self, tmp_path):
        text = build_vcf(
            [
                "1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=5000\tGT:DP\t0/0:30",
                "1\t6000\trs1\tA\tG,<NON_REF>\t50\tPASS\t.\tGT\t0/1",
            ],
            meta=['##ALT=<ID=NON_REF,Description="any possible alternative allele">'],
        )
        result = parse_sequencing_file(write_vcf(tmp_path, "a.g.vcf", text))
        assert result["provider"] == "gvcf"

    def test_gvcf_reference_block_is_not_expanded(self, tmp_path):
        text = build_vcf([
            "1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=5000\tGT:DP\t0/0:30",
            "1\t6000\trs1\tA\tG,<NON_REF>\t50\tPASS\t.\tGT\t0/1",
        ])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.g.vcf", text))
        # One call, not 4901 of them.
        assert result["snp_count"] == 1
        assert result["skipped"]["symbolic"] == 1
        assert any("reference block" in w for w in result["warnings"])

    def test_grch38_file_raises_build_mismatch(self, tmp_path):
        text = build_vcf([record()], contigs=GRCH38_CONTIGS)
        with pytest.raises(BuildMismatch):
            parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))

    def test_grch38_file_is_accepted_when_that_is_what_was_expected(self, tmp_path):
        text = build_vcf([record()], contigs=GRCH38_CONTIGS)
        result = parse_sequencing_file(
            write_vcf(tmp_path, "a.vcf", text), expected_build="GRCh38"
        )
        assert result["build"] == "GRCh38"

    def test_unknown_build_parses_but_warns_loudly(self, tmp_path):
        text = build_vcf([record()], contigs=[])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["build"] is None
        assert result["build_confidence"] == "none"
        assert any("could not be determined" in w for w in result["warnings"])

    def test_unknown_build_is_not_silently_defaulted(self, tmp_path):
        text = build_vcf([record()], contigs=[])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["build"] != REFERENCE_BUILD

    def test_contradictory_header_text_produces_a_warning(self, tmp_path):
        text = build_vcf([record()], meta=["##reference=/refs/hg38.fa"])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["build"] == "GRCh37"
        assert any("Contig lengths say" in w for w in result["warnings"])

    def test_conflicting_contig_lengths_are_reported(self, tmp_path):
        contigs = ["##contig=<ID=1,length=249250621>", "##contig=<ID=2,length=242193529>"]
        text = build_vcf([record()], contigs=contigs)
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["build"] is None
        assert result["build_confidence"] == "conflict"
        assert any("self-contradictory" in w for w in result["warnings"])

    def test_multi_sample_reports_sample_count_and_choice(self, tmp_path):
        text = build_vcf(
            ["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1\t0/1"],
            samples=("MUM", "DAD", "KID"),
        )
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["sample_count"] == 3
        assert result["sample"] == "MUM"

    def test_multi_sample_selection_by_name(self, tmp_path):
        text = build_vcf(
            ["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1\t0/1"],
            samples=("MUM", "DAD", "KID"),
        )
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text), sample="KID")
        assert result["sample"] == "KID"
        assert (result["snps"][0]["allele1"], result["snps"][0]["allele2"]) == ("A", "G")

    def test_multi_sample_selection_by_index(self, tmp_path):
        text = build_vcf(
            ["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1\t0/1"],
            samples=("MUM", "DAD", "KID"),
        )
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text), sample=1)
        assert result["sample"] == "DAD"

    def test_unknown_sample_name_raises(self, tmp_path):
        text = build_vcf(["1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1"],
                         samples=("MUM", "DAD"))
        with pytest.raises(ParseError, match="not in this file"):
            parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text), sample="GRAN")

    def test_sites_only_vcf_raises(self, tmp_path):
        lines = ["##fileformat=VCFv4.2"] + GRCH37_CONTIGS
        lines.append("\t".join(COLUMNS[:8]))
        lines.append("1\t100\trs1\tA\tG\t50\tPASS\t.")
        path = write_vcf(tmp_path, "sites.vcf", "\n".join(lines) + "\n")
        with pytest.raises(ParseError, match="no sample columns"):
            parse_sequencing_file(path)

    def test_empty_body_raises(self, tmp_path):
        path = write_vcf(tmp_path, "a.vcf", build_vcf([]))
        with pytest.raises(ParseError, match="No variant records"):
            parse_sequencing_file(path)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ParseError, match="File not found"):
            parse_sequencing_file(tmp_path / "nope.vcf")

    def test_non_vcf_text_raises(self, tmp_path):
        path = tmp_path / "array.txt"
        path.write_text("rsid\tchromosome\tposition\tgenotype\n", encoding="utf-8")
        with pytest.raises(ParseError, match="not a format this module reads"):
            parse_sequencing_file(path)

    def test_all_records_skipped_is_a_warning_not_an_error(self, tmp_path):
        text = build_vcf(["1\t100\trs1\tA\tATT\t50\tPASS\t.\tGT\t0/1"])
        result = parse_sequencing_file(write_vcf(tmp_path, "a.vcf", text))
        assert result["snp_count"] == 0
        assert result["skipped"]["indel"] == 1
        assert any("No annotatable SNVs" in w for w in result["warnings"])

    def test_bam_degrades_instead_of_raising(self, tmp_path):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        result = parse_sequencing_file(path)
        assert result["provider"] == "bam"
        assert result["snp_count"] == 0
        assert result["build"] is None
        assert any("targeted extraction" in w for w in result["warnings"])

    def test_cram_degrades_instead_of_raising(self, tmp_path):
        path = tmp_path / "a.cram"
        path.write_bytes(b"CRAM\x03\x00rubbish")
        result = parse_sequencing_file(path)
        assert result["provider"] == "cram"
        assert result["snps"] == []


# ---------------------------------------------------------------------------
# Liftover
# ---------------------------------------------------------------------------

class TestLiftoverUnavailable:
    def test_liftover_available_is_false_without_a_chain_file(self, isolated_home):
        assert liftover_available("GRCh38", "GRCh37") is False

    def test_find_chain_file_returns_none(self, isolated_home):
        assert find_chain_file("GRCh38", "GRCh37") is None

    def test_liftover_returns_the_unavailable_payload(self, isolated_home):
        payload = liftover([{"rsid": "rs1", "chromosome": "1", "position": 100}],
                           "GRCh38", "GRCh37")
        assert payload["available"] is False
        assert payload["not_attempted"] is True

    def test_unavailable_payload_has_the_external_shape(self, isolated_home):
        payload = liftover([], "GRCh38", "GRCh37")
        expected = {"available", "capability", "tool", "tool_id", "state",
                    "reason", "not_attempted", "results", "how_to_enable"}
        assert expected.issubset(set(payload))

    def test_unavailable_payload_returns_no_results(self, isolated_home):
        # The failure mode this guards against is handing back the input
        # coordinates relabelled as the destination build.
        payload = liftover([{"rsid": "rs1", "chromosome": "1", "position": 100}],
                           "GRCh38", "GRCh37")
        assert payload["results"] == []

    def test_unavailable_reason_states_coordinates_were_not_converted(self, isolated_home):
        payload = liftover([], "GRCh38", "GRCh37")
        assert "NOT been relabelled" in payload["reason"]

    def test_unavailable_payload_explains_how_to_enable(self, isolated_home):
        payload = liftover([], "GRCh38", "GRCh37")
        assert str(chain_dir()) == payload["how_to_enable"]["install_to"]

    def test_same_build_liftover_is_a_no_op(self, isolated_home):
        records = [{"rsid": "rs1", "chromosome": "1", "position": 100}]
        payload = liftover(records, "GRCh37", "GRCh37")
        assert payload["available"] is True
        assert payload["results"] == records


class TestChainParser:
    def test_chain_count(self):
        assert ChainFile.from_text(SIMPLE_CHAIN).chain_count == 2

    def test_first_block_start_maps(self):
        hit = ChainFile.from_text(SIMPLE_CHAIN).map_position("1", 101)
        assert (hit["chromosome"], hit["position"]) == ("1", 501)

    def test_position_inside_first_block_maps(self):
        hit = ChainFile.from_text(SIMPLE_CHAIN).map_position("1", 150)
        assert hit["position"] == 550

    def test_position_after_the_gap_maps_with_the_new_offset(self):
        hit = ChainFile.from_text(SIMPLE_CHAIN).map_position("1", 161)
        assert hit["position"] == 571

    def test_position_inside_a_gap_is_unmappable(self):
        assert ChainFile.from_text(SIMPLE_CHAIN).map_position("1", 155) is None

    def test_position_outside_the_chain_is_unmappable(self):
        assert ChainFile.from_text(SIMPLE_CHAIN).map_position("1", 900000) is None

    def test_unknown_contig_is_unmappable(self):
        assert ChainFile.from_text(SIMPLE_CHAIN).map_position("21", 101) is None

    def test_chr_prefixed_query_matches_a_bare_chain(self):
        hit = ChainFile.from_text(SIMPLE_CHAIN).map_position("chr1", 101)
        assert hit is not None

    def test_minus_strand_chain_flips_the_coordinate(self):
        hit = ChainFile.from_text(SIMPLE_CHAIN).map_position("2", 1001)
        assert hit["strand"] == "-"
        assert hit["position"] == 242193529 - 2000

    def test_malformed_chain_header_raises(self):
        with pytest.raises(ParseError, match="chain header"):
            ChainFile.from_text("chain 4900 chr1 249250621 +\n10\n")

    def test_missing_chain_file_raises(self, tmp_path):
        with pytest.raises(ParseError, match="Chain file not found"):
            ChainFile.load(tmp_path / "nope.chain")


class TestLiftoverWithChain:
    @pytest.fixture
    def chain_installed(self, isolated_home):
        target = chain_dir()
        target.mkdir(parents=True, exist_ok=True)
        (target / "GRCh37_to_GRCh38.chain").write_text(SIMPLE_CHAIN, encoding="utf-8")
        return target

    def test_liftover_available_is_true(self, chain_installed):
        assert liftover_available("GRCh37", "GRCh38") is True

    def test_gzipped_chain_file_is_found(self, isolated_home):
        target = chain_dir()
        target.mkdir(parents=True, exist_ok=True)
        with gzip.open(target / "hg19ToHg38.over.chain.gz", "wt", encoding="utf-8") as fh:
            fh.write(SIMPLE_CHAIN)
        assert find_chain_file("GRCh37", "GRCh38") is not None

    def test_records_are_mapped(self, chain_installed):
        records = [{"rsid": "rs1", "chromosome": "1", "position": 101,
                    "allele1": "A", "allele2": "G"}]
        payload = liftover(records, "GRCh37", "GRCh38")
        assert payload["available"] is True
        assert payload["results"][0]["position"] == 501
        assert payload["mapped_count"] == 1

    def test_unmappable_records_are_reported_not_passed_through(self, chain_installed):
        records = [{"rsid": "rs1", "chromosome": "1", "position": 155,
                    "allele1": "A", "allele2": "G"}]
        payload = liftover(records, "GRCh37", "GRCh38")
        assert payload["results"] == []
        assert payload["unmapped_count"] == 1
        assert "no GRCh37 to GRCh38 chain block" in payload["unmapped"][0]["reason"]

    def test_minus_strand_mapping_complements_the_alleles(self, chain_installed):
        records = [{"rsid": "rs2", "chromosome": "2", "position": 1001,
                    "allele1": "A", "allele2": "G"}]
        payload = liftover(records, "GRCh37", "GRCh38")
        out = payload["results"][0]
        assert out["strand_flipped"] is True
        assert (out["allele1"], out["allele2"]) == ("T", "C")

    def test_plus_strand_mapping_leaves_alleles_alone(self, chain_installed):
        records = [{"rsid": "rs1", "chromosome": "1", "position": 101,
                    "allele1": "A", "allele2": "G"}]
        out = liftover(records, "GRCh37", "GRCh38")["results"][0]
        assert (out["allele1"], out["allele2"]) == ("A", "G")
        assert out["strand_flipped"] is False

    def test_chain_file_path_is_reported(self, chain_installed):
        payload = liftover([], "GRCh37", "GRCh38")
        assert payload["chain_file"].endswith("GRCh37_to_GRCh38.chain")


# ---------------------------------------------------------------------------
# BED writer
# ---------------------------------------------------------------------------

class TestWritePositionsBed:
    def test_intervals_are_zero_based_half_open(self, tmp_path):
        dest = write_positions_bed(
            [{"chromosome": "1", "position": 100, "rsid": "rs1"}], tmp_path / "p.bed"
        )
        assert dest.read_text(encoding="utf-8") == "1\t99\t100\trs1\n"

    def test_rows_are_sorted_by_chromosome_then_position(self, tmp_path):
        dest = write_positions_bed(
            [("X", 5, "rsX"), ("2", 10, "rsB"), ("1", 30, "rsA"), ("1", 20, "rsC")],
            tmp_path / "p.bed",
        )
        names = [line.split("\t")[3] for line in dest.read_text(encoding="utf-8").splitlines()]
        assert names == ["rsC", "rsA", "rsB", "rsX"]

    def test_mitochondria_sort_last(self, tmp_path):
        dest = write_positions_bed(
            [("MT", 1, "rsM"), ("Y", 1, "rsY"), ("22", 1, "rs22")], tmp_path / "p.bed"
        )
        names = [line.split("\t")[3] for line in dest.read_text(encoding="utf-8").splitlines()]
        assert names == ["rs22", "rsY", "rsM"]

    def test_duplicate_positions_are_written_once(self, tmp_path):
        dest = write_positions_bed(
            [("1", 100, "rs1"), ("1", 100, "rs1")], tmp_path / "p.bed"
        )
        assert len(dest.read_text(encoding="utf-8").splitlines()) == 1

    def test_chrom_prefix_is_applied(self, tmp_path):
        dest = write_positions_bed(
            [("1", 100, "rs1")], tmp_path / "p.bed", chrom_prefix="chr"
        )
        assert dest.read_text(encoding="utf-8").startswith("chr1\t")

    def test_chr_prefixed_input_is_normalised(self, tmp_path):
        dest = write_positions_bed([("chr1", 100, "rs1")], tmp_path / "p.bed")
        assert dest.read_text(encoding="utf-8").startswith("1\t")

    def test_zero_position_is_rejected(self, tmp_path):
        with pytest.raises(ParseError, match="1-based"):
            write_positions_bed([("1", 0, "rs1")], tmp_path / "p.bed")

    def test_entry_without_a_chromosome_is_rejected(self, tmp_path):
        with pytest.raises(ParseError, match="no chromosome"):
            write_positions_bed([{"position": 100, "rsid": "rs1"}], tmp_path / "p.bed")

    def test_unreadable_entry_is_rejected(self, tmp_path):
        with pytest.raises(ParseError, match="Unreadable position"):
            write_positions_bed(["nonsense"], tmp_path / "p.bed")


# ---------------------------------------------------------------------------
# Pileup parsing and calling
# ---------------------------------------------------------------------------

class TestPileupParsing:
    def test_reference_matches_are_counted_as_the_reference_base(self):
        row = parse_mpileup_line("1\t100\tA\t4\t....\tIIII")
        assert row["counts"]["A"] == 4

    def test_reverse_strand_reference_matches_count_too(self):
        row = parse_mpileup_line("1\t100\tA\t4\t..,,\tIIII")
        assert row["counts"]["A"] == 4

    def test_mismatches_are_counted_by_base(self):
        row = parse_mpileup_line("1\t100\tA\t4\t..GG\tIIII")
        assert row["counts"]["A"] == 2 and row["counts"]["G"] == 2

    def test_lowercase_mismatches_count_the_same(self):
        row = parse_mpileup_line("1\t100\tA\t2\tgg\tII")
        assert row["counts"]["G"] == 2

    def test_read_start_mapping_quality_char_is_not_an_allele(self):
        # "^A." is a read start whose mapping quality happens to render as 'A'.
        row = parse_mpileup_line("1\t100\tC\t2\t^A.^].\tII")
        assert row["counts"]["A"] == 0
        assert row["counts"]["C"] == 2

    def test_read_end_marker_is_ignored(self):
        row = parse_mpileup_line("1\t100\tC\t2\t.$.\tII")
        assert row["counts"]["C"] == 2

    def test_insertion_sequence_is_not_counted(self):
        row = parse_mpileup_line("1\t100\tC\t2\t.+2AG.\tII")
        assert row["counts"]["A"] == 0 and row["counts"]["G"] == 0
        assert row["counts"]["C"] == 2

    def test_deletion_sequence_is_not_counted(self):
        row = parse_mpileup_line("1\t100\tC\t2\t.-2ag.\tII")
        assert row["counts"]["A"] == 0 and row["counts"]["G"] == 0

    def test_deleted_base_placeholder_is_not_an_allele(self):
        row = parse_mpileup_line("1\t100\tC\t3\t.*.\tIII")
        assert sum(row["counts"].values()) == 2

    def test_chromosome_is_normalised(self):
        row = parse_mpileup_line("chr1\t100\tA\t1\t.\tI")
        assert row["chromosome"] == "1"

    def test_blank_line_returns_none(self):
        assert parse_mpileup_line("\n") is None

    def test_short_row_raises(self):
        with pytest.raises(ParseError, match="mpileup"):
            parse_mpileup_line("1\t100\tA")


class TestCallGenotypeFromCounts:
    def test_thin_coverage_is_a_no_call(self):
        assert call_genotype_from_counts({"A": 3, "G": 0, "C": 0, "T": 0}) == ("N", "N")

    def test_homozygous_call(self):
        counts = {"A": 20, "G": 0, "C": 0, "T": 0}
        assert call_genotype_from_counts(counts) == ("A", "A")

    def test_heterozygous_call_is_returned_sorted(self):
        counts = {"A": 10, "G": 10, "C": 0, "T": 0}
        assert call_genotype_from_counts(counts) == ("A", "G")

    def test_low_fraction_noise_is_ignored(self):
        counts = {"A": 19, "G": 1, "C": 0, "T": 0}
        assert call_genotype_from_counts(counts) == ("A", "A")

    def test_minimum_depth_is_the_documented_constant(self):
        counts = {"A": MIN_PILEUP_DEPTH, "G": 0, "C": 0, "T": 0}
        assert call_genotype_from_counts(counts) == ("A", "A")

    def test_empty_counts_are_a_no_call(self):
        assert call_genotype_from_counts({}) == ("N", "N")


# ---------------------------------------------------------------------------
# Alignment extraction, degraded
# ---------------------------------------------------------------------------

class TestExtractPositions:
    def test_samtools_tool_id_is_recorded(self):
        assert SAMTOOLS_TOOL_ID == "samtools"

    def test_samtools_is_registered(self):
        # Updated by the v3.0 wiring pass, which added the entry this module
        # was written against. The predecessor of this test asserted the entry
        # was ABSENT, which was true while the registry had not been touched.
        # Keeping the assertion in place, inverted, means the registry and this
        # module cannot drift apart silently in either direction.
        entry = external.get(SAMTOOLS_TOOL_ID)
        assert entry is not None
        assert entry["capability"] == SAMTOOLS_CAPABILITY
        assert entry["spdx"] == "MIT"

    def test_samtools_registration_does_not_make_it_available(self, isolated_home):
        # Registered is not installed, and installed is not licence-accepted.
        # Every alignment path must still degrade on a machine without the
        # binary, which is most machines.
        assert external.is_available(SAMTOOLS_TOOL_ID) is False

    def test_extraction_degrades_without_samtools(self, tmp_path, isolated_home):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        payload = extract_positions(path, [("1", 100, "rs1")])
        assert payload["available"] is False
        assert payload["not_attempted"] is True
        assert payload["results"] == []

    def test_degraded_payload_explains_the_registry_gap(self, tmp_path, isolated_home):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        payload = extract_positions(path, [("1", 100, "rs1")])
        assert "not attempted" in payload["reason"]

    def test_wrong_build_raises_before_anything_runs(self, tmp_path):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        with pytest.raises(BuildMismatch):
            extract_positions(path, [("1", 100, "rs1")], build="GRCh38")

    def test_cram_without_a_reference_refuses(self, tmp_path, isolated_home):
        path = tmp_path / "a.cram"
        path.write_bytes(b"CRAM\x03\x00rubbish")
        payload = extract_positions(path, [("1", 100, "rs1")])
        assert payload["state"] == "reference_missing"
        assert payload["available"] is False

    def test_missing_alignment_raises(self, tmp_path):
        with pytest.raises(ParseError, match="Alignment file not found"):
            extract_positions(tmp_path / "nope.bam", [("1", 100, "rs1")])

    def test_non_alignment_input_raises(self, tmp_path):
        path = write_vcf(tmp_path, "a.vcf", build_vcf([record()]))
        with pytest.raises(ParseError, match="not a BAM or CRAM"):
            extract_positions(path, [("1", 100, "rs1")])

    def test_empty_position_list_raises(self, tmp_path):
        path = tmp_path / "a.bam"
        path.write_bytes(b"BAM\x01rubbish")
        with pytest.raises(ParseError, match="No positions"):
            extract_positions(path, [])
