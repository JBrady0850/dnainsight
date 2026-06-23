"""
parsers.py -- DNA file format detection and parsing.

Supported providers:
    - AncestryDNA (V1/V2 array, tab-delimited, allele1 + allele2 columns)
    - 23andMe (tab-delimited, single 'genotype' column)
    - MyHeritage (tab-delimited, similar to 23andMe)
    - FamilyTreeDNA (FTDNA) (tab-delimited, single 'RESULT' column)
    - LivingDNA (tab-delimited, single 'genotype' column)
    - Generic (attempts best-guess detection)

Returns a list of dicts: [{rsid, chromosome, position, allele1, allele2}, ...]
"""

import re
from pathlib import Path


class ParseError(Exception):
    pass


# ---------------------------------------------------------------------------
# Provider detection
# ---------------------------------------------------------------------------

def detect_provider(header_lines: list[str], column_line: str) -> str:
    """Identify DNA provider from comment block and column headers."""
    combined = " ".join(header_lines + [column_line]).lower()

    if "ancestrydna" in combined or "ancestry.com" in combined:
        return "ancestrydna"
    if "23andme" in combined:
        return "23andme"
    if "myheritage" in combined:
        return "myheritage"
    if "familytreedna" in combined or "ftdna" in combined or "familytree dna" in combined:
        return "ftdna"
    if "livingdna" in combined:
        return "livingdna"
    # Fallback: inspect column headers
    cols = [c.strip().lower() for c in column_line.split("\t")]
    if "allele1" in cols and "allele2" in cols:
        return "ancestrydna"
    if "genotype" in cols:
        return "23andme"
    if "result" in cols:
        return "ftdna"
    return "generic"


def _split_genotype(gt: str) -> tuple[str, str]:
    """Split a 1-2 char genotype string into (allele1, allele2)."""
    gt = gt.strip().upper()
    if len(gt) == 2:
        return gt[0], gt[1]
    if len(gt) == 1:
        return gt, gt
    if "--" in gt or gt in ("0", "00", "NN", "--"):
        return "N", "N"
    return gt[:1], gt[1:2] if len(gt) > 1 else "N"


# ---------------------------------------------------------------------------
# Format-specific parsers
# ---------------------------------------------------------------------------

def _parse_tsv_dual_allele(lines: list[str], col_map: dict) -> list[dict]:
    """Parse files with separate allele1 / allele2 columns (AncestryDNA style)."""
    results = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < max(col_map.values()) + 1:
            continue
        rsid = parts[col_map["rsid"]].strip()
        if not rsid.startswith("rs"):
            continue
        try:
            chrom = parts[col_map["chromosome"]].strip()
            pos   = int(parts[col_map["position"]].strip())
            a1    = parts[col_map["allele1"]].strip().upper()
            a2    = parts[col_map["allele2"]].strip().upper()
        except (ValueError, IndexError):
            continue
        if a1 in ("0", "") or a2 in ("0", ""):
            continue
        results.append({
            "rsid": rsid, "chromosome": chrom, "position": pos,
            "allele1": a1, "allele2": a2,
        })
    return results


def _parse_tsv_single_genotype(lines: list[str], col_map: dict, gt_col: str) -> list[dict]:
    """Parse files with a single merged genotype column (23andMe / MyHeritage style)."""
    results = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < max(col_map.values()) + 1:
            continue
        rsid = parts[col_map["rsid"]].strip()
        if not rsid.startswith("rs"):
            continue
        try:
            chrom = parts[col_map["chromosome"]].strip()
            pos   = int(parts[col_map["position"]].strip())
            gt    = parts[col_map[gt_col]].strip().upper()
        except (ValueError, IndexError):
            continue
        if gt in ("--", "00", "NN", "DD", "II", "DI"):
            continue
        a1, a2 = _split_genotype(gt)
        results.append({
            "rsid": rsid, "chromosome": chrom, "position": pos,
            "allele1": a1, "allele2": a2,
        })
    return results


def _col_index(header: list[str], *names: str) -> int | None:
    for name in names:
        for i, h in enumerate(header):
            if h.strip().lower() == name.lower():
                return i
    return None


# ---------------------------------------------------------------------------
# Public parse function
# ---------------------------------------------------------------------------

def parse_dna_file(filepath: str) -> dict:
    """
    Parse any supported DNA file.

    Returns:
        {
            "provider": str,
            "snp_count": int,
            "snps": [{"rsid", "chromosome", "position", "allele1", "allele2"}, ...]
        }
    """
    path = Path(filepath)
    if not path.exists():
        raise ParseError(f"File not found: {filepath}")

    content = path.read_text(encoding="utf-8", errors="replace")
    lines   = content.splitlines()

    # Separate comment block from data
    header_lines = [l[1:].strip() for l in lines if l.startswith("#")]
    data_lines   = [l for l in lines if not l.startswith("#") and l.strip()]

    if not data_lines:
        raise ParseError("No data found in file. Verify the file is an uncompressed raw DNA export.")

    # First non-comment line is the column header
    col_header = data_lines[0]
    data_rows  = data_lines[1:]

    provider = detect_provider(header_lines, col_header)
    cols = [c.strip().lower() for c in col_header.split("\t")]

    # --- AncestryDNA ---
    if provider == "ancestrydna":
        col_map = {
            "rsid":       _col_index(cols, "rsid") or 0,
            "chromosome": _col_index(cols, "chromosome") or 1,
            "position":   _col_index(cols, "position") or 2,
            "allele1":    _col_index(cols, "allele1") or 3,
            "allele2":    _col_index(cols, "allele2") or 4,
        }
        snps = _parse_tsv_dual_allele(data_rows, col_map)

    # --- 23andMe / MyHeritage / LivingDNA / Generic with genotype col ---
    elif provider in ("23andme", "myheritage", "livingdna"):
        gt_col = "genotype"
        col_map = {
            "rsid":       _col_index(cols, "rsid", "snpid", "snp id") or 0,
            "chromosome": _col_index(cols, "chromosome", "chr") or 1,
            "position":   _col_index(cols, "position", "coordinate") or 2,
            "genotype":   _col_index(cols, "genotype", "alleles") or 3,
        }
        snps = _parse_tsv_single_genotype(data_rows, col_map, gt_col)

    # --- FamilyTreeDNA ---
    elif provider == "ftdna":
        gt_col = "result"
        col_map = {
            "rsid":       _col_index(cols, "rsid", "snpid") or 0,
            "chromosome": _col_index(cols, "chromosome", "chr") or 1,
            "position":   _col_index(cols, "position") or 2,
            "result":     _col_index(cols, "result", "genotype") or 3,
        }
        snps = _parse_tsv_single_genotype(data_rows, col_map, gt_col)

    # --- Generic fallback ---
    else:
        if "allele1" in cols or "allele 1" in cols:
            col_map = {
                "rsid":       _col_index(cols, "rsid", "snpid") or 0,
                "chromosome": _col_index(cols, "chromosome", "chr") or 1,
                "position":   _col_index(cols, "position") or 2,
                "allele1":    _col_index(cols, "allele1", "allele 1") or 3,
                "allele2":    _col_index(cols, "allele2", "allele 2") or 4,
            }
            snps = _parse_tsv_dual_allele(data_rows, col_map)
        else:
            gt_col = "genotype"
            col_map = {
                "rsid":       _col_index(cols, "rsid", "snpid") or 0,
                "chromosome": _col_index(cols, "chromosome", "chr") or 1,
                "position":   _col_index(cols, "position") or 2,
                "genotype":   _col_index(cols, "genotype", "result", "alleles") or 3,
            }
            snps = _parse_tsv_single_genotype(data_rows, col_map, gt_col)

    if not snps:
        raise ParseError(
            f"Parsed 0 valid SNPs from {path.name}. "
            "Ensure the file is an uncompressed raw DNA export, not a compressed (.zip/.gz) file."
        )

    return {"provider": provider, "snp_count": len(snps), "snps": snps}
