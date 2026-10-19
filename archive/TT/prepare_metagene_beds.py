#!/usr/bin/env python3
"""Prepare unique gene-level TSS BEDs split by significant fold-change direction.

Result row IDs are RefSeq transcript accessions. The provided RefSeq GTF maps
significant rows to gene symbols. Each gene is included once, using its MANE
gene interval/TSS when available; remaining genes use the outermost RefSeq
transcript span from the same hg38 GTF.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pyBigWig


ROOT = Path("/mnt/d/work/pby")
PADJ_CUTOFF = 0.05
LOG2FC_CUTOFF: float | None = None
RESULTS = {
    "MED12": ROOT / "TT/med12/results.txt",
    "SPT6": ROOT / "TT/spt6/results.txt",
}
GTF = Path("/mnt/d/work/reference/hg38/refseq_sort.gtf")
MANE_BED = Path("/mnt/d/work/reference/hg38/MANE_gene.bed")
FAI = Path("/mnt/d/work/reference/hg38/genome.fa.fai")
BIGWIG = ROOT / "MNase/CTRL.bw"
OUT_DIRS = {
    "MED12": ROOT / "TT/med12",
    "SPT6": ROOT / "TT/spt6",
}

TX_RE = re.compile(r'(?:^|;)\s*transcript_id "([^"]+)"')
GENE_ID_RE = re.compile(r'(?:^|;)\s*gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'(?:^|;)\s*gene_name "([^"]+)"')


def read_significant_ids(path: Path) -> list[tuple[str, float, float]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if "padj" not in header or "log2FoldChange" not in header:
            raise ValueError(f"Expected padj and log2FoldChange columns in {path}")
        padj_idx = header.index("padj")
        lfc_idx = header.index("log2FoldChange")
        significant: list[tuple[str, float, float]] = []
        seen: set[str] = set()
        for row in reader:
            if not row:
                continue
            transcript_id = row[0].strip()
            if not transcript_id:
                continue
            try:
                padj = float(row[padj_idx])
                log2fc = float(row[lfc_idx])
            except (ValueError, IndexError):
                continue
            if (
                math.isfinite(padj)
                and padj < PADJ_CUTOFF
                and math.isfinite(log2fc)
                and (LOG2FC_CUTOFF is None or abs(log2fc) >= LOG2FC_CUTOFF)
            ):
                if transcript_id in seen:
                    raise ValueError(f"Duplicate result ID {transcript_id} in {path}")
                seen.add(transcript_id)
                significant.append((transcript_id, padj, log2fc))
    return significant


def read_fai(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                sizes[fields[0]] = int(fields[1])
    return sizes


def attribute(pattern: re.Pattern[str], text: str) -> str | None:
    match = pattern.search(text)
    return match.group(1) if match else None


def read_mane_bed(path: Path) -> dict[str, set[tuple[str, int, int, str]]]:
    genes: dict[str, set[tuple[str, int, int, str]]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start, end, gene, _score, strand = fields[:6]
            genes[gene].add((chrom, int(start), int(end), strand))
    return genes


def main() -> None:
    significant = {name: read_significant_ids(path) for name, path in RESULTS.items()}
    all_ids = {transcript_id for rows in significant.values() for transcript_id, _, _ in rows}
    all_bases = {transcript_id.rsplit(".", 1)[0] for transcript_id in all_ids}

    bw = pyBigWig.open(str(BIGWIG))
    if bw is None:
        raise RuntimeError(f"Could not open BigWig: {BIGWIG}")
    bw_chroms = bw.chroms()
    bw.close()
    fai_sizes = read_fai(FAI)
    mismatches = [
        (chrom, bw_chroms[chrom], fai_sizes.get(chrom))
        for chrom in bw_chroms
        if chrom in fai_sizes and bw_chroms[chrom] != fai_sizes[chrom]
    ]
    if mismatches:
        raise ValueError(f"BigWig and hg38 FASTA chromosome lengths differ: {mismatches[:5]}")
    common_chroms = set(bw_chroms) & set(fai_sizes)
    if len(common_chroms) < 20:
        raise ValueError(f"Too few BigWig/hg38 chromosome names agree ({len(common_chroms)})")

    exact_gene: dict[str, str] = {}
    genes_by_base: dict[str, set[str]] = defaultdict(set)
    gene_models: dict[tuple[str, str, str], list[int]] = {}
    with GTF.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            transcript_id = attribute(TX_RE, attrs)
            if transcript_id is None:
                continue
            gene_id = attribute(GENE_ID_RE, attrs) or ""
            gene_name = attribute(GENE_NAME_RE, attrs) or gene_id
            if not gene_name:
                continue
            base_id = transcript_id.rsplit(".", 1)[0]
            if transcript_id in all_ids:
                old = exact_gene.get(transcript_id)
                if old is not None and old != gene_name:
                    raise ValueError(f"Conflicting gene mapping for {transcript_id}: {old} vs {gene_name}")
                exact_gene[transcript_id] = gene_name
            if base_id in all_bases:
                genes_by_base[base_id].add(gene_name)

            chrom, strand = fields[0], fields[6]
            if chrom not in common_chroms or strand not in {"+", "-"}:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            key = (gene_name, chrom, strand)
            span = gene_models.get(key)
            if span is None:
                gene_models[key] = [start, end]
            else:
                span[0] = min(span[0], start)
                span[1] = max(span[1], end)

    mane_genes = read_mane_bed(MANE_BED)

    gene_best_result: dict[str, dict[str, tuple[float, float, str]]] = {}
    unmapped_ids: dict[str, list[str]] = {}
    for name, rows in significant.items():
        selected: dict[str, tuple[float, float, str]] = {}
        directions_by_gene: dict[str, set[str]] = defaultdict(set)
        missing: list[str] = []
        for transcript_id, padj, log2fc in rows:
            gene_name = exact_gene.get(transcript_id)
            if gene_name is None:
                candidates = genes_by_base.get(transcript_id.rsplit(".", 1)[0], set())
                if len(candidates) == 1:
                    gene_name = next(iter(candidates))
            if gene_name is None:
                missing.append(transcript_id)
                continue
            directions_by_gene[gene_name].add("UP" if log2fc > 0 else "DOWN")
            previous = selected.get(gene_name)
            if previous is None or padj < previous[0]:
                selected[gene_name] = (padj, log2fc, transcript_id)
        gene_best_result[name] = selected
        unmapped_ids[name] = missing

    models_by_gene: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    for (gene_name, chrom, strand), span in gene_models.items():
        models_by_gene[gene_name].append((chrom, span[0], span[1], strand))

    for name, selected in gene_best_result.items():
        mane_count = 0
        refseq_fallback_count = 0
        unavailable: list[str] = []
        bed_rows: dict[str, list[str]] = {"ALL": [], "UP": [], "DOWN": []}
        for gene_name in sorted(selected):
            mane_records = mane_genes.get(gene_name, set())
            record: tuple[str, int, int, str] | None = None
            if len(mane_records) == 1:
                candidate = next(iter(mane_records))
                if candidate[0] in common_chroms:
                    record = candidate
                    mane_count += 1
            if record is None:
                model_records = models_by_gene.get(gene_name, [])
                if len(model_records) == 1:
                    record = model_records[0]
                    refseq_fallback_count += 1
            if record is None:
                unavailable.append(gene_name)
                continue
            chrom, start, end, strand = record
            if start < 0 or end <= start or end > bw_chroms[chrom]:
                raise ValueError(f"Invalid gene BED interval for {gene_name}: {record}")
            bed_row = f"{chrom}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n"
            bed_rows["ALL"].append(bed_row)
            direction = "UP" if selected[gene_name][1] > 0 else "DOWN"
            bed_rows[direction].append(bed_row)
        output_paths = {
            "ALL": OUT_DIRS[name] / f"{name.lower()}_PADJ0.05_TSS.bed",
            "UP": OUT_DIRS[name] / f"{name.lower()}_PADJ0.05_UP_TSS.bed",
            "DOWN": OUT_DIRS[name] / f"{name.lower()}_PADJ0.05_DOWN_TSS.bed",
        }
        for group, out_path in output_paths.items():
            with out_path.open("w", encoding="utf-8", newline="") as out:
                out.writelines(bed_rows[group])
        mixed_direction = sum(len(directions) > 1 for gene, directions in directions_by_gene.items())
        rows_written = len(bed_rows["ALL"])
        filter_label = f"padj_lt_{PADJ_CUTOFF}"
        if LOG2FC_CUTOFF is not None:
            filter_label += f"_and_abs_log2FC_ge_{LOG2FC_CUTOFF}"
        print(
            f"{name}: result_rows_{filter_label}={len(significant[name])} "
            f"mapped_unique_genes={len(selected)} unique_gene_TSS_BED_rows={rows_written} "
            f"UP_genes={len(bed_rows['UP'])} DOWN_genes={len(bed_rows['DOWN'])} "
            f"genes_with_mixed_transcript_directions={mixed_direction} "
            f"MANE_TSS={mane_count} RefSeq_GTF_fallback={refseq_fallback_count} "
            f"unmapped_result_IDs={len(unmapped_ids[name])} unavailable_gene_TSS={len(unavailable)} "
            f"UP_BED={output_paths['UP']} DOWN_BED={output_paths['DOWN']}"
        )
        if unmapped_ids[name]:
            print(f"{name} unmapped result IDs: {', '.join(unmapped_ids[name])}")
        if unavailable:
            print(f"{name} genes without a unique hg38 TSS (first 20): {', '.join(unavailable[:20])}")
        if not bed_rows["UP"] or not bed_rows["DOWN"]:
            raise ValueError(f"One fold-change direction has no gene TSS intervals for {name}")

    print(
        f"CTRL BigWig chromosome entries={len(bw_chroms)}; "
        f"matching hg38 FASTA entries={len(common_chroms)}; "
        f"length mismatches={len(mismatches)}"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
