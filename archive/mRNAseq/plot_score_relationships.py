#!/usr/bin/env python3
"""Create PDF figures for the H2AZ/RNA-seq composite score.

The score is defined as

    SCORE = z(H2AZ signal) + (MED12OE log2FC + SPT6OE log2FC) / 2

H2AZ rows are mapped from RefSeq transcript IDs to mouse gene symbols using
the UCSC mm39 refGene table. If a gene has multiple mapped transcripts, its
H2AZ signal is averaged before z-scoring so that genes are not weighted by
isoform count.

The script intentionally uses ReportLab directly so the resulting PDF is a
vector PDF and does not depend on a system R/Python plotting installation.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas


ROOT = Path(r"D:\work\pby\mRNAseq")
H2AZ_FILE = ROOT / "h2az" / "mousetss500_h2az.txt"
REFSEQ_FILE = ROOT / "h2az" / "refGene_mm39.txt.gz"
MED12_FILE = ROOT / "med12_oe" / "results.txt"
SPT6_FILE = ROOT / "spt6_oe" / "results.txt"
OUT_DATA = ROOT / "score_relationships_data.csv"
OUT_PAIRWISE = ROOT / "score_relationships_pairwise.pdf"
OUT_3D = ROOT / "score_relationships_3d.pdf"


def read_refseq_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 13:
                # UCSC refGene: field 1 is name (RefSeq transcript), field
                # 12 is name2 (gene symbol). The table has a leading bin.
                tx = fields[1].split(".", 1)[0]
                gene = fields[12].strip()
                if tx and gene:
                    mapping.setdefault(tx, gene)
    return mapping


def load_data() -> tuple[pd.DataFrame, dict[str, int]]:
    ref_map = read_refseq_map(REFSEQ_FILE)

    h2az = pd.read_csv(
        H2AZ_FILE,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "transcript", "bed_score", "strand", "h2az_signal", "h2az_value"],
    )
    h2az["tx_key"] = h2az["transcript"].astype(str).str.split(".").str[0]
    h2az["gene"] = h2az["tx_key"].map(ref_map)
    h2az["h2az_signal"] = pd.to_numeric(h2az["h2az_signal"], errors="coerce")
    h2az = h2az.dropna(subset=["gene", "h2az_signal"])
    # Avoid giving genes with multiple RefSeq isoforms disproportionate weight.
    h2az_gene = h2az.groupby("gene", as_index=False)["h2az_signal"].mean()
    h2az_gene["h2az_zscore"] = (
        h2az_gene["h2az_signal"] - h2az_gene["h2az_signal"].mean()
    ) / h2az_gene["h2az_signal"].std(ddof=1)

    med12 = pd.read_csv(MED12_FILE, index_col=0)
    spt6 = pd.read_csv(SPT6_FILE, index_col=0)
    med12_fc = pd.to_numeric(med12["log2FoldChange"], errors="coerce").rename("med12_log2fc")
    spt6_fc = pd.to_numeric(spt6["log2FoldChange"], errors="coerce").rename("spt6_log2fc")

    data = h2az_gene.set_index("gene").join(med12_fc, how="inner").join(spt6_fc, how="inner")
    data = data.dropna(subset=["h2az_signal", "h2az_zscore", "med12_log2fc", "spt6_log2fc"]).copy()
    data["rna_mean_log2fc"] = (data["med12_log2fc"] + data["spt6_log2fc"]) / 2.0
    data["score"] = data["h2az_zscore"] + data["rna_mean_log2fc"]
    data = data.reset_index().rename(columns={"gene": "gene_symbol"})
    data = data.sort_values("score", ascending=False).reset_index(drop=True)

    counts = {
        "h2az_rows": int(len(h2az)),
        "h2az_mapped_rows": int(h2az["gene"].notna().sum()),
        "h2az_genes": int(len(h2az_gene)),
        "merged_genes": int(len(data)),
        "med12_non_na": int(med12_fc.notna().sum()),
        "spt6_non_na": int(spt6_fc.notna().sum()),
    }
    return data, counts


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    den = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    return float(np.dot(x, y) / den) if den else float("nan")


def rank_values(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    # Average ranks for ties.
    vals, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    for i, n in enumerate(counts):
        if n > 1:
            ranks[inverse == i] = ranks[inverse == i].mean()
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rank_values(np.asarray(x)), rank_values(np.asarray(y)))


def fmt(v: float, digits: int = 2) -> str:
    if not np.isfinite(v):
        return "NA"
    return f"{v:.{digits}f}"


def color_for_score(score: float, limit: float) -> colors.Color:
    # Blue -> near-white -> red, centered on zero.
    t = max(-1.0, min(1.0, float(score) / limit if limit else 0.0))
    if t < 0:
        u = t + 1.0
        r, g, b = 0.16 + 0.82 * u, 0.39 + 0.55 * u, 0.78 + 0.17 * u
    else:
        u = t
        r, g, b = 0.98, 0.95 - 0.62 * u, 0.95 - 0.66 * u
    return colors.Color(r, g, b)


def safe_alpha(c: canvas.Canvas, alpha: float) -> None:
    try:
        c.setFillAlpha(alpha)
        c.setStrokeAlpha(alpha)
    except Exception:
        pass


def draw_text(c: canvas.Canvas, x: float, y: float, text: str, size: float = 8, color=colors.HexColor("#252525"), bold: bool = False) -> None:
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, str(text))


def draw_center(c: canvas.Canvas, x: float, y: float, text: str, size: float = 8, color=colors.HexColor("#252525"), bold: bool = False) -> None:
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawCentredString(x, y, str(text))


def draw_right(c: canvas.Canvas, x: float, y: float, text: str, size: float = 8, color=colors.HexColor("#252525"), bold: bool = False) -> None:
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawRightString(x, y, str(text))


def draw_page_title(c: canvas.Canvas, title: str, subtitle: str, width: float, height: float) -> None:
    draw_text(c, 0.45 * inch, height - 0.38 * inch, title, 15, bold=True)
    draw_text(c, 0.45 * inch, height - 0.60 * inch, subtitle, 8.3, colors.HexColor("#555555"))


def axis_ticks(lo: float, hi: float, n: int = 5) -> np.ndarray:
    return np.linspace(lo, hi, n)


def map_value(v: np.ndarray | float, lo: float, hi: float, a: float, b: float):
    if hi <= lo:
        return np.full_like(np.asarray(v, dtype=float), (a + b) / 2.0)
    return a + (np.asarray(v, dtype=float) - lo) / (hi - lo) * (b - a)


def draw_scatter(c: canvas.Canvas, df: pd.DataFrame, xcol: str, ycol: str, xlabel: str, ylabel: str, title: str, x: float, y: float, w: float, h: float, score_limit: float, color_by: str = "score") -> None:
    left, right = x + 0.55 * inch, x + w - 0.16 * inch
    bottom, top = y + 0.43 * inch, y + h - 0.35 * inch
    xv = df[xcol].to_numpy(float)
    yv = df[ycol].to_numpy(float)
    # Use robust limits for legibility, while clipping points at the frame.
    xlo, xhi = np.nanpercentile(xv, [0.5, 99.5])
    ylo, yhi = np.nanpercentile(yv, [0.5, 99.5])
    if xlo == xhi:
        xlo, xhi = xlo - 1, xhi + 1
    if ylo == yhi:
        ylo, yhi = ylo - 1, yhi + 1
    c.setStrokeColor(colors.HexColor("#B7B7B7"))
    c.setLineWidth(0.65)
    c.rect(left, bottom, right - left, top - bottom, stroke=1, fill=0)
    for tick in axis_ticks(xlo, xhi):
        px = float(map_value(tick, xlo, xhi, left, right))
        c.setStrokeColor(colors.HexColor("#E4E4E4"))
        c.line(px, bottom, px, top)
        draw_center(c, px, bottom - 0.16 * inch, fmt(tick, 1), 7, colors.HexColor("#555555"))
    for tick in axis_ticks(ylo, yhi):
        py = float(map_value(tick, ylo, yhi, bottom, top))
        c.setStrokeColor(colors.HexColor("#E4E4E4"))
        c.line(left, py, right, py)
        c.setFillColor(colors.HexColor("#555555"))
        c.setFont("Helvetica", 7)
        c.drawRightString(left - 0.08 * inch, py - 2, fmt(tick, 1))
    # Reference lines at zero when visible.
    c.setStrokeColor(colors.HexColor("#9B9B9B"))
    c.setLineWidth(0.7)
    if xlo <= 0 <= xhi:
        px = float(map_value(0, xlo, xhi, left, right))
        c.line(px, bottom, px, top)
    if ylo <= 0 <= yhi:
        py = float(map_value(0, ylo, yhi, bottom, top))
        c.line(left, py, right, py)
    px = map_value(xv, xlo, xhi, left, right)
    py = map_value(yv, ylo, yhi, bottom, top)
    px = np.clip(px, left, right)
    py = np.clip(py, bottom, top)
    c.saveState()
    safe_alpha(c, 0.43)
    for i in np.argsort(py):
        col = color_for_score(float(df.iloc[i][color_by]), score_limit)
        c.setFillColor(col)
        c.setStrokeColor(col)
        c.circle(float(px[i]), float(py[i]), 1.55, stroke=0, fill=1)
    c.restoreState()
    draw_center(c, (left + right) / 2, y + 0.12 * inch, xlabel, 8.5, bold=True)
    c.saveState()
    c.translate(x + 0.12 * inch, (bottom + top) / 2)
    c.rotate(90)
    draw_center(c, 0, 0, ylabel, 8.5, bold=True)
    c.restoreState()
    draw_text(c, x + 0.15 * inch, y + h - 0.19 * inch, title, 9.5, bold=True)
    r = pearson(xv, yv)
    rs = spearman(xv, yv)
    draw_right(c, right, y + h - 0.19 * inch, f"r={fmt(r)}; rho={fmt(rs)}", 7.0, colors.HexColor("#555555"))


def draw_histogram(c: canvas.Canvas, values: np.ndarray, xlabel: str, title: str, x: float, y: float, w: float, h: float, score_limit: float) -> None:
    left, right = x + 0.55 * inch, x + w - 0.16 * inch
    bottom, top = y + 0.43 * inch, y + h - 0.35 * inch
    counts, edges = np.histogram(values, bins=32)
    lo, hi = float(edges[0]), float(edges[-1])
    ymax = max(int(counts.max()), 1)
    c.setStrokeColor(colors.HexColor("#B7B7B7"))
    c.rect(left, bottom, right - left, top - bottom, stroke=1, fill=0)
    barw = (right - left) / len(counts)
    for i, count in enumerate(counts):
        bh = (top - bottom) * count / ymax
        mid = (edges[i] + edges[i + 1]) / 2
        c.setFillColor(color_for_score(float(mid), score_limit))
        c.rect(left + i * barw + 0.25, bottom, max(0.1, barw - 0.5), bh, stroke=0, fill=1)
    for tick in axis_ticks(lo, hi):
        px = float(map_value(tick, lo, hi, left, right))
        c.setStrokeColor(colors.HexColor("#E4E4E4"))
        c.line(px, bottom, px, top)
        draw_center(c, px, bottom - 0.16 * inch, fmt(tick, 1), 7, colors.HexColor("#555555"))
    for tick in axis_ticks(0, ymax):
        py = float(map_value(tick, 0, ymax, bottom, top))
        draw_text(c, left - 0.08 * inch - stringWidth(str(int(tick)), "Helvetica", 7), py - 2, str(int(tick)), 7, colors.HexColor("#555555"))
    draw_center(c, (left + right) / 2, y + 0.12 * inch, xlabel, 8.5, bold=True)
    draw_text(c, x + 0.15 * inch, y + h - 0.19 * inch, title, 9.5, bold=True)
    draw_text(c, right - 1.48 * inch, y + h - 0.19 * inch, f"n={len(values):,}; median={fmt(float(np.median(values)))}", 7.0, colors.HexColor("#555555"))


def draw_formula_panel(c: canvas.Canvas, counts: dict[str, int], x: float, y: float, w: float, h: float) -> None:
    draw_text(c, x + 0.15 * inch, y + h - 0.19 * inch, "Score construction", 9.5, bold=True)
    draw_text(c, x + 0.15 * inch, y + h - 0.48 * inch, "SCORE = z(H2AZ signal) + (MED12OE log2FC + SPT6OE log2FC) / 2", 10.2, colors.HexColor("#1C4D78"), bold=True)
    # Three component boxes and a score box.
    box_y = y + 0.54 * inch
    box_w = 1.35 * inch
    box_h = 0.5 * inch
    starts = [x + 0.16 * inch, x + 1.72 * inch, x + 3.28 * inch]
    labels = ["H2AZ signal\n z-score", "MED12OE\n log2FC", "SPT6OE\n log2FC"]
    fills = [colors.HexColor("#DCECF8"), colors.HexColor("#FDE2D0"), colors.HexColor("#E6EFD8")]
    for i, sx in enumerate(starts):
        c.setFillColor(fills[i])
        c.setStrokeColor(colors.HexColor("#9A9A9A"))
        c.roundRect(sx, box_y, box_w, box_h, 5, stroke=1, fill=1)
        a, b = labels[i].split("\n")
        draw_center(c, sx + box_w / 2, box_y + 0.30 * inch, a, 8.5, bold=True)
        draw_center(c, sx + box_w / 2, box_y + 0.12 * inch, b, 8, colors.HexColor("#555555"))
    score_x = x + w - 1.55 * inch
    c.setFillColor(colors.HexColor("#F5E6F5"))
    c.setStrokeColor(colors.HexColor("#8C6D8C"))
    c.roundRect(score_x, box_y, 1.32 * inch, box_h, 5, stroke=1, fill=1)
    draw_center(c, score_x + 0.66 * inch, box_y + 0.30 * inch, "COMPOSITE", 8.5, bold=True)
    draw_center(c, score_x + 0.66 * inch, box_y + 0.12 * inch, "SCORE", 8, colors.HexColor("#555555"))
    c.setStrokeColor(colors.HexColor("#757575"))
    c.setLineWidth(1)
    for sx in starts:
        c.line(sx + box_w, box_y + box_h / 2, score_x - 0.11 * inch, box_y + box_h / 2)
        c.line(score_x - 0.11 * inch, box_y + box_h / 2, score_x - 0.18 * inch, box_y + box_h / 2 + 0.05 * inch)
        c.line(score_x - 0.11 * inch, box_y + box_h / 2, score_x - 0.18 * inch, box_y + box_h / 2 - 0.05 * inch)
    draw_text(c, x + 0.16 * inch, y + 0.21 * inch, f"Analysis unit: gene; H2AZ isoforms averaged per gene; merged genes: {counts['merged_genes']:,}", 7.3, colors.HexColor("#555555"))


def draw_score_legend(c: canvas.Canvas, x: float, y: float, score_limit: float) -> None:
    draw_text(c, x, y + 0.12 * inch, "SCORE", 7.5, bold=True)
    n = 80
    for i in range(n):
        v = -score_limit + (2 * score_limit) * (i + 0.5) / n
        c.setFillColor(color_for_score(v, score_limit))
        c.rect(x + 0.44 * inch + i * 1.1, y + 0.13 * inch, 1.2, 0.11 * inch, stroke=0, fill=1)
    draw_text(c, x + 0.44 * inch, y - 0.02 * inch, f"low {fmt(-score_limit, 1)}", 6.7, colors.HexColor("#555555"))
    draw_center(c, x + 0.44 * inch + n * 0.55, y - 0.02 * inch, "0", 6.7, colors.HexColor("#555555"))
    draw_text(c, x + 0.44 * inch + n * 1.1 - 26, y - 0.02 * inch, f"{fmt(score_limit, 1)} high", 6.7, colors.HexColor("#555555"))


def draw_pairwise_pdf(df: pd.DataFrame, counts: dict[str, int]) -> None:
    width, height = 11 * inch, 8.5 * inch
    c = canvas.Canvas(str(OUT_PAIRWISE), pagesize=(width, height))
    c.setTitle("H2AZ MED12OE SPT6OE composite score relationships")
    draw_page_title(c, "Composite SCORE: calculation and pairwise relationships", "Gene-level data; z-score of mean H2AZ signal + mean(MED12OE, SPT6OE) log2 fold-change", width, height)
    score_limit = max(float(np.max(np.abs(df["score"]))), 1.0)
    # Top formula / workflow strip.
    draw_formula_panel(c, counts, 0.45 * inch, 5.30 * inch, 10.05 * inch, 1.65 * inch)
    # Three plots below.
    draw_scatter(c, df, "h2az_zscore", "rna_mean_log2fc", "H2AZ signal z-score", "Mean RNA-seq log2FC", "A. H2AZ vs mean RNA response", 0.45 * inch, 2.55 * inch, 3.28 * inch, 2.55 * inch, score_limit)
    draw_scatter(c, df, "med12_log2fc", "spt6_log2fc", "MED12OE log2FC", "SPT6OE log2FC", "B. MED12OE vs SPT6OE", 3.86 * inch, 2.55 * inch, 3.28 * inch, 2.55 * inch, score_limit)
    draw_histogram(c, df["score"].to_numpy(float), "SCORE", "C. SCORE distribution", 7.27 * inch, 2.55 * inch, 3.23 * inch, 2.55 * inch, score_limit)
    draw_score_legend(c, 0.6 * inch, 1.87 * inch, score_limit)
    # Footer with representative extremes and exact counts.
    high = df.nlargest(5, "score")["gene_symbol"].tolist()
    low = df.nsmallest(5, "score")["gene_symbol"].tolist()
    draw_text(c, 0.6 * inch, 1.48 * inch, "Highest SCORE genes: " + ", ".join(high), 7.4, colors.HexColor("#444444"))
    draw_text(c, 0.6 * inch, 1.24 * inch, "Lowest SCORE genes: " + ", ".join(low), 7.4, colors.HexColor("#444444"))
    draw_text(c, 0.6 * inch, 0.86 * inch, f"Input H2AZ rows: {counts['h2az_rows']:,}; mapped rows: {counts['h2az_mapped_rows']:,}; unique H2AZ genes: {counts['h2az_genes']:,}; complete merged genes: {counts['merged_genes']:,}", 7.2, colors.HexColor("#555555"))
    c.showPage()
    c.save()


def axis_point(origin: tuple[float, float], vx: tuple[float, float], vy: tuple[float, float], vz: tuple[float, float], nx: float, ny: float, nz: float) -> tuple[float, float]:
    return (origin[0] + nx * vx[0] + ny * vy[0] + nz * vz[0], origin[1] + nx * vx[1] + ny * vy[1] + nz * vz[1])


def draw_3d_pdf(df: pd.DataFrame) -> None:
    width, height = 11 * inch, 8.5 * inch
    c = canvas.Canvas(str(OUT_3D), pagesize=(width, height))
    c.setTitle("3D projection of H2AZ MED12OE SPT6OE components")
    draw_page_title(c, "3D projection of the three score components", "Each point is a gene; point color is the derived SCORE. Axes are independently rescaled for the perspective view.", width, height)
    origin = (2.0 * inch, 1.14 * inch)
    vx = (4.55 * inch, 0.0)
    vy = (1.62 * inch, 1.02 * inch)
    vz = (0.0, 4.35 * inch)
    # Normalize each component to its displayed axis range.
    cols = ["med12_log2fc", "spt6_log2fc", "h2az_zscore"]
    norms = {}
    for col in cols:
        lo, hi = np.nanpercentile(df[col].to_numpy(float), [0.5, 99.5])
        if hi <= lo:
            lo, hi = lo - 1, hi + 1
        norms[col] = (float(lo), float(hi))
    # Light floor and side-plane guides.
    c.setStrokeColor(colors.HexColor("#D8D8D8"))
    c.setLineWidth(0.55)
    for frac in np.linspace(0, 1, 5):
        p1 = axis_point(origin, vx, vy, vz, frac, 0, 0)
        p2 = axis_point(origin, vx, vy, vz, frac, 1, 0)
        c.line(*p1, *p2)
        p3 = axis_point(origin, vx, vy, vz, 0, frac, 0)
        p4 = axis_point(origin, vx, vy, vz, 1, frac, 0)
        c.line(*p3, *p4)
        p5 = axis_point(origin, vx, vy, vz, 0, 0, frac)
        p6 = axis_point(origin, vx, vy, vz, 1, 0, frac)
        c.line(*p5, *p6)
        p7 = axis_point(origin, vx, vy, vz, 0, 1, frac)
        p8 = axis_point(origin, vx, vy, vz, 1, 1, frac)
        c.line(*p7, *p8)
    # Axes.
    c.setStrokeColor(colors.HexColor("#555555"))
    c.setLineWidth(1.1)
    for vec in [vx, vy, vz]:
        c.line(origin[0], origin[1], origin[0] + vec[0], origin[1] + vec[1])
    axis_specs = [(vx, "MED12OE log2FC", norms["med12_log2fc"]), (vy, "SPT6OE log2FC", norms["spt6_log2fc"]), (vz, "H2AZ z-score", norms["h2az_zscore"])]
    for vec, label, (lo, hi) in axis_specs:
        for frac in np.linspace(0, 1, 5):
            px = origin[0] + frac * vec[0]
            py = origin[1] + frac * vec[1]
            # A small tick approximately normal to the axis.
            if abs(vec[0]) > abs(vec[1]):
                c.line(px, py - 3, px, py + 3)
            else:
                c.line(px - 3, py, px + 3, py)
            draw_center(c, px + (5 if vec is vy else 0), py - (13 if vec is not vz else 4), fmt(lo + frac * (hi - lo), 1), 6.6, colors.HexColor("#555555"))
        lx = origin[0] + vec[0] * 0.63
        ly = origin[1] + vec[1] * 0.63
        if vec is vz:
            draw_text(c, lx + 8, ly, label, 9.5, bold=True)
        else:
            draw_center(c, lx, ly - 0.28 * inch, label, 9.5, bold=True)
    score_limit = max(float(np.max(np.abs(df["score"]))), 1.0)
    # Draw far points first so foreground points remain visible.
    records = []
    for row in df.itertuples(index=False):
        nx = (row.med12_log2fc - norms["med12_log2fc"][0]) / (norms["med12_log2fc"][1] - norms["med12_log2fc"][0])
        ny = (row.spt6_log2fc - norms["spt6_log2fc"][0]) / (norms["spt6_log2fc"][1] - norms["spt6_log2fc"][0])
        nz = (row.h2az_zscore - norms["h2az_zscore"][0]) / (norms["h2az_zscore"][1] - norms["h2az_zscore"][0])
        nx, ny, nz = [max(0.0, min(1.0, float(v))) for v in (nx, ny, nz)]
        px, py = axis_point(origin, vx, vy, vz, nx, ny, nz)
        records.append((py, px, row.score))
    c.saveState()
    safe_alpha(c, 0.58)
    for _, px, score in sorted(records):
        col = color_for_score(float(score), score_limit)
        c.setFillColor(col)
        c.setStrokeColor(col)
        c.circle(px, _, 1.45, stroke=0, fill=1)
    c.restoreState()
    # Right-side explanatory legend.
    lx, ly = 8.55 * inch, 5.65 * inch
    draw_text(c, lx, ly, "How to read the projection", 10.5, bold=True)
    draw_text(c, lx, ly - 0.28 * inch, "x-axis: MED12OE response", 8.5)
    draw_text(c, lx, ly - 0.48 * inch, "diagonal axis: SPT6OE response", 8.5)
    draw_text(c, lx, ly - 0.68 * inch, "vertical axis: H2AZ signal z-score", 8.5)
    draw_text(c, lx, ly - 0.95 * inch, "The SCORE is the additive combination", 8.5, colors.HexColor("#555555"))
    draw_text(c, lx, ly - 1.13 * inch, "of the vertical component and the", 8.5, colors.HexColor("#555555"))
    draw_text(c, lx, ly - 1.31 * inch, "mean of the two RNA-seq components.", 8.5, colors.HexColor("#555555"))
    draw_score_legend(c, lx, ly - 1.82 * inch, score_limit)
    draw_text(c, lx, ly - 2.32 * inch, f"n = {len(df):,} complete genes", 8.5, colors.HexColor("#555555"))
    draw_text(c, lx, ly - 2.52 * inch, "Note: independently rescaled axes", 7.7, colors.HexColor("#777777"))
    draw_text(c, lx, ly - 2.69 * inch, "improve 3D readability; score values", 7.7, colors.HexColor("#777777"))
    draw_text(c, lx, ly - 2.86 * inch, "are calculated from unscaled inputs.", 7.7, colors.HexColor("#777777"))
    c.showPage()
    c.save()


def main() -> None:
    data, counts = load_data()
    data.to_csv(OUT_DATA, index=False)
    draw_pairwise_pdf(data, counts)
    draw_3d_pdf(data)
    print(f"wrote {OUT_DATA}")
    print(f"wrote {OUT_PAIRWISE}")
    print(f"wrote {OUT_3D}")
    print(f"complete genes: {len(data)}")
    print(f"score range: {data['score'].min():.3f} to {data['score'].max():.3f}")
    print(f"top genes: {', '.join(data.head(5)['gene_symbol'])}")
    print(f"bottom genes: {', '.join(data.tail(5)['gene_symbol'])}")


if __name__ == "__main__":
    main()
