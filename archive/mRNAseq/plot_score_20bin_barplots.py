#!/usr/bin/env python3
"""Create two aligned barplots summarizing 20 SCORE-ranked gene bins."""

from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(r"D:\work\pby\mRNAseq")
DATA = ROOT / "score_relationships_data.csv"
OUTPUT_SCORE = ROOT / "score_20bins_mean_score_optimized.pdf"
OUTPUT_COMPONENTS = ROOT / "score_20bins_three_components_optimized.pdf"

NAVY = colors.HexColor("#27364A")
SLATE = colors.HexColor("#596779")
GRID = colors.HexColor("#E5E9EE")
BLUE = colors.HexColor("#2F69B3")
GREEN = colors.HexColor("#27834D")
ORANGE = colors.HexColor("#E79618")
LOW = colors.HexColor("#3D5A91")
MID = colors.HexColor("#9AA3AC")
HIGH = colors.HexColor("#D96352")
LOW_BAND = colors.HexColor("#EEF3F9")
HIGH_BAND = colors.HexColor("#FCF1EC")


def text(c, x, y, value, size=8, col=NAVY, bold=False):
    c.setFillColor(col)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawString(x, y, value)


def center(c, x, y, value, size=8, col=NAVY, bold=False):
    c.setFillColor(col)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawCentredString(x, y, value)


def right(c, x, y, value, size=8, col=NAVY, bold=False):
    c.setFillColor(col)
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.drawRightString(x, y, value)


def lerp(a, b, t):
    return colors.Color(a.red + (b.red - a.red) * t, a.green + (b.green - a.green) * t, a.blue + (b.blue - a.blue) * t)


def score_color(value, limit):
    t = max(-1.0, min(1.0, float(value) / max(limit, 1e-9)))
    return lerp(LOW, MID, t + 1.0) if t < 0 else lerp(MID, HIGH, t)


def nice_ticks(lo, hi, count=6):
    raw = max((hi - lo) / max(count - 1, 1), 1e-9)
    power = 10 ** np.floor(np.log10(raw))
    scaled = raw / power
    step = (1 if scaled <= 1 else 2 if scaled <= 2 else 2.5 if scaled <= 2.5 else 5 if scaled <= 5 else 10) * power
    start = np.floor(lo / step) * step
    end = np.ceil(hi / step) * step
    return np.arange(start, end + step * 0.5, step)


def y_map(value, lo, hi, y0, height):
    return y0 + (value - lo) / (hi - lo) * height


def draw_axes(c, x0, y0, width, height, ylo, yhi, ylabel, ticks):
    c.setStrokeColor(GRID)
    c.setLineWidth(0.45)
    for tick in ticks:
        if ylo <= tick <= yhi:
            yy = y_map(tick, ylo, yhi, y0, height)
            c.line(x0, yy, x0 + width, yy)
            right(c, x0 - 7, yy - 2.5, f"{tick:g}", 7.0, SLATE)
    c.setStrokeColor(NAVY)
    c.setLineWidth(0.9)
    c.line(x0, y0, x0, y0 + height)
    zero_y = y_map(0, ylo, yhi, y0, height)
    if y0 <= zero_y <= y0 + height:
        c.setLineWidth(1.1)
        c.line(x0, zero_y, x0 + width, zero_y)
    c.saveState()
    c.translate(x0 - 42, y0 + height / 2)
    c.rotate(90)
    center(c, 0, 0, ylabel, 8.2, NAVY, True)
    c.restoreState()


def build_summary():
    df = pd.read_csv(DATA).dropna(subset=["score", "h2az_zscore", "med12_log2fc", "spt6_log2fc"]).copy()
    df = df.sort_values("score", ascending=True).reset_index(drop=True)
    df["score_bin"] = np.floor(np.arange(len(df)) * 20 / len(df)).astype(int) + 1
    summary = (
        df.groupby("score_bin", as_index=False)
        .agg(
            n=("score", "size"),
            mean_score=("score", "mean"),
            mean_h2az=("h2az_zscore", "mean"),
            mean_med12=("med12_log2fc", "mean"),
            mean_spt6=("spt6_log2fc", "mean"),
        )
    )
    return df, summary


def draw_score_panel(c, summary, x0, y0, width, height):
    values = summary["mean_score"].to_numpy(float)
    ylo = min(float(values.min()) * 1.14, -0.3)
    yhi = max(float(values.max()) * 1.14, 0.3)
    ticks = nice_ticks(ylo, yhi, 6)
    ylo, yhi = float(ticks.min()), float(ticks.max())
    draw_axes(c, x0, y0, width, height, ylo, yhi, "Mean SCORE", ticks)
    step = width / 20.0
    bar_w = step * 0.84
    limit = max(abs(values.min()), abs(values.max()))
    zero_y = y_map(0, ylo, yhi, y0, height)
    for idx, row in summary.iterrows():
        x = x0 + idx * step + (step - bar_w) / 2
        yy = y_map(row.mean_score, ylo, yhi, y0, height)
        c.setFillColor(score_color(row.mean_score, limit))
        c.rect(x, min(zero_y, yy), bar_w, abs(yy - zero_y), stroke=0, fill=1)
        center(c, x + bar_w / 2, y0 - 13, str(int(row.score_bin)), 6.8, SLATE)
    center(c, x0 + width / 2, y0 - 28, "SCORE bin (low to high)", 8.0, NAVY, True)


def draw_component_panel(c, summary, x0, y0, width, height):
    h2az_values = summary["mean_h2az"].to_numpy(float)
    log2fc_values = summary[["mean_med12", "mean_spt6"]].to_numpy(float)

    # Use symmetric, independent scales so both axes share the same visual zero.
    left_limit = max(0.5, np.ceil(np.max(np.abs(h2az_values)) * 2) / 2)
    right_limit = max(1.0, np.ceil(np.max(np.abs(log2fc_values))))
    left_ticks = np.linspace(-left_limit, left_limit, 7)
    right_ticks = np.linspace(-right_limit, right_limit, 7)

    c.setStrokeColor(GRID)
    c.setLineWidth(0.45)
    for tick in left_ticks:
        yy = y_map(tick, -left_limit, left_limit, y0, height)
        c.line(x0, yy, x0 + width, yy)
        right(c, x0 - 7, yy - 2.5, f"{tick:g}", 7.0, GREEN)
    for tick in right_ticks:
        yy = y_map(tick, -right_limit, right_limit, y0, height)
        text(c, x0 + width + 7, yy - 2.5, f"{tick:g}", 7.0, SLATE)

    c.setStrokeColor(NAVY)
    c.setLineWidth(0.9)
    c.line(x0, y0, x0, y0 + height)
    c.line(x0 + width, y0, x0 + width, y0 + height)
    zero_y = y_map(0, -left_limit, left_limit, y0, height)
    c.setLineWidth(1.1)
    c.line(x0, zero_y, x0 + width, zero_y)

    c.saveState()
    c.translate(x0 - 47, y0 + height / 2)
    c.rotate(90)
    center(c, 0, 0, "Mean H2AZ z-score", 8.2, GREEN, True)
    c.restoreState()
    c.saveState()
    c.translate(x0 + width + 49, y0 + height / 2)
    c.rotate(-90)
    center(c, 0, 0, "Mean log2FC", 8.2, NAVY, True)
    c.restoreState()

    step = width / 20.0
    series = [
        ("mean_h2az", GREEN, -left_limit, left_limit),
        ("mean_med12", BLUE, -right_limit, right_limit),
        ("mean_spt6", ORANGE, -right_limit, right_limit),
    ]
    group_w = step * 0.92
    bar_w = group_w / 3.0
    for idx, row in summary.iterrows():
        group_x = x0 + idx * step + (step - group_w) / 2
        for offset, (col, color, ylo, yhi) in enumerate(series):
            value = float(row[col])
            yy = y_map(value, ylo, yhi, y0, height)
            c.setFillColor(color)
            c.rect(group_x + offset * bar_w, min(zero_y, yy), bar_w - 0.55, abs(yy - zero_y), stroke=0, fill=1)
        center(c, x0 + idx * step + step / 2, y0 - 13, str(int(row.score_bin)), 6.8, SLATE)
    center(c, x0 + width / 2, y0 - 28, "SCORE bin (low to high)", 8.0, NAVY, True)

    # Compact vector legend.
    legend_y = y0 + height + 8
    legend_x = x0 + width - 320
    for label, color, dx in (("H2AZ z-score (left)", GREEN, 0), ("MED12OE (right)", BLUE, 125), ("SUPT6H (right)", ORANGE, 235)):
        c.setFillColor(color)
        c.rect(legend_x + dx, legend_y, 9, 9, stroke=0, fill=1)
        text(c, legend_x + dx + 14, legend_y + 1, label, 7.4, SLATE)


def main():
    df, summary = build_summary()
    page_w, page_h = 11 * inch, 6.1 * inch
    plot_x = 0.82 * inch
    plot_w = 9.65 * inch
    min_n, max_n = int(summary.n.min()), int(summary.n.max())

    c = canvas.Canvas(str(OUTPUT_SCORE), pagesize=(page_w, page_h))
    c.setTitle("Mean SCORE across 20 ranked gene bins")
    text(c, 0.55 * inch, page_h - 0.42 * inch, "SCORE across 20 ranked gene bins", 15.5, NAVY, True)
    text(c, 0.55 * inch, page_h - 0.67 * inch, "SCORE = z(H2AZ) + [MED12OE log2FC + SUPT6H log2FC] / 2", 8.6, SLATE)
    draw_score_panel(c, summary, plot_x, 0.95 * inch, plot_w, 3.90 * inch)
    text(c, 0.55 * inch, 0.30 * inch, f"20 equal-size bins; n = {len(df):,} genes.", 7.3, SLATE)
    right(c, page_w - 0.55 * inch, 0.30 * inch, "Low SCORE -> high SCORE; suitable for ranked GSEA.", 7.3, SLATE)
    c.showPage()
    c.save()

    c = canvas.Canvas(str(OUTPUT_COMPONENTS), pagesize=(page_w, page_h))
    c.setTitle("Three SCORE components across 20 ranked gene bins")
    text(c, 0.55 * inch, page_h - 0.42 * inch, "Three components across the same SCORE bins", 15.5, NAVY, True)
    text(c, 0.55 * inch, page_h - 0.67 * inch, "SCORE = z(H2AZ) + [MED12OE log2FC + SUPT6H log2FC] / 2", 8.6, SLATE)
    draw_component_panel(c, summary, 0.92 * inch, 0.95 * inch, 9.05 * inch, 3.90 * inch)
    text(c, 0.55 * inch, 0.30 * inch, f"20 equal-size bins; n = {len(df):,} genes.", 7.3, SLATE)
    right(c, page_w - 0.55 * inch, 0.30 * inch, "Concordant low/high components support SCORE-ranked GSEA.", 7.3, SLATE)
    c.showPage()
    c.save()

    print(f"wrote {OUTPUT_SCORE}")
    print(f"wrote {OUTPUT_COMPONENTS}")


if __name__ == "__main__":
    main()
