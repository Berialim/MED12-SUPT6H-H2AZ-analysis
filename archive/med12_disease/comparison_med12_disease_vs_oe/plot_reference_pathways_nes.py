"""Plot the pathway set supplied by the user for MED12 disease vs OE GSEA."""

from pathlib import Path
import math

import pandas as pd
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(r"D:\work\pby")
OUT = ROOT / "med12_disease" / "comparison_med12_disease_vs_oe"
INPUTS = {
    "Disease": ROOT / "med12_disease",
    "OE": ROOT / "mRNAseq" / "med12_oe",
}

# [K] and [H] indicate the KEGG and Hallmark collections in the reference image.
# The first two items are metabolism pathways; the remaining items are cancer-related.
REFERENCE = [
    ("KEGG", "Biosynthesis of unsaturated fatty acids", "Biosynthesis of unsaturated fatty acids [K]", "Metabolism"),
    ("KEGG", "Oxidative phosphorylation", "Oxidative phosphorylation [K]", "Metabolism"),
    ("HALLMARK", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT [H]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_P53_PATHWAY", "p53 [H]", "Cancer-related"),
    ("KEGG", "Wnt signaling pathway", "Wnt signaling pathway [K]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "WNT / beta-catenin signaling [H]", "Cancer-related"),
    ("KEGG", "TGF-beta signaling pathway", "TGF-beta signaling pathway [K]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_TGF_BETA_SIGNALING", "TGF-beta signaling [H]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_KRAS_SIGNALING_DN", "KRAS signaling DN [H]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "PI3K-AKT-mTOR signaling [H]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNFA / NF-kB signaling [H]", "Cancer-related"),
    ("KEGG", "Pathways in cancer", "Pathways in cancer [K]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_MYC_TARGETS_V1", "MYC targets V1 [H]", "Cancer-related"),
    ("HALLMARK", "HALLMARK_MYC_TARGETS_V2", "MYC targets V2 [H]", "Cancer-related"),
    ("KEGG", "Transcriptional misregulation in cancer", "Transcriptional misregulation in cancer [K]", "Cancer-related"),
    ("KEGG", "MAPK signaling pathway", "MAPK signaling pathway [K]", "Cancer-related"),
    ("KEGG", "Ras signaling pathway", "RAS signaling pathway [K]", "Cancer-related"),
    ("KEGG", "Melanoma", "Melanoma [K]", "Cancer-related"),
]

BLUE = HexColor("#2878B5")
YELLOW = HexColor("#E2B51B")
INK = HexColor("#26364A")
GREY = HexColor("#6B7788")
GRID = HexColor("#D8DEE7")
WHITE = HexColor("#FFFFFF")


def register_fonts():
    fonts = Path(r"C:\Windows\Fonts")
    pdfmetrics.registerFont(TTFont("Arial", str(fonts / "arial.ttf")))
    pdfmetrics.registerFont(TTFont("Arial-Bold", str(fonts / "arialbd.ttf")))


def load_report(base: Path, collection: str) -> pd.DataFrame:
    report = base / f"gsea_results_{collection}" / "gseapy.gene_set.prerank.report.csv"
    df = pd.read_csv(report)
    return df.set_index("Term")[["NES", "FDR q-val"]]


def assemble_data() -> pd.DataFrame:
    disease_reports = {name: load_report(INPUTS["Disease"], name) for name in ("KEGG", "HALLMARK")}
    oe_reports = {name: load_report(INPUTS["OE"], name) for name in ("KEGG", "HALLMARK")}
    rows = []
    for collection, term, display, category in REFERENCE:
        disease = disease_reports[collection].loc[term]
        oe = oe_reports[collection].loc[term]
        rows.append({
            "collection": collection,
            "pathway": term,
            "display_name": display,
            "category": category,
            "disease_NES": float(disease["NES"]),
            "disease_FDR_q": float(disease["FDR q-val"]),
            "oe_NES": float(oe["NES"]),
            "oe_FDR_q": float(oe["FDR q-val"]),
        })
    return pd.DataFrame(rows)


def draw_bar(c, x0, x1, y, height, color, hatched=False):
    x, width = min(x0, x1), abs(x1 - x0)
    c.saveState()
    if hatched:
        c.setFillColor(WHITE)
        c.setStrokeColor(color)
        c.setLineWidth(0.75)
        c.rect(x, y, width, height, fill=1, stroke=1)
        c.setStrokeColor(color)
        c.setLineWidth(0.65)
        path = c.beginPath()
        path.rect(x, y, width, height)
        c.clipPath(path, stroke=0, fill=0)
        step = 5
        start = x - height
        while start < x + width:
            c.line(start, y, start + height, y + height)
            start += step
    else:
        c.setFillColor(color)
        c.setStrokeColor(color)
        c.rect(x, y, width, height, fill=1, stroke=0)
    c.restoreState()


def make_plot(data: pd.DataFrame, output: Path):
    register_fonts()
    page_w, page_h = 950, 820
    c = canvas.Canvas(str(output), pagesize=(page_w, page_h))
    c.setTitle("Reference pathways NES: MED12 disease vs OE")
    c.setFillColor(INK)
    c.setFont("Arial-Bold", 19)
    c.drawString(42, page_h - 45, "Reference pathways: NES in MED12 disease and OE")
    c.setFillColor(GREY)
    c.setFont("Arial", 9)
    c.drawString(42, page_h - 66, "All pathways from the supplied list; sorted as listed. [K] KEGG, [H] Hallmark.")
    c.setStrokeColor(GRID)
    c.line(42, page_h - 80, page_w - 42, page_h - 80)

    left, right = 382, page_w - 52
    plot_w = right - left
    low = min(data["disease_NES"].min(), data["oe_NES"].min()) - 0.28
    high = max(data["disease_NES"].max(), data["oe_NES"].max()) + 0.28
    x_pos = lambda value: left + (value - low) / (high - low) * plot_w
    zero = x_pos(0)
    top, row_gap, bar_h = page_h - 112, 29, 8
    c.setFillColor(GREY)
    c.setFont("Arial-Bold", 8)
    c.drawString(left, page_h - 92, "NES")

    tick_step = 0.5
    tick = math.ceil(low / tick_step) * tick_step
    while tick <= high + 1e-8:
        x = x_pos(tick)
        c.setStrokeColor(INK if abs(tick) < 1e-8 else GRID)
        c.setLineWidth(0.9 if abs(tick) < 1e-8 else 0.45)
        c.line(x, top - len(data) * row_gap - 1, x, top + 10)
        c.setFillColor(GREY)
        c.setFont("Arial", 7.5)
        c.drawCentredString(x, top + 15, f"{tick:g}")
        tick += tick_step

    for i, row in enumerate(data.itertuples(index=False)):
        center = top - i * row_gap
        color = BLUE if row.category == "Cancer-related" else YELLOW
        c.setFillColor(INK)
        c.setFont("Arial", 7.8)
        c.drawString(42, center - 2, row.display_name)
        for value, y, hatched in (
            (row.disease_NES, center + 2, False),
            (row.oe_NES, center - 10, True),
        ):
            x1 = x_pos(value)
            draw_bar(c, zero, x1, y, bar_h, color, hatched)
            c.setFillColor(INK)
            c.setFont("Arial", 7.2)
            label_x = x1 + 4 if value >= 0 else x1 - 4
            label = f"{value:.2f}"
            if value >= 0:
                c.drawString(label_x, y + 1, label)
            else:
                c.drawRightString(label_x, y + 1, label)

    legend_y = 65
    legend_items = [
        (52, BLUE, False, "Cancer-related"),
        (190, YELLOW, False, "Metabolism"),
        (340, BLUE, False, "Disease (solid)"),
        (488, BLUE, True, "OE (hatched)"),
    ]
    for x, color, hatched, label in legend_items:
        draw_bar(c, x, x + 18, legend_y, 10, color, hatched)
        c.setFillColor(INK)
        c.setFont("Arial", 8)
        c.drawString(x + 24, legend_y + 1, label)

    c.setFillColor(GREY)
    c.setFont("Arial", 7)
    c.drawCentredString(
        page_w / 2,
        39,
        "Color = pathway class (cancer-related blue; metabolism yellow); solid = disease; hatched = OE.",
    )
    c.drawCentredString(
        page_w / 2,
        27,
        "NES values are from each original report; cross-run differences are descriptive. FDR q-values are in the CSV.",
    )
    c.save()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = assemble_data()
    csv_path = OUT / "reference_pathways_NES_disease_vs_OE.csv"
    pdf_path = OUT / "reference_pathways_NES_disease_vs_OE.pdf"
    data.to_csv(csv_path, index=False, float_format="%.8g")
    make_plot(data, pdf_path)
    print(f"Pathways: {len(data)} (cancer {sum(data.category == 'Cancer-related')}, metabolism {sum(data.category == 'Metabolism')})")
    print(f"PDF: {pdf_path}")
    print(f"CSV: {csv_path}")


if __name__ == "__main__":
    main()
