from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CDK8_FILE = ROOT.parent / "KO_MED12_CDK8_APEX2" / "volcannodata.csv"
RNAPII_FILE = ROOT.parent / "KO_MED12_PRB1_APEX2" / "volcannodata.csv"
OUTPUT_PDF = ROOT / "med_core_cdk8_vs_rnapii_heatmap.pdf"
OUTPUT_TABLE = ROOT / "med_core_cdk8_vs_rnapii_log2fc.csv"

MED_CORE = ["MED1", "MED4", "MED6", "MED14", "MED15", "MED17", "MED19", "MED24"]
MIN_MAX_ABS_LOG2FC = 0.5

# Arial TrueType text remains editable in Illustrator.
ARIAL = Path(r"C:\Windows\Fonts\arial.ttf")
ARIAL_BOLD = Path(r"C:\Windows\Fonts\arialbd.ttf")
font_manager.fontManager.addfont(str(ARIAL))
font_manager.fontManager.addfont(str(ARIAL_BOLD))
plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def read_dataset(path: Path, prefix: str) -> pd.DataFrame:
    data = pd.read_csv(path)
    data["name"] = data["name"].astype(str).str.upper().str.strip()
    data["padj"] = pd.to_numeric(data["padj"], errors="coerce")
    data["log2FoldChange"] = pd.to_numeric(data["log2FoldChange"], errors="coerce")
    data = data.dropna(subset=["name"]).drop_duplicates("name", keep="first")
    data = data.set_index("name")
    out = data.reindex(MED_CORE)[["padj", "log2FoldChange"]].copy()
    out.columns = [f"{prefix}_{column}" for column in out.columns]
    return out


def main() -> None:
    cdk8 = read_dataset(CDK8_FILE, "CDK8")
    rnapii = read_dataset(RNAPII_FILE, "RNAPII")
    table = cdk8.join(rnapii)
    table["max_abs_log2FC"] = table[["CDK8_log2FoldChange", "RNAPII_log2FoldChange"]].abs().max(axis=1)
    table = table.loc[table["max_abs_log2FC"] >= MIN_MAX_ABS_LOG2FC].copy()
    genes = table.index.tolist()
    table["CDK8_significant_down"] = (
        (table["CDK8_padj"] < 0.05) & (table["CDK8_log2FoldChange"] < 0)
    )
    table["RNAPII_significant_down"] = (
        (table["RNAPII_padj"] < 0.05) & (table["RNAPII_log2FoldChange"] < 0)
    )
    table.to_csv(OUTPUT_TABLE, index_label="gene")

    values = table[["CDK8_log2FoldChange", "RNAPII_log2FoldChange"]].to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    limit = max(1.0, float(np.ceil(np.max(np.abs(finite_values)) * 10) / 10))

    fig, ax = plt.subplots(figsize=(5.3, 5.2), constrained_layout=True)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#D9D9D9")
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    # Opaque vector rectangles avoid PDF image masks in Illustrator.
    for row in range(len(genes)):
        for col in range(2):
            value = values[row, col]
            color = "#D9D9D9" if not np.isfinite(value) else cmap(norm(value))
            ax.add_patch(Rectangle((col - 0.5, row - 0.5), 1, 1,
                                   facecolor=color, edgecolor="white", linewidth=1.2))
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    image = ScalarMappable(norm=norm, cmap=cmap)
    ax.set_xticks([0, 1], ["CDK8-APEX2", "RNAPII-APEX2"])
    ax.set_yticks(np.arange(len(genes)), genes)
    ax.tick_params(axis="both", length=0, labelsize=10)
    ax.set_xlabel("* padj < 0.05; NA = not detected", fontsize=8, labelpad=18, color="#555555")
    ax.set_title("MED-CORE change after MED12-KD", pad=12, fontsize=13, fontweight="bold")

    for row, gene in enumerate(genes):
        for col, dataset in enumerate(["CDK8", "RNAPII"]):
            value = table.loc[gene, f"{dataset}_log2FoldChange"]
            padj = table.loc[gene, f"{dataset}_padj"]
            if pd.isna(value):
                text = "NA"
                color = "#555555"
            else:
                star = "*" if pd.notna(padj) and padj < 0.05 else ""
                text = f"{value:.2f}{star}"
                color = "white" if abs(value) >= limit * 0.45 else "black"
            ax.text(col, row, text, ha="center", va="center", color=color, fontsize=10)

    for edge in range(len(genes) + 1):
        ax.axhline(edge - 0.5, color="white", linewidth=1.2)
    for edge in range(3):
        ax.axvline(edge - 0.5, color="white", linewidth=1.2)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.04)
    colorbar.set_label("log2 fold change (+IAA vs -IAA)", rotation=90, labelpad=12)
    colorbar.ax.tick_params(labelsize=9)
    # Keep the color scale vector-based as well.
    colorbar.solids.set_rasterized(False)
    fig.savefig(OUTPUT_PDF, format="pdf")
    plt.close(fig)

    print(f"Wrote {OUTPUT_PDF}")
    print(f"Wrote {OUTPUT_TABLE}")


if __name__ == "__main__":
    main()
