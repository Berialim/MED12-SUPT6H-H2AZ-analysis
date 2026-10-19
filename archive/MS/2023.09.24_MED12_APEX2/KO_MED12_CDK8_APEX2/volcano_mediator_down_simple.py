from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "volcannodata.csv"
OUTPUT = ROOT / "volcano_significant_down_mediator_core_supt6h.pdf"

# Set Arial before plotting; TrueType text remains editable in Illustrator.
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


# Same MED-CORE definition used in the preceding GSEA script.
MEDIATOR_CORE = {"MED1", "MED4", "MED6", "MED14", "MED15", "MED17", "MED19", "MED24"}


def main():
    diff = pd.read_csv(INPUT)
    diff["name"] = diff["name"].astype(str).str.upper().str.strip()
    diff["padj"] = pd.to_numeric(diff["padj"], errors="coerce")
    diff["log2FoldChange"] = pd.to_numeric(diff["log2FoldChange"], errors="coerce")
    diff = diff.dropna(subset=["name", "padj", "log2FoldChange"]).copy()
    diff["minusLog10Padj"] = -np.log10(diff["padj"].clip(lower=np.finfo(float).tiny))

    # Previous differential-analysis rule: padj < 0.05 and log2FC < 0.
    significant_down = (diff["padj"] < 0.05) & (diff["log2FoldChange"] < 0)

    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.scatter(
        diff.loc[~significant_down, "log2FoldChange"],
        diff.loc[~significant_down, "minusLog10Padj"],
        color="#BDBDBD",
        s=22,
        alpha=0.65,
        linewidths=0,
    )
    ax.scatter(
        diff.loc[significant_down, "log2FoldChange"],
        diff.loc[significant_down, "minusLog10Padj"],
        color="#377EB8",
        s=24,
        alpha=0.85,
        linewidths=0,
        label="Significantly down",
    )

    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.axhline(-np.log10(0.05), color="#AAAAAA", linewidth=0.7, linestyle="--")

    # Label significantly down-regulated MED-CORE genes and SUPT6H.
    down_genes = set(diff.loc[significant_down, "name"])
    label_genes = sorted((MEDIATOR_CORE & down_genes) | ({"SUPT6H"} & down_genes))
    offsets = {
        "MED6": (8, 8),
        "MED14": (8, -18),
        "MED24": (8, 8),
        "SUPT6H": (-8, -18),
    }
    for gene in label_genes:
        row = diff.loc[diff["name"] == gene].iloc[0]
        dx, dy = offsets.get(gene, (8, 8))
        ax.annotate(
            gene,
            xy=(row["log2FoldChange"], row["minusLog10Padj"]),
            xytext=(dx, dy),
            textcoords="offset points",
            color="#377EB8",
            fontsize=10,
            ha="right" if gene == "SUPT6H" else "left",
            arrowprops={"arrowstyle": "-", "color": "#777777", "linewidth": 0.6},
        )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(padj)")
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(OUTPUT, format="pdf")
    plt.close(fig)

    print(f"Significantly down-regulated genes: {int(significant_down.sum())}")
    print(f"Labels: {', '.join(label_genes)}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
