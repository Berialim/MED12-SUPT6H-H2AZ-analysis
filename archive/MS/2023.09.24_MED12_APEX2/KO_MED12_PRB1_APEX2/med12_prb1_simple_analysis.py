from pathlib import Path
import sys

# Use the same local scientific Python packages as the previous analysis.
ROOT = Path(__file__).resolve().parent
GSEAPY_RUNTIME = ROOT.parents[2] / ".deps" / "gseapy_runtime"
sys.path.insert(0, str(GSEAPY_RUNTIME))

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
import gseapy as gp
import numpy as np
import pandas as pd


INPUT = ROOT / "volcannodata.csv"
OUTPUT_COMPLEX = ROOT / "gseapy_mediator_complex.pdf"
OUTPUT_CORE = ROOT / "gseapy_mediator_core.pdf"
OUTPUT_VOLCANO = ROOT / "volcano_significant_down_mediator_core_supt6h.pdf"
OUTPUT_GSEA_TABLE = ROOT / "gseapy_mediator_results.csv"
OUTPUT_DOWN_TABLE = ROOT / "volcano_significant_down_genes.csv"
OUTPUT_CORE_TABLE = ROOT / "mediator_core_status.csv"

MEDIATOR_GENES = [
    "MED1", "MED4", "MED6", "MED12", "MED13", "MED14", "MED15",
    "MED17", "MED19", "MED24", "CDK8", "CDK19", "RBM14", "THRAP3",
]
MEDIATOR_CORE = ["MED1", "MED4", "MED6", "MED14", "MED15", "MED17", "MED19", "MED24"]

# Set Arial before any plotting; TrueType text remains editable in Illustrator.
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


def prepare_data():
    diff = pd.read_csv(INPUT)
    diff["name"] = diff["name"].astype(str).str.upper().str.strip()
    diff["padj"] = pd.to_numeric(diff["padj"], errors="coerce")
    diff["log2FoldChange"] = pd.to_numeric(diff["log2FoldChange"], errors="coerce")
    diff = diff.dropna(subset=["name", "padj", "log2FoldChange"]).copy()
    diff["minusLog10Padj"] = -np.log10(diff["padj"].clip(lower=np.finfo(float).tiny))
    diff["significant_down"] = (diff["padj"] < 0.05) & (diff["log2FoldChange"] < 0)
    return diff


def run_gsea(rnk, genes, term, output):
    detected_genes = sorted(set(genes) & set(rnk["name"]))
    if len(detected_genes) < 3:
        print(f"Skip {term}: only {len(detected_genes)} genes detected")
        return None

    # Keep the previous minimum of 5 when possible; use 3 for the smaller core set.
    min_size = min(5, len(detected_genes))
    result = gp.prerank(
        rnk=rnk[["name", "log2FoldChange"]],
        gene_sets={term: genes},
        min_size=min_size,
        max_size=5000,
        permutation_num=1000,
        seed=20260928,
        threads=1,
        outdir=None,
        verbose=False,
    )

    row = result.res2d.iloc[0]
    detail = result.results[term]
    gp.gseaplot(
        term=term,
        hits=detail["hits"],
        nes=float(row["NES"]),
        pval=float(row["NOM p-val"]),
        fdr=float(row["FDR q-val"]),
        RES=detail["RES"],
        rank_metric=rnk["log2FoldChange"].values,
        ofname=str(output),
        figsize=(6, 4),
    )
    plt.close("all")

    table_row = row.to_dict()
    table_row.update({
        "Term": term,
        "Detected_genes": ";".join(detected_genes),
        "Detected_gene_count": len(detected_genes),
        "Min_size": min_size,
    })
    print(
        f"{term}: detected={len(detected_genes)}, "
        f"NES={float(row['NES']):.4f}, "
        f"NOM_p={float(row['NOM p-val']):.4g}, "
        f"FDR_q={float(row['FDR q-val']):.4g}"
    )
    print(f"Wrote {output}")
    return table_row


def draw_volcano(diff):
    significant_down = diff["significant_down"]
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.scatter(
        diff.loc[~significant_down, "log2FoldChange"],
        diff.loc[~significant_down, "minusLog10Padj"],
        color="#BDBDBD", s=22, alpha=0.65, linewidths=0,
    )
    ax.scatter(
        diff.loc[significant_down, "log2FoldChange"],
        diff.loc[significant_down, "minusLog10Padj"],
        color="#377EB8", s=24, alpha=0.85, linewidths=0,
        label="Significantly down",
    )
    ax.axvline(0, color="#777777", linewidth=0.7)
    ax.axhline(-np.log10(0.05), color="#AAAAAA", linewidth=0.7, linestyle="--")

    # Label significant down-regulated MED-CORE genes and SUPT6H.
    down_genes = set(diff.loc[significant_down, "name"])
    label_genes = sorted((set(MEDIATOR_CORE) & down_genes) | ({"SUPT6H"} & down_genes))
    offsets = {"SUPT6H": (-8, -18), "MED6": (8, 8), "MED14": (8, -18), "MED24": (8, 8)}
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
    fig.savefig(OUTPUT_VOLCANO, format="pdf")
    plt.close(fig)
    print(f"Significantly down-regulated genes: {int(significant_down.sum())}")
    print(f"Volcano labels: {', '.join(label_genes) if label_genes else 'none'}")
    print(f"Wrote {OUTPUT_VOLCANO}")


def main():
    diff = prepare_data()
    draw_volcano(diff)

    rnk = diff[["name", "log2FoldChange"]].drop_duplicates("name", keep="first")
    rnk = rnk.sort_values(["log2FoldChange", "name"], ascending=[False, True]).reset_index(drop=True)
    rnk["log2FoldChange"] += 1e-12 * (len(rnk) - rnk.index) / len(rnk)

    gsea_rows = []
    for genes, term, output in [
        (MEDIATOR_GENES, "MEDIATOR-COMPLEX", OUTPUT_COMPLEX),
        (MEDIATOR_CORE, "MEDIATOR-CORE", OUTPUT_CORE),
    ]:
        row = run_gsea(rnk, genes, term, output)
        if row is not None:
            gsea_rows.append(row)
    if gsea_rows:
        pd.DataFrame(gsea_rows).to_csv(OUTPUT_GSEA_TABLE, index=False)

    down_columns = ["name", "ID", "padj", "log2FoldChange", "minusLog10Padj"]
    diff.loc[diff["significant_down"], down_columns].sort_values("padj").to_csv(
        OUTPUT_DOWN_TABLE, index=False
    )

    core_columns = ["name", "padj", "log2FoldChange", "minusLog10Padj", "significant_down"]
    diff.loc[diff["name"].isin(MEDIATOR_CORE), core_columns].sort_values("name").to_csv(
        OUTPUT_CORE_TABLE, index=False
    )
    print(f"Wrote {OUTPUT_GSEA_TABLE}")
    print(f"Wrote {OUTPUT_DOWN_TABLE}")
    print(f"Wrote {OUTPUT_CORE_TABLE}")


if __name__ == "__main__":
    main()
