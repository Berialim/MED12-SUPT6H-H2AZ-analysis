from pathlib import Path
import sys

# Local gseapy installation used for this project.
ROOT = Path(__file__).resolve().parent
GSEAPY_RUNTIME = ROOT.parents[2] / ".deps" / "gseapy_runtime"
sys.path.insert(0, str(GSEAPY_RUNTIME))

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
import pandas as pd
import gseapy as gp


INPUT = ROOT / "volcannodata.csv"
OUTPUT_COMPLEX = ROOT / "gseapy_mediator_complex.pdf"
OUTPUT_CORE = ROOT / "gseapy_mediator_core.pdf"

# Mediator complex genes represented in this APEX experiment.
MEDIATOR_GENES = [
    "MED1", "MED4", "MED6", "MED12", "MED13", "MED14", "MED15",
    "MED17", "MED19", "MED24", "CDK8", "CDK19", "RBM14", "THRAP3",
]

CORE_GENES = ["MED1", "MED4", "MED6", "MED14", "MED15", "MED17", "MED19", "MED24"]


# Use embedded Arial TrueType text in the PDF.
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


def run_gsea(rnk, genes, term, output):
    result = gp.prerank(
        rnk=rnk,
        gene_sets={term: genes},
        min_size=5,
        max_size=5000,
        permutation_num=1000,
        seed=20260928,
        threads=1,
        outdir=None,
        verbose=False,
    )

    # Native gseapy GSEA plot with gseapy's NES, P value and FDR parameters.
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
    print(f"{term}: NES={float(row['NES']):.4f}, NOM_p={float(row['NOM p-val']):.4g}, FDR_q={float(row['FDR q-val']):.4g}")
    print(f"Wrote {output}")


def main():
    diff = pd.read_csv(INPUT)
    rnk = diff[["name", "log2FoldChange"]].dropna().copy()
    rnk["name"] = rnk["name"].astype(str).str.upper()
    rnk["log2FoldChange"] = pd.to_numeric(rnk["log2FoldChange"], errors="coerce")
    rnk = rnk.dropna().drop_duplicates("name", keep="first")
    rnk = rnk.sort_values(["log2FoldChange", "name"], ascending=[False, True]).reset_index(drop=True)

    # Break rounded-value ties deterministically without changing the ranking.
    rnk["log2FoldChange"] += 1e-12 * (len(rnk) - rnk.index) / len(rnk)

    run_gsea(rnk, MEDIATOR_GENES, "MEDIATOR-COMPLEX", OUTPUT_COMPLEX)
    run_gsea(rnk, CORE_GENES, "MEDIATOR-CORE", OUTPUT_CORE)


if __name__ == "__main__":
    main()
