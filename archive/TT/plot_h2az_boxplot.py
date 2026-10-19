"""Plot log2FoldChange by H2A.Z density tertiles (Q1, Q2, Q3).

Install: python -m pip install pandas matplotlib seaborn scipy
Usage:   python plot_h2az_boxplot.py med12/results.txt
         python plot_h2az_boxplot.py results.txt --density density.bed -o plot.pdf
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="CSV or tab-separated results.txt; first column is the gene ID")
    parser.add_argument("--density", type=Path,
                        default=Path(__file__).resolve().with_name("output_with_density.bed"),
                        help="BED: column 4 = gene ID, column 7 = H2A.Z density")
    parser.add_argument("-o", "--output", type=Path, help="Output PDF path")
    parser.add_argument("--title", help="Optional panel title")
    args = parser.parse_args()

    # Require real Arial rather than silently substituting another font.
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
    except ValueError:
        parser.error("Arial is not installed. Install Arial before running this script.")
    plt.rcParams.update({"font.family": "Arial", "font.size": 12,
                         "pdf.fonttype": 42, "ps.fonttype": 42})

    results = pd.read_csv(args.results, sep=None, engine="python", index_col=0)
    if "log2FoldChange" not in results:
        parser.error("results.txt must contain a log2FoldChange column.")
    bed = pd.read_csv(args.density, sep="\t", header=None, comment="#")
    if bed.shape[1] < 7:
        parser.error("The density BED must have at least 7 columns.")
    density = bed.set_index(3)[6].rename("density")
    results.index = results.index.astype(str)
    density.index = density.index.astype(str)
    if not results.index.is_unique or not density.index.is_unique:
        parser.error("Gene IDs must be unique in both input files.")
    data = pd.concat([results["log2FoldChange"], density], axis=1)
    data = data.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 3:
        parser.error("Fewer than 3 genes have both a finite log2FoldChange and density.")

    order = ["Q1", "Q2", "Q3"]
    try:
        data["group"], edges = pd.qcut(data["density"], q=3, labels=order, retbins=True)
    except ValueError as exc:
        parser.error(f"Cannot form three density groups (too many tied values): {exc}")
    groups = [data.loc[data["group"] == group, "log2FoldChange"] for group in order]
    if any(group.empty for group in groups):
        parser.error("At least one density group is empty.")

    fig, ax = plt.subplots(figsize=(4.8, 6.5))
    sns.boxplot(data=data, x="group", y="log2FoldChange", order=order,
                hue="group", hue_order=order, dodge=False,
                palette=["#FA7F6F", "#8C8A99", "#1E95D4"],
                showfliers=False, width=0.65, ax=ax)
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    ax.set_xlabel("H2A.Z density tertile")
    ax.set_ylabel("log2 fold change")
    ax.set_title(args.title or args.results.parent.name, pad=90)

    # Pairwise two-sided Mann-Whitney U tests, without multiple-test correction.
    # Use axes coordinates so hidden outliers do not expand the displayed range.
    transform = ax.get_xaxis_transform()
    for (i, j), height in zip([(0, 1), (1, 2), (0, 2)], [1.02, 1.10, 1.18]):
        p = mannwhitneyu(groups[i], groups[j], alternative="two-sided").pvalue
        ax.plot([i, i, j, j], [height, height + .015, height + .015, height],
                color="black", linewidth=1, transform=transform, clip_on=False)
        ax.text((i + j) / 2, height + .02, f"M.W.U. p = {p:.3e}",
                ha="center", va="bottom", transform=transform, fontsize=10)
        print(f"{order[i]} vs {order[j]}: p = {p:.6g}")
    ax.set_xticks(range(3), [f"{name}\n(n={len(values):,})" for name, values in zip(order, groups)])
    fig.subplots_adjust(left=.23, right=.96, bottom=.16, top=.70)
    output = args.output or args.results.with_name(f"{args.results.stem}_h2az_boxplot.pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print("Density boundaries:", edges.tolist())
    print("Saved:", output.resolve())


if __name__ == "__main__":
    main()
