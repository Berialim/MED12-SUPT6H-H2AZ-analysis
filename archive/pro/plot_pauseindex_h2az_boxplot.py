"""Plot pause-index and log2(KD/CTRL) boxplots by H2A.Z density tertiles.

Dependencies:
    python -m pip install pandas matplotlib scipy

Example:
    python plot_pauseindex_h2az_boxplot.py pauseindex_MED12.txt
    python plot_pauseindex_h2az_boxplot.py pauseindex_SPT6.txt
    python plot_pauseindex_h2az_boxplot.py puaseindex_CDK8.txt --ctrl-pause 9_pause --kd-pause 12_pause
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


COLORS = ["#FA7F6F", "#9994A5", "#1E95D4"]
GROUPS = ["Q1", "Q2", "Q3"]
PAIRS = [(0, 1), (1, 2), (0, 2)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pauseindex", type=Path,
                        help="Pause-index TSV with a gene column and one or more *_pause columns")
    parser.add_argument("--density", type=Path,
                        default=Path(__file__).resolve().parent.parent / "TT" / "output_with_density.bed",
                        help="H2A.Z density BED; defaults to ../TT/output_with_density.bed")
    parser.add_argument("--density-column", type=int, default=7,
                        help="One-based density column in the BED file (default: 7)")
    parser.add_argument("--ctrl-pause", help="Name of the control pause-index column")
    parser.add_argument("--kd-pause", help="Name of the KD pause-index column")
    parser.add_argument("-o", "--output", type=Path,
                        help="Output PDF; defaults beside the pauseindex input")
    args = parser.parse_args()

    # Use Arial throughout and preserve editable vector text in the PDF.
    try:
        from matplotlib import font_manager
        windows_arial = Path("/mnt/c/Windows/Fonts/arial.ttf")
        if windows_arial.exists():
            font_manager.fontManager.addfont(str(windows_arial))
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["pdf.fonttype"] = 42
        plt.rcParams["ps.fonttype"] = 42
        font_manager.findfont("Arial", fallback_to_default=False)
    except ValueError:
        parser.error("Arial was not found. Install Arial, then run this script again.")
    plt.rcParams["font.size"] = 11

    pause = pd.read_csv(args.pauseindex, sep="\t")
    if "gene" not in pause.columns:
        parser.error("The pauseindex table must contain a 'gene' column.")
    pause_cols = [column for column in pause.columns
                  if str(column).endswith("_pause")]
    if not pause_cols:
        parser.error("No pause-index columns ending in '_pause' were found.")
    if pause["gene"].duplicated().any():
        parser.error("The 'gene' column contains duplicate IDs.")

    density_col = args.density_column - 1
    if density_col < 0:
        parser.error("--density-column is one-based and must be at least 1.")
    bed = pd.read_csv(args.density, sep="\t", header=None, comment="#")
    if bed.shape[1] <= max(3, density_col):
        parser.error(f"Density file has {bed.shape[1]} columns; requested column {args.density_column}.")
    density = bed.iloc[:, [3, density_col]].copy()
    density.columns = ["gene", "density"]
    if density["gene"].duplicated().any():
        parser.error("The density file contains duplicate gene IDs in column 4.")

    # Infer obvious CTRL/KD names; require explicit names when columns are numeric labels.
    ctrl_matches = [column for column in pause_cols
                    if str(column).lower() in {"ctrl_pause", "control_pause"}]
    kd_matches = [column for column in pause_cols
                  if "kd" in str(column).lower() and "ctrl" not in str(column).lower()]
    ctrl_pause = args.ctrl_pause or (ctrl_matches[0] if len(ctrl_matches) == 1 else None)
    kd_pause = args.kd_pause or (kd_matches[0] if len(kd_matches) == 1 else None)
    if ctrl_pause not in pause_cols or kd_pause not in pause_cols:
        parser.error(
            "Could not unambiguously identify the KD and CTRL columns. "
            f"Available pause columns: {', '.join(map(str, pause_cols))}. "
            "Set --ctrl-pause and --kd-pause explicitly."
        )

    data = pause[["gene", *pause_cols]].merge(density, on="gene", how="inner")
    for column in [*pause_cols, "density"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=["density"])
    if len(data) < 3:
        parser.error("Fewer than 3 genes have usable pause-index and density values.")

    positive_pause = (data[kd_pause] > 0) & (data[ctrl_pause] > 0)
    data["log2FC_KD_over_CTRL"] = np.nan
    data.loc[positive_pause, "log2FC_KD_over_CTRL"] = np.log2(
        data.loc[positive_pause, kd_pause] / data.loc[positive_pause, ctrl_pause]
    )

    # Assign equal-sized thirds using the same density cut points for every panel.
    try:
        bins = pd.qcut(data["density"], q=3, labels=False,
                       retbins=True, duplicates="drop")
        group_index, edges = bins
    except ValueError as exc:
        parser.error(f"Could not divide density into three groups: {exc}")
    if len(edges) != 4:
        parser.error("Density values have too many ties to form three distinct groups.")
    data["density_group"] = pd.Categorical(
        [GROUPS[int(index)] for index in group_index], categories=GROUPS, ordered=True)

    plot_specs = [(column, column, "Pause index") for column in pause_cols]
    plot_specs.append(("log2FC_KD_over_CTRL", f"log2(KD/CTRL): {kd_pause} / {ctrl_pause}",
                       "log2(KD/CTRL pause index)"))
    fig, axes = plt.subplots(1, len(plot_specs),
                             figsize=(4.2 * len(plot_specs), 6.2), squeeze=False)
    for ax, (column, title, ylabel) in zip(axes[0], plot_specs):
        values = [data.loc[data["density_group"] == group, column].dropna().to_numpy()
                  for group in GROUPS]
        if any(len(group_values) == 0 for group_values in values):
            parser.error(f"{column} has no usable observations in at least one density group.")
        boxes = ax.boxplot(values, tick_labels=GROUPS, patch_artist=True,
                           showfliers=False, widths=0.62)
        for box, color in zip(boxes["boxes"], COLORS):
            box.set_facecolor(color)
            box.set_alpha(0.85)
        for item in [*boxes["whiskers"], *boxes["caps"], *boxes["medians"]]:
            item.set_color("#333333")

        ax.set_title(title, y=1.33, pad=0)
        ax.set_xlabel("H2A.Z density tertile")
        ax.set_ylabel(ylabel)
        ax.set_xticks(range(1, 4),
                      [f"{group}\n(n={len(group_values):,})"
                       for group, group_values in zip(GROUPS, values)])

        # Pairwise two-sided Mann-Whitney U tests, without multiple-test correction.
        transform = ax.get_xaxis_transform()
        for (i, j), height in zip(PAIRS, [1.02, 1.10, 1.18]):
            p_value = mannwhitneyu(values[i], values[j], alternative="two-sided").pvalue
            ax.plot([i + 1, i + 1, j + 1, j + 1],
                    [height, height + 0.015, height + 0.015, height],
                    color="black", linewidth=0.9, transform=transform, clip_on=False)
            ax.text((i + j + 2) / 2, height + 0.018,
                    f"M.W.U. p = {p_value:.2e}", ha="center", va="bottom",
                    fontsize=9, transform=transform)

    fig.suptitle("Pause index and log2 KD/CTRL by H2A.Z density tertiles",
                 y=0.99, fontsize=13)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.70, wspace=0.42)
    output = args.output or args.pauseindex.with_name(
        f"{args.pauseindex.stem}_h2az_boxplot.pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Matched genes: {len(data):,}")
    print(f"Density boundaries (column {args.density_column}): {edges.tolist()}")
    print(f"log2FC calculated as log2({kd_pause} / {ctrl_pause}) for "
          f"{int(data['log2FC_KD_over_CTRL'].notna().sum()):,} genes")
    print(f"Saved PDF: {output.resolve()}")


if __name__ == "__main__":
    main()
