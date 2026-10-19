#!/usr/bin/env python3
"""Analyze shared MED12/SUPT6H DE groups and H2AZ differences."""

from math import erfc, sqrt
from pathlib import Path

import numpy as np
import pandas as pd

from plot_score_relationships import load_data


ROOT = Path(r"D:\work\pby\mRNAseq")
MED12_FILE = ROOT / "med12_oe" / "results.txt"
SPT6_FILE = ROOT / "spt6_oe" / "results.txt"


def rankdata(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0 + 1.0
        start = end
    return ranks


def mann_whitney_p(a, b):
    """Two-sided normal-approximation Mann-Whitney p-value with ties."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    combined = np.concatenate([a, b])
    ranks = rankdata(combined)
    n1, n2 = len(a), len(b)
    u = ranks[:n1].sum() - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    _, counts = np.unique(combined, return_counts=True)
    tie_term = np.sum(counts**3 - counts)
    n = n1 + n2
    variance = n1 * n2 / 12.0 * (n + 1.0 - tie_term / (n * (n - 1.0)))
    z = (u - mean_u - 0.5 * np.sign(u - mean_u)) / sqrt(max(variance, 1e-12))
    p = erfc(abs(z) / sqrt(2.0))
    return float(u), float(z), float(p)


def hedges_g(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = len(a), len(b)
    pooled = sqrt(((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) / max(n1 + n2 - 2, 1))
    d = (np.mean(a) - np.mean(b)) / pooled if pooled else float("nan")
    correction = 1.0 - 3.0 / max(4.0 * (n1 + n2) - 9.0, 1.0)
    return float(d * correction)


def summarize_group(df, mask):
    raw = df.loc[mask, "h2az_signal"].to_numpy(float)
    z = df.loc[mask, "h2az_zscore"].to_numpy(float)
    log = np.log2(raw + 1.0)
    return {
        "n": int(len(raw)),
        "h2az_signal_mean": float(np.mean(raw)),
        "h2az_signal_median": float(np.median(raw)),
        "h2az_zscore_mean": float(np.mean(z)),
        "h2az_zscore_median": float(np.median(z)),
        "log2_h2az_mean": float(np.mean(log)),
    }


def analyze(effect_cutoff=0.0):
    data, _ = load_data()
    med = pd.read_csv(MED12_FILE, index_col=0)
    spt = pd.read_csv(SPT6_FILE, index_col=0)
    de = pd.DataFrame(index=med.index.union(spt.index))
    de["med12_log2fc"] = pd.to_numeric(med["log2FoldChange"], errors="coerce")
    de["med12_padj"] = pd.to_numeric(med["padj"], errors="coerce")
    de["spt6_log2fc"] = pd.to_numeric(spt["log2FoldChange"], errors="coerce")
    de["spt6_padj"] = pd.to_numeric(spt["padj"], errors="coerce")
    h2az_data = data.set_index("gene_symbol")[["h2az_signal", "h2az_zscore"]]
    merged = h2az_data.join(de, how="inner").dropna(subset=["h2az_signal", "h2az_zscore", "med12_log2fc", "spt6_log2fc", "med12_padj", "spt6_padj"])
    sig = (merged["med12_padj"] < 0.05) & (merged["spt6_padj"] < 0.05)
    cutoff = float(effect_cutoff)
    up = sig & (merged["med12_log2fc"] >= cutoff) & (merged["spt6_log2fc"] >= cutoff)
    down = sig & (merged["med12_log2fc"] <= -cutoff) & (merged["spt6_log2fc"] <= -cutoff)
    result = {"effect_cutoff": cutoff, "n_complete": int(len(merged)), "groups": {"both_up": summarize_group(merged, up), "both_down": summarize_group(merged, down)}}
    a = np.log2(merged.loc[up, "h2az_signal"].to_numpy(float) + 1.0)
    b = np.log2(merged.loc[down, "h2az_signal"].to_numpy(float) + 1.0)
    if len(a) >= 2 and len(b) >= 2:
        _, z, p = mann_whitney_p(a, b)
        result["comparison_log2_h2az"] = {"mean_difference_up_minus_down": float(np.mean(a) - np.mean(b)), "median_difference_up_minus_down": float(np.median(a) - np.median(b)), "mann_whitney_z": z, "mann_whitney_p": p, "hedges_g": hedges_g(a, b)}
    return result


if __name__ == "__main__":
    for cutoff in (0.0, 0.5):
        result = analyze(cutoff)
        print(f"effect_cutoff={cutoff}")
        print(result)
