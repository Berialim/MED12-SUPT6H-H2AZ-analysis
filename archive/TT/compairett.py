#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 4 21:02:50 2025
@author: berial
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from pathlib import Path
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

os.chdir("/Users/berial/ljy/lailab2025/pby_med12/TT/")

plt.rcParams.update({
    'pdf.fonttype': 42,
    'font.family':  'Arial',
    'font.size':    12,
})


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.1f}%\n({v:d})'.format(p=pct, v=val)
    return my_autopct


def save_pdf(fig, path):
    fig.savefig(path, format='pdf', bbox_inches='tight',
                metadata={'Creator': ''})
    plt.close(fig)


def sort_pair(col1, col2):
    """确保 cdk8 相关列始终作为横坐标 (col_x)"""
    if 'cdk8' in col2.lower() and 'cdk8' not in col1.lower():
        return col2, col1
    return col1, col2


def main():
    current_dir = Path('.')
    folders = [f for f in current_dir.iterdir() if f.is_dir()]

    all_data = {}
    for folder in folders:
        results_path = folder / 'results.txt'
        if results_path.exists():
            df = pd.read_csv(results_path, index_col=0)
            folder_name = folder.name
            all_data[f"Log2FC_{folder_name}KD"] = df['log2FoldChange']
            all_data[f"padj_{folder_name}KD"]   = df['padj']

    merged_df = pd.DataFrame(all_data)

    bed_path = current_dir / 'output_with_density.bed'
    if bed_path.exists():
        bed_df = pd.read_csv(bed_path, sep='\t', header=None)
        density_series = bed_df.set_index(3)[6].rename('density')
        merged_df = pd.concat([merged_df, density_series], axis=1).dropna(how='all')
        print("Merged DataFrame with density shape:", merged_df.shape)
        print(merged_df.head())
    else:
        print("output_with_density.bed not found!")
        return

    logfc_cols  = [col for col in merged_df.columns if col.startswith('Log2FC_')]
    padj_cols   = [col for col in merged_df.columns if col.startswith('padj_')]
    density_col = 'density' if 'density' in merged_df.columns else None

    if len(logfc_cols) < 1:
        print("No Log2FC data found.")
        return

    # ── 生成所有 pair，cdk8 始终为 col_x ──────────────────────────
    pairs_logfc = []
    for ii in range(len(logfc_cols)):
        for jj in range(ii + 1, len(logfc_cols)):
            cx, cy = sort_pair(logfc_cols[ii], logfc_cols[jj])
            pairs_logfc.append((cx, cy))

    DIRECTION_COLORS = {
        'Up-Up':    'red',
        'Up-Down':  'green',
        'Down-Up':  'orange',
        'Down-Down':'blue',
    }

    # ════════════════════════════════════════════════════════════
    # 1. Boxplot by density quartile + cumulative line plot
    # ════════════════════════════════════════════════════════════
    from statannotations.Annotator import Annotator
    from matplotlib import colors as mcolors

    n_logfc = len(logfc_cols)

    if n_logfc > 0 and density_col:
        cmap_box = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", ["#FA7F6F", "#1E95D4"])
        norm_box   = mcolors.Normalize(vmin=0, vmax=2)
        color_list_box = [cmap_box(norm_box(i)) for i in range(3)]

        fig_box, axes_box = plt.subplots(1, n_logfc, figsize=(3 * n_logfc, 6))
        if n_logfc == 1:
            axes_box = [axes_box]

        for idx, logfc_col in enumerate(logfc_cols):
            ax = axes_box[idx]
            data = merged_df[[logfc_col, density_col]].dropna().copy()
            if len(data) < 3:
                ax.text(0.5, 0.5, 'Insufficient data',
                        ha='center', va='center', transform=ax.transAxes)
                continue
            data['quartile'] = pd.qcut(
                data[density_col], q=3,
                labels=['Q1', 'Q2', 'Q3'], duplicates='drop')

            order = ['Q1', 'Q2', 'Q3']
            sns.boxplot(data=data, x='quartile', y=logfc_col, ax=ax,
                        showfliers=False, order=order,
                        palette={k: v for k, v in zip(order, color_list_box)})
            ax.set_xlabel('H2A.Z Density Quartile')
            ax.set_ylabel(logfc_col)
            ax.set_title(logfc_col)

            pairs_ann = [('Q1', 'Q2'), ('Q1', 'Q3'), ('Q2', 'Q3')]
            annotator = Annotator(ax, pairs_ann, data=data,
                                  x='quartile', y=logfc_col, order=order)
            annotator.configure(test='Mann-Whitney',
                                text_format='full', loc='outside')
            annotator.apply_and_annotate()

            sm = plt.cm.ScalarMappable(cmap=cmap_box, norm=norm_box)
            sm.set_array([])
            plt.colorbar(sm, ax=ax, pad=0.02).set_label('Density Quartile')

        fig_box.suptitle('Boxplots of Log2FC by Density Quartiles', fontsize=14)
        fig_box.tight_layout()
        save_pdf(fig_box, 'boxplots_by_density_quartiles.pdf')

        # cumulative line plots
        for logfc_col in logfc_cols:
            data = merged_df[[logfc_col, density_col]].dropna().copy()
            data['quartile'] = pd.qcut(
                data[density_col], q=3,
                labels=['Q1', 'Q2', 'Q3'], duplicates='drop')

            fig_cdf, ax_cdf = plt.subplots(figsize=(8, 6))
            for i, q in enumerate(['Q1', 'Q2', 'Q3']):
                subset = data[data['quartile'] == q][logfc_col].dropna()
                if subset.empty:
                    continue
                values, base = np.histogram(subset, bins=200, density=True)
                cumulative = np.cumsum(values) / np.sum(values)
                ax_cdf.plot(base[:-1], cumulative, label=q,
                            color=color_list_box[i], linewidth=2)

            ax_cdf.set_xlabel(logfc_col)
            ax_cdf.set_ylabel('Cumulative Density')
            ax_cdf.set_title(f'CDF of {logfc_col} by Density Quartile')
            ax_cdf.legend(title='Quartile')
            sm = plt.cm.ScalarMappable(cmap=cmap_box, norm=norm_box)
            sm.set_array([])
            plt.colorbar(sm, ax=ax_cdf, pad=0.02).set_label('Density Quartile')
            fig_cdf.tight_layout()
            safe = logfc_col.replace('Log2FC_', '').replace('KD', '')
            save_pdf(fig_cdf, f'cdf_{safe}.pdf')

    # ════════════════════════════════════════════════════════════
    # 2. Correlation heatmaps by quartile
    # ════════════════════════════════════════════════════════════
    if n_logfc > 1 and density_col:
        full_data = merged_df[logfc_cols + ['density']].dropna(how='any').copy()
        if len(full_data) >= 3:
            full_data['quartile'] = pd.qcut(
                full_data['density'], q=3,
                labels=['Q1', 'Q2', 'Q3'], duplicates='drop')

            fig_corr, axes_corr = plt.subplots(1, 3, figsize=(10, 5))
            for idx, q_label in enumerate(['Q1', 'Q2', 'Q3']):
                sub = full_data[full_data['quartile'] == q_label][logfc_cols]
                if len(sub) < 2:
                    axes_corr[idx].text(0.5, 0.5, f'Insufficient data for {q_label}',
                                        ha='center', va='center',
                                        transform=axes_corr[idx].transAxes)
                    continue
                corr_matrix = sub.corr(method='pearson')
                sns.heatmap(corr_matrix, annot=True, cmap='coolwarm',
                            center=0, ax=axes_corr[idx],
                            fmt='.3f', square=True)
                axes_corr[idx].set_title(
                    f'Correlation – {q_label}\n(n={len(sub)} genes)')

            fig_corr.suptitle(
                'Pearson Correlation Matrices by Quartile', fontsize=14)
            fig_corr.tight_layout()
            save_pdf(fig_corr, 'correlation_matrices_by_quartiles.pdf')

        # ── per-quartile regplot + pie ────────────────────────────────
        full_data_padj = merged_df[
            logfc_cols + padj_cols + ['density']].dropna(how='any').copy()
        full_data_padj['quartile'] = pd.qcut(
            full_data_padj['density'], q=3,
            labels=['Q1', 'Q2', 'Q3'], duplicates='drop')

        stacked_records = []

        for q_label in ['Q1', 'Q2', 'Q3']:
            sub = full_data_padj[
                full_data_padj['quartile'] == q_label][
                logfc_cols + padj_cols].dropna(how='any')
            if len(sub) < 2:
                continue

            for col_x, col_y in pairs_logfc:           # ← 使用排序后的 pair
                padj_x = col_x.replace('Log2FC_', 'padj_')
                padj_y = col_y.replace('Log2FC_', 'padj_')

                x = sub[col_x];  y = sub[col_y]
                r, pvalue = stats.pearsonr(x, y)

                both = (sub[padj_x] < 0.05) & (sub[padj_y] < 0.05)
                non_b = ~both
                bdf   = sub[both]

                uu = (bdf[col_x] > 0) & (bdf[col_y] > 0)
                ud = (bdf[col_x] > 0) & (bdf[col_y] < 0)
                du = (bdf[col_x] < 0) & (bdf[col_y] > 0)
                dd = (bdf[col_x] < 0) & (bdf[col_y] < 0)

                counts = {'Up-Up': uu.sum(), 'Up-Down': ud.sum(),
                          'Down-Up': du.sum(), 'Down-Down': dd.sum()}
                total  = len(bdf)

                safe_x = col_x.replace('Log2FC_', '').replace('KD', '')
                safe_y = col_y.replace('Log2FC_', '').replace('KD', '')

                stacked_records.append({
                    'quartile': q_label,
                    'pair':     f'{safe_x}_vs_{safe_y}',
                    **counts,
                    'total_both_sig': total,
                })

                print(f"\n[{q_label}] {col_x}(x) vs {col_y}(y): "
                      f"double sig={total}, r={r:.3f}")

                # combined regplot + pie
                fig, (ax_reg, ax_pie) = plt.subplots(
                    1, 2, figsize=(10, 5),
                    gridspec_kw={'width_ratios': [2, 1]})

                if non_b.sum() > 0:
                    sns.scatterplot(x=x[non_b], y=y[non_b],
                                    color='lightblue', alpha=0.5,
                                    ax=ax_reg, label='Non-double sig')

                for d_key, d_mask in zip(
                        ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
                        [uu, ud, du, dd]):
                    if d_mask.sum() > 0:
                        idx_b = bdf.index[d_mask]
                        mf    = sub.index.isin(idx_b)
                        sns.scatterplot(
                            x=x[mf], y=y[mf],
                            color=DIRECTION_COLORS[d_key], alpha=0.85,
                            ax=ax_reg,
                            label=f'{d_key} (n={d_mask.sum()})')

                sns.regplot(x=x, y=y, scatter=False, color='black',
                            ax=ax_reg, line_kws={'lw': 1.5})
                ax_reg.annotate(
                    f'Up-Up: {uu.sum()}\nDown-Down: {dd.sum()}',
                    xy=(0.05, 0.92), xycoords='axes fraction', fontsize=10,
                    bbox=dict(boxstyle='round,pad=0.3',
                              fc='lightyellow', ec='gray', alpha=0.8))
                ax_reg.set_xlim(-2, 2);  ax_reg.set_ylim(-2, 2)
                ax_reg.axhline(0, color='gray', lw=0.8, ls='--')
                ax_reg.axvline(0, color='gray', lw=0.8, ls='--')
                ax_reg.set_xlabel(col_x)        # cdk8 在 x 轴
                ax_reg.set_ylabel(col_y)
                ax_reg.set_title(
                    f'{col_x} vs {col_y}\nR={r:.3f}, P={pvalue:.3e}')
                ax_reg.legend(loc='lower right', fontsize=8)

                if total > 0:
                    sizes = [counts[k] for k in
                             ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down']]
                    ax_pie.pie(
                        sizes,
                        labels=['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
                        colors=[DIRECTION_COLORS[k] for k in
                                ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down']],
                        autopct=make_autopct(sizes), startangle=90)
                    ax_pie.set_title(
                        f'Direction (double sig)\nTotal: {total}')
                else:
                    ax_pie.text(0.5, 0.5, 'No Double Sig',
                                ha='center', va='center',
                                transform=ax_pie.transAxes)

                fig.suptitle(
                    f'{q_label}: {safe_x}(x) vs {safe_y}(y)  |  '
                    f'Up-Up={uu.sum()}, Down-Down={dd.sum()}',
                    fontsize=14)
                fig.tight_layout()
                save_pdf(fig,
                    f'combined_regplot_pie_{safe_x}_{safe_y}_{q_label.lower()}.pdf')

        # ════════════════════════════════════════════════════════
        # 3. Stacked barplot: percentage, Q1/Q2/Q3 per pair
        # ════════════════════════════════════════════════════════
        if stacked_records:
            stacked_df = pd.DataFrame(stacked_records)
            stacked_df.to_csv('quartile_direction_counts.csv', index=False)
            print("\n→ Saved: quartile_direction_counts.csv")

            for pair in stacked_df['pair'].unique():
                pdf_raw = stacked_df[stacked_df['pair'] == pair].set_index('quartile')
                direction_cols = ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down']
                pdf_raw = pdf_raw[direction_cols + ['total_both_sig']].reindex(
                    ['Q1', 'Q2', 'Q3'])

                counts_df = pdf_raw[direction_cols]
                totals    = pdf_raw['total_both_sig'].fillna(0).astype(int)
                pdf_pct   = counts_df.div(counts_df.sum(axis=1), axis=0) * 100

                fig_bar, ax_bar = plt.subplots(figsize=(7, 5))
                bottom = np.zeros(3)
                x_pos  = np.arange(3)

                for d_key in direction_cols:
                    vals = pdf_pct[d_key].fillna(0).values.astype(float)
                    raw  = counts_df[d_key].fillna(0).values.astype(int)
                    bars = ax_bar.bar(
                        x_pos, vals, bottom=bottom,
                        color=DIRECTION_COLORS[d_key],
                        label=d_key, edgecolor='white', linewidth=0.5)
                    for bar, pct, n, bot in zip(bars, vals, raw, bottom):
                        if pct > 3:
                            ax_bar.text(
                                bar.get_x() + bar.get_width() / 2,
                                bot + pct / 2,
                                f'{pct:.1f}%\n(n={n})',
                                ha='center', va='center',
                                fontsize=8, color='white', fontweight='bold')
                    bottom += vals

                ax_bar.set_ylim(0, 108)
                ax_bar.set_xticks(x_pos)
                ax_bar.set_xticklabels(
                    [f'{q}\n(total={totals[q]})' for q in ['Q1', 'Q2', 'Q3']])
                ax_bar.set_xlabel('H2A.Z Density Quartile')
                ax_bar.set_ylabel('Percentage (%)')
                ax_bar.set_title(f'Direction percentage by quartile\n{pair}')
                ax_bar.legend(title='Direction',
                              bbox_to_anchor=(1.02, 1), loc='upper left')
                fig_bar.tight_layout()
                save_pdf(fig_bar, f'stacked_barplot_{pair}.pdf')
                print(f"→ Saved: stacked_barplot_{pair}.pdf")

    # ════════════════════════════════════════════════════════════
    # 4. Pairwise regplot grid (upper triangle)
    # ════════════════════════════════════════════════════════════
    plot_cols = logfc_cols + ([density_col] if density_col else [])
    n_plot = len(plot_cols)
    if n_plot < 2:
        print("Need at least 2 plot columns for pairwise.")
        return

    fig_pw, axes_pw = plt.subplots(n_plot, n_plot, figsize=(20, 16))
    axes_pw = axes_pw if n_plot > 1 else axes_pw.reshape(1, 1)

    for i in range(n_plot):
        for j in range(n_plot):
            if i >= j:
                axes_pw[i, j].set_visible(False)
                continue

            raw_x = plot_cols[i];  raw_y = plot_cols[j]
            # 上三角格子也遵守 cdk8-as-x 原则
            if (raw_x.startswith('Log2FC_') and raw_y.startswith('Log2FC_')):
                x_col, y_col = sort_pair(raw_x, raw_y)
            else:
                x_col, y_col = raw_x, raw_y

            common = merged_df.dropna(subset=[x_col, y_col]).index
            if len(common) < 2:
                axes_pw[i, j].text(0.5, 0.5, 'Insufficient data',
                                   ha='center', va='center',
                                   transform=axes_pw[i, j].transAxes)
                continue

            x = merged_df.loc[common, x_col]
            y = merged_df.loc[common, y_col]
            r, pvalue = stats.pearsonr(x, y)
            ax = axes_pw[i, j]

            is_both_logfc = (x_col.startswith('Log2FC_') and
                             y_col.startswith('Log2FC_'))
            if is_both_logfc:
                px = x_col.replace('Log2FC_', 'padj_')
                py = y_col.replace('Log2FC_', 'padj_')
                if px in merged_df.columns and py in merged_df.columns:
                    padj_x = merged_df.loc[common, px]
                    padj_y = merged_df.loc[common, py]
                    both  = (padj_x < 0.05) & (padj_y < 0.05)
                    non_b = ~both
                    bdf2  = merged_df.loc[common][[x_col, y_col]][both]
                    uu2 = (bdf2[x_col] > 0) & (bdf2[y_col] > 0)
                    ud2 = (bdf2[x_col] > 0) & (bdf2[y_col] < 0)
                    du2 = (bdf2[x_col] < 0) & (bdf2[y_col] > 0)
                    dd2 = (bdf2[x_col] < 0) & (bdf2[y_col] < 0)

                    if non_b.sum() > 0:
                        sns.scatterplot(x=x[non_b], y=y[non_b],
                                        color='lightblue', alpha=0.5, ax=ax)
                    for d_key, d_mask in zip(
                            ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
                            [uu2, ud2, du2, dd2]):
                        if d_mask.sum() > 0:
                            idx_b = bdf2.index[d_mask]
                            mf = merged_df.loc[common].index.isin(idx_b)
                            sns.scatterplot(
                                x=x[mf], y=y[mf],
                                color=DIRECTION_COLORS[d_key],
                                alpha=0.8, ax=ax)
                    sns.regplot(x=x, y=y, scatter=False,
                                color='black', ax=ax, line_kws={'lw': 1.5})
                    ax.set_xlim(-2, 2);  ax.set_ylim(-2, 2)

                    legend_elements = [
                        Line2D([0],[0], marker='o', color='w',
                               label='Non-double sig',
                               markerfacecolor='lightblue', markersize=8),
                    ] + [
                        Line2D([0],[0], marker='o', color='w', label=k,
                               markerfacecolor=DIRECTION_COLORS[k], markersize=8)
                        for k in ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down']
                    ]
                    ax.legend(handles=legend_elements,
                              loc='upper right', fontsize=6)
                else:
                    sns.regplot(x=x, y=y, ax=ax, scatter_kws={'alpha': 0.5})
                    ax.set_xlim(-2, 2);  ax.set_ylim(-2, 2)
            else:
                sns.regplot(x=x, y=y, ax=ax, scatter_kws={'alpha': 0.5})

            ax.set_xlabel(x_col)    # cdk8 在 x 轴
            ax.set_ylabel(y_col)
            ax.set_title(f'{x_col} vs {y_col}\nR={r:.3f}, P={pvalue:.3e}')

    fig_pw.suptitle('Pairwise Regression Plots (including density)', fontsize=16)
    fig_pw.tight_layout()
    save_pdf(fig_pw, 'pairwise_regression_plots.pdf')

    # ════════════════════════════════════════════════════════════
    # 5. All-genes combined regplot + pie + summary CSV
    # ════════════════════════════════════════════════════════════
    if n_logfc >= 2:
        summary_records = []

        for col_x, col_y in pairs_logfc:            # ← cdk8 已在 x
            padj_x = col_x.replace('Log2FC_', 'padj_')
            padj_y = col_y.replace('Log2FC_', 'padj_')

            data_pair = merged_df.dropna(
                subset=[col_x, col_y, padj_x, padj_y]).copy()
            x = data_pair[col_x];  y = data_pair[col_y]
            r, pvalue = stats.pearsonr(x, y)

            both  = (data_pair[padj_x] < 0.05) & (data_pair[padj_y] < 0.05)
            non_b = ~both
            bdf   = data_pair[both]

            uu = (bdf[col_x] > 0) & (bdf[col_y] > 0)
            ud = (bdf[col_x] > 0) & (bdf[col_y] < 0)
            du = (bdf[col_x] < 0) & (bdf[col_y] > 0)
            dd = (bdf[col_x] < 0) & (bdf[col_y] < 0)

            uu_c, ud_c = uu.sum(), ud.sum()
            du_c, dd_c = du.sum(), dd.sum()
            total = len(bdf)

            safe_x = col_x.replace('Log2FC_', '').replace('KD', '')
            safe_y = col_y.replace('Log2FC_', '').replace('KD', '')

            print(f"\n[All genes] {col_x}(x) vs {col_y}(y): double sig={total}")

            # gene-level CSV
            direction_label = pd.Series('Other', index=bdf.index)
            direction_label[uu] = 'Up-Up'
            direction_label[ud] = 'Up-Down'
            direction_label[du] = 'Down-Up'
            direction_label[dd] = 'Down-Down'

            out_df = bdf[[col_x, padj_x, col_y, padj_y]].copy()
            out_df['direction'] = direction_label
            out_df.index.name = 'gene'
            csv_name = f'allgenes_bothsig_{safe_x}_{safe_y}.csv'
            out_df.to_csv(csv_name)
            print(f"  → Saved: {csv_name}")

            summary_records.append({
                'pair':          f'{safe_x}_vs_{safe_y}',
                'col_x':         col_x,   'col_y':        col_y,
                'n_total':       len(data_pair),
                'n_both_sig':    total,
                'n_Up_Up':       uu_c,    'n_Up_Down':    ud_c,
                'n_Down_Up':     du_c,    'n_Down_Down':  dd_c,
                'pct_Up_Up':     round(uu_c/total*100,1) if total else 0,
                'pct_Up_Down':   round(ud_c/total*100,1) if total else 0,
                'pct_Down_Up':   round(du_c/total*100,1) if total else 0,
                'pct_Down_Down': round(dd_c/total*100,1) if total else 0,
                'pearson_r':     round(r, 4),
                'pearson_p':     pvalue,
            })

            # combined regplot + pie
            fig_all, (ax_reg, ax_pie) = plt.subplots(
                1, 2, figsize=(15, 6),
                gridspec_kw={'width_ratios': [2, 1]})

            if non_b.sum() > 0:
                sns.scatterplot(x=x[non_b], y=y[non_b],
                                color='lightblue', alpha=0.5,
                                ax=ax_reg, label='Non-double sig')

            for d_key, d_mask, d_cnt in zip(
                    ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
                    [uu, ud, du, dd],
                    [uu_c, ud_c, du_c, dd_c]):
                if d_mask.sum() > 0:
                    idx_b = bdf.index[d_mask]
                    mf    = data_pair.index.isin(idx_b)
                    sns.scatterplot(
                        x=x[mf], y=y[mf],
                        color=DIRECTION_COLORS[d_key], alpha=0.85,
                        ax=ax_reg, label=f'{d_key} (n={d_cnt})')

            sns.regplot(x=x, y=y, scatter=False, color='black',
                        ax=ax_reg, line_kws={'lw': 1.5})
            ax_reg.annotate(
                f'Up-Up: {uu_c}\nDown-Down: {dd_c}',
                xy=(0.05, 0.92), xycoords='axes fraction', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.3',
                          fc='lightyellow', ec='gray', alpha=0.8))
            ax_reg.set_xlim(-2, 2);  ax_reg.set_ylim(-2, 2)
            ax_reg.axhline(0, color='gray', lw=0.8, ls='--')
            ax_reg.axvline(0, color='gray', lw=0.8, ls='--')
            ax_reg.set_xlabel(col_x)        # cdk8 在 x 轴
            ax_reg.set_ylabel(col_y)
            ax_reg.set_title(
                f'{col_x} vs {col_y}\nR={r:.3f}, P={pvalue:.3e}  '
                f'(n={len(data_pair)})')
            ax_reg.legend(loc='lower right', fontsize=8)

            if total > 0:
                sizes = [uu_c, ud_c, du_c, dd_c]
                ax_pie.pie(
                    sizes,
                    labels=['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
                    colors=[DIRECTION_COLORS[k] for k in
                            ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down']],
                    autopct=make_autopct(sizes), startangle=90)
                ax_pie.set_title(f'Direction (double sig)\nTotal: {total}')
            else:
                ax_pie.text(0.5, 0.5, 'No Double Sig',
                            ha='center', va='center',
                            transform=ax_pie.transAxes)

            fig_all.suptitle(
                f'All genes: {safe_x}(x) vs {safe_y}(y)  |  '
                f'Up-Up={uu_c}, Down-Down={dd_c}', fontsize=14)
            fig_all.tight_layout()
            save_pdf(fig_all, f'allgenes_regplot_pie_{safe_x}_{safe_y}.pdf')

        if summary_records:
            summary_df = pd.DataFrame(summary_records)
            summary_df.to_csv('allgenes_direction_summary.csv', index=False)
            print("\n→ Saved: allgenes_direction_summary.csv")
            print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
    
    
    
# %%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from pathlib import Path
import numpy as np
from matplotlib.lines import Line2D

os.chdir("/Users/berial/ljy/lailab2025/pby_med12/TT/")

plt.rcParams.update({
    'pdf.fonttype': 42,
    'font.family': 'Arial',
    'font.size': 12,
})

# 防止大数据绘图卡顿
plt.rcParams['agg.path.chunksize'] = 10000


def save_pdf(fig, path):
    fig.savefig(path, format='pdf', bbox_inches='tight',
                metadata={'Creator': ''})
    plt.close(fig)


def sort_pair(col1, col2):
    if 'cdk8' in col2.lower() and 'cdk8' not in col1.lower():
        return col2, col1
    return col1, col2


def main():
    current_dir = Path('.')
    folders = [f for f in current_dir.iterdir() if f.is_dir()]

    all_data = {}
    for folder in folders:
        results_path = folder / 'results.txt'
        if results_path.exists():
            df = pd.read_csv(results_path, index_col=0)
            name = folder.name
            all_data[f"Log2FC_{name}KD"] = df['log2FoldChange']
            all_data[f"padj_{name}KD"] = df['padj']

    merged_df = pd.DataFrame(all_data)

    bed_df = pd.read_csv("output_with_density.bed", sep="\t", header=None)
    density_series = bed_df.set_index(3)[6].rename('density')
    merged_df = pd.concat([merged_df, density_series], axis=1)

    logfc_cols = [c for c in merged_df.columns if c.startswith("Log2FC_")]
    padj_cols = [c for c in merged_df.columns if c.startswith("padj_")]

    # ===============================
    # 核心绘图函数（统一风格）
    # ===============================
    def plot_reg(ax, x, y, both_mask, colors, title, xlabel, ylabel):
        non_b = ~both_mask

        # 非双显著
        ax.scatter(
            x[non_b], y[non_b],
            color='lightblue', alpha=0.4,
            s=8, rasterized=True
        )

        # 双显著分类
        bdf = pd.DataFrame({'x': x, 'y': y})[both_mask]

        uu = (bdf['x'] > 0) & (bdf['y'] > 0)
        ud = (bdf['x'] > 0) & (bdf['y'] < 0)
        du = (bdf['x'] < 0) & (bdf['y'] > 0)
        dd = (bdf['x'] < 0) & (bdf['y'] < 0)

        for label, mask in zip(
            ['Up-Up', 'Up-Down', 'Down-Up', 'Down-Down'],
            [uu, ud, du, dd]
        ):
            if mask.sum() > 0:
                sub = bdf[mask]
                ax.scatter(
                    sub['x'], sub['y'],
                    color=colors[label],
                    s=10, alpha=0.85,
                    rasterized=True
                )

        # 回归线
        sns.regplot(x=x, y=y, scatter=False,
                    color='black', ax=ax,
                    line_kws={'lw': 1.2})

        # 坐标轴
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_aspect('equal', adjustable='box')

        ax.axhline(0, color='gray', lw=0.8, ls='--')
        ax.axvline(0, color='gray', lw=0.8, ls='--')

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        sns.despine(ax=ax)

    # ===============================
    # 颜色定义
    # ===============================
    COLORS = {
        'Up-Up': 'red',
        'Up-Down': 'green',
        'Down-Up': 'orange',
        'Down-Down': 'blue'
    }

    # ===============================
    # 主分析：All genes regplot
    # ===============================
    for i in range(len(logfc_cols)):
        for j in range(i + 1, len(logfc_cols)):

            col_x, col_y = sort_pair(logfc_cols[i], logfc_cols[j])
            padj_x = col_x.replace('Log2FC_', 'padj_')
            padj_y = col_y.replace('Log2FC_', 'padj_')

            df = merged_df.dropna(subset=[col_x, col_y, padj_x, padj_y])

            if len(df) < 10:
                continue

            x = df[col_x]
            y = df[col_y]

            both = (df[padj_x] < 0.05) & (df[padj_y] < 0.05)

            r, p = stats.pearsonr(x, y)

            fig, ax = plt.subplots(figsize=(4, 4))

            plot_reg(
                ax=ax,
                x=x,
                y=y,
                both_mask=both,
                colors=COLORS,
                title=f"{col_x} vs {col_y}\nR={r:.3f}, P={p:.2e}",
                xlabel=col_x,
                ylabel=col_y
            )

            safe_x = col_x.replace('Log2FC_', '')
            safe_y = col_y.replace('Log2FC_', '')

            save_pdf(fig, f"regplot_{safe_x}_{safe_y}.pdf")

    # ===============================
    # Pairwise grid（高级版）
    # ===============================
    plot_cols = logfc_cols + ['density']
    n = len(plot_cols)

    fig, axes = plt.subplots(n, n, figsize=(18, 16))

    for i in range(n):
        for j in range(n):
            if i >= j:
                axes[i, j].set_visible(False)
                continue

            x_col = plot_cols[i]
            y_col = plot_cols[j]

            df = merged_df.dropna(subset=[x_col, y_col])

            if len(df) < 10:
                continue

            x = df[x_col]
            y = df[y_col]

            ax = axes[i, j]

            if x_col.startswith("Log2FC_") and y_col.startswith("Log2FC_"):
                padj_x = x_col.replace('Log2FC_', 'padj_')
                padj_y = y_col.replace('Log2FC_', 'padj_')

                both = (df[padj_x] < 0.05) & (df[padj_y] < 0.05)

                plot_reg(
                    ax, x, y, both, COLORS,
                    f"{x_col} vs {y_col}",
                    x_col, y_col
                )
            else:
                ax.scatter(x, y, s=8, alpha=0.4, rasterized=True)
                sns.regplot(x=x, y=y, scatter=False, ax=ax)

            ax.set_aspect('equal', adjustable='box')

    fig.tight_layout()
    save_pdf(fig, "pairwise_regression.pdf")


if __name__ == "__main__":
    main()