# -*- coding: utf-8 -*-
"""
This script performs Gene Set Enrichment Analysis (GSEA) on a differential
expression results file using the gseapy library.

It now ranks genes using a combined metric of log2FoldChange and padj
to ensure a more robust analysis. It then performs GSEA against three
different gene set databases: KEGG, GO, and Reactome.

Dependencies:
- pandas: For data manipulation.
- gseapy: For GSEA analysis.

To install the dependencies, run:
pip install pandas gseapy
"""

import pandas as pd
import gseapy as gp
import os
import sys
import numpy as np
from gseapy import dotplot
import matplotlib.pyplot as plt
from gseapy import Msigdb
msig = Msigdb
# mouse hallmark gene sets
gmt = msig.get_gmt(category='mh.all', dbver="2023.1.Mm")
# --- Configuration ---
# Set the name of your input file
INPUT_FILE = "results.txt"
# Specify the gene set libraries to use for analysis.
# These are stable and supported versions for mouse.
GENE_SET_LIBRARIES = {
    'KEGG': 'KEGG_2021_Human',
    'HALLMARK': '/Users/berial/reference/hg38/h.all.v2025.1.Hs.symbols.gmt' 
}

# --- Main Script ---
def run_gsea():
    """
    Loads data, preprocesses, and runs GSEA analysis for multiple gene sets.
    """
    print(f"Loading data from '{INPUT_FILE}'...")
    if not os.path.exists(INPUT_FILE):
        print(f"Error: The file '{INPUT_FILE}' was not found.")

    try:
        # Load the data, using the first column as the index (gene names)
        df = pd.read_csv(INPUT_FILE, index_col=0)
    except Exception as e:
        print(f"Error reading the file: {e}")

    print("Data loaded successfully.")
    print("Preprocessing data...")

    # The file has a header like '""', so we rename the index to 'Gene_Name' for clarity.
    df.index.name = 'Gene_Name'

    # Filter out rows with NA values in the columns required for ranking.
    # This is crucial for a robust analysis.
    df_cleaned = df.dropna(subset=['log2FoldChange', 'padj'])

    if df_cleaned.empty:
        print("Error: No valid data found in 'log2FoldChange' or 'padj' columns after cleaning.")

    # We need to handle padj values that are exactly 0, which can lead to
    # infinite values when taking the -log10. We replace them with a very small number.
    df_cleaned.loc[:, 'padj'] = df_cleaned['padj'].replace(0, np.finfo(float).eps)

    # Create a new ranking metric that combines log2FoldChange and padj.
    # The new rank is -log10(padj) * sign(log2FoldChange).
    # This gives a high positive score to significantly upregulated genes
    # and a high negative score to significantly downregulated genes.
    # df_cleaned = df_cleaned[df_cleaned['padj'] < 0.05]
    df_cleaned.loc[:, 'combined_rank'] = df_cleaned['log2FoldChange']

    # Rank the genes based on this new combined metric in descending order.
    # 替换 padj 为 0 的值，避免 -log10(0)
    df_cleaned.loc[:, 'padj'] = df_cleaned['padj'].replace(0, np.finfo(float).eps)
    
    # 新的 ranking metric
    df_cleaned.loc[:, 'ranking_score'] = df_cleaned['log2FoldChange']  # * (-np.log10(df_cleaned['padj']))
    
    # 按 ranking_score 排序
    ranked_genes = df_cleaned.sort_values(by='ranking_score', ascending=False)
    
    print("Data preprocessed and ranked by 'combined_rank' (log2FC and padj).")
    
    # Iterate through each gene set library and perform GSEA
    for name, library in GENE_SET_LIBRARIES.items():
        output_dir = f'gsea_results_{name}'
        print(f"\n--- Running GSEA for '{name}' using library '{library}' ---")
        
        # Perform the GSEA analysis using the prerank function.
        prerank_results = gp.prerank(
            rnk=ranked_genes['ranking_score'],
            gene_sets=library,
            outdir=output_dir,
            permutation_num=1000,
            min_size=5,
            max_size=1000,
            threads=4,
            no_plot=False
        )
        # 绘制 PI3K-Akt pathway富集图
        terms = prerank_results.res2d['Term'].tolist()
        pi3k_terms = ["Oxidative phosphorylation", "Peroxisome", "ECM-receptor interaction", "Hippo signaling pathway"]# [t for t in terms if 'PI3K' in t]
        if pi3k_terms:
            prerank_results.plot(terms=pi3k_terms, ofname=os.path.join(output_dir, "PI3K_Akt_GSEA_plot.pdf"))
            print("✅ PI3K-Akt GSEA plot generated.")
        else:
            print("⚠️ PI3K-Akt pathway not found in results.")
        dotplot(prerank_results.res2d, ofname=os.path.join(output_dir, "GSEA_dotplot.pdf"), 
             column="FDR q-val",
             title='KEGG_2016',
             cmap=plt.cm.viridis,
             size=6, # adjust dot size
             figsize=(6,7), cutoff=0.25, show_ring=False)
        prerank_results.plot(terms[:5], show_ranking=True, ofname=os.path.join(output_dir, "GSEA_top5plot.pdf"))
        print(f"GSEA for '{name}' complete.")
        print(f"Results are saved in the '{output_dir}' directory.")
        print(f"Analysis summary for '{name}':")
        print(prerank_results.res2d.head())
        
    return prerank_results

if __name__ == "__main__":
    prerank_results = run_gsea()
# %%

# Load the data, using the first column as the index (gene names)
df = pd.read_csv(INPUT_FILE, index_col=0)


print("Data loaded successfully.")
print("Preprocessing data...")

# The file has a header like '""', so we rename the index to 'Gene_Name' for clarity.
df.index.name = 'Gene_Name'

# Filter out rows with NA values in the columns required for ranking.
# This is crucial for a robust analysis.
df_cleaned = df.dropna(subset=['log2FoldChange', 'padj'])
# %%

# Define significance thresholds (common defaults; adjust if needed)
FOLD_CHANGE_THRESHOLD = 1  # For up/down: >0 for up, <0 for down
PADJ_THRESHOLD = 0.05    # Adjusted p-value < 0.05 for significance

# Screen out significantly upregulated genes
upregulated = df_cleaned[
    (df_cleaned['log2FoldChange'] > FOLD_CHANGE_THRESHOLD) & 
    (df_cleaned['padj'] < PADJ_THRESHOLD)
]
up_genes = set(upregulated.index.tolist())  # Extract gene names as set for fast lookup

# Screen out significantly downregulated genes
downregulated = df_cleaned[
    (df_cleaned['log2FoldChange'] < -FOLD_CHANGE_THRESHOLD) & 
    (df_cleaned['padj'] < PADJ_THRESHOLD)
]
down_genes = set(downregulated.index.tolist())  # Extract gene names as set for fast lookup

print(f"Found {len(up_genes)} significantly upregulated genes.")
print(f"Found {len(down_genes)} significantly downregulated genes.")

# Load the gene BED file (from previous step; update path if needed)
BED_FILE = '/Users/berial/reference/hg38/MANE_gene.bed'
bed = pd.read_csv(BED_FILE, sep='\t', header=None)

# Rename columns for clarity (standard BED12)
bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 
               'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

# Filter BED for upregulated genes
up_bed = bed[bed['name'].isin(up_genes)].copy()
up_bed_file = 'upregulated_genes.bed'
up_bed.to_csv(up_bed_file, sep='\t', header=False, index=False)
print(f"Upregulated BED file saved to: {up_bed_file}")

# Filter BED for downregulated genes
down_bed = bed[bed['name'].isin(down_genes)].copy()
down_bed_file = 'downregulated_genes.bed'
down_bed.to_csv(down_bed_file, sep='\t', header=False, index=False)
print(f"Downregulated BED file saved to: {down_bed_file}")

# Optional: Save gene lists to CSV for reference
upregulated.to_csv('/Users/berial/reference/hg38/upregulated_genes.csv')
downregulated.to_csv('/Users/berial/reference/hg38/downregulated_genes.csv')
print("Gene lists saved to CSV files.")

# %%
import pyBigWig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# --- 1. 定义文件路径 ---
# 请确保这些文件存在于您的运行环境中
BIGWIG_FILE = "SRR14237124h2azchip.bw"
DOWN_BED_FILE = "downregulated_genes.bed"
UP_BED_FILE = "upregulated_genes.bed"
OUTPUT_PLOT_FILE = "ratio_boxplot.png"

# --- 2. 辅助函数：BED 文件解析 ---

def parse_bed_genes(filepath):
    """
    解析 BED 文件，确定每个基因的 TSS 坐标。
    假设 BED 文件至少包含 6 列 (chrom, start, end, name, score, strand)。
    TSS 定义：+ 链取 start 坐标 (0-based)，- 链取 end 坐标 (1-based)。
    """
    genes = []
    try:
        # 读取 BED 文件，假设为制表符或空格分隔
        df = pd.read_csv(filepath, sep='\s+', header=None, comment='#',
                         usecols=[0, 1, 2, 5], names=['chrom', 'start', 'end', 'strand'])
    except Exception as e:
        print(f"无法读取或解析文件 {filepath}: {e}")
        return []

    for _, row in df.iterrows():
        chrom, start, end, strand = row['chrom'], row['start'], row['end'], row['strand']

        # 确定 TSS 坐标 (0-based)
        if strand == '+':
            tss = start
        elif strand == '-':
            # BED start是0-based, end是1-based。-链的TSS是end。
            # 为了保持0-based (半开区间)，我们使用 end-1 作为TSS的实际位置。
            # 但在计算区域时，pyBigWig的stats是0-based的，start inclusive, end exclusive。
            # 通常将TSS定为基因体开始的碱基。
            # 对于 - 链，TSS通常是 end 坐标 (1-based) 之前的那个位置。
            tss = end - 1
        else:
            # 忽略没有明确链信息的基因
            continue

        genes.append({
            'chrom': chrom,
            'tss': tss,
            'strand': strand,
            'original_start': start,
            'original_end': end
        })

    print(f"文件 {filepath} 成功解析了 {len(genes)} 个基因。")
    return genes

# --- 3. 核心函数：信号提取和比值计算 ---

def calculate_gene_signal_ratios(bw_handle, genes):
    """
    对一组基因，计算指定区域的平均信号和比值。
    """
    results = []

    for gene in genes:
        chrom = gene['chrom']
        tss = gene['tss']
        strand = gene['strand']

        # 区域 1: TSS 上下游 1kb (TSS ± 1000 bp)
        # 始终以 TSS 为中心，取基因组坐标范围 [TSS - 1000, TSS + 1000]
        tss_1k_start = max(0, tss - 1000)
        tss_1k_end = tss + 1000

        # 区域 2: TSS 下游 200bp (相对于转录方向)
        if strand == '+':
            # +链: [TSS, TSS + 200]
            tss_200_start = tss
            tss_200_end = tss + 200
        elif strand == '-':
            # -链: [TSS - 200, TSS]
            tss_200_start = max(0, tss - 200)
            tss_200_end = tss
        else:
            # 忽略无链信息的基因
            continue

        try:
            # 获取 TSS ± 1kb 区域的平均信号
            signal_1k = bw_handle.stats(
                chrom, tss_1k_start, tss_1k_end, type="mean", numpy=True
            )[0]

            # 获取 TSS 下游 200bp 区域的平均信号
            signal_200 = bw_handle.stats(
                chrom, tss_200_start, tss_200_end, type="mean", numpy=True
            )[0]

            # 检查信号值是否有效 (非 NaN)
            if np.isnan(signal_1k) or np.isnan(signal_200):
                ratio = np.nan
            else:
                # 避免除以零：如果分母为零，比值也设为 NaN 或零
                if signal_1k == 0:
                    ratio = np.nan # 无法计算比值
                else:
                    ratio = signal_200 # / signal_1k

            results.append({
                'gene_set': 'N/A', # 占位符，稍后填充
                'signal_1k': signal_1k,
                'signal_200': signal_200,
                'ratio': ratio
            })

        except Exception as e:
            # 捕获 BigWig 读取错误 (如染色体名不匹配)
            # print(f"警告: 无法读取 {chrom}:{tss} 附近的信号: {e}")
            pass # 静默忽略无法读取的基因

    # 过滤掉无法计算比值的基因 (ratio 为 NaN)
    df = pd.DataFrame(results).dropna(subset=['ratio'])
    return df

# --- 4. 主执行流程 ---

def main():
    if not os.path.exists(BIGWIG_FILE):
        print(f"错误: BigWig 文件 {BIGWIG_FILE} 不存在。请检查路径。")
        return

    # 1. 初始化 BigWig 文件句柄
    try:
        bw = pyBigWig.open(BIGWIG_FILE)
    except Exception as e:
        print(f"错误: 无法打开 BigWig 文件 {BIGWIG_FILE}。请确保文件未损坏且安装了 pyBigWig。")
        return

    # 2. 解析 BED 文件
    down_genes = parse_bed_genes(DOWN_BED_FILE)
    up_genes = parse_bed_genes(UP_BED_FILE)

    if not down_genes or not up_genes:
        print("错误: 至少一个基因列表为空。请检查 BED 文件内容。")
        return

    # 3. 计算信号比值
    print("开始计算 downregulated 基因的信号比值...")
    df_down = calculate_gene_signal_ratios(bw, down_genes)
    df_down['gene_set'] = 'Downregulated'

    print("开始计算 upregulated 基因的信号比值...")
    df_up = calculate_gene_signal_ratios(bw, up_genes)
    df_up['gene_set'] = 'Upregulated'

    bw.close() # 关闭 BigWig 文件句柄

    # 4. 合并数据
    final_df = pd.concat([df_down, df_up], ignore_index=True)

    if final_df.empty:
        print("错误: 计算后没有有效的信号比值数据。")
        return

    # 5. 统计检验 (Mann-Whitney U 检验, 即 Wilcoxon Rank-Sum Test)
    # 适用于非正态分布的数据比较
    data_down = final_df[final_df['gene_set'] == 'Downregulated']['ratio']
    data_up = final_df[final_df['gene_set'] == 'Upregulated']['ratio']

    # 确保每个组都有足够的数据进行检验
    if len(data_down) < 2 or len(data_up) < 2:
        print("警告: 至少有一个基因集的数据点过少，无法进行统计检验。")
        p_value = 1.0
        p_label = "N/A"
    else:
        # 使用 Mann-Whitney U 检验
        stat, p_value = stats.mannwhitneyu(data_down, data_up, alternative='two-sided')

        if p_value < 0.0001:
            p_label = f"P < 0.0001"
        else:
            p_label = f"P = {p_value:.4f}"

    # 6. 绘图 (Boxplot)
    plt.figure(figsize=(6, 8))
    sns.boxplot(x='gene_set', y='ratio', data=final_df,
                order=['Downregulated', 'Upregulated'],
                palette=['#7EC8E3', '#F28585'],  # 蓝色系和红色系
                linewidth=1.5, fliersize=3)
    
    # 添加散点图层 (Stripplot) 以显示数据点分布
    sns.stripplot(x='gene_set', y='ratio', data=final_df,
                  order=['Downregulated', 'Upregulated'],
                  color='black', size=2, jitter=True, alpha=0.3)

    plt.title(f'{os.path.basename(BIGWIG_FILE).replace(".bw", "")}\nSignal Ratio (TSS 200bp / TSS 1kb)', fontsize=14)
    plt.xlabel('Gene Set', fontsize=12)
    plt.ylabel('Signal Ratio (TSS+200bp Avg / TSS±1kb Avg)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    # 标注 P 值
    # 计算连接线的 y 坐标
    y_max = final_df['ratio'].max() * 1.05
    h = final_df['ratio'].max() * 0.05
    
    # 确保 P 值标签在图的合适位置
    if not np.isnan(y_max) and y_max > 0 and len(data_down) >= 2 and len(data_up) >= 2:
        # 绘制横线
        plt.plot([0, 1], [y_max, y_max], lw=1.5, c='k')
        # 标注 P 值
        plt.text(0.5, y_max + h/2, p_label, ha='center', va='bottom', fontsize=12,
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_FILE, dpi=300)
    print(f"\n成功生成 Boxplot 图像并保存到 {OUTPUT_PLOT_FILE}")
    print(f"Downregulated 组平均比值: {data_down.mean():.4f} (N={len(data_down)})")
    print(f"Upregulated 组平均比值: {data_up.mean():.4f} (N={len(data_up)})")
    print(f"统计检验结果: P 值 = {p_value:.4g}")

if __name__ == '__main__':
    main()
