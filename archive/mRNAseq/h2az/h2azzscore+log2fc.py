#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 15:23:27 2025

@author: berial
"""

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
import os 
import pandas as pd
os.chdir("/Users/berial/ljy/lailab2025/pby_med12/mRNAseq/h2az/")
import gseapy as gp
import sys
import numpy as np
from gseapy import dotplot
import matplotlib.pyplot as plt
# --- Configuration ---
# Set the name of your input file
# Specify the gene set libraries to use for analysis.
# These are stable and supported versions for mouse.
GENE_SET_LIBRARIES = {
    'KEGG': 'KEGG_2019_Mouse',
    'HALLMARK': '/Users/berial/reference/mm39//mh.all.v2025.1.Mm.symbols.gmt' 
}


print("Data loaded successfully.")
print("Preprocessing data...")

info = pd.read_table("/Users/berial/reference/mm39/RefSeq_all.txt")
info = info[['name', 'name2']]

# The file has a header like '""', so we rename the index to 'Gene_Name' for clarity.
df = pd.read_table("mousetss500_h2az.txt", header=None, usecols=[3, 6])
df.columns=['Gene_Name','log2FoldChange']
df2 = pd.read_table("mousetss500_input.txt", header=None, usecols=[3, 6])
df2.columns=['Gene_Name','log2FoldChange']
info_unique = info.drop_duplicates(subset=['name'], keep='first')

name_to_gene_map = info_unique.set_index('name')['name2']
df['Gene_Name'] = df['Gene_Name'].map(name_to_gene_map).combine_first(df['Gene_Name'])

df_cleaned = df.dropna(subset=['log2FoldChange'])

# %%
from scipy.stats import zscore


df3 = pd.read_csv("../results.txt", index_col=0)
df3["Gene_Name"] = df3.index
df3_cleaned = df3.dropna(subset=['log2FoldChange'])
df3_cleaned = df3_cleaned[["Gene_Name", "log2FoldChange"]]
df_cleaned  = df_cleaned.merge(df3_cleaned, on="Gene_Name")
df_cleaned["zscore_log2FoldChange_x"] = zscore(df_cleaned["log2FoldChange_x"])

df_cleaned["log2FoldChange"]  = df_cleaned["zscore_log2FoldChange_x"] + df_cleaned['log2FoldChange_y']
# %%






df_cleaned.loc[:, 'combined_rank'] = df_cleaned['log2FoldChange']

# Rank the genes based on this new combined metric in descending order.
ranked_genes = df_cleaned.sort_values(by='combined_rank', ascending=False)

print("Data preprocessed and ranked by 'combined_rank' (log2FC and padj).")
ranked_series = ranked_genes.set_index('Gene_Name')['combined_rank']

# Iterate through each gene set library and perform GSEA
for name, library in GENE_SET_LIBRARIES.items():
    output_dir = f'gsea_results_{name}'
    print(f"\n--- Running GSEA for '{name}' using library '{library}' ---")
    
    # Perform the GSEA analysis using the prerank function.
    prerank_results = gp.prerank(
        rnk=ranked_series,
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
    pi3k_terms = [t for t in terms if 'Chemical carcinogenesis' in t]
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
    
    
# %%
    
os.chdir("/Users/berial/ljy/lailab2025/pby_med12/mRNAseq/h2az/gsea_results_KEGG/")
df = pd.read_csv("gseapy.gene_set.prerank.report.csv")
df

os.mkdir("gsea_bed")

df = df[["Term", "Lead_genes"]]

bed = pd.read_table("~/reference/mm39/longest_nm_gene_transcripts.bed", header=None)
bed = bed[~bed[0].str.contains("_", na=False)]
bed

def sep(term, df=df, bed=bed):
    gene_series = df[df['Term'] == term]['Lead_genes'].str.split(";").tolist()[0]
    bed = bed[bed[3].isin(gene_series)]
    bed.to_csv(f"gsea_bed/{term}.bed", sep="\t", header=False, index=False)

sep("Valine, leucine and isoleucine degradation")
sep("Peroxisome")
sep("Oxidative phosphorylation")
sep("Propanoate metabolism")
sep("ECM-receptor interaction")

sep("Hippo signaling pathway")