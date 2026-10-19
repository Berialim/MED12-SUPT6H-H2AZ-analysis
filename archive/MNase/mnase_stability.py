"""MNase positional comparison converted from the source notebook.

Usage: python mnase_stability.py /path/to/working_directory
Requires med12tss_to_plus1_peak.bed and spt6tss_to_plus1_peak.bed.
Preserves the original analysis, including the known two-bin BED export.
See docs/KNOWN_ISSUES.md before using exported groups.
"""
import os 
import sys
import pandas as pd
import seaborn as sns 
# Run against a working copy of the two input BED files.
if len(sys.argv) != 2:
    raise SystemExit("Usage: python mnase_stability.py /path/to/working_directory")
os.chdir(sys.argv[1])
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # 使 PDF 中字体可编辑
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"

med12 = pd.read_table("med12tss_to_plus1_peak.bed", header=None)
med12_bed = med12.copy()
med12_bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
med12['length'] = med12[2] - med12[1]
med12.rename(columns={3:'name'}, inplace=True)
med12 = med12[['name', 'length']]

spt6 = pd.read_table("spt6tss_to_plus1_peak.bed", header=None)
spt6_bed = spt6.copy()
spt6_bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
spt6['length'] = spt6[2] - spt6[1]
spt6.rename(columns={3:'name'}, inplace=True)
spt6 = spt6[['name', 'length']]

df = med12.merge(spt6, on='name')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import gaussian_kde

# 示例数据（可替换为 df['length_x'], df['length_y']）
x = df['length_x']
y = df['length_y']

# 计算相关性
corr, pval = spearmanr(x, y)

# 计算密度用于颜色
xy = np.vstack([x, y])
z = gaussian_kde(xy)(xy)

# 按密度排序（高密度画在上层）
idx = z.argsort()
x_plot, y_plot, z_plot = x.values[idx], y.values[idx], z[idx]

# 绘图
plt.figure(figsize=(6, 5))
sc = plt.scatter(x_plot, y_plot, c=z_plot, cmap='viridis', s=10)
plt.colorbar(sc, label='Density')
plt.xlabel('length_x')
plt.ylabel('length_y')
plt.title(f"Spearman r = {corr:.2f}, p = {pval:.2e}")
# sns.regplot(x=x, y=y, scatter=False, color='red', line_kws={'linewidth': 1.2})
plt.tight_layout()
plt.savefig("p1tss_distance_densityscatterplot.pdf")
plt.show()


# Q1 (smallest difference) is stable; Q4 (largest difference) is unstable.
df['length_difference'] = (df['length_x'] - df['length_y']).abs()
df['category'] = pd.qcut(
    df['length_difference'],
    q=4,
    labels=['Q1: Smallest Diff', 'Q2: Small Diff', 'Q3: Medium Diff', 'Q4: Largest Diff']
)

print("\nDataFrame with difference quantiles (first 5 rows):")
print(df.head())
print("-" * 30)
print("\nCounts by category:")
print(df['category'].value_counts().sort_index())


# Set plot style
sns.set_style("whitegrid")
plt.figure(figsize=(5, 4))
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # Keep PDF text editable
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"

# Q1 is stable (orange), Q4 is unstable (blue), and Q2/Q3 are transition colors.
scatterplot = sns.scatterplot(
    data=df,
    x='length_x',
    y='length_y',
    hue='category',
    hue_order=['Q1: Smallest Diff', 'Q2: Small Diff', 'Q3: Medium Diff', 'Q4: Largest Diff'],
    palette={
        'Q1: Smallest Diff': '#F39800',
        'Q2: Small Diff': '#F9C876',
        'Q3: Medium Diff': '#7CAAC6',
        'Q4: Largest Diff': '#2171a9',
    },
    s=10,
    alpha=1,
    edgecolor=None
)

# Add the diagonal y=x as a perfect-match reference.
min_val = min(df['length_x'].min(), df['length_y'].min())
max_val = max(df['length_x'].max(), df['length_y'].max())
plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.8, label='y = x (Perfect Match)')

plt.title('Distribution of Lengths by Difference Quantile Category')
plt.xlabel('Length X')
plt.ylabel('Length Y')
plt.legend(title='Difference Category')
plt.grid(False)
plt.tight_layout()
plt.savefig("p1tss_distance_scatterplot.pdf")
plt.show()


import pandas as pd


df['quartile'] = pd.qcut(df['length_difference'], 2, labels=['Q1', 'Q2']) # , 'Q3', 'Q4'])


# 3. 合并 df 与 med12_bed，按 name 匹配
merged = pd.merge(med12_bed, df[['name', 'quartile']], on='name', how='inner')

# 4. 拆分并输出为 4 个 BED 文件
for q in ['Q1', 'Q2']: #, 'Q3', 'Q4']:
    subset = merged[merged['quartile'] == q]
    bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    subset[bed_cols].to_csv(f"med12_{q}.bed", sep='\t', header=False, index=False)
    print(f"Saved {len(subset)} entries to med12_{q}.bed")


import pandas as pd


df['quartile'] = pd.qcut(df['length_difference'], 2, labels=['Q1', 'Q2']) # , 'Q3', 'Q4'])


# 3. 合并 df 与 spt6_bed，按 name 匹配
merged = pd.merge(spt6_bed, df[['name', 'quartile']], on='name', how='inner')

# 4. 拆分并输出为 4 个 BED 文件
for q in ['Q1', 'Q2']: #, 'Q3', 'Q4']:
    subset = merged[merged['quartile'] == q]
    bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    subset[bed_cols].to_csv(f"spt6_{q}.bed", sep='\t', header=False, index=False)
    print(f"Saved {len(subset)} entries to spt6_{q}.bed")
