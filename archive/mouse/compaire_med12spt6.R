library(tidyverse)
homo_df = read.csv("~/ljy/lailab2025/mart_export.txt", sep="\t")
homo_df = homo_df %>% filter(!is.na(Mouse.gene.name))%>% filter(Mouse.gene.name != "") %>% filter(Gene.name != "") %>% distinct(Mouse.gene.name, .keep_all = TRUE) %>% distinct(Gene.name, .keep_all = TRUE)
mouse_gene_list = homo_df$Mouse.gene.name
fc=1.5
pval_threshold <- 0.05
logfc_threshold <- log2(2)
significance_fc <- 2
significance_p <- 0.001

## Deseq2 for RNA-seq analysis
options(download.file.method = 'libcurl')
options(url.method = 'libcurl')
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(VennDiagram)
library(clusterProfiler)
library(enrichplot)
library("RColorBrewer")

loading_data = function(treatment){
  control = "CTRL"
  data <- read.table('featureCounts.txt', header = T, quote = '\t', skip = 1)
  data = data %>% filter(Geneid %in% mouse_gene_list)
  database <- read.csv('database.csv', header = T, check.names = FALSE)
  ## Arrange the database, same as the range of featureCounts
  database_rep <- database %>% group_by(condition) %>% 
    mutate(su = length(name)) %>% 
    mutate(rep = 1:unique(su)) %>% 
    dplyr::select(-su) %>% 
    mutate(condition_rep = paste(condition, "rep", rep, sep = "_")) %>% 
    mutate(name = condition_rep)
  sampleNames <- as.factor(database_rep$condition_rep)
  countData <- as.matrix(data[, 7:ncol(data)])
  colnames(countData) <- sampleNames
  rownames(countData) <- data$Geneid
  rownames(database) <- sampleNames
  database$condition <- as.factor(database$condition)
  ## DESeq2
  dds <- DESeqDataSetFromMatrix(countData, colData = database, design = ~condition)
  x <- data$Length
  mcols(dds)$basepairs <- x
  dds <- DESeq(dds)
  dds <- dds[ rowMaxs(counts(dds))>1,]
  
  # Variance stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  vst <- assay(vsd)
  
  # Sample distance heatmap
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- database_rep$condition_rep
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  p <- pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
  ggsave("distance.pdf", plot = p)
  condition = as.vector(unique(database$condition))
  res <- results(dds, contrast = c('condition', treatment, "CTRL"))
  normalized_counts <- counts(dds, normalized = TRUE)
  filtered_cols <- colnames(normalized_counts)[grep(paste0("^", control, "|^", treatment), colnames(normalized_counts))]
  normalized_counts <- normalized_counts[, filtered_cols]
  return(list(res, dds, vst, normalized_counts))
}

setwd("~/ljy/lailab2025/pby_med12/mouse/spt6_oe//")
spt6 = loading_data("OESPT6") 
spt6.res_raw = spt6[1]
spt6.res = spt6.res_raw  %>% as.data.frame() %>% rownames_to_column('gene')
spt6.sig.up = spt6.res %>% filter(log2FoldChange >logfc_threshold, padj<0.05)
spt6.sig.down = spt6.res %>% filter(log2FoldChange < -logfc_threshold, padj<0.05)
setwd("~/ljy/lailab2025/pby_med12/mouse/med12_oe//")
med12 = loading_data("OEMED12") 
med12.res_raw = med12[1]
med12.res = med12.res_raw  %>% as.data.frame() %>% rownames_to_column('gene')
med12.sig.up = med12.res %>% filter(log2FoldChange >logfc_threshold, padj<0.05)
med12.sig.down = med12.res %>% filter(log2FoldChange < -logfc_threshold, padj<0.05)
setwd("~/ljy/lailab2025/pby_med12/mouse/")
intersect(med12.sig.down$gene, spt6.sig.down$gene)

## venn plot 
library(eulerr)
library(tidyverse)
venn <- function(ls, name){
  a <- plot(euler(ls, shape = "ellipse"), quantities = TRUE, main = name)
  ggsave(paste0(name, ".pdf"), plot = a, width = 5, height = 5)
}
venn(list("SPT6.sig.up"=spt6.sig.up$gene,  "MED12.sig.up"=med12.sig.up$gene), "sig.up.venn.pdf")
venn(list("SPT6.sig.down"=spt6.sig.down$gene,  "MED12.sig.down"=med12.sig.down$gene), "sig.down.venn.pdf")

intersect_up_gene = intersect(spt6.sig.up$gene, med12.sig.up$gene)
intersect_down_gene = intersect(spt6.sig.down$gene, med12.sig.down$gene)

# trans into human gene
library(pheatmap)

# Function to map mouse genes to human genes
map_mouse_to_human <- function(mouse_genes, homo_df) {
  # Create a data frame with mouse_genes in original order
  mouse_df <- data.frame(Mouse.gene.name = mouse_genes, order = seq_along(mouse_genes))
  
  # Merge with homo_df to get corresponding human genes, preserving order
  mapped <- merge(mouse_df, homo_df, by = "Mouse.gene.name", all.x = TRUE)
  
  # Sort by original order and extract human gene names
  mapped_genes <- mapped %>% 
    arrange(order) %>% 
    pull(Gene.name)
  
  return(mapped_genes)
}

# Generate heatmap
plot_heatmap <- function(intersect_genes, normalized_counts, output_file) {
  # Filter normalized counts to include only intersected genes
  heatmap_data <- normalized_counts[rownames(normalized_counts) %in% intersect_genes, ]
  
  # Scale the data for better visualization
  heatmap_data <- t(scale(t(heatmap_data)))
  
  # Plot heatmap
  p <- pheatmap(heatmap_data, 
                color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), 
                cluster_rows = TRUE, 
                cluster_cols = FALSE, 
                show_rownames = FALSE,  scale = "row",
                main = "Heatmap of Intersected Genes")
  ggsave(output_file, plot = p, width = 8, height = 6)
}

info <- read.csv("~/reference/hg38/RefSeq_all.txt",  sep = "\t") %>% filter(!str_detect(chrom, "_"))
info_id <- info %>% dplyr::select(ID = name, name2) 
colnames(info_id) = c('trans', "gene")
# Map intersected up and down genes to human genes
intersect_up_human <- map_mouse_to_human(intersect_up_gene, homo_df)
intersect_down_human <- map_mouse_to_human(intersect_down_gene, homo_df)




plot_volcano <- function(res, output_file, logfc_cutoff = 0, pval_cutoff = 0.05) {
  # Add significance labels
  res <- res %>%
    mutate(Significance = case_when(
      padj < pval_cutoff & log2FoldChange > logfc_cutoff ~ "Upregulated",
      padj < pval_cutoff & log2FoldChange < -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  
  # Count significant genes
  up_count <- sum(res$Significance == "Upregulated")
  down_count <- sum(res$Significance == "Downregulated")
  
  # Generate volcano plot
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    # geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", linewidth = 0.6) +
    theme_minimal() +
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      subtitle = paste("Upregulated:", up_count, " | Downregulated:", down_count),
    ) 
  
  # Save plot
  ggsave(output_file, plot = p, width = 8, height = 6)
}
# Plot volcano for SPT6
plot_volcano(spt6.res, "spt6_volcano_annotated.pdf")
# Plot volcano for MED12
plot_volcano(med12.res, "med12_volcano_annotated.pdf")



library(dplyr)
library(ggplot2)

# 设置 fold change 阈值
fc_threshold <- log2(2)

# 提取 log2 fold changes 和 gene 名称
spt6_logfc <- spt6.res_raw %>% as.data.frame() %>%
  rownames_to_column("gene") %>% 
  dplyr::select(gene, log2FoldChange, padj, baseMean)

med12_logfc <- med12.res_raw %>% as.data.frame() %>%
  rownames_to_column("gene") %>% 
  dplyr::select(gene, log2FoldChange, padj)

# 合并
correlation_data <- inner_join(spt6_logfc, med12_logfc, by = "gene", suffix = c("_spt6", "_med12"))

# 分类颜色（根据 fold change）
correlation_data <- correlation_data %>%
  mutate(
    color_group = case_when(
      log2FoldChange_spt6 > fc_threshold & log2FoldChange_med12 > fc_threshold & padj_med12 < 0.05 & padj_spt6 < 0.05 ~ "Up-Up",
      log2FoldChange_spt6 < -fc_threshold & log2FoldChange_med12 < -fc_threshold & padj_med12 < 0.05 & padj_spt6 < 0.05 ~ "Down-Down",
      TRUE ~ "Other"
    )
  )

# 相关性检验
cor_test <- cor.test(correlation_data$log2FoldChange_spt6, correlation_data$log2FoldChange_med12, use = "complete.obs")
cor_value <- cor_test$estimate
p_value <- cor_test$p.value
formatted_p_value <- ifelse(p_value == 0, "< 2.2e-16", sprintf("%.3f", p_value))

# 绘图
# 统计数量
up_up_count <- sum(correlation_data$color_group == "Up-Up")
down_down_count <- sum(correlation_data$color_group == "Down-Down")

correlation_data$color_group <- factor(correlation_data$color_group, 
                                       levels = c("Other", "Up-Up", "Down-Down"))

p <- ggplot(correlation_data, aes(x = log2FoldChange_spt6, 
                                  y = log2FoldChange_med12, 
                                  color = color_group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  scale_color_manual(values = c("Up-Up" = "red", 
                                "Down-Down" = "blue", 
                                "Other" = "gray")) +
  labs(
    title = "Correlation of Log2 Fold Changes",
    subtitle = paste0("Pearson correlation: ", round(cor_value, 2),
                      "\nP-value: ", formatted_p_value,
                      "\nUp-Up: ", up_up_count,
                      " genes, Down-Down: ", down_down_count, " genes"),
    x = "Log2 Fold Change (SPT6)",
    y = "Log2 Fold Change (MED12)",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )



p

# Save plot
ggsave("log2fc_correlation_scatter_with_pvalue.pdf", plot = p, width = 8, height = 6)
# Display plot
print(p)



###################### Loading H2AZ signal close to TSS
h2az = read_table("output_with_density.bed",col_names = FALSE)
minfo <- read.csv("~/reference/mm39/RefSeq_all.txt",  sep = "\t") %>% filter(!str_detect(chrom, "_"))
minfo_id <- minfo %>% dplyr::select(ID = name, name2) 
# write all intersect up/down transcripts name 
intersect_up_trans = minfo_id %>% filter(name2 %in% intersect_up_gene) 
intersect_up_trans = intersect_up_trans$ID
intersect_down_trans = minfo_id %>% filter(name2 %in% intersect_down_gene) 
intersect_down_trans = intersect_down_trans$ID
bed = read.table("tsstoh2azp1.bed")
write.table(bed %>% filter(V4 %in% intersect_down_trans), file = 'intersectdown_tsstop1.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(bed %>% filter(V4 %in% intersect_up_trans), file = 'intersectup_tsstop1.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
upupgene = correlation_data %>% filter(color_group == 'Up-Up')
upupgene = minfo_id %>% filter(name2 %in% upupgene$gene)
write.csv(file='upup.csv',upupgene)
downdowngene = correlation_data %>% filter(color_group == 'Down-Down')
downdowngene = minfo_id %>% filter(name2 %in% downdowngene$gene)
write.csv(file='downdown.csv',downdowngene)


write.table(bed %>% filter(V4 %in% upupgene$ID), file = 'upup_tsstop1.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(bed %>% filter(V4 %in% downdowngene$ID), file = 'downdown_tsstop1.bed', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

colnames(minfo_id) = c('trans', "gene")
h2az$length = h2az$X3 - h2az$X2
h2az = h2az %>% dplyr::select(trans=X4, h2az=X7, length) %>% inner_join(minfo_id, by='trans')

correlation_data = correlation_data %>% inner_join(h2az, by="gene")
library(ggpubr)

library(ggpubr)


# 创建 h2az_group 分组（分为三组：Low, Medium, High）
correlation_data <- correlation_data %>%
  mutate(h2az_group = cut(h2az,
                          breaks = quantile(h2az, probs = seq(0, 1, length.out = 4), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = c("Low", "Medium", "High")))
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)

boxplot_by_group <- function(data, gene_column, group_column = "h2az") {
  gene_column <- rlang::ensym(gene_column)
  group_column <- rlang::ensym(group_column)
  
  # 分组变量
  data <- data %>%
    mutate(group = cut(!!group_column,
                       breaks = quantile(!!group_column, probs = seq(0, 1, length.out = 4), na.rm = TRUE),
                       include.lowest = TRUE)) %>%
    filter(!is.na(!!gene_column))
  
  # 计算 ymin/ymax：所有组中 boxplot.stats 的 min/max 非异常值
  group_stats <- data %>%
    group_by(group) %>%
    summarise(
      ymin = boxplot.stats(!!gene_column)$stats[1],
      ymax = boxplot.stats(!!gene_column)$stats[5]
    )
  
  ymin <- min(group_stats$ymin, na.rm = TRUE)
  ymax <- max(group_stats$ymax, na.rm = TRUE)
  
  # 自动生成比较组合
  comparisons <- combn(levels(data$group), 2, simplify = FALSE)
  
  # 统计检验
  stat.test <- data %>%
    wilcox_test(as.formula(paste0(rlang::as_string(gene_column), " ~ group"))) %>%
    adjust_pvalue() %>%
    add_significance("p.adj") %>%
    add_x_position(x = "group")
  
  # 设置标注位置在 boxplot 范围之上
  stat.test$y.position <- seq(ymax * 1.05, by = (ymax - ymin) * 0.15, length.out = nrow(stat.test))
  
  # 绘图
  p <- ggplot(data, aes(x = group, y = !!gene_column)) +
    geom_boxplot(fill = "lightblue", outlier.shape = NA) +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0) +
    labs(
      title = paste0("log2FoldChange (", rlang::as_string(gene_column), ") vs ", rlang::as_string(group_column), " Expression Group"),
      x = paste(rlang::as_string(group_column), "Expression Group"),
      y = paste("log2FoldChange (", rlang::as_string(gene_column), ")")
    ) +
    coord_cartesian(ylim = c(ymin, max(stat.test$y.position) * 1.05)) +
    theme_minimal(base_size = 14)
  
  return(p)
}
pdf('boxplot.pdf')
boxplot_by_group(correlation_data %>% filter(log2FoldChange_spt6 >fc | log2FoldChange_spt6 < -fc), gene_column = log2FoldChange_spt6)
boxplot_by_group(correlation_data %>% filter(log2FoldChange_med12 >0 & log2FoldChange_med12 >0), gene_column = log2FoldChange_med12)
boxplot_by_group(correlation_data %>% mutate(lfc = log2FoldChange_spt6 + log2FoldChange_med12 ) %>% filter(padj_med12 <0.05 & padj_spt6 <0.05) , gene_column = lfc)
boxplot_by_group(correlation_data %>% mutate(p=-log10(padj_spt6)) , gene_column = p)
boxplot_by_group(correlation_data %>% mutate(expression=baseMean), gene_column = expression)


boxplot_by_group(correlation_data %>% filter(log2FoldChange_med12 <0 & log2FoldChange_med12 <0), gene_column = log2FoldChange_spt6)
boxplot_by_group(correlation_data %>% filter(log2FoldChange_med12 <0 & log2FoldChange_med12 <0), gene_column = log2FoldChange_med12)


correlation_data <- correlation_data %>%
  mutate(h2az_group = cut(log2FoldChange_med12,
                          breaks = quantile(log2FoldChange_med12, probs = seq(0, 1, length.out = 4), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = c("Low", "Medium", "High")))

ggplot(correlation_data, aes(x = h2az_group, y = h2az)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "log2FoldChange (Spt6) vs H2AZ Quartiles",
       x = "H2AZ Expression Quartile",
       y = "log2FoldChange (Spt6)") +
  theme_minimal()

correlation_data$change = ifelse(correlation_data$gene %in% med12.sig.down$gene, "down", ifelse(correlation_data$gene %in% med12.sig.up$gene, "up", 'others'))

ggplot(correlation_data, aes(x = change, y = h2az)) +
  geom_boxplot(fill = "lightblue") + theme_minimal()

library(ggpubr)

# 绘图 + 添加 p-value 标注
ggplot(correlation_data, aes(x = change, y = h2az)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("down", "others"), c("up", "down")),
                     label = "p.format")  # 或用 "p.format" 显示具体数值



intersect_up_gene = intersect(spt6.sig.up$gene, med12.sig.up$gene)
intersect_down_gene = intersect(spt6.sig.down$gene, med12.sig.down$gene)


library(ggplot2)
library(ggpubr)

# 假设 spt6.sig.up / med12.sig.up / spt6.sig.down / med12.sig.down 都包含 gene 列

# 分类基因
spt6_up <- setdiff(spt6.sig.up$gene, med12.sig.up$gene)
med12_up <- setdiff(med12.sig.up$gene, spt6.sig.up$gene)
both_up <- intersect(spt6.sig.up$gene, med12.sig.up$gene)

spt6_down <- setdiff(spt6.sig.down$gene, med12.sig.down$gene)
med12_down <- setdiff(med12.sig.down$gene, spt6.sig.down$gene)
both_down <- intersect(spt6.sig.down$gene, med12.sig.down$gene)

# 创建新分类列
correlation_data$change <- "others"
correlation_data$change[correlation_data$gene %in% spt6_up] <- "spt6_up"
correlation_data$change[correlation_data$gene %in% med12_up] <- "med12_up"
correlation_data$change[correlation_data$gene %in% both_up] <- "both_up"
correlation_data$change[correlation_data$gene %in% spt6_down] <- "spt6_down"
correlation_data$change[correlation_data$gene %in% med12_down] <- "med12_down"
correlation_data$change[correlation_data$gene %in% both_down] <- "both_down"

# 因子顺序（可按你需要调整）
correlation_data$change <- factor(
  correlation_data$change,
  levels = c("spt6_down", "med12_down", "both_down", "others", "spt6_up", "med12_up", "both_up")
)

# 设置比较组（例如下面是和“others”比较）
my_comparisons <- list(
  c("spt6_down", "both_down"),
  c("med12_down", "both_down"),
  c("both_down", "both_down"),
  c("spt6_up", "both_down"),
  c("med12_up", "both_down"),
  c("both_up", "both_down")
)

# 绘图
ggplot(correlation_data %>% filter(change != 'others'), aes(x = change, y = h2az)) +
  geom_boxplot(aes(fill = change), outlier.shape = NA, alpha = 0.6) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "spt6_down" = "#d73027",
    "med12_down" = "#fc8d59",
    "both_down" = "#fee090",
    "others" = "grey80",
    "spt6_up" = "#91bfdb",
    "med12_up" = "#4575b4",
    "both_up" = "#313695"
  )) +
  labs(title = "H2A.Z Signal by Gene Regulation Group", x = NULL, y = "H2A.Z")


fc=log2(2)
correlation_data <- correlation_data %>%
  mutate(change = ifelse(log2FoldChange_spt6 > fc & log2FoldChange_med12 > fc, 
                         'up', 
                         ifelse(log2FoldChange_spt6 < -fc & log2FoldChange_med12 < -fc, 
                                'down', 
                                'others')))


ggplot(correlation_data %>% filter(change != 'others'), 
       aes(x = change, y = h2az)) +
  geom_boxplot(aes(fill = change), outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(
    comparisons = list(c("up", "down")), 
    method = "wilcox.test", 
    label = "p.format"
  ) +
  scale_fill_manual(
    values = c("up" = "red", "down" = "blue")
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Gene Change Group", y = "H2AZ Expression")

dev.off()

library(dplyr)

# 筛选条件：H2AZ 组为 High，且两个log2FC都大于0
gene_up_high_h2az <- correlation_data %>%
  filter(h2az_group == "High",
         log2FoldChange_spt6 > fc,
         log2FoldChange_med12 > fc) %>%
  pull(gene)
library(org.Mm.eg.db)
library(clusterProfiler)

# 转换 Symbol 到 Entrez ID
gene_entrez <- bitr(gene_up_high_h2az,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)


# for down
gene_up_high_h2az <- correlation_data %>%
  filter( # h2az_group == "High",
         log2FoldChange_spt6 < -fc,
         log2FoldChange_med12 < -fc) %>%
  pull(gene)
library(org.Mm.eg.db)
library(clusterProfiler)

# 转换 Symbol 到 Entrez ID
gene_entrez <- bitr(gene_up_high_h2az,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)


library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# 提取 ENTREZ ID 列表
entrez_ids <- gene_entrez$ENTREZID

# ======== GO BP 分析 ========
go_result <- enrichGO(gene         = entrez_ids,
                      OrgDb        = org.Mm.eg.db,
                      keyType      = "ENTREZID",
                      ont          = "ALL",         # Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

# GO 可视化
barplot(go_result, showCategory = 5, title = "GO Biological Process")
dotplot(go_result, showCategory = 5, title = "GO Biological Process")


# ======== KEGG 分析 ========
# KEGG 要求的物种代码：小鼠是 "mmu"
kegg_result <- enrichKEGG(gene         = entrez_ids,
                          organism     = 'mmu',       # Mus musculus
                          keyType      = "kegg",       # 虽然参数存在，实际上 enrichKEGG 默认用 Entrez ID
                          pvalueCutoff = 0.05)

# KEGG 可视化
barplot(kegg_result, showCategory = 5, title = "KEGG Pathway")
dotplot(kegg_result, showCategory = 5, title = "KEGG Pathway")

