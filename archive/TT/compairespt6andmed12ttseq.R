############# COMPAIR KD TT-SEQ AND OE MRNA-SEQ
library(tidyverse)
library(ggpubr)
homo_df = read.csv("~/ljy/lailab2025/mart_export.txt", sep="\t")
homo_df = homo_df %>% filter(!is.na(Mouse.gene.name))%>% filter(Mouse.gene.name != "") %>% filter(Gene.name != "") %>% distinct(Mouse.gene.name, .keep_all = TRUE) %>% distinct(Gene.name, .keep_all = TRUE)
mouse_gene_list = homo_df$Mouse.gene.name

pval_threshold <- 0.05
logfc_threshold <- 0
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
  data <- read.table("featureCounts.txt", header = T, quote = '\t', skip = 1)
  countData <-as.matrix(data[,7:ncol(data)])
  database <- read.csv('database.csv',header = T,check.names = FALSE)
  database_rep <- database %>% group_by(condition) %>% mutate(su = length(name)) %>% mutate(rep = 1:unique(su)) %>% 
    dplyr::select(-su) %>% mutate(condition_rep = paste(condition,"rep",rep,sep = "_")) %>% mutate(name = condition_rep)
  sampleNmaes <- as.factor(database_rep$condition_rep)
  colnames(countData) <- sampleNmaes
  rownames(countData) <- data$Geneid
  database$condition <- as.factor(database$condition)
  # database$condition <- relevel(database$condition,ref=control)
  ##DESeq2
  dds <- DESeqDataSetFromMatrix(countData, colData = database, design = ~condition)
  x <- data$Length
  mcols(dds)$basepairs <- x
  # size factor = 1/scale factor
  sf_df = read.csv("scalingfactors.txt", sep = "\t")
  sf_df$name = str_replace_all(sf_df$sample, ".bam", "") %>% as.numeric()
  sf_df = left_join(database, sf_df, by="name")
  sizefactor = 1 / sf_df$scalingFactor
  print("sizefactor")
  # normalized with spike-in's size factor
  dds = estimateSizeFactors(dds)
  old_size = sizeFactors(dds)
  # further normalize new size factor to the same level of old size factor (same mean)
  sizefactor = sizefactor * (mean(old_size) / mean(sizefactor))
  sizeFactors(dds) = sizefactor
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  # dds = DESeq(dds)
  dds <- dds[ rowMaxs(counts(dds))> 100,]
  
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
  # ggsave("distance.pdf", plot = p)
  condition = as.vector(unique(database$condition))
  res <- results(dds, contrast = c('condition', treatment, "CTRL"))
  normalized_counts <- counts(dds, normalized = TRUE)
  filtered_cols <- colnames(normalized_counts)[grep(paste0("^", control, "|^", treatment), colnames(normalized_counts))]
  normalized_counts <- normalized_counts[, filtered_cols]
  return(list(res, dds, vst, normalized_counts))
}
fc=0
padj_cf=0.05
setwd("~/ljy/lailab2025/pby_med12/TT/spt6/")
spt6 = loading_data("SPT6KD2h") 
spt6.res_raw = spt6[1]
spt6.res = spt6.res_raw  %>% as.data.frame() %>% rownames_to_column('gene')
spt6.sig.up = spt6.res %>% filter(log2FoldChange >fc, padj<padj_cf)
spt6.sig.down = spt6.res %>% filter(log2FoldChange < -0, padj<padj_cf)
setwd("~/ljy/lailab2025/pby_med12/TT/med12/")
med12 = loading_data("MED12KD") 
med12.res_raw = med12[1]
med12.res = med12.res_raw  %>% as.data.frame() %>% rownames_to_column('gene')
med12.sig.up = med12.res %>% filter(log2FoldChange >fc, padj<padj_cf)
med12.sig.down = med12.res %>% filter(log2FoldChange < -1, padj<padj_cf)

intersect(med12.sig.down$gene, spt6.sig.down$gene)





info <- read.csv("~/reference/hg38/RefSeq_all.txt",  sep = "\t") %>% filter(!str_detect(chrom, "_"))
info_id <- info %>% dplyr::select(ID = name, name2) 


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



# Extract log2 fold changes and gene names
spt6_logfc <- spt6.res_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::select(gene, log2FoldChange)
med12_logfc <- med12.res_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::select(gene, log2FoldChange)
# Merge the two datasets by gene
correlation_data <- inner_join(spt6_logfc, med12_logfc, by = "gene", suffix = c("_spt6", "_med12"))
# Perform correlation test
cor_test <- cor.test(correlation_data$log2FoldChange_spt6, correlation_data$log2FoldChange_med12, use = "complete.obs")
cor_value <- cor_test$estimate
p_value <- cor_test$p.value
# Format p-value for display
formatted_p_value <- ifelse(p_value == 0, "< 2.2e-16", sprintf("%.3f", p_value))
# Plot scatter plot
p <- ggplot(correlation_data, aes(x = log2FoldChange_spt6, y = log2FoldChange_med12)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  labs(
    title = "Correlation of Log2 Fold Changes",
    subtitle = paste(
      "Pearson correlation:", round(cor_value, 2), 
      "\nP-value:", formatted_p_value
    ),
    x = "Log2 Fold Change (SPT6)",
    y = "Log2 Fold Change (MED12)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# Save plot
ggsave("log2fc_correlation_scatter_with_pvalue.pdf", plot = p, width = 8, height = 6)
# Display plot
print(p)

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
write.csv(info_id %>% filter(name2 %in% intersect_down_gene), file = "intersectdown.csv",row.names = FALSE)
colnames(info_id) = c('trans', "gene")
############## RPO-seq pauseindex 
med12.pauseindex = read.csv("~/ljy/lailab2025/pby_med12/pro/pauseindex_MED12.txt", sep="\t")
colnames(med12.pauseindex)[1] <- "trans"
med12.pauseindex = med12.pauseindex %>% inner_join(info_id, by="trans")
spt6.pauseindex = read.csv("~/ljy/lailab2025/pby_med12/pro/pauseindex_SPT6.txt", sep="\t")
colnames(spt6.pauseindex)[1] <- "trans"
spt6.pauseindex = spt6.pauseindex %>% inner_join(info_id, by="trans")

h2azpuseindex = read.csv("~/ljy/lailab2025/pby_med12/pro/h2az/h2azpauseindex.tsv", sep="\t")
colnames(h2azpuseindex)[1] <- "trans"
h2azpuseindex = h2azpuseindex %>% inner_join(info_id, by="trans")

med12.pauseindex$regulate = ifelse(med12.pauseindex$gene %in% intersect_down_gene, "intersect_down", "others")
ggplot(med12.pauseindex, aes(x=regulate, y=CTRL_pause, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='bottom') + theme_bw()
ggplot(med12.pauseindex, aes(x=regulate, y=MED12KD_4h_pause, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='bottom') + theme_bw()
med12.pauseindex$log2FC_MED12KD = log2(med12.pauseindex$MED12KD_4h_pause + 1) - log2(med12.pauseindex$CTRL_pause + 1)
ggplot(med12.pauseindex, aes(x=regulate, y=log2FC_MED12KD, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='top') + theme_bw()

spt6.pauseindex$regulate = ifelse(spt6.pauseindex$gene %in% intersect_down_gene, "intersect_down", "others")
ggplot(spt6.pauseindex, aes(x=regulate, y=CTRL_pause, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='bottom') + theme_bw()
ggplot(spt6.pauseindex, aes(x=regulate, y=SPT6KD_4h_pause, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='bottom') + theme_bw()
spt6.pauseindex$log2FC_SPT6KD = log2(spt6.pauseindex$SPT6KD_4h_pause + 1) - log2(spt6.pauseindex$CTRL_pause + 1)
ggplot(spt6.pauseindex, aes(x=regulate, y=log2FC_SPT6KD, fill = regulate)) + geom_boxplot(outliers = FALSE) + 
  stat_compare_means(label.y.npc='top') + theme_bw()

h2azpuseindex$log2FC_h2azkd = log2(h2azpuseindex$H2AZ_24h_pause + 1) - log2(h2azpuseindex$H2AZ_0h_pause + 1)



############################################# H2AZ quantile 

library(ggplot2)
library(dplyr)
library(ggpubr)

n <- 3  # 设置 quantile 数量
alpha_vals <- seq(0.7, 1, length.out = n)
names(alpha_vals) <- as.character(1:n)
mnase <- read.csv("~/ljy/lailab2025/pby_med12/h2az/output_with_density.bed", sep = "\t", header = FALSE)
colnames(mnase)[4] <- "trans"
colnames(mnase)[7] <- "signal"

mnase <- mnase %>% inner_join(info_id, by = "trans")
mnase$regulate <- ifelse(mnase$gene %in% intersect_down_gene, "intersect_down", "others")
mnase_keep <- mnase

pdf("H2AZ_signal_comparison.pdf", width = 3, height = 5)

# Boxplot: signal by regulation category
ggplot(mnase, aes(x = regulate, y = signal, fill = regulate)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(label.y.npc = 'bottom') +
  theme_classic() + ggtitle("H2A.Z signal in intersect_down vs others") +
  theme(legend.position = "none")

# SPT6 DESeq2 vs H2A.Z
mnase <- mnase_keep %>% dplyr::select(signal, gene, regulate) %>% inner_join(spt6.res, by = 'gene') %>%
  na.omit() %>%
  mutate(q = ntile(signal, n)) %>%
  mutate(q = as.character(q))
paired <- combn(as.character(1:n), 2, simplify = FALSE)
paired = paired[1: length(paired)-1]
ggplot(mnase, aes(x = q, y = log2(baseMean), fill = "#FA7F6F", alpha = q)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") + ylab("Log2 Expression Level") +
  ggtitle("SPT6: Expression vs H2A.Z") +
  scale_alpha_manual(values = alpha_vals) +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

ggplot(mnase, aes(x = q, y = log2FoldChange, fill = "#FA7F6F", alpha = q)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("SPT6: Log2 Fold Change vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +  scale_alpha_manual(values = alpha_vals) +
  theme_classic() + theme(legend.position = "none")

# MED12 DESeq2 vs H2A.Z
mnase <- mnase_keep %>% dplyr::select(signal, gene, regulate) %>% inner_join(med12.res, by = 'gene') %>%
  na.omit() %>%
  mutate(q = ntile(signal, n)) %>%
  mutate(q = as.character(q))

ggplot(mnase, aes(x = q, y = log2(baseMean), fill = "#82B0D2", alpha = q)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") + ylab("Log2 Expression Level") +
  ggtitle("MED12: Expression vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +scale_alpha_manual(values = alpha_vals) +
  theme_classic() + theme(legend.position = "none")

ggplot(mnase, aes(x = q, y = log2FoldChange, fill = "#82B0D2", alpha = q)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("MED12: Log2 Fold Change vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +scale_alpha_manual(values = alpha_vals) +
  theme_classic() + theme(legend.position = "none")

# SPT6 Pause Index vs H2A.Z
mnase <- mnase_keep %>% dplyr::select(signal, gene, regulate) %>% inner_join(spt6.pauseindex, by = 'gene') %>%
  na.omit() %>%
  mutate(q = ntile(signal, n)) %>%
  mutate(q = as.character(q))

ggplot(mnase, aes(x = q, y = log2FC_SPT6KD, fill = "#FA7F6F", alpha = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("SPT6: Pause Index FC vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

ggplot(mnase, aes(x = q, y = log2(CTRL_pause), fill = "#FA7F6F", alpha = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("SPT6: CTRL Pause Index vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

# MED12 Pause Index vs H2A.Z
mnase <- mnase_keep %>% dplyr::select(signal, gene, regulate) %>% inner_join(med12.pauseindex, by = 'gene') %>%
  na.omit() %>%
  mutate(q = ntile(signal, n)) %>%
  mutate(q = as.character(q))

ggplot(mnase, aes(x = q, y = log2FC_MED12KD, fill = "#82B0D2", alpha = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("MED12: Pause Index FC vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

b = "#82B0D2"
print(ggplot(mnase, aes(x = q, y = log2(CTRL_pause), alpha  = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA, aes(fill = "#82B0D2")) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("MED12: CTRL Pause Index vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")
)
print(
  ggplot(mnase, aes(x = log2(CTRL_pause), colour = q)) +
    stat_ecdf(geom = "step", size = 1) + xlim(-1,6)+
    xlab("log2 CTRL Pause Index") +
    ylab("Cumulative fraction") +
    ggtitle("Cumulative distribution of CTRL Pause Index by H2A.Z quantiles") +
    theme_classic()
)


# h2az Pause Index vs H2A.Z
mnase <- mnase_keep %>% dplyr::select(signal, gene, regulate) %>% inner_join(h2azpuseindex, by = 'gene') %>%
  na.omit() %>%
  mutate(q = ntile(signal, n)) %>%
  mutate(q = as.character(q))

ggplot(mnase, aes(x = q, y = log2FC_h2azkd, fill = "#FA7F6F", alpha = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("h2az: Pause Index FC vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

ggplot(mnase, aes(x = q, y = log2(H2AZ_0h_pause), fill = "#FA7F6F", alpha = q)) +scale_alpha_manual(values = alpha_vals) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Quantile of H2A.Z strength") +
  ggtitle("h2az: CTRL Pause Index vs H2A.Z") +
  stat_compare_means(label.y.npc = 'bottom', comparisons = paired, method = "wilcox.test") +
  theme_classic() + theme(legend.position = "none")

dev.off()


# ######################################### GO for co down regulate
a=mnase %>% filter(q=="1")
x = intersect(a$gene, intersect_down_gene)
library(clusterProfiler)
library(org.Hs.eg.db)
# GO富集分析（BP为例，也可以做CC, MF）
ego <- enrichGO(gene         = intersect_down_gene,
                OrgDb        = org.Hs.eg.db,
                keyType      = "SYMBOL",
                ont = "ALL",
                readable      = TRUE)

# 查看结果前几行
head(ego)

# 绘图：柱状图、气泡图
barplot(ego, showCategory = 10)
dotplot(ego, showCategory = 10)

