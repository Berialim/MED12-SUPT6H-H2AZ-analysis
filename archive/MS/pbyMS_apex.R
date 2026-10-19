# 加载必要的包
library(tidyverse)
library(org.Hs.eg.db)
library(DEP)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ComplexHeatmap)
library(readxl)
library(ggrepel)

# ---------------------------
# 参数设置及数据读取
# ---------------------------
organism <- "org.Hs.eg.db"
setwd("~/ljy/lailab2025/pby_med12/MS/2023.09.24_MED12_APEX2/KO_MED12_CDK8_APEX2/")
file_path <- "data analysis-protein-MED12-CDK8-20230923.xlsx"  # 输入文件路径
# setwd("~/ljy/lailab2025/pby_med12/MS/2023.09.24_MED12_APEX2/KO_MED12_PRB1_APEX2/")
# file_path <- "data analysis-protein-MED12-PRB1-20230923.xlsx"  # 输入文件路径
sheet_name <- 1  # 选择第一个 sheet

# 读取数据（跳过第一行）
df_input <- readxl::read_excel(file_path, sheet = sheet_name, skip = 1)


colnames(df_input) = gsub("\\+", ".plus", colnames(df_input))
colnames(df_input) = gsub("\\-IAA", ".CTRL", colnames(df_input))
colnames(df_input) = gsub("\\-", ".", colnames(df_input))
colnames(df_input) = str_to_upper( colnames(df_input))

# 自动识别所有包含 "Intensity" 的列
intensity_columns <- grep("LFQ", colnames(df_input), value = TRUE)
df_input <- df_input %>%
  mutate(across(all_of(intensity_columns), ~ 2^as.numeric(.x)))
df_input[] <- lapply(df_input, function(x) ifelse(is.nan(x), NA, x))

# 获取重复样品（具有相同前缀，后缀为_1和_2的列），并提取样品名称前缀作为处理类型
samples <- grep("1$", intensity_columns, value = TRUE)
treatments <- gsub("1$", "", samples)

# ---------------------------
# 数据预处理
# ---------------------------
df_input <- df_input %>%
  as.data.frame() %>%
  dplyr::select(ENSG, all_of(intensity_columns)) %>%
  # 如果一个细胞中有多个 ENGS，用";"分割后展开，每个基因只保留一个
  separate_rows(ENSG, sep = ";") 
# %>% distinct(ENSG, .keep_all = TRUE)
# 基因符号转换：使用 bitr 将ENSEMBL转换为基因符号
gene <- bitr(df_input$ENSG, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = organism, drop = FALSE)
colnames(gene) <- c('ENSG', 'Gene')
gene = gene %>% distinct(ENSG, .keep_all = TRUE)
# 合并数据，确保 intensity 列为数值型，并去除全部为 NA 的行；按各样品均值降序排列并对重复基因保留一个
df_input <- df_input %>%
  inner_join(gene, by = "ENSG") %>%
  filter(!if_all(all_of(intensity_columns), ~ is.na(.))) %>%
  arrange(desc(rowMeans(across(all_of(intensity_columns)), na.rm = TRUE)))
# %>% distinct(Gene, .keep_all = TRUE)

gene_list <- df_input$Gene
ego <- enrichGO(gene          = gene_list,
                OrgDb         = organism,
                keyType       = "SYMBOL",
                ont           = "CC",           # 这里选择生物过程，可根据需求选择 MF 或 CC
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                maxGSSize = 50000,
                readable      = TRUE)
# 查看富集结果
ego_results <- as.data.frame(ego)

# 2. 选取描述中包含“transcription”的 GO term，得到相关基因
transcription_terms <- ego_results %>%
  filter(Description == "nucleus")
genes_trans <- unique(unlist(str_split(transcription_terms$geneID, "/")))
df_input = df_input %>% filter(Gene %in% genes_trans)
# ---------------------------
# 这里将以 sample 为列计算每个样品非 NA 的蛋白数
protein_counts <- sapply(df_input %>% dplyr::select(all_of(intensity_columns)), function(x) sum(!is.na(x)))
print(protein_counts)

# 可视化蛋白数（条形图示例）
protein_counts_df <- data.frame(Sample = names(protein_counts), Count = protein_counts)
ggplot(protein_counts_df, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  labs(title = "Peptide number", x = "Sample", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------
# 构建用于 DEP 分析的数据对象
# ---------------------------
# 将 df_input 中的 intensity 数据转换为矩阵，并将 Gene 作为行名
# 对于 df_input 数据框，根据 Gene 列生成唯一的基因名称（如果有重复则添加后缀）
df_input <- df_input %>%
  group_by(Gene) %>%
  mutate(Gene_unique = if(n() > 1) paste0(Gene, "_", row_number()) else Gene) %>%
  ungroup()

# 去除 'Intensity_' 前缀
sample_names <- gsub("^LFQ INTENSITY ", "", colnames(df_input))

# 提取 condition 和 replicate
condition <- sub(".\\d+$", "", sample_names)        # 去掉结尾的 -数字，就是 condition
replicate <- as.integer(sub(".*.(\\d+)$", "\\1", sample_names))  # 提取结尾的数字作为 replicate

# 构建数据框
exp_design <- data.frame(
  label = sample_names,
  condition = factor(condition),
  replicate = replicate
)
rownames(exp_design) <- colnames(df_input)
exp_design = exp_design %>% filter(!is.na(replicate))

LFQ_columns <- grep("LFQ", colnames(df_input)) # get LFQ column numbers
df_input = make_unique(df_input,"Gene","Gene_unique")

pdf("quality.pdf")
data_se <- make_se(df_input, LFQ_columns, exp_design)

data_se_parsed <- make_se_parse(df_input, LFQ_columns)
plot_frequency(data_se)



data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
plot_coverage(data_filt)

data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)
plot_missval(data_filt)
plot_detect(data_filt)

proteins_MNAR <- get_df_long(data_filt) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(data_filt) %in% proteins_MNAR

#Perform a mixed imputation
# data_imp <- impute(
#   data_filt,
#   fun = "mixed",
#   randna = !MNAR, # we have to define MAR which is the opposite of MNAR
#   mar = "knn", # imputation function for MAR
#   mnar = "MinProb") # imputation function for MNAR
# data_imp <- impute(data_norm, fun = "knn", rowmax = 0.9)
# data_imp <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.05)

plot_imputation(data_norm, data_imp)

exp_design$condition <- gsub(" ", ".", exp_design$condition)
ctrl_label <- as.character(unique(exp_design$condition[grepl("CTRL", exp_design$condition)])[1])
data_diff <- test_diff(data_imp, type = "control", control = ctrl_label)
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(2))
plot_pca(dep, x = 1, y = 2, n = 1000, point_size = 4)
dev.off()

# 获取整体差异分析结果
diff_results <- get_results(dep)
head(diff_results)

# ---------------------------
# 火山图绘制及标注 top20 基因
# 获取差异分析结果
diff_results <- get_results(dep) %>% as.data.frame()

# 自动识别 padj 和 log2 fold change 列名
padj_col <- grep("_p.val$", colnames(diff_results), value = TRUE)
lfc_col <- grep("_ratio$", colnames(diff_results), value = TRUE)

# 安全检查
if (length(padj_col) != 1 | length(lfc_col) != 1) {
  stop("找不到唯一的 padj 或 log2 fold change 列，请检查列名")
}

# 重命名并准备绘图数据
volcano_data <- diff_results %>%
  rename(
    padj = all_of(padj_col),
    log2FoldChange = all_of(lfc_col)
  ) %>%
  mutate(
    minusLog10Padj = -log10(padj),
    significant = ifelse( padj < 0.05 & (log2FoldChange) < 0, "Yes", "No") # padj < 0.05 &
  )

# 按 padj 排序选前20个显著基因用于标注
top20 <- volcano_data %>%
  filter(significant == "Yes") %>%
  arrange(padj) %>%
  slice_head(n = 5)

library(ggrepel)
# 绘制火山图
ggplot(volcano_data, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "blue")) +
  geom_text_repel(
    data = top20,
    aes(label = name),
    size = 6,
    color = "blue",
    max.overlaps = Inf,  # 避免过滤标签
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  theme_classic() +
  labs(
    title = paste0("Volcano n=", nrow(volcano_data), " sig.up=", nrow(volcano_data %>% filter(significant=="Yes"))),
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  )
ggsave("volcanno.pdf")

write.csv(volcano_data, file = "volcannodata.csv")
write.csv(volcano_data %>% filter(significant == 'Yes'), file = "sigdonw.csv")
write.csv(volcano_data %>% filter(padj < 0.05 & log2FoldChange > 0), file = "sigup.csv")

med12_CDK8 = read.csv("../KO_MED12_CDK8_APEX2/CDK8_MED12KDvolcannodata.csv")

cdk8up = med12_CDK8 %>% filter(significant == "Yes")
rpb1up = volcano_data %>% filter(significant == "Yes")
intersect(cdk8up$name , rpb1up$name)


df = inner_join(med12_CDK8, volcano_data, by = "name")
df = df %>% dplyr::select(name, log2FC_CDK8 = log2FoldChange.x, log2FC_PRB1 = log2FoldChange.y)
ggplot(df, aes(x = log2FC_CDK8, y = log2FC_PRB1)) +
  geom_point(color="grey45") +
  ggrepel::geom_text_repel(aes(label = ifelse(name == "SUPT6H", "SUPT6H", "")), color="blue", size=6) + theme_classic()
ggsave("log2FC_CDK8_PRB1.pdf")

# # ---------------------------
# # 针对转录相关基因的差异分析
# # ---------------------------
# # 1. 利用 clusterProfiler 对所有差异蛋白进行 GO 分析，筛选描述中包含“transcription”的 GO term
# gene_list <- df_input$Gene
# ego <- enrichGO(gene          = gene_list,
#                 OrgDb         = organism,
#                 keyType       = "SYMBOL",
#                 ont           = "BP",           # 这里选择生物过程，可根据需求选择 MF 或 CC
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 1,
#                 maxGSSize = 3000, 
#                 readable      = TRUE)
# # 查看富集结果
# ego_results <- as.data.frame(ego)
# print(ego_results)
# 
# # 2. 选取描述中包含“transcription”的 GO term，得到相关基因
# transcription_terms <- ego_results %>% 
#   filter(str_detect(Description, regex("transcription", ignore_case = TRUE)))
# # 提取所有转录相关 GO term 中出现的基因（geneID 字段中多个基因以 "/" 分隔）
# genes_trans <- unique(unlist(str_split(transcription_terms$geneID, "/")))
# print(genes_trans)
# 
# # 3. 从整体差异分析结果中筛选转录相关基因
# diff_results_trans <- volcano_data[volcano_data$name %in% genes_trans, ]
# print(head(diff_results_trans))
# # 按 padj 排序选前20个显著基因用于标注
# top20 <- diff_results_trans %>%
#   filter(significant == "Yes") %>%
#   arrange(padj) %>%
#   slice_head(n = 1)
# 
# library(ggrepel)
# # 绘制火山图
# ggplot(diff_results_trans, aes(x = log2FoldChange, y = minusLog10Padj)) +
#   geom_point(aes(color = significant), alpha = 0.6) +
#   scale_color_manual(values = c("No" = "grey", "Yes" = "blue")) +
#   geom_text_repel(
#     data = top20,
#     aes(label = name),
#     size = 6,
#     color = "blue",
#     max.overlaps = Inf,  # 避免过滤标签
#     box.padding = 0.3,
#     point.padding = 0.2,
#     segment.color = "grey50"
#   ) +
#   theme_classic() +
#   labs(
#     title = "Volcano Plot",
#     x = "log2 Fold Change",
#     y = "-log10(adjusted p-value)"
#   )



# ---------------------------
# Volcano plot with point size = -log10(padj)
# ---------------------------
volcano_data <- diff_results %>%
  rename(
    padj = all_of(padj_col),
    log2FoldChange = all_of(lfc_col)
  ) %>%
  mutate(
    minusLog10Padj = -log10(padj),
    significant = ifelse(padj < 0.05 & (log2FoldChange) < 0, "Yes", "No")  
  )

# Top 20
top20 <- volcano_data %>%
  filter(significant == "Yes") %>%
  arrange(padj) %>%
  slice_head(n = 6)

ggplot(volcano_data, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(aes(color = significant, size = minusLog10Padj), alpha = 0.8) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "blue")) +
  scale_size_continuous(range = c(0.5, 4)) +
  geom_text_repel(
    data = top20,
    aes(label = name),
    size = 5,
    color = "blue",
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  theme_classic() +
  labs(
    title = paste0("Volcano n=", nrow(volcano_data), " sig=", nrow(volcano_data %>% filter(significant=="Yes"))),
    x = "log2 Fold Change (normalized)",
    y = "-log10(p-value)"
  )
ggsave("volcanno.pdf", width = 6,height = 6)

med12_CDK8 <- read.csv("../KO_MED12_CDK8_APEX2/volcannodata.csv")
cdk8up <- med12_CDK8 %>% filter(significant == "Yes")
rpb1up <- volcano_data %>% filter(significant == "Yes")
intersect(cdk8up$name, rpb1up$name)

# 假设两个数据集都有padj
df <- inner_join(med12_CDK8, volcano_data, by = "name") %>%
  mutate(
    log2FC_CDK8 = log2FoldChange.x,
    log2FC_PRB1 = log2FoldChange.y,
    min_neglog10padj = pmin(-log10(padj.x), -log10(padj.y), na.rm = TRUE),
    highlight = (log2FC_CDK8 < -1 & log2FC_PRB1 < -1 & min_neglog10padj > -log10(0.05)) # 高于标准的点
  )

ggplot(df, aes(x = log2FC_CDK8, y = log2FC_PRB1)) +
  geom_point(aes(size = min_neglog10padj, color = highlight), alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "blue")) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(highlight, name, "")),
    color = "blue",
    size = 4,
    max.overlaps = 20
  ) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(size = "Minimal -log10(padj)", color = "Above threshold")

ggsave("log2FC_CDK8_PRB1_highlight.pdf", width = 7, height = 5)




# ---------------------------
# Row-scaled heatmap for selected genes
# ---------------------------
# 提取 normalized matrix
library(SummarizedExperiment)
library(gtools) 
mat_norm <- assay(data_imp)

# 添加基因名
rownames(data_imp) <- rowData(data_imp)$name

# 筛选目标基因
selected_genes <- grep("^SUPT|^MED|^CDK|^POLR2|^NELF|^PAF", rownames(mat_norm), value = TRUE)
selected_genes <- selected_genes[!grepl("\\.|SUPT16H|CDK5RAP1|CDK5RAP3", selected_genes)]

# 子集并排序
mat_sel <- mat_norm[selected_genes, ]
mat_sel <- mat_sel[mixedorder(rownames(mat_sel)), ]

# Row Z-score
# mat_scaled <- t(scale(t(mat_sel)))

# 绘制热图
pdf("selected_genes_heatmap.pdf", width = 6, height = 10)
Heatmap(mat_sel,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_names_side = "left",
        column_names_rot = 45,
        show_row_names = TRUE)
dev.off()

