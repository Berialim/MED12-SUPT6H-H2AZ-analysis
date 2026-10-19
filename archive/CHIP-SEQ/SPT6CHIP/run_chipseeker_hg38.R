#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(ggplot2)
})

input_file <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP/CTRL_summits.bed"
out_dir <- dirname(input_file)
prefix <- "CTRL_summits_hg38"

annotation_tsv <- file.path(out_dir, paste0(prefix, "_annotation.tsv"))
category_tsv <- file.path(out_dir, paste0(prefix, "_annotation_category_counts.tsv"))
promoter_bed <- file.path(out_dir, paste0(prefix, "_promoter_pm3kb.bed"))
pie_pdf <- file.path(out_dir, paste0(prefix, "_promoter_fraction_pie_Arial.pdf"))
summary_txt <- file.path(out_dir, paste0(prefix, "_analysis_summary.txt"))

if (!file.exists(input_file)) {
  stop("Input peak file does not exist: ", input_file)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ChIPseeker's promoter annotation here is explicitly defined as TSS +/- 3 kb.
tss_region <- c(-3000L, 3000L)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

raw <- read.delim(
  input_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (ncol(raw) < 3L) {
  stop("The input must contain at least 3 BED columns.")
}
if (any(raw[[2L]] < 0L) || any(raw[[3L]] <= raw[[2L]])) {
  stop("Invalid BED coordinates detected.")
}

peak_gr <- ChIPseeker::readPeakFile(input_file, header = FALSE)
if (length(peak_gr) != nrow(raw)) {
  stop("The number of peaks read by ChIPseeker does not match the input row count.")
}
mcols(peak_gr)$input_id <- seq_len(length(peak_gr))

peak_annotation <- ChIPseeker::annotatePeak(
  peak_gr,
  TxDb = txdb,
  tssRegion = tss_region,
  verbose = FALSE
)
annotation_df <- as.data.frame(peak_annotation)
if (!all(c("annotation", "input_id") %in% colnames(annotation_df))) {
  stop("ChIPseeker output did not preserve the input peak identifier.")
}

# ChIPseeker drops peaks on contigs absent from the hg38 knownGene TxDb.
# Reinsert those input peaks as explicit unannotated records so percentages
# use all input peaks rather than only the TxDb-mappable subset.
annotation_df$input_id <- as.integer(annotation_df$input_id)
missing_id <- setdiff(seq_len(nrow(raw)), annotation_df$input_id)
if (length(missing_id) > 0L) {
  missing_df <- as.data.frame(matrix(NA, nrow = length(missing_id), ncol = ncol(annotation_df)), stringsAsFactors = FALSE)
  colnames(missing_df) <- colnames(annotation_df)
  if ("seqnames" %in% colnames(missing_df)) missing_df$seqnames <- as.character(raw[[1L]][missing_id])
  if ("start" %in% colnames(missing_df)) missing_df$start <- as.integer(raw[[2L]][missing_id]) + 1L
  if ("end" %in% colnames(missing_df)) missing_df$end <- as.integer(raw[[3L]][missing_id])
  if ("width" %in% colnames(missing_df)) missing_df$width <- as.integer(raw[[3L]][missing_id] - raw[[2L]][missing_id])
  if ("strand" %in% colnames(missing_df)) missing_df$strand <- "*"
  missing_df$input_id <- missing_id
  missing_df$annotation <- "Unannotated (contig absent from hg38 TxDb)"
  annotation_df <- rbind(annotation_df, missing_df)
}
annotation_df <- annotation_df[order(annotation_df$input_id), , drop = FALSE]

if (nrow(annotation_df) != nrow(raw)) {
  stop("The completed annotation table does not contain one row per input peak.")
}

write.table(
  annotation_df,
  file = annotation_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

if (!"annotation" %in% colnames(annotation_df)) {
  stop("ChIPseeker output does not contain an annotation column.")
}
annotation <- as.character(annotation_df$annotation)
annotation[is.na(annotation) | annotation == ""] <- "Unannotated"
is_promoter <- grepl("^Promoter", annotation)

collapse_annotation <- function(x) {
  out <- rep("Other", length(x))
  out[grepl("^Promoter", x)] <- "Promoter (TSS ±3kb)"
  out[grepl("^5' UTR", x)] <- "5' UTR"
  out[grepl("^3' UTR", x)] <- "3' UTR"
  out[grepl("^Exon", x)] <- "Exon"
  out[grepl("^Intron", x)] <- "Intron"
  out[grepl("^Downstream", x)] <- "Downstream"
  out[grepl("^Distal Intergenic", x)] <- "Distal Intergenic"
  out[grepl("^Unannotated", x)] <- "Unannotated"
  out
}
broad_annotation <- collapse_annotation(annotation)

category_counts <- as.data.frame(table(broad_annotation), stringsAsFactors = FALSE)
colnames(category_counts) <- c("annotation", "peak_count")
category_counts <- category_counts[order(category_counts$peak_count, decreasing = TRUE), , drop = FALSE]
category_counts$percent <- 100 * category_counts$peak_count / sum(category_counts$peak_count)
category_counts$is_promoter <- category_counts$annotation == "Promoter (TSS ±3kb)"
write.table(
  category_counts,
  file = category_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

# BED is kept in the original 0-based, half-open coordinate convention.
peak_name <- if (ncol(raw) >= 4L) as.character(raw[[4L]]) else paste0("peak_", seq_len(nrow(raw)))
peak_score <- if (ncol(raw) >= 5L) as.character(raw[[5L]]) else rep(".", nrow(raw))
bed_df <- data.frame(
  chrom = as.character(raw[[1L]]),
  start = as.integer(raw[[2L]]),
  end = as.integer(raw[[3L]]),
  name = peak_name,
  score = peak_score,
  stringsAsFactors = FALSE
)
promoter_bed_df <- bed_df[is_promoter, , drop = FALSE]
write.table(
  promoter_bed_df,
  file = promoter_bed,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = "."
)

promoter_count <- sum(is_promoter)
total_count <- length(is_promoter)
promoter_percent <- 100 * promoter_count / total_count
annotated_count <- sum(!grepl("^Unannotated", annotation))
unannotated_count <- total_count - annotated_count

pie_category <- broad_annotation
pie_counts <- as.data.frame(table(pie_category), stringsAsFactors = FALSE)
colnames(pie_counts) <- c("annotation", "peak_count")
pie_counts <- pie_counts[order(pie_counts$peak_count, decreasing = TRUE), , drop = FALSE]
pie_counts$is_promoter <- pie_counts$annotation == "Promoter (TSS ±3kb)"
pie_counts$percent <- 100 * pie_counts$peak_count / sum(pie_counts$peak_count)

other_categories <- pie_counts$annotation[!pie_counts$is_promoter]
other_colors <- if (length(other_categories) > 0L) {
  grDevices::hcl.colors(length(other_categories), palette = "Dark 3")
} else {
  character(0)
}
fill_values <- c(setNames("#D73027", pie_counts$annotation[pie_counts$is_promoter]),
                 setNames(other_colors, other_categories))
pie_counts$label <- ifelse(
  pie_counts$percent >= 2,
  paste0(pie_counts$annotation, "\n", pie_counts$peak_count, " (", sprintf("%.1f", pie_counts$percent), "%)"),
  ""
)

pie_plot <- ggplot(pie_counts, aes(x = "", y = peak_count, fill = annotation)) +
  geom_col(width = 1, color = "white", linewidth = 0.25) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    family = "Arial",
    size = 3.0,
    lineheight = 0.9,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 1,
    y = 0,
    label = paste0("Promoter\n", promoter_count, " / ", total_count, " (", sprintf("%.1f", promoter_percent), "%)"),
    family = "Arial",
    fontface = "bold",
    size = 5.0,
    lineheight = 0.95
  ) +
  labs(
    title = "ChIPseeker peak annotation (hg38)",
    subtitle = paste0("Promoter = TSS ± 3 kb; input: ", basename(input_file)),
    fill = "Genomic annotation"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(family = "Arial", size = 10, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(family = "Arial", face = "bold", size = 10),
    legend.text = element_text(family = "Arial", size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

grDevices::pdf(pie_pdf, width = 9, height = 7, family = "Arial", useDingbats = FALSE)
print(pie_plot)
grDevices::dev.off()

summary_lines <- c(
  paste0("Input: ", input_file),
  paste0("Genome annotation: hg38, ", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  paste0("ChIPseeker version: ", as.character(packageVersion("ChIPseeker"))),
  "Gene symbol mapping: not applied; TxDb geneId/transcriptId fields retained in the annotation table.",
  "PDF device: grDevices::pdf(family = 'Arial')",
  paste0("Promoter definition: TSS ", tss_region[1L], " to +", tss_region[2L], " bp"),
  paste0("Total peaks: ", total_count),
  paste0("Peaks annotated against hg38 TxDb: ", annotated_count),
  paste0("Peaks without matching hg38 TxDb contig: ", unannotated_count),
  paste0("Promoter peaks: ", promoter_count),
  paste0("Promoter fraction (%): ", sprintf("%.4f", promoter_percent)),
  paste0("Annotation table: ", annotation_tsv),
  paste0("Category counts: ", category_tsv),
  paste0("Promoter BED: ", promoter_bed),
  paste0("Promoter pie PDF: ", pie_pdf),
  "",
  "Session info:",
  capture.output(sessionInfo())
)
writeLines(summary_lines, summary_txt)

cat("Completed ChIPseeker hg38 annotation.\n")
cat("Total peaks:", total_count, "\n")
cat("Promoter peaks:", promoter_count, "\n")
cat("Promoter fraction:", sprintf("%.4f%%", promoter_percent), "\n")
cat("Promoter BED:", promoter_bed, "\n")
cat("Pie PDF:", pie_pdf, "\n")
