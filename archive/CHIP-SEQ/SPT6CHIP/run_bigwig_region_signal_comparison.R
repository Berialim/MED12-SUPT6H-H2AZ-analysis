#!/usr/bin/env Rscript

## Compare CTRL and MED12addIAA BigWig signal over ChIPseeker hg38
## annotation categories.  Each one-base summit is expanded to a
## 501-bp window (summit +/- 250 bp) before signal extraction.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

input_file <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP/CTRL_summits.bed"
annotation_file <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP/CTRL_summits_hg38_annotation.tsv"
ctrl_bw_path <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP/CTRL.bw"
treatment_bw_path <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP/MED12addIAA.bw"
out_dir <- "/mnt/d/work/pby/CHIP-SEQ/SPT6CHIP"

prefix <- "CTRL_vs_MED12addIAA_pm250bp"
peak_signal_tsv <- file.path(out_dir, paste0(prefix, "_peak_signal.tsv"))
category_stats_tsv <- file.path(out_dir, paste0(prefix, "_category_statistics.tsv"))
windows_bed <- file.path(out_dir, paste0(prefix, "_windows.bed"))
plot_pdf <- file.path(out_dir, paste0(prefix, "_category_boxplot.pdf"))
summary_txt <- file.path(out_dir, paste0(prefix, "_analysis_summary.txt"))

required_files <- c(input_file, annotation_file, ctrl_bw_path, treatment_bw_path)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required file(s):\n",
    paste(required_files[!file.exists(required_files)], collapse = "\n")
  )
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Read the original BED as 0-based, half-open coordinates.
raw <- read.delim(
  input_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)
if (ncol(raw) < 3L) {
  stop("The input peak file must contain at least three BED columns.")
}
if (any(is.na(raw[[2L]]) | is.na(raw[[3L]]) |
        raw[[2L]] < 0L | raw[[3L]] <= raw[[2L]])) {
  stop("Invalid BED coordinates detected.")
}

peak_id <- seq_len(nrow(raw))
chrom <- as.character(raw[[1L]])
start0 <- as.integer(raw[[2L]])
end0 <- as.integer(raw[[3L]])

## The input is a one-base summit BED.  In BED coordinates, a +/-250-bp
## window around [start0, start0+1) is [start0-250, start0+251).
window_start0 <- pmax(0L, start0 - 250L)
window_end0 <- end0 + 250L
windows <- GRanges(
  seqnames = chrom,
  ranges = IRanges(start = window_start0 + 1L, end = window_end0)
)
names(windows) <- paste0("peak_", peak_id)

## Reuse the completed ChIPseeker hg38 annotation from the preceding step.
## It contains one row per original peak and an input_id column preserving
## the BED order.
annotation_df <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)
if (!all(c("input_id", "annotation") %in% colnames(annotation_df))) {
  stop("The ChIPseeker annotation table must contain input_id and annotation columns.")
}
annotation_df$input_id <- as.integer(annotation_df$input_id)
annotation_df <- annotation_df[order(annotation_df$input_id), , drop = FALSE]
if (nrow(annotation_df) != nrow(raw) ||
    !identical(annotation_df$input_id, peak_id)) {
  stop("Annotation table does not contain exactly one row per input peak in BED order.")
}

collapse_annotation <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unannotated"
  out <- rep("Other", length(x))
  out[grepl("^Promoter", x)] <- "Promoter (TSS +/-3kb)"
  out[grepl("^5' UTR", x)] <- "5' UTR"
  out[grepl("^3' UTR", x)] <- "3' UTR"
  out[grepl("^Exon", x)] <- "Exon"
  out[grepl("^Intron", x)] <- "Intron"
  out[grepl("^Downstream", x)] <- "Downstream"
  out[grepl("^Distal Intergenic", x)] <- "Distal Intergenic"
  out[grepl("^Unannotated", x)] <- "Unannotated"
  out
}

broad_category <- collapse_annotation(annotation_df$annotation)

## Extract one mean BigWig value per window.  Ranges on contigs absent from
## a track remain NA and are subsequently treated as zero signal.
extract_bigwig_mean <- function(path, ranges) {
  bw <- BigWigFile(path)
  on.exit(try(close(bw), silent = TRUE), add = TRUE)

  bw_seqinfo <- seqinfo(bw)
  bw_chroms <- seqlevels(bw_seqinfo)
  query_chrom <- as.character(seqnames(ranges))
  signal <- rep(NA_real_, length(ranges))

  for (this_chrom in intersect(unique(query_chrom), bw_chroms)) {
    idx <- which(query_chrom == this_chrom)
    ans <- rtracklayer::summary(
      bw,
      which = ranges[idx],
      size = 1L,
      type = "mean"
    )
    if (is.matrix(ans) || is.data.frame(ans)) {
      ans <- ans[, 1L]
    } else if (inherits(ans, "GRangesList")) {
      ans <- mcols(unlist(ans, use.names = FALSE))[["score"]]
    }
    signal[idx] <- as.numeric(ans)
  }

  list(signal = signal, seqinfo = bw_seqinfo, chroms = bw_chroms)
}

ctrl_result <- extract_bigwig_mean(ctrl_bw_path, windows)
treatment_result <- extract_bigwig_mean(treatment_bw_path, windows)

ctrl_na_before_zero <- is.na(ctrl_result$signal)
treatment_na_before_zero <- is.na(treatment_result$signal)
ctrl_signal <- ctrl_result$signal
treatment_signal <- treatment_result$signal
ctrl_signal[ctrl_na_before_zero] <- 0
treatment_signal[treatment_na_before_zero] <- 0

## A small pseudocount is used only for log2 fold change and the plot scale.
## The Wilcoxon tests use the extracted signal values directly.
pseudocount <- 1e-06
log2fc_treatment_vs_ctrl <- log2(
  (treatment_signal + pseudocount) / (ctrl_signal + pseudocount)
)

peak_signal <- data.frame(
  peak_id = peak_id,
  chrom = chrom,
  summit_start0 = start0,
  summit_end0 = end0,
  window_start0 = window_start0,
  window_end0 = window_end0,
  window_width = window_end0 - window_start0,
  original_name = if (ncol(raw) >= 4L) as.character(raw[[4L]]) else paste0("peak_", peak_id),
  original_score = if (ncol(raw) >= 5L) as.character(raw[[5L]]) else ".",
  chipseeker_annotation = as.character(annotation_df$annotation),
  region_category = broad_category,
  ctrl_signal_mean = ctrl_signal,
  med12addIAA_signal_mean = treatment_signal,
  log2FC_med12addIAA_vs_ctrl = log2fc_treatment_vs_ctrl,
  med12addIAA_lower_than_ctrl = treatment_signal < ctrl_signal,
  ctrl_signal_was_NA_set_to_zero = ctrl_na_before_zero,
  med12addIAA_signal_was_NA_set_to_zero = treatment_na_before_zero,
  stringsAsFactors = FALSE
)
write.table(
  peak_signal,
  file = peak_signal_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

## Export the exact windows used for signal extraction as BED.
window_bed_df <- data.frame(
  chrom = chrom,
  start0 = window_start0,
  end0 = window_end0,
  name = paste0("peak_", peak_id),
  score = ".",
  region_category = broad_category,
  stringsAsFactors = FALSE
)
write.table(
  window_bed_df,
  file = windows_bed,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = "."
)

safe_wilcox <- function(x, y, alternative = "two.sided") {
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  if (all(is.finite(x)) && all(is.finite(y))) {
    return(suppressWarnings(
      tryCatch(
        wilcox.test(x, y, paired = TRUE, alternative = alternative, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    ))
  }
  NA_real_
}

category_levels <- unique(broad_category)
category_stats <- do.call(rbind, lapply(category_levels, function(cat) {
  idx <- which(broad_category == cat)
  x <- ctrl_signal[idx]
  y <- treatment_signal[idx]
  logfc <- log2fc_treatment_vs_ctrl[idx]
  data.frame(
    region_category = cat,
    peak_count = length(idx),
    ctrl_median = median(x),
    med12addIAA_median = median(y),
    ctrl_mean = mean(x),
    med12addIAA_mean = mean(y),
    median_log2FC_med12addIAA_vs_ctrl = median(logfc),
    mean_log2FC_med12addIAA_vs_ctrl = mean(logfc),
    median_signal_delta_med12addIAA_minus_ctrl = median(y - x),
    lower_peak_count = sum(y < x),
    lower_peak_fraction = mean(y < x),
    p_value_two_sided_paired_wilcoxon = safe_wilcox(x, y, "two.sided"),
    p_value_decrease_ctrl_greater_than_treatment = safe_wilcox(x, y, "greater"),
    stringsAsFactors = FALSE
  )
}))

category_stats$FDR_two_sided <- p.adjust(
  category_stats$p_value_two_sided_paired_wilcoxon,
  method = "BH"
)
category_stats$FDR_decrease <- p.adjust(
  category_stats$p_value_decrease_ctrl_greater_than_treatment,
  method = "BH"
)
category_stats$decrease_rank <- rank(
  category_stats$median_log2FC_med12addIAA_vs_ctrl,
  ties.method = "first",
  na.last = "keep"
)
category_stats <- category_stats[order(
  category_stats$median_log2FC_med12addIAA_vs_ctrl,
  na.last = TRUE
), , drop = FALSE]
write.table(
  category_stats,
  file = category_stats_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

## Base-R PDF boxplot.  The plotted values are log2(mean signal + 1) for
## readability; the displayed p values are paired two-sided Wilcoxon tests
## on the original mean signal values within each category.
plot_stats <- category_stats
plot_categories <- as.character(plot_stats$region_category)
box_values <- list()
box_at <- numeric(0)
box_col <- character(0)
category_mid <- numeric(nrow(plot_stats))
category_max <- numeric(nrow(plot_stats))

for (i in seq_len(nrow(plot_stats))) {
  cat <- plot_categories[i]
  idx <- which(broad_category == cat)
  ctrl_plot <- log2(ctrl_signal[idx] + 1)
  treatment_plot <- log2(treatment_signal[idx] + 1)
  box_values[[2L * i - 1L]] <- ctrl_plot
  box_values[[2L * i]] <- treatment_plot
  box_at <- c(box_at, 2L * i - 1L, 2L * i)
  box_col <- c(box_col, "#4C78A8", "#E45756")
  category_mid[i] <- 2L * i - 0.5
  category_max[i] <- max(c(ctrl_plot, treatment_plot), na.rm = TRUE)
}

all_plot_values <- unlist(box_values, use.names = FALSE)
y_range <- range(all_plot_values[is.finite(all_plot_values)])
if (!all(is.finite(y_range))) {
  y_range <- c(0, 1)
}
if (diff(y_range) == 0) {
  y_range <- y_range + c(-1, 1)
}
y_pad <- max(0.2, diff(y_range) * 0.08)
p_y <- category_max + y_pad
y_limits <- c(y_range[1L] - diff(y_range) * 0.05, max(p_y) + y_pad * 1.5)

grDevices::pdf(plot_pdf, width = 15, height = 9, useDingbats = FALSE)
op <- par(no.readonly = TRUE)
par(mar = c(16, 5, 5, 2) + 0.1)
boxplot(
  box_values,
  at = box_at,
  col = box_col,
  border = "grey25",
  outline = FALSE,
  xaxt = "n",
  ylim = y_limits,
  ylab = "log2(mean BigWig signal + 1)",
  xlab = "",
  main = "CTRL vs MED12addIAA signal over summit +/-250 bp",
  sub = ""
)
mtext(
  "Paired two-sided Wilcoxon p values; one BigWig per condition (exploratory peak-level comparison)",
  side = 3,
  line = 0.5,
  cex = 0.78
)
axis(
  1,
  at = box_at,
  labels = rep(c("CTRL", "MED12addIAA"), nrow(plot_stats)),
  las = 2,
  cex.axis = 0.68,
  line = 0
)
axis(
  1,
  at = category_mid,
  labels = plot_categories,
  tick = FALSE,
  las = 2,
  cex.axis = 0.72,
  line = 8
)
for (i in seq_len(nrow(plot_stats))) {
  p <- plot_stats$p_value_two_sided_paired_wilcoxon[i]
  fdr <- plot_stats$FDR_two_sided[i]
  p_label <- if (is.na(p)) {
    "p=NA"
  } else {
    paste0("p=", format.pval(p, digits = 2, eps = 1e-3),
           "\nFDR=", format.pval(fdr, digits = 2, eps = 1e-3))
  }
  text(category_mid[i], p_y[i], labels = p_label, cex = 0.68, xpd = NA)
}
legend(
  "topright",
  legend = c("CTRL", "MED12addIAA"),
  fill = c("#4C78A8", "#E45756"),
  border = "grey25",
  bty = "n",
  cex = 0.85
)
par(op)
grDevices::dev.off()

top_decreased <- category_stats[
  order(category_stats$median_log2FC_med12addIAA_vs_ctrl, na.last = TRUE),
  , drop = FALSE
]

summary_lines <- c(
  "CTRL vs MED12addIAA BigWig signal comparison",
  paste0("Input peak file: ", input_file),
  paste0("CTRL BigWig: ", ctrl_bw_path),
  paste0("Treatment BigWig: ", treatment_bw_path),
  "Genome build stated by user: hg38",
  "Annotation source: existing ChIPseeker hg38 annotation table",
  "Promoter category: ChIPseeker promoter annotations collapsed as TSS +/-3 kb",
  "Signal window: each BED summit expanded to +/-250 bp; window width is 501 bp unless clipped at chromosome start",
  "Signal statistic: mean BigWig value across each window; NA/no-track coverage values were set to 0",
  "P value: paired two-sided Wilcoxon signed-rank test within each region category, matching the same peak in CTRL and treatment",
  "FDR: Benjamini-Hochberg adjustment across region categories",
  "Interpretation: one BigWig per condition means these are exploratory peak-level comparisons, not replicate-aware differential binding statistics",
  paste0("Total peaks: ", nrow(peak_signal)),
  paste0("CTRL NA values replaced by zero: ", sum(ctrl_na_before_zero)),
  paste0("MED12addIAA NA values replaced by zero: ", sum(treatment_na_before_zero)),
  paste0("Peak-level signal table: ", peak_signal_tsv),
  paste0("Category statistics: ", category_stats_tsv),
  paste0("Analysis windows BED: ", windows_bed),
  paste0("Boxplot PDF: ", plot_pdf),
  "",
  "Categories ordered by median log2FC (most decreased first):",
  capture.output(print(top_decreased, row.names = FALSE)),
  "",
  "Session info:",
  capture.output(sessionInfo())
)
writeLines(summary_lines, summary_txt)

cat("Completed CTRL vs MED12addIAA BigWig comparison.\n")
cat("Peaks:", nrow(peak_signal), "\n")
cat("Category statistics:", category_stats_tsv, "\n")
cat("Boxplot PDF:", plot_pdf, "\n")
