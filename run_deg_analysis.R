# Load libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(GEOquery)
  library(DESeq2)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(pheatmap)
  library(ggrepel)
  library(BiocParallel)
})

# Parse CLI arguments ----------------------------------------------------
option_list <- list(
  make_option(c("--counts"), type = "character", metavar = "FILE",
              help = "Path to the merged counts CSV produced by merge_transcripts.py [required]"),
  make_option(c("--gse"), type = "character", metavar = "ACCESSION",
              help = "GEO Series accession (e.g. GSE80336) used to fetch sample metadata [required]"),
  make_option(c("--condition-field"), type = "character", default = "title",
              metavar = "COLUMN",
              help = paste("Column in GEO phenoData whose values identify sample groups.",
                           "Trailing '-<number>' suffixes are stripped automatically,",
                           "so 'control-2' becomes 'control'.",
                           "[default: title]")),
  make_option(c("--out"), type = "character", default = "deg_results.csv",
              metavar = "FILE",
              help = "Output CSV for ordered DEG results [default: deg_results.csv]"),
  make_option(c("--plots-dir"), type = "character", default = "data/plots",
              metavar = "DIR",
              help = "Directory for output PNG plots [default: data/plots]"),
  make_option(c("--workers"), type = "integer", default = 1,
              metavar = "N",
              help = paste("BiocParallel worker count for the DESeq2 dispersion fit.",
                           "DESeq2 is CPU-only (no GPU port exists); this gives a",
                           "modest multi-core speedup on large datasets only.",
                           "[default: 1]"))
)

opt <- parse_args(
  OptionParser(option_list = option_list),
  convert_hyphens_to_underscores = TRUE
)

if (is.null(opt$counts) || is.null(opt$gse)) {
  stop("--counts and --gse are required. Run with --help for usage.", call. = FALSE)
}

dir.create(opt$plots_dir, recursive = TRUE, showWarnings = FALSE)
cat("Writing plots to:", opt$plots_dir, "\n")

# Register the parallel backend used by DESeq2's parametric fit step.
# MulticoreParam(1) is equivalent to serial and is the safe default; setting
# workers > 1 forks N processes for the dispersion-fit and Wald-test loops.
if (opt$workers > 1) {
  register(MulticoreParam(workers = opt$workers))
  cat("BiocParallel: registered MulticoreParam with", opt$workers, "workers.\n")
}

# Step I: Preparing count data -------------------------------------------
counts_data <- read.csv(opt$counts, stringsAsFactors = FALSE, check.names = FALSE)
rownames(counts_data) <- counts_data$Ensembl_ID

# Drop the 4 annotation columns by name (safer than positional slicing)
meta_cols   <- c("Ensembl_ID", "GeneSymbol", "Biotype", "Chromosome")
counts_data <- counts_data[, !(names(counts_data) %in% meta_cols), drop = FALSE]

# Load sample metadata from GEO
gse     <- getGEO(opt$gse)
coldata <- pData(gse[[1]])

# Align counts columns to coldata rows using GSM IDs as the shared key
common_gsm <- intersect(colnames(counts_data), rownames(coldata))
if (length(common_gsm) == 0) {
  stop(paste(
    "No GSM IDs are shared between the counts matrix and GEO metadata for", opt$gse,
    "\nEnsure merge_transcripts.py produced GSM accession IDs as column names."
  ), call. = FALSE)
}
counts_data <- counts_data[, common_gsm, drop = FALSE]
coldata     <- coldata[common_gsm, , drop = FALSE]

# Derive condition labels: strip trailing '-<digits>' (e.g. "control-2" → "control")
condition_raw    <- coldata[[opt$condition_field]]
coldata$condition <- factor(gsub("-\\d+$", "", condition_raw))
cat("Condition levels:", paste(levels(coldata$condition), collapse = ", "), "\n")

stopifnot(all(colnames(counts_data) == rownames(coldata)))

# Step II: Constructing DESeqDataSet -------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData   = coldata,
  design    = ~ condition
)

# Pre-filtering: remove genes with zero total counts
dds <- dds[rowSums(counts(dds)) >= 1, ]

# Reference level: first alphabetically (typically "control")
dds$condition <- relevel(dds$condition, ref = levels(dds$condition)[1])
cat("Reference condition:", levels(dds$condition)[1], "\n")

# Step III: Running DESeq2 -----------------------------------------------
# parallel = TRUE forwards work to the BiocParallel backend registered above;
# it's a no-op when workers == 1 (default).
dds <- DESeq(dds, parallel = (opt$workers > 1))

normalized_counts <- counts(dds, normalized = TRUE)

# Default contrast: treatment vs reference (derived from factor levels)
res <- results(dds, alpha = 0.05)
summary(res)
cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")

resOrdered <- res[order(res$padj), ]

# LFC-shrunken results for visualization (second coefficient = treatment vs ref)
coef_name <- resultsNames(dds)[2]
cat("Using coefficient for LFC shrinkage:", coef_name, "\n")
resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# Step IV: Visualization -------------------------------------------------

# MA Plot (raw)
png(file.path(opt$plots_dir, "MAPlot.png"), width = 800, height = 600, res = 120)
plotMA(res, ylim = c(-10, 10), cex = 0.7, main = "MA Plot")
abline(h = c(-1, 1), col = "red", lwd = 3)
dev.off()

# MA Plot (LFC-shrunken)
png(file.path(opt$plots_dir, "resMAPlot.png"), width = 800, height = 600, res = 120)
plotMA(resLFC, ylim = c(-10, 10), cex = 0.7, main = "MA Plot (LFC shrinkage)")
abline(h = c(-1, 1), col = "red", lwd = 3)
dev.off()

# Dispersion plot
png(file.path(opt$plots_dir, "DispersionPlot.png"), width = 800, height = 600, res = 120)
plotDispEsts(dds, main = "Dispersion Plot")
dev.off()

# PCA plot
rld  <- vst(dds, blind = FALSE)
PCAA <- plotPCA(rld, intgroup = "condition")
PCAA <- PCAA +
  geom_text_repel(aes(label = name), size = 2.5) +
  ggtitle("PCA Plot") +
  theme_bw()
ggsave(file.path(opt$plots_dir, "PCAPlot.png"), plot = PCAA,
       width = 10, height = 7, dpi = 120)

# Volcano plot
res_df <- as.data.frame(res)
res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(res_df),
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)
volcano <- EnhancedVolcano(
  res_df,
  x        = "log2FoldChange",
  y        = "padj",
  lab      = res_df$symbol,
  pCutoff  = 0.1,
  FCcutoff = 1,
  ylim     = c(0, 2.5),
  title    = "Volcano plot of Differentially Expressed Genes",
  subtitle = bquote(italic("in human dorsal striatum"))
)
ggsave(file.path(opt$plots_dir, "VolcanoPlot.png"), plot = volcano,
       width = 12, height = 8, dpi = 120)

# Heatmap of sample-to-sample distances
sample_dists      <- dist(t(assay(rld)))
sample_dist_mat   <- as.matrix(sample_dists)
heat_colors       <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(file.path(opt$plots_dir, "HeatmapPairwisePlot.png"),
    width = 1200, height = 1000, res = 150)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         color = heat_colors,
         main  = "Heatmap of sample-to-sample distances")
dev.off()

# Heatmap of top DEGs (up to 2 000 rows)
n_top     <- min(2000, nrow(resOrdered))
top_genes <- na.omit(rownames(resOrdered)[seq_len(n_top)])
log_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
anno_col   <- data.frame(condition = coldata$condition, row.names = rownames(coldata))
png(file.path(opt$plots_dir, "HeatmapDEGPlot.png"),
    width = 1200, height = 1000, res = 150)
pheatmap(assay(rld)[top_genes, ],
         cluster_rows   = TRUE,
         show_rownames  = FALSE,
         show_colnames  = TRUE,
         cluster_cols   = TRUE,
         fontsize_row   = 8,
         annotation_col = anno_col,
         color          = log_colors,
         main           = "Heatmap of significant DEGs")
dev.off()

cat("All plots saved to:", opt$plots_dir, "\n")

# Step V: Write DEG results CSV ------------------------------------------
resOrdered_df            <- as.data.frame(resOrdered)
resOrdered_df$Ensembl_ID <- rownames(resOrdered_df)
resOrdered_df$symbol     <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(resOrdered_df),
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)
# Put identifiers first, then statistics
stat_cols     <- setdiff(names(resOrdered_df), c("Ensembl_ID", "symbol"))
resOrdered_df <- resOrdered_df[, c("Ensembl_ID", "symbol", stat_cols)]
write.csv(resOrdered_df, opt$out, row.names = FALSE)
cat("DEG results written to:", opt$out, "\n")
