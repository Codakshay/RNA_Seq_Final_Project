# Load libraries ---------------------------------------------------------
library(org.Hs.eg.db)
library(tidyverse)
library(GEOquery)
library(airway)
library(DESeq2)
library(EnhancedVolcano)
library(RColorBrewer)
library(ggplotify)
library(pheatmap)
library(ggrepel)

# Step I: Preparing count data --------------------------------------------

# Load RNA-Seq data
counts_data <- read.csv("counts.csv", stringsAsFactors = FALSE)
rownames(counts_data) <- counts_data$Ensembl_ID  # Set Ensemble ID as rownames
counts_data <- counts_data[, -c(1:4)]  # Drop metadata columns

# Load metadata from GEO
gse <- getGEO("GSE80336")
metadata <- gse[["GSE80336_series_matrix.txt.gz"]]@phenoData
coldata <- metadata@data  # Extract metadata

# Preprocess metadata
colnames(counts_data) <- coldata$geo_accession
coldata$title <- gsub("control-.*", "control", coldata$title)
coldata$title <- gsub("bipolar-.*", "bipolar", coldata$title)
colnames(coldata)[1] <- "condition"
coldata <- coldata[, -2]
keeps <- c("condition", "age (years):ch1", "Sex:ch1", 
           "postmortem interval (hours):ch1", "rin:ch1")
coldata <- coldata[keeps]

# Sanity check for column names
stopifnot(all(colnames(counts_data) %in% rownames(coldata)))
stopifnot(all(colnames(counts_data) == rownames(coldata)))

# Step II: Constructing DESeqDataSet --------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = coldata,
  design = ~ condition
)

# Pre-filtering low-count genes
dds <- dds[rowSums(counts(dds)) >= 1, ]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "control")

# Step III: Running DESeq procedure ---------------------------------------

dds <- DESeq(dds)

# Save normalized read counts
normalized_counts <- counts(dds, normalized = TRUE)

# Extract and summarize results
res <- results(dds, contrast = c("condition", "control", "bipolar"), alpha = 0.05)
summary(res)

# Significant genes with adjusted p-value < 0.05
cat("Significant genes:", sum(res$padj < 0.05, na.rm = TRUE), "\n")

# Order results by smallest adjusted p-value
resOrdered <- res[order(res$padj), ]

# Step IV: Visualization --------------------------------------------------

# MA Plot
plotMA(res, ylim = c(-10, 10), cex = 0.7)
abline(h = c(-1, 1), col = "red", lwd = 3)

# Shrink log fold changes for visualization
resLFC <- lfcShrink(dds, coef = "condition_bipolar_vs_control", type = "apeglm")
plotMA(resLFC, ylim = c(-10, 10), cex = 0.7)
abline(h = c(-1, 1), col = "red", lwd = 3)

# Dispersion plot
plotDispEsts(dds, main = "Dispersion Plot")

# Principal Component Analysis (PCA)
rld <- vst(dds, blind = FALSE)
PCAA <- plotPCA(rld, intgroup = "condition")
PCAA + geom_text(aes(label = name), size = 2.5) + ggtitle("PCA Plot")

# Volcano plot
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res.df),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)
EnhancedVolcano(res.df,
                x = "log2FoldChange",
                y = "padj",
                lab = res.df$symbol,
                pCutoff = 0.1,
                FCcutoff = 1,
                ylim = c(0,2.5),
                title = "Volcano plot of Differentially Expressed Genes",
                subtitle = bquote(italic('in human dorsal striatum'))
                )

# display up-regulated genes
up_regulated <- na.omit(res.df[res.df$padj < 0.1 & res.df$log2FoldChange > 1.0,])
up_regulated

# display down-regulated genes
down_regulated <- na.omit(res.df[res.df$padj < 0.1 & res.df$log2FoldChange < -1.0,])
down_regulated


# Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colors,
         main = 'Heatmap of sample-to-sample distances')

# Heatmap of significant DEGs
top_genes <- rownames(resOrdered)[1:2000]
logcolors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
pheatmap(assay(rld)[top_genes, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         annotation_col = coldata["condition"],
         color = logcolors,
         main = 'Heatmap of significant DEGs')
