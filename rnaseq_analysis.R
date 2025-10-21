#!/usr/bin/env Rscript
# RNA-seq Differential Expression Analysis with DESeq2 and edgeR
# 
# This script performs differential expression analysis on gene count data
# from STAR alignment output.
#
# Usage: Rscript rnaseq_analysis.R [count_matrix] [metadata] [output_dir]
#
# Arguments:
#   count_matrix: Path to gene count matrix TSV file
#   metadata: Path to sample metadata TSV file  
#   output_dir: Output directory for results (default: "results")

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(EnhancedVolcano)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript rnaseq_analysis.R <count_matrix> <metadata> [output_dir]")
}

count_file <- args[1]
metadata_file <- args[2]
output_dir <- if (length(args) >= 3) args[3] else "results"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("RNA-seq Differential Expression Analysis\n")
cat("=======================================\n")
cat("Count matrix:", count_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Read data
cat("Reading data...\n")
count_matrix <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure sample names match
common_samples <- intersect(colnames(count_matrix), metadata$sample)
count_matrix <- count_matrix[, common_samples]
metadata <- metadata[metadata$sample %in% common_samples, ]

cat("Samples:", ncol(count_matrix), "\n")
cat("Genes:", nrow(count_matrix), "\n\n")

# Quality control plots
cat("Generating quality control plots...\n")

# Total read counts per sample
pdf(file.path(output_dir, "qc_total_counts.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
barplot(colSums(count_matrix), las = 2, main = "Total Read Counts per Sample")
boxplot(colSums(count_matrix), main = "Distribution of Total Counts")
dev.off()

# Gene count distribution
pdf(file.path(output_dir, "qc_gene_counts.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
hist(rowSums(count_matrix), breaks = 50, main = "Total Counts per Gene", xlab = "Total Counts")
hist(log10(rowSums(count_matrix) + 1), breaks = 50, main = "Log10 Total Counts per Gene", xlab = "Log10(Total Counts + 1)")
dev.off()

# DESeq2 Analysis
cat("Running DESeq2 analysis...\n")

# Prepare DESeq2 data
if ("treatment" %in% colnames(metadata)) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = metadata,
                               design = ~ treatment)
} else if ("group" %in% colnames(metadata)) {
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = metadata,
                               design = ~ group)
} else {
  stop("No suitable grouping variable found in metadata (need 'treatment' or 'group')")
}

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# Get results
if ("treatment" %in% colnames(metadata)) {
  treatments <- unique(metadata$treatment)
  if (length(treatments) == 2) {
    res <- results(dds, contrast = c("treatment", treatments[2], treatments[1]))
  } else {
    res <- results(dds)
  }
} else {
  groups <- unique(metadata$group)
  if (length(groups) == 2) {
    res <- results(dds, contrast = c("group", groups[2], groups[1]))
  } else {
    res <- results(dds)
  }
}

# Add gene names
res$gene_id <- rownames(res)

# Save DESeq2 results
write.table(res, file.path(output_dir, "deseq2_results.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary of significant genes
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]
cat("Significant genes (padj < 0.05):", nrow(sig_genes), "\n")

# MA plot
pdf(file.path(output_dir, "deseq2_ma_plot.pdf"), width = 8, height = 6)
plotMA(res, main = "DESeq2 MA Plot")
dev.off()

# Volcano plot
pdf(file.path(output_dir, "deseq2_volcano_plot.pdf"), width = 10, height = 8)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DESeq2 Volcano Plot',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 3.0)
dev.off()

# Heatmap of top variable genes
cat("Generating heatmap of top variable genes...\n")
topVarGenes <- head(order(rowVars(assay(dds)), decreasing = TRUE), 50)
mat <- assay(dds)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[, c("treatment", "group")])

pdf(file.path(output_dir, "heatmap_top_variable_genes.pdf"), width = 10, height = 8)
pheatmap(mat, annotation_col = df, show_rownames = FALSE, 
         main = "Top 50 Most Variable Genes")
dev.off()

# edgeR Analysis
cat("Running edgeR analysis...\n")

# Create DGEList object
if ("treatment" %in% colnames(metadata)) {
  group <- metadata$treatment
} else {
  group <- metadata$group
}

dge <- DGEList(counts = count_matrix, group = group)

# Filter low count genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("Genes after edgeR filtering:", nrow(dge), "\n")

# Normalize
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)

# Fit model
fit <- glmQLFit(dge, design)

# Test for differential expression
qlf <- glmQLFTest(fit, coef = 2)

# Get results
edger_results <- topTags(qlf, n = nrow(dge))$table
edger_results$gene_id <- rownames(edger_results)

# Save edgeR results
write.table(edger_results, file.path(output_dir, "edger_results.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary
cat("edgeR significant genes (FDR < 0.05):", sum(edger_results$FDR < 0.05, na.rm = TRUE), "\n")

# Comparison plot
pdf(file.path(output_dir, "deseq2_vs_edger.pdf"), width = 8, height = 6)
common_genes <- intersect(res$gene_id, edger_results$gene_id)
res_common <- res[res$gene_id %in% common_genes, ]
edger_common <- edger_results[edger_results$gene_id %in% common_genes, ]

plot(res_common$log2FoldChange, edger_common$logFC,
     xlab = "DESeq2 log2FoldChange", ylab = "edgeR logFC",
     main = "DESeq2 vs edgeR Comparison",
     pch = 16, col = alpha("black", 0.6))
abline(0, 1, col = "red", lty = 2)
dev.off()

# Generate summary report
cat("Generating summary report...\n")
sink(file.path(output_dir, "analysis_summary.txt"))
cat("RNA-seq Differential Expression Analysis Summary\n")
cat("===============================================\n\n")
cat("Input files:\n")
cat("  Count matrix:", count_file, "\n")
cat("  Metadata:", metadata_file, "\n\n")
cat("Data summary:\n")
cat("  Total samples:", ncol(count_matrix), "\n")
cat("  Total genes:", nrow(count_matrix), "\n")
cat("  Genes after filtering (DESeq2):", nrow(dds), "\n")
cat("  Genes after filtering (edgeR):", nrow(dge), "\n\n")
cat("DESeq2 results:\n")
cat("  Significant genes (padj < 0.05):", nrow(sig_genes), "\n")
cat("  Upregulated genes:", sum(sig_genes$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated genes:", sum(sig_genes$log2FoldChange < 0, na.rm = TRUE), "\n\n")
cat("edgeR results:\n")
cat("  Significant genes (FDR < 0.05):", sum(edger_results$FDR < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated genes:", sum(edger_results$FDR < 0.05 & edger_results$logFC > 0, na.rm = TRUE), "\n")
cat("  Downregulated genes:", sum(edger_results$FDR < 0.05 & edger_results$logFC < 0, na.rm = TRUE), "\n\n")
cat("Output files:\n")
cat("  DESeq2 results: deseq2_results.tsv\n")
cat("  edgeR results: edger_results.tsv\n")
cat("  Quality control plots: qc_*.pdf\n")
cat("  Analysis plots: *_plot.pdf\n")
cat("  Heatmap: heatmap_top_variable_genes.pdf\n")
sink()

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("Check analysis_summary.txt for detailed results.\n")
