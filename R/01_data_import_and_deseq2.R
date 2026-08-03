# Differential Gene Expression Analysis with DESeq2

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -----------------------------
# 1. Set input and output paths
# -----------------------------
count_file <- "results/counts/combined_gene_count_matrix.tsv"
metadata_file <- "data/metadata/sample_info.tsv"

out_dir <- "results/deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Read count matrix and metadata
# -----------------------------
counts <- read.delim(count_file,
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

metadata <- read.delim(metadata_file,
                       header = TRUE,
                       check.names = FALSE)

# -----------------------------
# 3. Clean sample names if needed
# -----------------------------
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- gsub(".sorted.bam", "", colnames(counts), fixed = TRUE)
colnames(counts) <- gsub(".bam", "", colnames(counts), fixed = TRUE)
colnames(counts) <- gsub("_forward", "", colnames(counts), fixed = TRUE)
colnames(counts) <- gsub("_reverse", "", colnames(counts), fixed = TRUE)

# -----------------------------
# 4. Match metadata to count matrix columns
# -----------------------------
metadata <- metadata[metadata$SampleID %in% colnames(counts), ]

if (nrow(metadata) == 0) {
  stop("No matching SampleID values found between metadata and count matrix.")
}

rownames(metadata) <- metadata$SampleID
counts <- counts[, rownames(metadata)]

if (!all(colnames(counts) == rownames(metadata))) {
  stop("Count matrix columns and metadata rows are not in the same order.")
}

# -----------------------------
# 5. Prepare condition factor
# -----------------------------
metadata$Condition <- factor(metadata$Condition)

# Set control/reference condition
metadata$Condition <- relevel(metadata$Condition, ref = "NonInfected_Day0")

# -----------------------------
# 6. Ensure counts are integers
# -----------------------------
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# -----------------------------
# 7. Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Condition)

# -----------------------------
# 8. Pre-filter very low-count genes
# -----------------------------
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# -----------------------------
# 9. Run DESeq2
# -----------------------------
dds <- DESeq(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts,
          file = file.path(out_dir, "normalized_counts_DESeq2.csv"))

# -----------------------------
# 10. VST-normalised count matrix
# -----------------------------
vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)

write.csv(vst_counts,
          file = file.path(out_dir, "vst_normalized_count_matrix.csv"))

# -----------------------------
# 11. PCA plot
# -----------------------------
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = Condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave(file.path(out_dir, "PCA_plot_vst.png"), plot = p, width = 7, height = 5)

# -----------------------------
# 12. Sample distance heatmap
# -----------------------------
sample_dists <- dist(t(vst_counts))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst_counts)
colnames(sample_dist_matrix) <- colnames(vst_counts)

png(file.path(out_dir, "sample_distance_heatmap.png"),
    width = 900, height = 800)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample Distance Heatmap")
dev.off()

# -----------------------------
# 13. Function to export DESeq2 results
# -----------------------------
export_deseq_results <- function(dds, contrast_condition, reference_condition) {
  
  contrast_name <- paste0(contrast_condition, "_vs_", reference_condition)
  
  res <- results(dds,
                 contrast = c("Condition", contrast_condition, reference_condition),
                 alpha = 0.05)
  
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]
  res_df <- res_df[order(res_df$padj), ]
  
  # Significant DEGs: adjusted p-value < 0.05 and absolute log2FC >= 1
  sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
  up <- subset(sig, log2FoldChange > 0)
  down <- subset(sig, log2FoldChange < 0)
  
  write.csv(res_df,
            file = file.path(out_dir, paste0("DESeq2_all_results_", contrast_name, ".csv")),
            row.names = FALSE)
  
  write.csv(sig,
            file = file.path(out_dir, paste0("DESeq2_significant_DEGs_", contrast_name, ".csv")),
            row.names = FALSE)
  
  summary_row <- data.frame(
    Contrast = contrast_name,
    Total_genes_tested = nrow(res_df),
    Significant_DEGs = nrow(sig),
    Upregulated = nrow(up),
    Downregulated = nrow(down)
  )
  
  return(summary_row)
}

# -----------------------------
# 14. Run comparisons
# -----------------------------
summary_day4 <- export_deseq_results(dds,
                                     contrast_condition = "InfluenzaA_Day4",
                                     reference_condition = "NonInfected_Day0")

summary_day8 <- export_deseq_results(dds,
                                     contrast_condition = "InfluenzaA_Day8",
                                     reference_condition = "NonInfected_Day0")

deg_summary <- rbind(summary_day4, summary_day8)

write.csv(deg_summary,
          file = file.path(out_dir, "DEG_summary.csv"),
          row.names = FALSE)

print(deg_summary)

# -----------------------------
# 15. Save session information
# -----------------------------
sink(file.path(out_dir, "R_sessionInfo.txt"))
sessionInfo()
sink()

cat("DESeq2 analysis completed successfully.\n")
cat("Results saved in:", out_dir, "\n")