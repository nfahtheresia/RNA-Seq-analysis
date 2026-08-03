# Visualization, GO/KEGG Enrichment, and Report Preparation

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

# -----------------------------
# 1. Paths
# -----------------------------
count_file <- "results/counts/combined_gene_count_matrix.tsv"
metadata_file <- "data/metadata/sample_info.tsv"

out_dir <- "results/visualization"
fig_dir <- file.path(out_dir, "figures")
enrich_dir <- file.path(out_dir, "enrichment")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

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

# Clean sample names
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- gsub(".sorted.bam", "", colnames(counts), fixed = TRUE)
colnames(counts) <- gsub(".bam", "", colnames(counts), fixed = TRUE)

# Match count matrix with metadata
metadata <- metadata[metadata$SampleID %in% colnames(counts), ]
rownames(metadata) <- metadata$SampleID
counts <- counts[, rownames(metadata)]

if (!all(colnames(counts) == rownames(metadata))) {
  stop("Sample names in count matrix and metadata do not match.")
}

metadata$Condition <- factor(metadata$Condition)
metadata$Condition <- relevel(metadata$Condition, ref = "NonInfected_Day0")

counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# -----------------------------
# 3. Run DESeq2 again for reproducibility
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

write.csv(vst_mat,
          file = file.path(out_dir, "vst_normalized_count_matrix_visualization.csv"))

# -----------------------------
# 4. PCA plot
# -----------------------------
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot of VST-normalised Counts")

ggsave(file.path(fig_dir, "PCA_plot.png"), pca_plot, width = 7, height = 5, dpi = 300)

# -----------------------------
# 5. Sample distance heatmap
# -----------------------------
sample_dists <- dist(t(vst_mat))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst_mat)
colnames(sample_dist_matrix) <- colnames(vst_mat)

png(file.path(fig_dir, "sample_distance_heatmap.png"), width = 900, height = 800)
pheatmap(sample_dist_matrix,
         main = "Sample Distance Heatmap",
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists)
dev.off()

# -----------------------------
# 6. Helper function for DEG plots and enrichment
# -----------------------------
analyse_contrast <- function(contrast_condition, reference_condition = "NonInfected_Day0") {
  
  contrast_name <- paste0(contrast_condition, "_vs_", reference_condition)
  message("Analysing: ", contrast_name)
  
  res <- results(dds,
                 contrast = c("Condition", contrast_condition, reference_condition),
                 alpha = 0.05)
  
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]
  res_df <- res_df[order(res_df$padj), ]
  
  res_df$Significance <- "Not significant"
  res_df$Significance[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        res_df$log2FoldChange >= 1] <- "Upregulated"
  res_df$Significance[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        res_df$log2FoldChange <= -1] <- "Downregulated"
  
  sig_df <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
  
  write.csv(res_df,
            file.path(out_dir, paste0("DESeq2_all_results_", contrast_name, ".csv")),
            row.names = FALSE)
  
  write.csv(sig_df,
            file.path(out_dir, paste0("DESeq2_significant_DEGs_", contrast_name, ".csv")),
            row.names = FALSE)
  
  # -----------------------------
  # Volcano plot
  # -----------------------------
  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    ggtitle(paste("Volcano Plot:", contrast_name)) +
    xlab("log2 Fold Change") +
    ylab("-log10 Adjusted p-value")
  
  ggsave(file.path(fig_dir, paste0("volcano_", contrast_name, ".png")),
         volcano, width = 7, height = 5, dpi = 300)
  
  # -----------------------------
  # Top DEG heatmap
  # -----------------------------
  top_genes <- sig_df %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(Geneid)
  
  if (length(top_genes) >= 2) {
    heat_mat <- vst_mat[top_genes, , drop = FALSE]
    heat_mat <- heat_mat - rowMeans(heat_mat)
    
    annotation_col <- data.frame(Condition = metadata$Condition)
    rownames(annotation_col) <- rownames(metadata)
    
    png(file.path(fig_dir, paste0("top_DEG_heatmap_", contrast_name, ".png")),
        width = 1000, height = 1000)
    pheatmap(heat_mat,
             annotation_col = annotation_col,
             show_rownames = TRUE,
             fontsize_row = 6,
             main = paste("Top 30 DEGs:", contrast_name))
    dev.off()
  }
  
  # -----------------------------
  # Gene ID conversion: Ensembl to Entrez
  # -----------------------------
  sig_genes <- sig_df$Geneid
  sig_genes_clean <- gsub("\\..*$", "", sig_genes)
  
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = sig_genes_clean,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  
  entrez_ids <- unique(na.omit(entrez_ids))
  
  # -----------------------------
  # GO Biological Process enrichment
  # -----------------------------
  go_df <- data.frame()
  kegg_df <- data.frame()
  
  if (length(entrez_ids) >= 5) {
    
    ego <- enrichGO(gene = entrez_ids,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
    
    go_df <- as.data.frame(ego)
    
    write.csv(go_df,
              file.path(enrich_dir, paste0("GO_BP_enrichment_", contrast_name, ".csv")),
              row.names = FALSE)
    
    if (nrow(go_df) > 0) {
      png(file.path(fig_dir, paste0("GO_BP_dotplot_", contrast_name, ".png")),
          width = 1000, height = 800)
      print(dotplot(ego, showCategory = 10) +
              ggtitle(paste("GO Biological Process:", contrast_name)))
      dev.off()
    }
    
    # -----------------------------
    # KEGG enrichment
    # Mouse KEGG organism code = mmu
    # -----------------------------
    ekegg <- enrichKEGG(gene = entrez_ids,
                        organism = "mmu",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH")
    
    kegg_df <- as.data.frame(ekegg)
    
    write.csv(kegg_df,
              file.path(enrich_dir, paste0("KEGG_enrichment_", contrast_name, ".csv")),
              row.names = FALSE)
    
    if (nrow(kegg_df) > 0) {
      png(file.path(fig_dir, paste0("KEGG_dotplot_", contrast_name, ".png")),
          width = 1000, height = 800)
      print(dotplot(ekegg, showCategory = 10) +
              ggtitle(paste("KEGG Pathway Enrichment:", contrast_name)))
      dev.off()
    }
  }
  
  summary_row <- data.frame(
    Contrast = contrast_name,
    Total_genes_tested = nrow(res_df),
    Significant_DEGs = nrow(sig_df),
    Upregulated = sum(sig_df$log2FoldChange > 0),
    Downregulated = sum(sig_df$log2FoldChange < 0),
    Top_5_GO_BP_terms = ifelse(nrow(go_df) > 0,
                               paste(head(go_df$Description, 5), collapse = "; "),
                               "No significant GO BP terms"),
    Top_3_KEGG_pathways = ifelse(nrow(kegg_df) > 0,
                                 paste(head(kegg_df$Description, 3), collapse = "; "),
                                 "No significant KEGG pathways")
  )
  
  return(summary_row)
}

# -----------------------------
# 7. Run analyses for both contrasts
# -----------------------------
summary_day4 <- analyse_contrast("InfluenzaA_Day4")
summary_day8 <- analyse_contrast("InfluenzaA_Day8")

visualization_summary <- bind_rows(summary_day4, summary_day8) 

write.csv(visualization_summary,
          file.path(out_dir, "visualization_DEG_GO_KEGG_summary.csv"),
          row.names = FALSE)

print(visualization_summary)

# -----------------------------
# 8. Save session information
# -----------------------------
sink(file.path(out_dir, "visualization_sessionInfo.txt"))
sessionInfo()
sink()

message("visualization and pathway analysis completed.")
message("Outputs saved in: ", out_dir)
