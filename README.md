
# RNA-Seq Differential Expression and Pathway Analysis Project

## Name

Theresia Ngonda Nfah

## Program Title

Bioinformatics Training: RNA-seq Data Analysis Project

## Project Topic

Differential gene expression analysis of mouse cerebellum during influenza A infection.

## Introduction

RNA sequencing provides a genome-wide method for measuring gene expression under different biological conditions. This project analysed publicly available mouse RNA-seq data to identify genes whose expression changed following Influenza A infection. The analysis compared non-infected Day 0 samples with Influenza A-infected samples at Day 4 and Day 8.

The purpose of the project was to identify differentially expressed genes, visualise expression patterns, and interpret biological meaning using Gene Ontology Biological Process and KEGG pathway enrichment analysis.

## Methods Summary

Six RNA-seq samples from NCBI GEO/SRA series GSE96870 were used. The workflow included raw read download, quality control using FastQC and MultiQC, read trimming using fastp, alignment to the Mus musculus GRCm39 reference genome using HISAT2, read counting using featureCounts, differential gene expression analysis using DESeq2, and pathway enrichment analysis using clusterProfiler.

The DESeq2 model used the design formula:

`~ Condition`

The following comparisons were performed:

- InfluenzaA_Day4 vs NonInfected_Day0
- InfluenzaA_Day8 vs NonInfected_Day0

Significant differentially expressed genes were defined as genes with adjusted p-value < 0.05 and absolute log2 fold change ≥ 1.

## Differential Gene Expression Results

| Contrast | Total genes tested | Significant DEGs | Upregulated | Downregulated |
|---|---:|---:|---:|---:|
| InfluenzaA_Day4_vs_NonInfected_Day0 | 29722 | 21 | 19 | 2 |
| InfluenzaA_Day8_vs_NonInfected_Day0 | 29722 | 291 | 174 | 117 |

The Day 8 infected group showed more differentially expressed genes than the Day 4 infected group, suggesting a stronger transcriptional response at Day 8.

## Visualisation Results

The following plots were generated:

- PCA plot
- Sample distance heatmap
- Volcano plots
- Top DEG heatmaps
- GO Biological Process dotplots
- KEGG pathway dotplots

### PCA Plot

![PCA Plot](results/visualization_pathway/figures/PCA_plot.png)

### Sample Distance Heatmap

![Sample Distance Heatmap](results/visualization_pathway/figures/sample_distance_heatmap.png)

### Volcano Plot: InfluenzaA Day 4 vs NonInfected Day 0

![Volcano Day 4](results/visualization_pathway/figures/volcano_InfluenzaA_Day4_vs_NonInfected_Day0.png)

### Volcano Plot: InfluenzaA Day 8 vs NonInfected Day 0

![Volcano Day 8](results/visualization_pathway/figures/volcano_InfluenzaA_Day8_vs_NonInfected_Day0.png)

## GO Biological Process Enrichment

The top 5 enriched GO Biological Process terms are reported in:

`results/visualization_pathway/visualization_DEG_GO_KEGG_summary.csv`

Full GO enrichment result files are available in:

- `results/visualization_pathway/enrichment/GO_BP_enrichment_InfluenzaA_Day4_vs_NonInfected_Day0.csv`
- `results/visualization_pathway/enrichment/GO_BP_enrichment_InfluenzaA_Day8_vs_NonInfected_Day0.csv`

## KEGG Pathway Enrichment

The top 3 enriched KEGG pathways are reported in:

`results/visualization_pathway/visualization_DEG_GO_KEGG_summary.csv`

Full KEGG enrichment result files are available in:

- `results/visualization_pathway/enrichment/KEGG_enrichment_InfluenzaA_Day4_vs_NonInfected_Day0.csv`
- `results/visualization_pathway/enrichment/KEGG_enrichment_InfluenzaA_Day8_vs_NonInfected_Day0.csv`

## Brief Discussion

The analysis identified changes in gene expression associated with Influenza A infection in mouse cerebellum samples. The number of significant differentially expressed genes increased from Day 4 to Day 8, indicating that infection-associated transcriptional changes were more pronounced at the later time point.

The PCA plot and sample distance heatmap were used to assess sample-level relationships. Volcano plots showed the distribution of upregulated and downregulated genes, while top DEG heatmaps showed expression patterns among the most significant genes.

GO Biological Process and KEGG pathway enrichment analyses provided functional interpretation of the significant genes. These enriched biological processes and pathways may reflect host responses to Influenza A infection, including immune-related and infection-response mechanisms.

## Conclusion

This project successfully completed a reproducible RNA-seq workflow from count matrix generation to differential expression analysis, visualisation, pathway enrichment, and final reporting. The repository contains scripts, result tables, figures, enrichment outputs, and final report files.

## Key Output Files

- `scripts/10_visualization_pathway_analysis.R`
- `results/visualization_pathway/visualization_DEG_GO_KEGG_summary.csv`
- `results/visualization_pathway/vst_normalized_count_matrix_visualization.csv`
- `results/visualization_pathway/figures/`
- `results/visualization_pathway/enrichment/`
- `reports/final_report.Rmd`
- `docs/final_report.html`

## References

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.

Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284–287.

The Gene Ontology Consortium. (2023). The Gene Ontology knowledgebase in 2023. Genetics, 224(1), iyad031.

Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M., & Ishiguro-Watanabe, M. (2023). KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research, 51(D1), D587–D592.
