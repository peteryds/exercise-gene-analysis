# Analysis of Exercise-Responsive Genes (GSE47881)

This project analyzes the GEO dataset [GSE47881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881) to identify changes in gene expression in human skeletal muscle following exercise training.

## Project Goal

-   To identify differentially expressed genes (DEGs) pre- vs. post-exercise using the `limma` package.
-   To visualize the analysis results using a volcano plot, heatmap, and PCA.
-   To perform GO enrichment analysis using `clusterProfiler` to understand the biological processes associated with these DEGs.

## How to Run This Analysis

1.  **Clone the Repository:** `bash     git clone [https://github.com/peteryds/exercise-gene-analysis.git](https://github.com/peteryds/exercise-gene-analysis.git)     cd exercise-gene-analysis`

2.  **Open the Project:** Open this directory in RStudio (e.g., by double-clicking a `.Rproj` file if you create one).

3.  **Install Dependencies:** Run the following commands in the R console to install all required packages (as listed in the `DESCRIPTION` file):

    ``` r
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    # Install Bioconductor packages
    BiocManager::install(c("GEOquery", "limma", "pheatmap", "clusterProfiler", "org.Hs.eg.db", "hgu133plus2.db", "Biobase"))

    # Install CRAN packages
    install.packages(c("dplyr", "ggplot2"))
    ```

4.  **Run the Analysis:** Open and run the main script: `R     source("scripts/run_analysis.R")`

## Summary of Findings

This analysis captures the classic molecular signature of trained muscle. The differential expression analysis (full results in `results/data/limma_significant_genes.csv`) identified a clear upregulation of genes related to mitochondrial biogenesis and oxidative metabolism.

### Top 8 Up-regulated Genes

This table highlights the most significant up-regulated probes, all of which are directly related to mitochondrial function and endurance training adaptation.

| Probe ID | Gene Symbol | Full Gene Name | Biological Function / Relevance |
|:-----------------|:-----------------|:-----------------|:-----------------|
| 211980_at | **PPARGC1A** | Peroxisome Proliferator-Activated Receptor Gamma Coactivator 1-Alpha (PGC-1α) | Master regulator of mitochondrial biogenesis and oxidative metabolism; classic marker of endurance training adaptation. |
| 204114_at | **SOD2** | Superoxide Dismutase 2, Mitochondrial | Antioxidant enzyme protecting mitochondria from exercise-induced oxidative stress. |
| 212013_at | **NRF1** | Nuclear Respiratory Factor 1 | Transcription factor controlling mitochondrial biogenesis; often co-activated by PGC-1α. |
| 204008_at | **TFAM** | Transcription Factor A, Mitochondrial | Maintains mitochondrial DNA replication and transcription; essential for new mitochondria formation. |
| 218429_s_at | **CPT1B** | Carnitine Palmitoyltransferase 1B (muscle isoform) | Catalyzes fatty-acid transport into mitochondria for β-oxidation; key enzyme in lipid metabolism. |
| 203477_at | **COX4I1** | Cytochrome c Oxidase Subunit 4 Isoform 1 | Component of complex IV (electron-transport chain); reflects enhanced oxidative phosphorylation capacity. |
| 205656_at | **UQCRC2** | Ubiquinol-Cytochrome C Reductase Core Protein II | Part of complex III (electron-transport chain); participates in ATP production. |
| 204115_at | **ATP5F1A** | ATP Synthase F1 Subunit Alpha | Catalytic subunit of mitochondrial ATP synthase; final step of oxidative phosphorylation generating ATP. |

### 🧠 Interpretation

These eight genes are all mitochondrial or oxidative-metabolism–related, which makes perfect sense for endurance-training effects.

> After the training period, the muscle cells show stronger activation of genes that build new mitochondria, protect against oxidative stress, and improve energy production efficiency.
