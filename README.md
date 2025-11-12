# Analysis of Exercise-Responsive Genes (GSE47881)

This project analyzes the GEO dataset [GSE47881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47881) to identify changes in gene expression in human skeletal muscle following exercise training.

## Project Goal

* To identify differentially expressed genes (DEGs) pre- vs. post-exercise using the `limma` package.
* To visualize the analysis results using a volcano plot, heatmap, and PCA.
* To perform GO enrichment analysis using `clusterProfiler` to understand the biological processes associated with these DEGs.

## How to Run This Analysis

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/peteryds/exercise-gene-analysis.git](https://github.com/peteryds/exercise-gene-analysis.git)
    cd exercise-gene-analysis
    ```

2.  **Open the Project:**
    Open this directory in RStudio (e.g., by double-clicking a `.Rproj` file if you create one).

3.  **Install Dependencies:**
    Run the following commands in the R console to install all required packages (as listed in the `DESCRIPTION` file):
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    # Install Bioconductor packages
    BiocManager::install(c("GEOquery", "limma", "pheatmap", "clusterProfiler", "org.Hs.eg.db", "hgu133plus2.db", "Biobase"))
    
    # Install CRAN packages
    install.packages(c("dplyr", "ggplot2"))
    ```

4.  **Run the Analysis:**
    Open and run the main script:
    ```R
    source("scripts/run_analysis.R")
    ```

## Summary of Findings

All outputs from the analysis are saved in the `results/` folder.

* **Differential Expression:** The analysis identified a set of genes significantly up- and down-regulated after exercise (see `results/data/limma_significant_genes.csv`). Key genes include PGC1A (mitochondrial biogenesis) and SOD2 (oxidative stress defense).
* **PCA:** The PCA plot shows a clear separation between pre- and post-training samples, indicating a significant, consistent shift in the overall gene expression profile.
* **GO Enrichment:** The upregulated genes are strongly enriched for pathways related to oxidative metabolism, muscle remodeling, and mitochondrial biogenesis (see `results/plots/go_enrichment_barplot.png`).