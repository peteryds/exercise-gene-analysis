# =============================================================
# Project: Exercise-Responsive Gene Analysis (GSE47881)
# Main analysis pipeline (with violin plot instead of heatmap)
# =============================================================

# 1️⃣ Load Functions and Libraries ----------------------------------
# Load custom functions
source("R/01_data_loading.R")
source("R/02_analysis.R")
source("R/03_plotting.R")

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(reshape2)
library(Biobase)

message("All functions and libraries loaded.")

# 2️⃣ Parameters ------------------------------------------------------

GSE_ID <- "GSE47881"
P_CUTOFF <- 0.05
FC_CUTOFF <- 1
TOP_N_GENES <- 50

dir.create("results/data", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# 3️⃣ Pipeline --------------------------------------------------------

# Load and prep data
data_list <- load_and_prep_gse(GSE_ID)
expr_matrix <- data_list$expr
group_factor <- data_list$group

# Differential expression
deg_results <- run_limma_analysis(expr_matrix, group_factor)
sig_genes <- filter_significant_genes(deg_results, P_CUTOFF, FC_CUTOFF)

write.csv(deg_results, "results/data/limma_all_genes.csv")
write.csv(sig_genes, "results/data/limma_significant_genes.csv")

# Volcano plot
p_volcano <- plot_volcano(deg_results, P_CUTOFF, FC_CUTOFF)
ggsave("results/plots/volcano_plot.png", p_volcano, width = 8, height = 6)

# New: Violin plot ---------------------------------------------------
if (nrow(sig_genes) >= 2) 
  {
  top_probes <- head(rownames(sig_genes), TOP_N_GENES)
  p_violin <- plot_violin_degs( expr_matrix[top_probes, ], group_factor,
    title = sprintf("Expression Distribution (Top %d DEGs)", min(TOP_N_GENES, nrow(sig_genes))))
  ggsave("results/plots/violin_top_degs.png",p_violin,width = 10,height = 8)
  
} else {
  message(" Not enough DEGs for violin plot.")
}

# PCA – all genes
p_pca_all <- plot_pca(expr_matrix, group_factor, "PCA of All Genes")
ggsave("results/plots/pca_all_genes.png", p_pca_all, width = 7, height = 5)

# PCA – significant genes
if (nrow(sig_genes) >= 2) {
  p_pca_sig <- plot_pca(expr_matrix[rownames(sig_genes), ], group_factor, "PCA of DEGs")
  ggsave("results/plots/pca_significant_genes.png", p_pca_sig, width = 7, height = 5)
}

# GO analysis
if (nrow(sig_genes) >= 2) {
  ego <- run_go_enrichment(
    probe_ids = rownames(sig_genes),
    annotation_db = "hgu133plus2.db",
    org_db = "org.Hs.eg.db"
  )
  if (!is.null(ego)) {
    saveRDS(ego, "results/data/go_enrichment_result.rds")
    write.csv(as.data.frame(ego), "results/data/go_enrichment_result.csv")
    
    p_go <- plot_go_barplot(
      ego,
      n_categories = 15,
      title = "GO Biological Processes Enriched by Exercise"
    )
    if (!is.null(p_go)) {
      ggsave("results/plots/go_enrichment_barplot.png", p_go, width = 10, height = 8)
    }
  }
}

message("Analysis complete! Results saved in 'results/'.")

