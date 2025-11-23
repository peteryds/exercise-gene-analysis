# =============================================================
# Project: Exercise-Responsive Gene Analysis (GSE47881)
# Main analysis pipeline
# =============================================================

# 1️⃣ Load Functions and Libraries ----------------------------------
# Load all custom functions from the R/ directory
source("R/01_data_loading.R")
source("R/02_analysis.R")
source("R/03_plotting.R")

# Load all packages (declared in DESCRIPTION)
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

# 2️⃣ Set Parameters and Paths --------------------------------------

GSE_ID <- "GSE47881"
P_CUTOFF <- 0.05
FC_CUTOFF <- 1
TOP_N_GENES <- 50

# Create output directories
dir.create("results/data", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# 3️⃣ Run Analysis Pipeline -------------------------------------------

# 3.1 Load and Prepare Data
data_list <- load_and_prep_gse(GSE_ID)
expr_matrix <- data_list$expr
group_factor <- data_list$group

# 3.2 Differential Expression Analysis
deg_results <- run_limma_analysis(expr_matrix, group_factor)
sig_genes <- filter_significant_genes(deg_results, P_CUTOFF, FC_CUTOFF)

# Save DE results
write.csv(deg_results, "results/data/limma_all_genes.csv")
write.csv(sig_genes, "results/data/limma_significant_genes.csv")

# 3.3 Plot Volcano
p_volcano <- plot_volcano(deg_results, P_CUTOFF, FC_CUTOFF)
ggsave("results/plots/volcano_plot.png", p_volcano, width = 8, height = 6)

# 3.4 Plot Heatmap (Top 50 significant genes)
top_probes <- head(rownames(sig_genes), 50)
plot_heatmap(
  expr_matrix[top_probes, ],
  title = "Heatmap of Top 50 DEGs (Post vs Pre)",
  filename = "results/plots/heatmap_top50_degs.png"
)

# 3.5 Plot PCA (All Genes)
p_pca_all <- plot_pca(
  expr_matrix, 
  group_factor, 
  "PCA of All Genes (Pre vs Post)"
)
ggsave("results/plots/pca_all_genes.png", p_pca_all, width = 7, height = 5)

# 3.6 Plot PCA (Significant Genes)
p_pca_sig <- plot_pca(
  expr_matrix[rownames(sig_genes), ], 
  group_factor, 
  "PCA of Significant Genes (Pre vs Post)"
)
ggsave("results/plots/pca_significant_genes.png", p_pca_sig, width = 7, height = 5)

# 3.7 GO Enrichment Analysis
if (nrow(sig_genes) >= 2) {
  ego <- run_go_enrichment(
    probe_ids = rownames(sig_genes),
    annotation_db = "hgu133plus2.db",
    org_db = "org.Hs.eg.db"
  )

  # Save and plot GO results
  if (!is.null(ego)) {
    saveRDS(ego, "results/data/go_enrichment_result.rds")
    write.csv(as.data.frame(ego), "results/data/go_enrichment_result.csv")
    
    # Plot GO Barplot
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

# New: Violin plot ---------------------------------------------------
if (nrow(sig_genes) >= 2) 
  {
  top_probes <- head(rownames(sig_genes), TOP_N_GENES)
  p_violin <- plot_violin_degs( expr_matrix[top_probes, ], group_factor,
    title = sprintf("Expression Distribution (Top %d DEGs)", min(TOP_N_GENES, nrow(sig_genes))))
  ggsave("results/plots/violin_top_degs.png",p_violin,width = 10,height = 8)
  
} else {
  message("Not enough DEGs for violin plot.")
}

# 4️⃣ Completion ----------------------------------------------------
message("Analysis complete! All results saved to the 'results/' directory.")

# =============================================================
# End of Script
# =============================================================