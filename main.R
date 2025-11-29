# ============================================================
# Gene Expression Analysis: Resistance Training & Age
# Main Execution Script (Final Updated Version)
# ============================================================

# 1. Load Environment & Functions
source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")  # process_gene_data, clean_and_normalize_data
source("R/models.R")   # run_limma_screening, run_limma_interaction
source("R/eda.R")      # plotting functions
source("R/facetplot.R")  # facet_scatter_interaction, facet_violin_interaction

# Setup libraries
setup_environment()

# ------------------------------------------------------------
# 2. Setup Directories
# ------------------------------------------------------------
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)

dir_raw <- file.path(output_dir, "QC_1_Raw")
if (!dir.exists(dir_raw)) dir.create(dir_raw)

dir_clean <- file.path(output_dir, "QC_2_Cleaned")
if (!dir.exists(dir_clean)) dir.create(dir_clean)

dir_plots <- file.path(output_dir, "Plots_Interaction")
if (!dir.exists(dir_plots)) dir.create(dir_plots)

dir_facet <- file.path(output_dir, "Plots_Facet")
if (!dir.exists(dir_facet)) dir.create(dir_facet)


# ============================================================
# 3. Data Loading & Phase 1 QC (Raw Data)
# ============================================================
message("\n=== STEP 1: Loading Raw Data ===")
eset_raw <- get_geo_data(gse_id = "GSE47881")

message("Generating Phase 1 QC Plots (Raw Data)...")

# PCA (Raw)
p_pca_raw <- plot_pca(eset_raw, title = "PCA: Raw Data (Before QC)")
ggsave(file.path(dir_raw, "QC_Raw_PCA.png"), p_pca_raw, width = 8, height = 6)

# Heatmap (Raw)
png(file.path(dir_raw, "QC_Raw_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_raw)
dev.off()

# Density (Raw)
p_dens_raw <- plot_density(eset_raw) +
  labs(subtitle = "Raw Data: Note Linear Scale & Outlier Spike")
ggsave(file.path(dir_raw, "QC_Raw_Density.png"), p_dens_raw, width = 8, height = 6)


# ============================================================
# 4. Data Cleaning & Normalization
# ============================================================
message("\n=== STEP 2: Data Cleaning & Normalization ===")

eset_clean <- clean_and_normalize_data(eset_raw, outliers_to_remove = NULL)


# ============================================================
# 5. Phase 2 QC (Cleaned Data)
# ============================================================
message("Generating Phase 2 QC Plots (Cleaned Data)...")

# PCA (Clean)
p_pca_clean <- plot_pca(eset_clean, title = "PCA: Cleaned & Log2 Transformed")
ggsave(file.path(dir_clean, "QC_Clean_PCA.png"), p_pca_clean, width = 8, height = 6)

# Heatmap (Clean)
png(file.path(dir_clean, "QC_Clean_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_clean)
dev.off()

# Density (Clean)
p_dens_clean <- plot_density(eset_clean) +
  labs(subtitle = "Cleaned Data: Log2 Scale & Normalized")
ggsave(file.path(dir_clean, "QC_Clean_Density.png"), p_dens_clean, width = 8, height = 6)


# ============================================================
# 6. Statistical Analysis (Limma)
# ============================================================
message("\n=== STEP 3: Running Limma Models ===")

# A. Main Effect (Pre vs Post)
limma_res <- run_limma_screening(eset_clean, p_cutoff = 0.05)

# B. Interaction Effect (Age * Time)
limma_int <- run_limma_interaction(eset_clean, p_cutoff = 0.05)

# Save Results
write.csv(limma_res$full_results, file.path(output_dir, "Results_Main_Effect.csv"))
write.csv(limma_res$sig_genes_df, file.path(output_dir, "Results_Main_Effect_Sig.csv"))
write.csv(limma_int$full_results, file.path(output_dir, "Results_Interaction_Age.csv"))
write.csv(limma_int$sig_genes_df, file.path(output_dir, "Results_Interaction_Age_Sig.csv"))


# ============================================================
# 7. Visualization of Results
# ============================================================
message("\n=== STEP 4: Visualizing Significant Findings ===")

# 7.1 Volcano (Main Effect)
volcano_plot <- plot_volcano(limma_res$full_results, p_cutoff = 0.05)
ggsave(file.path(output_dir, "Volcano_Main_Effect.png"), volcano_plot, width = 8, height = 6)


# 7.2 Individual Plots for Top Interaction Genes
final_df_viz <- process_gene_data(eset_clean)

# Select top genes
top_genes <- rownames(limma_int$sig_genes_df)
if (length(top_genes) == 0) {
  message("No significant interaction genes found (FDR < 0.05). Plotting top 10 by P-value.")
  top_genes <- rownames(head(limma_int$full_results, 10))
} else {
  top_genes <- head(top_genes, 20)
}

message(paste("Generating single-gene plots for", length(top_genes), "interaction candidates..."))

for (gene in top_genes) {
  
  # Scatter (Age vs Change)
  p_scatter <- plot_gene_age_scatter(final_df_viz, gene)
  if (!is.null(p_scatter)) {
    p_scatter <- p_scatter + labs(subtitle = "Selected via Limma Interaction (Age * Time)")
    ggsave(file.path(dir_plots, paste0("Scatter_", gene, ".png")),
           p_scatter, width = 5, height = 4)
  }
  
  # Violin (Pre vs Post)
  p_violin <- plot_gene_violin(eset_clean, gene)
  if (!is.null(p_violin)) {
    ggsave(file.path(dir_plots, paste0("Violin_", gene, ".png")),
           p_violin, width = 5, height = 4)
  }
}


# ============================================================
# 7.3 Combined FACET Plots
# ============================================================
message("Generating FACET plots for top interaction genes...")

# FACET SCATTER
p_facet_scatter <- facet_scatter_interaction(eset_clean, top_genes)
ggsave(
  file.path(dir_facet, "Facet_Scatter_TopGenes.png"),
  plot = p_facet_scatter, width = 14, height = 10
)

# FACET VIOLIN
p_facet_violin <- facet_violin_interaction(eset_clean, top_genes)
ggsave(
  file.path(dir_facet, "Facet_Violin_TopGenes.png"),
  plot = p_facet_violin, width = 14, height = 10
)


# ============================================================
# FINISH
# ============================================================
message("\n[DONE] Pipeline Finished Successfully!")
message(paste("Check output directory:", output_dir))

