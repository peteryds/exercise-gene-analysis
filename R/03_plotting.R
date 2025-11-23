#' @title Plot Volcano Plot
#' @param deg_table Data frame from 'run_limma_analysis'.
#' @param p_cutoff Adjusted p-value cutoff for highlighting.
#' @param fc_cutoff Absolute log Fold Change cutoff for highlighting.
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_volcano <- function(deg_table, p_cutoff = 0.05, fc_cutoff = 1) {
  deg_table$threshold <- ifelse(
    deg_table$adj.P.Val < p_cutoff & abs(deg_table$logFC) > fc_cutoff,
    "Significant",
    "Not Significant"
  )
  
  p <- ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = threshold), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey80")) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    labs(title = "Volcano Plot: Post vs Pre Training",
         x = "log2 Fold Change",
         y = "-log10 Adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' @title Plot Heatmap
#' @param expr_subset The expression matrix subset to plot.
#' @param title The main title for the heatmap.
#' @param filename The path to save the plot (e.g., "results/plots/heatmap.png").
#' @import pheatmap
#' @export
plot_heatmap <- function(expr_subset, title, filename) {
  pheatmap::pheatmap(
    expr_subset,
    scale = "row",
    show_rownames = TRUE,
    show_colnames = FALSE, # Usually too many samples to show
    main = title,
    filename = filename
  )
}

#' @title Plot PCA
#' @param expr Expression matrix (genes as rows, samples as columns) for PCA.
#' @param group A factor variable for coloring the points.
#' @param title The main title for the plot.
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats prcomp
#' @export
plot_pca <- function(expr, group, title) {
  expr_t <- t(expr)
  pca <- stats::prcomp(expr_t, scale. = TRUE)
  pca_df <- data.frame(pca$x, group = group)
  
  # Calculate percent variance
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  p <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = title,
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    )
  
  return(p)
}

#' @title Plot GO Barplot
#' @param ego An 'enrichResult' object from 'run_go_enrichment'.
#' @param n_categories Number of categories to show.
#' @param title The main title for the plot.
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_go_barplot <- function(ego, n_categories = 10, title = "GO Biological Processes") {
  if (is.null(ego) || nrow(ego@result) == 0) {
    message("No GO results to plot.")
    return(NULL)
  }
  
  # clusterProfiler's barplot function returns a ggplot object
  p <- barplot(ego, showCategory = n_categories, title = title)
  
  return(p)
}

#' @title Plot Violin Plot for DEGs
#' @description Creates a faceted violin plot showing expression distributions 
#'              for the top differentially expressed genes.
#' @param expr_subset Expression matrix (genes × samples), usually top N DEGs.
#' @param group_factor Factor vector indicating sample condition (e.g., pre/post).
#' @param title Title for the plot.
#' @return A ggplot object with faceted violin plots.
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @export
plot_violin_degs <- function(expr_subset, group_factor, title = "Violin Plot of Top DEGs") {
  
  # Convert expression matrix to long tidy format
  df <- expr_subset %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    )
  

  # Match condition for each sample
  if (is.null(names(group_factor))) {
    stop("group_factor must be a named vector with sample IDs as names")
  }
  missing_samples <- setdiff(df$Sample, names(group_factor))
  if (length(missing_samples) > 0) {
    stop(
      paste0(
        "Some samples in expr_subset do not have corresponding entries in group_factor: ",
        paste(missing_samples, collapse = ", ")
      )
    )
  }
  df$Condition <- group_factor[df$Sample]
  df$Condition <- factor(df$Condition,
                         levels = levels(group_factor))
 
  
  # Violin plot with jitter & facet per gene
  p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    scale_fill_manual(values = c(
      "pre-training" = "#1f78b4",
      "post-training" = "#e31a1c"
    )) +
    labs(
      title = title,
      x = "Condition",
      y = "Normalized Expression"
    ) +
    theme_minimal() +
    theme(
      strip.text.x = element_text(size = 7, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 10)
  
  return(p)
}


