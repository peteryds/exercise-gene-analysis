#' @title Run Limma Differential Expression Analysis
#' @param expr Expression matrix (genes as rows, samples as columns).
#' @param group A factor variable for the experimental groups, matching
#'              the column order of 'expr'.
#' @return A data frame with the 'topTable' results from limma.
#' @import limma
#' @export
run_limma_analysis <- function(expr, group) {
  message("Running limma differential expression analysis...")
  
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", "PostVsPre")
  
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  
  # Extract comparison of interest (Post vs Pre)
  deg_table <- limma::topTable(fit, coef = "PostVsPre", n = Inf)
  
  message("Limma analysis complete.")
  return(deg_table)
}

#' @title Filter Significant Genes
#' @param deg_table Data frame from 'run_limma_analysis'.
#' @param p_cutoff Adjusted p-value cutoff.
#' @param fc_cutoff Absolute log Fold Change cutoff.
#' @return A data frame containing only significant genes.
#' @export
filter_significant_genes <- function(deg_table, p_cutoff = 0.05, fc_cutoff = 1) {
  sig_genes <- deg_table[deg_table$adj.P.Val < p_cutoff & 
                           abs(deg_table$logFC) > fc_cutoff, ]
  
  message("Found ", nrow(sig_genes), " significant genes (adj.P.Val < ", p_cutoff, 
          " & |logFC| > ", fc_cutoff, ")")
  
  return(sig_genes)
}

#' @title Run GO Enrichment Analysis
#' @param probe_ids A character vector of significant Probe IDs 
#'                  (e.g., from rownames(sig_genes)).
#' @param annotation_db The annotation database package name (e.g., "hgu133plus2.db").
#' @param org_db The organism database package name (e.g., "org.Hs.eg.db").
#' @return An 'enrichResult' object from clusterProfiler.
#' @import clusterProfiler
#' @importFrom AnnotationDbi mapIds
#' @export
run_go_enrichment <- function(probe_ids, annotation_db = "hgu133plus2.db", org_db = "org.Hs.eg.db") {
  
  # Ensure annotation packages are available
  if (!requireNamespace(annotation_db, quietly = TRUE)) {
    stop("Please install annotation package: ", annotation_db)
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop("Please install organism package: ", org_db)
  }
  
  message("Converting Probe IDs to Entrez IDs...")
  entrez_ids <- AnnotationDbi::mapIds(get(annotation_db),
                                      keys = probe_ids,
                                      keytype = "PROBEID",
                                      column = "ENTREZID",
                                      multiVals = "first")
  
  entrez_ids <- na.omit(entrez_ids)
  message(length(entrez_ids), " Entrez IDs successfully mapped.")
  
  if (length(entrez_ids) == 0) {
    warning("No Entrez IDs were available for GO analysis.")
    return(NULL)
  }
  
  message("Running GO (Biological Process) enrichment analysis...")
  ego <- clusterProfiler::enrichGO(gene         = entrez_ids,
                                   OrgDb        = get(org_db),
                                   ont          = "BP",
                                   readable     = TRUE,
                                   pvalueCutoff = 0.05)
  
  return(ego)
}