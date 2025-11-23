#' @title Load and Prepare GSE Data
#' @description Downloads a specified GSE dataset from GEO and extracts
#'              the expression matrix and grouping information.
#' @param gse_id The GEO dataset ID to download (e.g., "GSE47881").
#' @return A list containing 'expr' (expression matrix), 
#'         'group' (grouping factor), and 'pdata' (phenotypic data).
#' @import GEOquery
#' @import Biobase
#' @export
load_and_prep_gse <- function(gse_id = "GSE47881") {
  message("Downloading data from GEO: ", gse_id)
  
  gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (length(gse) == 0) {
    stop("Could not find GSE dataset: ", gse_id)
  }
  gse <- gse[[1]]
  
  # Extract data
  expr <- Biobase::exprs(gse)
  pdata <- Biobase::pData(gse)
  
  message("Expression matrix dimensions: ", paste(dim(expr), collapse = " x "))
  
  # Define groups based on 'time:ch1'
  if (!"time:ch1" %in% colnames(pdata)) {
    stop("Column 'time:ch1' not found in pData.")
  }
  group <- factor(pdata$`time:ch1`,
                  levels = c("pre-training", "post-training"))
  
  # Attach sample names (GSM IDs)
  names(group) <- rownames(pdata)
  
  if (any(is.na(group))) {
    warning("NA values present in grouping variable.")
  }
  
  message("Group information:")
  print(table(group))
  
  return(list(expr = expr, group = group, pdata = pdata))
}