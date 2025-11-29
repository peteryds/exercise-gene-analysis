# ============================================================
# Facet Plot Functions for Interaction Genes
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(Biobase)

# ------------------------------------------------------------
# Convert ExpressionSet → long dataframe
# ------------------------------------------------------------
eset_to_long <- function(eset, genes) {
  
  # ---------------------------
  # 1. Expression → long format
  # ---------------------------
  expr <- as.data.frame(exprs(eset))
  expr$Gene <- rownames(expr)
  expr <- expr[expr$Gene %in% genes, ]
  
  long_expr <- expr %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    )
  
  # ---------------------------
  # 2. Phenotype table
  # ---------------------------
  pheno <- pData(eset)
  pheno$Sample <- rownames(pheno)
  
  # Clean sample IDs for merge
  long_expr$Sample <- long_expr$Sample %>% trimws() %>% toupper() %>% gsub("\\.CEL$", "", .)
  pheno$Sample     <- pheno$Sample %>% trimws() %>% toupper() %>% gsub("\\.CEL$", "", .)
  
  # ---------------------------
  # 3. JOIN + PROPER FIELD EXTRACTION
  # ---------------------------
  long_df <- long_expr %>%
    left_join(pheno, by = "Sample") %>%
    mutate(
      # time:ch1 (string)
      Time = ifelse(
        grepl("pre", characteristics_ch1.2, ignore.case = TRUE),
        "Pre", "Post"
      ),
      # extract numbers from "age: 25"
      Age = as.numeric(gsub(".*?(\\d+).*", "\\1", characteristics_ch1.1)),
      # extract subject ID (letters+numbers)
      Subject = gsub("patientid:?\\s*", "", characteristics_ch1)
    ) %>%
    dplyr::select(Gene, Sample, Expression, Time, Age, Subject)
  
  return(long_df)
}


# ------------------------------------------------------------
# FACET VIOLIN — Pre vs Post
# ------------------------------------------------------------
facet_violin_interaction <- function(eset_clean, genes) {
  
  df <- eset_to_long(eset_clean, genes)
  
  ggplot(df, aes(x = Time, y = Expression, fill = Time)) +
    geom_line(aes(group = Subject), color = "gray60", alpha = 0.4) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
    geom_violin(trim = FALSE, alpha = 0.35) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~ Gene, scales = "free_y") +
    theme_bw(base_size = 12) +
    labs(
      title = "Gene Expression: Pre vs Post",
      x = "Timepoint", y = "log2(Expression)"
    )
}

# ------------------------------------------------------------
# FACET SCATTER — Age vs Change
# ------------------------------------------------------------
facet_scatter_interaction <- function(eset_clean, genes) {
  
  df <- eset_to_long(eset_clean, genes)
  
 
  df_wide <- df %>%
    dplyr::select(Gene, Subject, Age, Time, Expression) %>%
    pivot_wider(names_from = Time, values_from = Expression) %>%
    mutate(Change = Post - Pre)
  
  ggplot(df_wide, aes(x = Age, y = Change)) +
    geom_point(alpha = 0.7, size = 2.1, color = "seagreen4") +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +    # <–– wider layout
    theme_bw(base_size = 12) +
    labs(
      title = "Age Interaction: Post–Pre Change",
      x = "Age (years)", y = "Change (Post - Pre)"
    )
}

