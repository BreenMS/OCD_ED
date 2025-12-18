# =============================================================================
# TITLE: Meta-Analysis of T21 Differential Expressed Gene (DEG) Summary Statistics
# AUTHOR: MSB
# DATE: 2025-07-01
#
# DESCRIPTION: 
#   This script performs a meta-analysis on DEG summary statistics 
#   from multiple datasets using the metafor package. It extracts 
#   relevant columns, filters incomplete data, runs REML-based 
#   meta-analysis per gene, and outputs a summary table with 
#   adjusted p-values (BH) and significance categories.
#
# INPUT: 
#   - 17 DEG summary statistic `.txt` files in tab-delimited format
# OUTPUT:
#   - summary_table: tibble of pooled logFC estimates with p-values
#   - failed_genes.txt: list of genes that failed meta-analysis
#
# REQUIREMENTS: (two R packages) tidyverse, metafor
# =============================================================================

# Load required packages
library(tidyverse)
library(metafor)

# ─────────────────────────────────────────────────────────────────────────────
# META-ANALYSIS FOR T21 DEG
# ─────────────────────────────────────────────────────────────────────────────

# 1. DEFINE FILE LIST AND COLUMNS TO EXTRACT

DEG_summary_files <- c(
  "GSE101942_Cortical_Neurons.txt", "GSE101942_IPSC.txt", "GSE128621_WBCs.txt",
  "GSE144857_hIPSC.txt", "GSE144857_GLU_Neurons.txt", "GSE185192_IPSC.txt",
  "GSE185192_NPC.txt", "GSE52249_IPSC.txt", "GSE55504_Fetal_Primary_Fibroblast.txt",
  "GSE55504_Primary_Fibrobast.txt", "GSE79842_Fibroblast.txt", "GSE84531_Monocyte.txt",
  "GSE84531_Tcell.txt", "PRJNA751601_NPC.txt", "PRJNA751601_Ventralized_NPC.txt",
  "DEG_DLPFC_final.txt", "DEG_HIPP_final.txt"
)

keep_columns <- c("gene.id", "gene.symbol", "Chr", "logFC", "t")

# 2. LOAD AND EXTRACT RELEVANT COLUMNS FROM EACH FILE
load_and_filter <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t")
  if (all(keep_columns %in% colnames(df))) {
    df[, keep_columns, drop = FALSE]
  } else {
    missing <- setdiff(keep_columns, colnames(df))
    warning(sprintf("File '%s' is missing columns: %s", file, paste(missing, collapse = ", ")))
    return(NULL)
  }
}

data_list <- lapply(DEG_summary_files, load_and_filter)
data_list <- Filter(Negate(is.null), data_list)  

# 3. COMBINE DATASETS AND CALCULATE ADDITIONAL STATS
combined_data <- bind_rows(data_list) %>%
  mutate(
    gene = gene.id,
    SE = abs(logFC / t)
  )


# 4. RUN META-ANALYSIS PER GENE USING metafor::rma
failed_genes <- character(0)

run_meta <- function(gene_symbol) {
  dat <- combined_data %>% filter(gene.symbol == gene_symbol)
  tryCatch({
    res <- rma(yi = logFC, sei = SE, data = dat, method = "REML")
    res$gene_id <- gene_symbol
    res$Chr <- unique(dat$Chr)
    return(res)
  }, error = function(e) {
    failed_genes <<- c(failed_genes, gene_symbol)
    return(NULL)
  })
}

results <- lapply(unique(combined_data$gene.symbol), run_meta)
results <- Filter(Negate(is.null), results)

if (length(failed_genes) > 0) {
  message(length(failed_genes), " genes failed meta-analysis.")
  writeLines(failed_genes, "failed_genes.txt")  
}

# 5. BUILD SUMMARY STATISTICS TABLE FROM META-ANLAYSIS RESULTS
summary_table <- lapply(results, function(res) {
  tibble(
    gene = res$gene_id,
    Chr = res$Chr,
    pooled_logFC = res$beta[1],
    se = res$se,
    zvalue = res$zval,
    pvalue = res$pval,
    ci.lb = res$ci.lb,
    ci.ub = res$ci.ub,
    I2 = res$I2
  )
}) %>%
  bind_rows() %>%
  mutate(
    adj_pvalue = p.adjust(pvalue, method = "BH"),
    direction = ifelse(pooled_logFC > 0, "T21", "Control"),
    Significance = ifelse(adj_pvalue < 0.05, "Significant", "Nonsignificant"),
    color = case_when(
      Significance == "Significant" & direction == "T21" ~ "T21",
      Significance == "Significant" & direction == "Control" ~ "Control",
      TRUE ~ "Nonsignificant"
    )
  ) %>%
  arrange(adj_pvalue)