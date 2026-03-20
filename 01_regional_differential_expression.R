# ============================================================
# 01_regional_differential_expression.R
#
# Purpose:
# Regional limma-voom differential expression analysis for
# DLPFC or caudate.
#
# Model:
# ~ 0 + Dx + Neuronal + RIN + Age + Sex + Ethnicity +
#     Perc.Intronic + Perc.rRNA
#
# Outputs:
# - library size table
# - voom matrix
# - limma result tables for:
#   * OCD_vs_Control
#   * ED_vs_Control
#   * OCD_vs_ED
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
})

source("utils/helpers.R")

# -----------------------------
# User input
# -----------------------------
region <- "DLPFC"  # set to "DLPFC" or "Caudate"

counts_file <- if (region == "DLPFC") {
  "data/raw/OCD_raw_DLPFC_noOut.txt"
} else {
  "data/raw/OCD_raw_caudate_noOut.txt"
}

meta_file <- if (region == "DLPFC") {
  "data/metadata/OCD_meta_DLPFC_noOut.txt"
} else {
  "data/metadata/OCD_meta_caudate_noOut.txt"
}

out_dir <- file.path("results", tolower(region))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
counts <- read_counts_matrix(counts_file)
meta <- read_metadata(meta_file)
meta <- prepare_metadata(meta)

stopifnot(ncol(counts) == nrow(meta))

# -----------------------------
# Expression filtering + voom
# -----------------------------
obj <- run_voom_filter(counts)
dge <- obj$dge
v <- obj$voom

write.table(
  dge$samples,
  file = file.path(out_dir, paste0(region, "_library_sizes.tsv")),
  sep = "\t", quote = FALSE, row.names = TRUE
)

write.table(
  v$E,
  file = file.path(out_dir, paste0(region, "_voom_matrix.tsv")),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

# -----------------------------
# Design and model fit
# -----------------------------
design <- build_design_matrix(meta)
contrast_matrix <- get_default_contrasts(design)

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

# -----------------------------
# Extract results
# -----------------------------
res_ocd <- standardize_limma_table(
  topTable(fit, coef = "OCD_vs_Control", number = Inf, sort.by = "P")
)
res_ed <- standardize_limma_table(
  topTable(fit, coef = "ED_vs_Control", number = Inf, sort.by = "P")
)
res_cmp <- standardize_limma_table(
  topTable(fit, coef = "OCD_vs_ED", number = Inf, sort.by = "P")
)

# -----------------------------
# Save results
# -----------------------------
write.table(
  res_ocd,
  file = file.path(out_dir, paste0(region, "_OCD_vs_Control.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_ed,
  file = file.path(out_dir, paste0(region, "_ED_vs_Control.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_cmp,
  file = file.path(out_dir, paste0(region, "_OCD_vs_ED.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Regional differential expression completed for ", region)