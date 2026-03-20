# ============================================================
# 02_paired_differential_expression.R
#
# Purpose:
# Paired limma-voom differential expression analysis using
# duplicateCorrelation to model within-subject dependence.
#
# Required metadata fields:
# - SubjectID
# - Dx, Sex, Age, RIN, Ethnicity
# - Excitatory, Inhibitory
# - Perc.Intronic, Perc.rRNA
#
# Model:
# ~ 0 + Dx + Neuronal + RIN + Age + Sex + Ethnicity +
#     Perc.Intronic + Perc.rRNA
# with block = SubjectID
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
counts_file <- "data/raw/OCD_paired_counts.txt"
meta_file   <- "data/metadata/OCD_paired_metadata.txt"
out_dir     <- "results/paired"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
counts <- read_counts_matrix(counts_file)
meta <- read_metadata(meta_file)
meta <- prepare_metadata(meta)

stopifnot("SubjectID" %in% colnames(meta))
meta$SubjectID <- factor(meta$SubjectID)

stopifnot(ncol(counts) == nrow(meta))

# -----------------------------
# Expression filtering
# -----------------------------
dge <- DGEList(counts = counts)
keep <- rowSums(cpm(dge) > 1) >= ceiling(ncol(counts) / 3)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# -----------------------------
# Design
# -----------------------------
design <- build_design_matrix(meta)
contrast_matrix <- get_default_contrasts(design)

# -----------------------------
# Voom + duplicateCorrelation
# -----------------------------
v <- voom(dge, design = design, plot = TRUE)

corfit <- duplicateCorrelation(
  v,
  design = design,
  block = meta$SubjectID
)

fit <- lmFit(
  v,
  design = design,
  block = meta$SubjectID,
  correlation = corfit$consensus
)

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
# Save outputs
# -----------------------------
write.table(
  res_ocd,
  file = file.path(out_dir, "paired_OCD_vs_Control.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_ed,
  file = file.path(out_dir, "paired_ED_vs_Control.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_cmp,
  file = file.path(out_dir, "paired_OCD_vs_ED.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(
  paste("Consensus correlation:", round(corfit$consensus, 4)),
  con = file.path(out_dir, "duplicateCorrelation_summary.txt")
)

message("Paired differential expression completed")