# ============================================================
# helpers.R
# Shared helper functions for transcriptomic analyses
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

read_counts_matrix <- function(file) {
  x <- read.delim(
    file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  as.matrix(x)
}

read_metadata <- function(file) {
  read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

standardize_limma_table <- function(tab, gene_col = NULL) {
  tab <- as.data.frame(tab)

  genes <- if (!is.null(gene_col) && gene_col %in% colnames(tab)) {
    tab[[gene_col]]
  } else {
    rownames(tab)
  }

  tibble(
    Gene = genes,
    logFC = if ("logFC" %in% colnames(tab)) as.numeric(tab$logFC) else NA_real_,
    AveExpr = if ("AveExpr" %in% colnames(tab)) as.numeric(tab$AveExpr) else NA_real_,
    t = if ("t" %in% colnames(tab)) as.numeric(tab$t) else NA_real_,
    P.Value = if ("P.Value" %in% colnames(tab)) as.numeric(tab$P.Value) else NA_real_,
    adj.P.Val = if ("adj.P.Val" %in% colnames(tab)) as.numeric(tab$adj.P.Val) else NA_real_,
    B = if ("B" %in% colnames(tab)) as.numeric(tab$B) else NA_real_
  )
}

prepare_metadata <- function(meta) {
  stopifnot(all(c("Dx", "Sex", "Age", "RIN", "Ethnicity", "Excitatory",
                  "Inhibitory", "Perc.Intronic", "Perc.rRNA") %in% colnames(meta)))

  meta$Dx <- factor(meta$Dx)
  meta$Sex <- factor(meta$Sex)
  meta$Ethnicity <- factor(meta$Ethnicity)
  meta$Neuronal <- meta$Excitatory + meta$Inhibitory
  meta
}

build_design_matrix <- function(meta) {
  design <- model.matrix(
    ~ 0 + Dx + Neuronal + RIN + Age + Sex + Ethnicity + Perc.Intronic + Perc.rRNA,
    data = meta
  )
  colnames(design) <- make.names(colnames(design))
  design
}

run_voom_filter <- function(counts) {
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(limma))

  dge <- DGEList(counts = counts)
  keep <- rowSums(cpm(dge) > 1) >= ceiling(ncol(counts) / 3)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design = NULL, plot = FALSE)

  list(dge = dge, voom = v, keep = keep)
}

get_default_contrasts <- function(design) {
  limma::makeContrasts(
    OCD_vs_Control = Dxcontrol - DxOCD,
    ED_vs_Control  = Dxcontrol - DxED,
    OCD_vs_ED      = DxOCD - DxED,
    levels = design
  )
}