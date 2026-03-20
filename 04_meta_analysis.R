# ============================================================
# 04_meta_analysis.R
#
# Purpose:
# Gene-wise random-effects meta-analysis across studies using
# metafor::rma with REML estimation.
#
# Expected input columns in each study file:
# - Gene
# - logFC
# - SE
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(metafor)
})

# -----------------------------
# User input
# -----------------------------
study_files <- c(
  "data/meta/d1.tsv",
  "data/meta/d2.tsv",
  "data/meta/d3.tsv"
)

study_names <- c("Study1", "Study2", "Study3")
out_file <- "results/meta/meta_analysis_results_3studies.tsv"
dir.create("results/meta", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper
# -----------------------------
read_study <- function(file, study_name) {
  dat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(all(c("Gene", "logFC", "SE") %in% colnames(dat)))

  dat %>%
    transmute(
      Gene = Gene,
      logFC = as.numeric(logFC),
      SE = as.numeric(SE),
      Study = study_name
    )
}

# -----------------------------
# Read and combine studies
# -----------------------------
study_list <- purrr::map2(study_files, study_names, read_study)
dat_all <- bind_rows(study_list)

# Keep genes present in all studies
complete_genes <- dat_all %>%
  count(Gene) %>%
  filter(n == length(study_files)) %>%
  pull(Gene)

dat_all <- dat_all %>%
  filter(Gene %in% complete_genes)

# -----------------------------
# Run meta-analysis per gene
# -----------------------------
meta_results <- dat_all %>%
  group_by(Gene) %>%
  group_modify(~{
    fit <- rma(yi = .x$logFC, sei = .x$SE, method = "REML")

    tibble(
      Effect = as.numeric(fit$beta),
      SE = fit$se,
      Zval = fit$zval,
      Pval = fit$pval,
      CI.LB = fit$ci.lb,
      CI.UB = fit$ci.ub,
      Tau2 = fit$tau2,
      QEp = fit$QEp
    )
  }) %>%
  ungroup() %>%
  mutate(
    AdjustedP_BH = p.adjust(Pval, method = "BH"),
    AdjustedP_Hochberg = p.adjust(Pval, method = "hochberg")
  ) %>%
  arrange(Pval)

write.table(
  meta_results,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Meta-analysis completed")