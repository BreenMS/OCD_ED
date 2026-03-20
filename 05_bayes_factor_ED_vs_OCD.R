# ============================================================
# 05_bayes_factor_ED_vs_OCD.R
#
# Purpose:
# Compute Wakefield approximate Bayes factors for the direct
# ED vs OCD comparison using limma summary statistics.
#
# Expected input:
# limma result tables with columns:
# - Gene
# - logFC
# - t   and/or SE
# - P.Value
# - adj.P.Val
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# -----------------------------
# User input
# -----------------------------
dlpfc_file   <- "results/dlpfc/DLPFC_OCD_vs_ED.tsv"
caudate_file <- "results/caudate/Caudate_OCD_vs_ED.tsv"

out_dir <- "results/bayes"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prior_sd <- 0.2
prior_var <- prior_sd^2

# -----------------------------
# Helpers
# -----------------------------
compute_se <- function(df) {
  if ("SE" %in% colnames(df)) {
    return(as.numeric(df$SE))
  }

  if (all(c("logFC", "t") %in% colnames(df))) {
    return(abs(as.numeric(df$logFC) / as.numeric(df$t)))
  }

  stop("Need either SE or both logFC and t columns.")
}

compute_wakefield <- function(beta_hat, V, W) {
  z2 <- (beta_hat^2) / V
  bf10 <- sqrt(V / (V + W)) * exp((z2 / 2) * (W / (V + W)))
  bf01 <- 1 / bf10
  tibble(bf10 = bf10, bf01 = bf01)
}

summarize_bf <- function(df, region_label) {
  tibble(
    region = region_label,
    n_genes = nrow(df),
    median_abs_logFC = median(abs(df$logFC), na.rm = TRUE),
    median_SE = median(df$SE, na.rm = TRUE),
    median_bf01 = median(df$bf01, na.rm = TRUE),
    median_bf10 = median(df$bf10, na.rm = TRUE),
    prop_bf01_gt1 = mean(df$bf01 > 1, na.rm = TRUE),
    prop_bf01_gt3 = mean(df$bf01 > 3, na.rm = TRUE),
    prop_bf01_gt10 = mean(df$bf01 > 10, na.rm = TRUE),
    prop_bf10_gt1 = mean(df$bf10 > 1, na.rm = TRUE),
    prop_bf10_gt3 = mean(df$bf10 > 3, na.rm = TRUE),
    prop_bf10_gt10 = mean(df$bf10 > 10, na.rm = TRUE)
  )
}

analyze_region <- function(file, region_label, prior_var, out_dir) {
  df <- fread(file, data.table = FALSE)

  stopifnot(all(c("Gene", "logFC") %in% colnames(df)))

  res <- tibble(
    gene = df$Gene,
    logFC = as.numeric(df$logFC),
    t = if ("t" %in% colnames(df)) as.numeric(df$t) else NA_real_,
    P.Value = if ("P.Value" %in% colnames(df)) as.numeric(df$P.Value) else NA_real_,
    adj.P.Val = if ("adj.P.Val" %in% colnames(df)) as.numeric(df$adj.P.Val) else NA_real_
  )

  res$SE <- compute_se(df)
  res$V <- res$SE^2

  res <- res %>%
    filter(!is.na(gene), !is.na(logFC), !is.na(SE), is.finite(SE), SE > 0)

  bf <- compute_wakefield(res$logFC, res$V, prior_var)

  res <- bind_cols(res, bf) %>%
    mutate(
      region = region_label,
      null_support = case_when(
        bf01 > 10 ~ "Strong null support",
        bf01 > 3  ~ "Moderate null support",
        bf01 > 1  ~ "Anecdotal null support",
        bf10 > 10 ~ "Strong alternative support",
        bf10 > 3  ~ "Moderate alternative support",
        bf10 > 1  ~ "Anecdotal alternative support",
        TRUE ~ "Equivocal"
      )
    )

  write_tsv(
    res,
    file.path(out_dir, paste0(region_label, "_BayesFactor_results.tsv"))
  )

  write_tsv(
    summarize_bf(res, region_label),
    file.path(out_dir, paste0(region_label, "_BayesFactor_summary.tsv"))
  )

  p1 <- ggplot(res, aes(x = bf01)) +
    geom_histogram(bins = 80) +
    geom_vline(xintercept = c(1, 3, 10), linetype = c("solid", "dashed", "dashed")) +
    scale_x_log10() +
    theme_bw() +
    labs(
      title = paste0(region_label, ": BF01 distribution"),
      x = "BF01 (log10 scale)",
      y = "Number of genes"
    )

  ggsave(
    file.path(out_dir, paste0(region_label, "_BF01_histogram.png")),
    p1, width = 7, height = 5, dpi = 300
  )

  p2 <- ggplot(res, aes(x = abs(logFC), y = bf01)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_hline(yintercept = c(1, 3, 10), linetype = c("solid", "dashed", "dashed")) +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = paste0(region_label, ": |logFC| vs BF01"),
      x = "|logFC|",
      y = "BF01 (log10 scale)"
    )

  ggsave(
    file.path(out_dir, paste0(region_label, "_abslogFC_vs_BF01.png")),
    p2, width = 7, height = 5, dpi = 300
  )

  res
}

# -----------------------------
# Run analyses
# -----------------------------
res_dlpfc <- analyze_region(dlpfc_file, "DLPFC", prior_var, out_dir)
res_caudate <- analyze_region(caudate_file, "Caudate", prior_var, out_dir)

summary_all <- bind_rows(
  summarize_bf(res_dlpfc, "DLPFC"),
  summarize_bf(res_caudate, "Caudate")
)

write_tsv(
  summary_all,
  file.path(out_dir, "BayesFactor_summary_all_regions.tsv")
)

message("Bayes factor analysis completed")