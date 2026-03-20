# ============================================================
# 03_variance_partition.R
#
# Purpose:
# Quantify the relative contribution of diagnosis and major
# biological/technical covariates to transcriptome-wide
# expression variation using variancePartition.
#
# Fixed effects:
# - RIN
# - Age
# - Neuronal
# - Perc.Intronic
# - Perc.rRNA
#
# Random effects:
# - Ethnicity
# - Dx
# - Sex
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(doParallel)
  library(reshape2)
  library(ggplot2)
})

source("utils/helpers.R")

# -----------------------------
# User input
# -----------------------------
region <- "DLPFC"  # "DLPFC" or "Caudate"

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

# -----------------------------
# Expression filtering + voom
# -----------------------------
obj <- run_voom_filter(counts)
v <- obj$voom

# -----------------------------
# Variance partition model
# -----------------------------
form <- ~ RIN + Age + Neuronal + Perc.Intronic + Perc.rRNA +
  (1|Ethnicity) + (1|Dx) + (1|Sex)

cl <- makeCluster(4)
registerDoParallel(cl)

vp <- fitExtractVarPartModel(v$E, form, meta)
vp <- sortCols(vp)

# -----------------------------
# Save standard plot
# -----------------------------
png(
  filename = file.path(out_dir, paste0(region, "_variancePartition.png")),
  width = 1800, height = 1200, res = 200
)
plotVarPart(vp, label.angle = 45)
dev.off()

# -----------------------------
# Save violin/box plot
# -----------------------------
vp_pct <- vp * 100
vp_long <- melt(vp_pct)

meds <- round(apply(vp_pct, 2, median), 2)
meds[meds == 0] <- "<0.1"
axis_labels <- paste0(colnames(vp_pct), "\n(", meds, "%)")

p <- ggplot(vp_long, aes(x = variable, y = value)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() +
  ylim(0, 100) +
  scale_x_discrete(labels = axis_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Variance explained (%)")

ggsave(
  filename = file.path(out_dir, paste0(region, "_variancePartition_violin.png")),
  plot = p,
  width = 9, height = 6, dpi = 300
)

write.table(
  vp_pct,
  file = file.path(out_dir, paste0(region, "_variancePartition_percent.tsv")),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
)

stopCluster(cl)

message("Variance partition completed for ", region)