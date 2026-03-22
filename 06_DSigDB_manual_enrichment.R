# ============================================================
# Manual DSigDB enrichment (no clusterProfiler)
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# -----------------------------
# 1. Inputs
# -----------------------------
deg_file <- "DLPFC_DEGs.txt"
deg_file <- "Caudate_DEGs.txt"
gmt_file <- "DSigDB.txt"

out_dir <- "DSigDB_manual_enrichment"
dir.create(out_dir, showWarnings = FALSE)

fdr_cutoff <- 0.01

# CNS keywords
cns_keywords <- c(
  "ssri","snri","antidepress","fluoxetine","sertraline","paroxetine",
  "citalopram","escitalopram","fluvoxamine","venlafaxine","duloxetine",
  "bupropion","mirtazapine",
  "antipsych","clozapine","haloperidol","risperidone","olanzapine",
  "quetiapine","aripiprazole",
  "lithium","valproate","lamotrigine",
  "benzodiazep","diazepam","lorazepam","alprazolam",
  "amphetamine","methylphenidate","modafinil"
)

# -----------------------------
# 2. Read DEG file
# -----------------------------
deg <- fread(deg_file)

deg$ENSEMBL <- sub("\\..*", "", deg$`Ensembl ID`)

# Map to SYMBOL
mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg$ENSEMBL,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

deg <- deg %>%
  left_join(mapping, by = c("ENSEMBL" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL))

# Keep best per gene
deg <- deg %>%
  arrange(adj.P.Val) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# Define DEG sets
deg_sig <- deg %>% filter(adj.P.Val < fdr_cutoff)

up_genes   <- deg_sig %>% filter(logFC > 0) %>% pull(SYMBOL)
down_genes <- deg_sig %>% filter(logFC < 0) %>% pull(SYMBOL)

background <- deg$SYMBOL

cat("DE genes:", length(deg_sig$SYMBOL), "\n")

# -----------------------------
# 3. Read GMT
# -----------------------------
gmt_lines <- readLines(gmt_file)

gmt_list <- lapply(gmt_lines, function(x) {
  parts <- strsplit(x, "\t")[[1]]
  term  <- parts[1]
  genes <- parts[-c(1,2)]
  list(term = term, genes = genes)
})

# -----------------------------
# 4. Filter CNS terms
# -----------------------------
pattern <- paste(cns_keywords, collapse = "|")

gmt_list_cns <- gmt_list[
  sapply(gmt_list, function(x)
    str_detect(tolower(x$term), pattern))
]

cat("CNS gene sets:", length(gmt_list_cns), "\n")

# -----------------------------
# 5. Enrichment function
# -----------------------------
run_enrichment <- function(gene_list, label) {

  results <- lapply(gmt_list_cns, function(gs) {

    geneset <- unique(gs$genes)

    overlap <- intersect(gene_list, geneset)

    k <- length(overlap)                    # overlap
    M <- length(geneset)                   # geneset size
    N <- length(background)                # universe
    n <- length(gene_list)                 # DEG set size

    if (k == 0) return(NULL)

    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)

    data.frame(
      term = gs$term,
      overlap = k,
      geneset_size = M,
      DEG_size = n,
      p_value = pval,
      overlap_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  res <- bind_rows(results)

  if (nrow(res) == 0) return(NULL)

  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res$contrast <- label

  res %>% arrange(p_adj)
}

# -----------------------------
# 6. Run enrichment
# -----------------------------
res_up   <- run_enrichment(up_genes, "UP")
res_down <- run_enrichment(down_genes, "DOWN")

all_res <- bind_rows(res_up, res_down)

write_tsv(all_res, file.path(out_dir, "DSigDB_CNS_Caudate_enrichment.tsv"))

# -----------------------------
# 7. Quick summary
# -----------------------------
summary <- all_res %>%
  group_by(contrast) %>%
  summarise(
    significant_terms = sum(p_adj < 0.05),
    .groups = "drop"
  )

print(summary)
write_tsv(summary, file.path(out_dir, "summary.tsv"))

cat("\nDone.\n")