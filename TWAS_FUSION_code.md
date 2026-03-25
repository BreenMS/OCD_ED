# TWAS Analysis Workflow (FUSION)

This repository documents a workflow for performing Transcriptome-Wide Association Studies (TWAS) using the FUSION framework, including custom weight generation from RNA-seq data and downstream association testing.

---

## References

- http://gusevlab.org/projects/fusion/
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4767558/

---

## 1. Installation

### FUSION
```bash
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
cd fusion_twas-master
```

### LD reference (1000 Genomes)
```bash
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2
```

### R packages
```r
install.packages(c("optparse", "RColorBrewer", "glmnet", "methods"))
library(devtools)
install_github("gabraham/plink2R/plink2R", repos = NULL)
```

### Additional requirements
- PLINK2 (in PATH)
- GCTA (robust version)
- GEMMA

---

## 2. Input GWAS Data

Required columns:
- SNP (rsID)
- A1 (effect allele)
- A2 (non-effect allele)
- Z (Z-score aligned to A1)

### Example: compute Z-scores
```r
library(data.table)

AN <- fread("AN.QC.Transformed")
AN$Z <- AN$BETA / AN$SE

AN.twas <- AN[, c(3,5,4,16,15)]
write.table(AN.twas, "AN.twas",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
```

---

## 3. Build Gene-Level PLINK Files

### Overview
- Normalize RNA-seq expression
- Adjust for covariates and hidden confounders (qSVs)
- Map SNPs within ±500kb of genes
- Generate gene-level PLINK files

### Example (core steps)
```r
load("rse_gene_modified.Rdata")
rse <- rse_gene

# Normalize expression
geneRpkm <- recount::getRPKM(rse_gene, length = "Length")
assays(rse)$raw_expr <- geneRpkm

# Model covariates
mod <- model.matrix(~Dx + AgeDeath + mitoRate + rRNA_rate +
                    totalAssignedGene + RIN + overallMapRate,
                    data = colData(rse))

# Clean expression
assays(rse)$clean_expr <- cleaningY(
    log2(assays(rse)$raw_expr + 1),
    mod,
    P = 2
)
```

---

## 4. Compute TWAS Weights

### Single gene
```bash
Rscript FUSION.compute_weights.R \
  --bfile path/to/gene \
  --tmp tmp_files/gene_X \
  --out out_files/gene_X \
  --PATH_plink /path/to/plink \
  --PATH_gcta /path/to/gcta \
  --PATH_gemma /path/to/gemma \
  --models top1,lasso \
  --hsq_p 1.0001 \
  --verbose 1 \
  --save_hsq
```

### Parallel execution
```bash
parallel -j 10 Rscript FUSION.compute_weights.R \
  --bfile /path/to/gene_{} \
  --tmp tmp_files/gene_{} \
  --out out_files/gene_{} \
  --PATH_plink /path/to/plink \
  --PATH_gcta /path/to/gcta \
  --PATH_gemma /path/to/gemma \
  --models top1,lasso \
  --hsq_p 1.0001 \
  --verbose 1 \
  --save_hsq ::: {1..17656}
```

---

## 5. Merge TWAS Weights

```r
rdat_files <- dir("out_files", ".wgt.RDat", full.names = TRUE)

write.table(rdat_files, "gene.list",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

system("Rscript fusion_twas-master/utils/FUSION.profile_wgt.R gene.list")
```

---

## 6. TWAS Association

```bash
Rscript FUSION.assoc_test.R \
  --sumstats AN.twas \
  --weights gene.pos \
  --weights_dir ./ \
  --ref_ld_chr LDREF/1000G.EUR. \
  --chr 22 \
  --out chr22.dat
```

### Parallel across chromosomes
```bash
parallel -j 10 Rscript FUSION.assoc_test.R \
  --sumstats AN.twas \
  --weights gene.pos \
  --weights_dir ./ \
  --ref_ld_chr LDREF/1000G.EUR. \
  --chr {} \
  --out chr{}.dat ::: {1..22}
```

---

## 7. Combine Results

```r
files <- list.files("AN_gwas", pattern = ".dat$")
dats <- lapply(files, function(x) {
  read.table(file.path("AN_gwas", x), header = TRUE)
})

twas <- do.call(rbind, dats)
save(twas, file = "twas_results.RData")
```

---

## 8. Visualization

```r
ggplot(twas, aes(x = BPcum, y = TWAS.Z)) +
  geom_point(aes(color = as.factor(CHR))) +
  theme_bw()
```

### Outputs
- Manhattan plot (PDF)
- Scatter plots
- TWAS results tables

---

## 9. Functional Enrichment

```r
library(clusterProfiler)
library(org.Hs.eg.db)

twas_sig <- twas[twas$fdr.p < 0.05,]

go <- enrichGO(
  gene = twas_sig$EntrezID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr"
)
```

---

## Notes

- SNP IDs must match LD reference (`.bim`)
- Only overlapping SNPs are used
- Expression is covariate-corrected
- Gene window: ±500kb
- Models used: top1, lasso

---

## Output Summary

| Step | Output |
|------|--------|
| Weights | `.wgt.RDat` files |
| Weight index | `gene.pos` |
| TWAS results | `chr*.dat` |
| Combined results | `.RData` |
| Plots | PDF |
| Enrichment | CSV |
