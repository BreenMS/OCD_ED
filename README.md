# Convergent molecular signatures across eating disorders and obsessive-compulsive disorder in the human brain

**Correspondence to: michael[dot] breen [at] mssm [dot] edu<br /> 
A-to-I Editing and Rare Disorders Lab: https://labs.icahn.mssm.edu/breenlab/

Eating disorders (ED) and obsessive-compulsive disorder (OCD) exhibit substantial clinical and genetic overlap, yet whether they share convergent molecular pathology within the human brain remains unclear. We performed large-scale transcriptomic profiling of the dorsolateral prefrontal cortex (DLPFC) and caudate in postmortem tissue from 86 controls, 57 individuals with ED, and 27 individuals with OCD. ED exhibited robust, region-specific transcriptional dysregulation, with 102 differentially expressed genes (DEGs) in the DLPFC and 222 DEGs in the caudate (FDR < 1%), consistently replicated in an independent cohort. In contrast, no single-cohort DEGs reached significance in OCD; however, meta-analysis across three independent datasets identified 57 OCD-associated DEGs in the caudate. Despite these differences in statistical power, transcriptome-wide effect sizes were strongly correlated between ED and OCD in both regions (DLPFC r = 0.67; caudate r = 0.75), indicating marked molecular convergence. Joint analysis of ED and OCD further amplified shared signal, identifying 233 DEGs in the DLPFC and 816 DEGs in the caudate, implicating coordinated disruption of GABAergic signaling, neuroendocrine regulation, metabolic pathways, and CHD8-associated networks. Integration with genetically regulated expression analyses identified five genes associated with ED risk and eight with OCD risk, including cross-disorder signals in WDR6, NCKIPSD, P4HTM, DALRD3, and SHISA5 that converge on common neuronal and metabolic processes. Together, these findings define a shared molecular architecture linking ED and OCD within cortico-striatal circuits, providing mechanistic insight into their comorbidity and a foundation for transdiagnostic target discovery.

![Abstract (1)](https://github.com/BreenMS/Trisomy21/blob/main/Figure_1.png)<br /> 

# WHERE IS THE RAW DATA?<br /> 

All original RNA-sequencing data are publicly at the National Center for Biotechnology Information Gene Expression Omnibus under the following accession number: [GSE262138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262138).  <br /> 
 <br /> <br /> 


# WHERE IS THE CODE? (it's all below)<br /> 
This repository contains the core transcriptomic analyses supporting the manuscript, including regional differential expression, paired differential expression, variance partitioning, cross-study meta-analysis, and Bayes factor evaluation of direct ED vs OCD contrasts.

## Core differential expression model

For all primary regional differential expression analyses, we used limma-voom linear modeling with adjustment for:

- estimated neuronal cell type proportions (`Excitatory + Inhibitory = Neuronal`)
- RNA integrity number (`RIN`)
- age
- sex
- ethnicity
- percentage of intronic reads (`Perc.Intronic`)
- percentage of rRNA reads (`Perc.rRNA`)

## Paired analyses

For paired analyses, we used the same fixed-effect covariates and modeled within-subject correlation using `duplicateCorrelation` in limma with `SubjectID` as the blocking variable.

## Repository structure

- `analysis/01_regional_differential_expression.R`: regional limma-voom analyses for DLPFC and caudate
- `analysis/02_paired_differential_expression.R`: paired limma-voom analyses using `duplicateCorrelation`
- `analysis/03_variance_partition.R`: transcriptome-wide variance decomposition
- `analysis/04_meta_analysis.R`: gene-wise random-effects meta-analysis using REML
- `analysis/05_bayes_factor_ED_vs_OCD.R`: Wakefield approximate Bayes factor analysis for direct ED vs OCD contrasts
- `utils/helpers.R`: shared helper functions

## Notes

Metadata are expected to contain the following fields where relevant:

- `Dx`
- `Sex`
- `Age`
- `RIN`
- `Ethnicity`
- `Excitatory`
- `Inhibitory`
- `Perc.Intronic`
- `Perc.rRNA`
- `SubjectID` (paired analyses only)

