# Trisomy 21 Drives ADARB1 Overexpression and Premature RNA Recoding in the Developing Fetal Brain

**Correspondence to: michael[dot] breen [at] mssm [dot] edu<br /> 
A-to-I Editing and Rare Disorders Lab: https://labs.icahn.mssm.edu/breenlab/

Understanding how chromosome 21 gene dosage contributes to neurodevelopmental and systemic phenotypes in trisomy 21 (T21) remains a fundamental challenge. We performed transcriptome-wide RNA sequencing on fetal cortical and hippocampal tissues from 20 T21 cases and 27 euploid controls collected between 13–22 weeks post-conception, a critical period for human brain development. Differential expression analysis revealed 572 dysregulated genes in the prefrontal cortex and 519 in the hippocampus (FDR < 5%), with significant enrichment for chromosome 21 genes. Functional enrichment analyses highlighted disruptions in neurodevelopmental, synaptic, and immune-related pathways. Among the most strongly dysregulated genes was ADARB1, a chromosome 21-encoded RNA editing enzyme, whose overexpression in T21 fetal brain was associated with increased adenosine-to-inosine (A-to-I) editing, including recoding sites in GRIA2 (p.R764G), GRIA3 (p.R775G), and GRIK2 (p.Y571C, p.Q621R). A meta-analysis incorporating nine independent transcriptomic datasets spanning early embryonic and progenitor cell types validated robust chromosome 21 dosage effects, including consistent ADARB1 overexpression. Extending these findings, a meta-analysis of A-to-I editing across datasets revealed widespread over-editing at 3′UTRs and at GRIA3 (p.R775G), a site critical for AMPA receptor desensitization. Together, these results implicate dysregulated RNA editing driven by ADARB1 overexpression as a post-transcriptional mechanism contributing to fetal neuropathology in T21 and provide a framework for understanding the broader molecular consequences of chromosome 21 dosage sensitivity during brain development.

![Abstract (1)](https://github.com/BreenMS/Trisomy21/blob/main/Figure_1.png)<br /> 

# WHERE IS THE DATA? (two locations below)<br /> 

1. All original RNA-sequencing data are publicly at the National Center for Biotechnology Information Gene Expression Omnibus under the following accession number: [GSE301886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301886).  <br /> 
2. Data can also be downloaded from our interactive [RShiny App](https://andyyang.shinyapps.io/Ts21-dashboard/). It also enables gene-by-gene inspection of between-study effect sizes and heterogeneity, generating meta-analytic forest plots for any gene of interest across the ten harmonized Trisomy 21 RNA-sequencing datasets. <br /> <br /> 


# WHERE IS THE CODE? (it's all below)<br /> 
Code falls into 2 core areas:

1. *Gene-level analyses*, including RNA-seq QC, mapping and counting via [a detailed and highly cited NextFlow pipeline](https://github.com/CommonMindConsortium/RAPiD-nf) <br />
Note that normalization [VOOM](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html), differential expression testing [limma](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html), network analysis [WGCNA](https://cran.r-project.org/web/packages/WGCNA/refman/WGCNA.html), gene set preservation analysis, and cell type deconvolution [bMIND](https://randel.github.io/MIND/) are done following these GitHub and publications while following our methods section. 

Additionally, gene-level meta-analyses were done with [provided R code](https://github.com/BreenMS/Trisomy21/blob/main/T21_meta_analysis.r)<br /> 

2. *RNA editing-level analyses*, including computint an Alu Editing Index (AEI) from a STAR mapped bam file [RNAEditingIndexer v1.0](https://github.com/a2iEditing/RNAEditingIndexer) and quantifying RNA editing sites from STAR mapped bam files using known editing sites [as previously described](https://github.com/BreenMS/Quantify-RNA-editing-from-bam-file). Further details provided below...<br /> <br /> 

Compute AEI from a STAR mapped bam file:<br />  <br /> 
We used already available software from the [RNAEditingIndexer GitHub account](https://github.com/a2iEditing/RNAEditingIndexer) to compute an AEI based on a mapped bam file. The method is describe in the original publication, [Nat. Methods (2019)](https://pubmed.ncbi.nlm.nih.gov/31636457/). Here, we provide an example of bash shell script that executes the AEI on one sample. Requirements and parameters are described in full in the bash script.  <br /> 
 
An example for computing AEI on human samples:
```ruby
RNAEditingIndex -d -f Aligned.sortedByCoord.out.bam -o .
--genes_expression ucscHg38GTExGeneExpression.bed.gz
--refseq ucscHg38RefSeqCurated.bed.gz
--snps ucscHg38CommonGenomicSNPs150.bed.gz
-gf ucscHg38Genome.fa
-rb ucscHg38Alu.bed.gz
--genome UserProvided  --paired_end --stranded
```
<br />  
<br /> 

Quantify RNA editing from STAR mapped bam files using a list of predefined list of sites (based on a predefined list of sites):<br />  <br /> 
It is often of interest to quantify RNA editing sites based on a user defined list of sites. Samtools mpileup has the functionality to execute this task. Here we provide two perl scripts that will achieve this task. The only requirement is installing a recent version of samtools. In the current study, we leveraged lists of known sites from these three resources: [REDIportal](https://academic.oup.com/nar/article/49/D1/D1012/5940507), [cellular and genetic drivers of RNA editing variability in the human brain](https://www.nature.com/articles/s41467-022-30531-0), and [an atlas of human recoding sites](https://www.nature.com/articles/s41467-022-28841-4). <br /> 


query_known_sites.pl= excute mpileup (samtools) to query a list of known editing sites.<br />
parse_pileup_query.pl = a requirement for query_known_sites.pl<br />  
Usage: perl query_known_sites.pl [A predefined list of known editing sites] [STAR mapped bam file] [Output file name]
```ruby
perl query_known_sites.pl CNS_A2G_events.txt SampleName.bam OutputFileName.txt
```
<br />  


