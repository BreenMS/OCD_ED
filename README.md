# Trisomy 21 Drives ADARB1 Overexpression and Premature RNA Recoding in the Developing Fetal Brain

**Correspondence to: michael[dot] breen [at] mssm [dot] edu<br /> 
A-to-I Editing and Rare Disorders Lab: https://labs.icahn.mssm.edu/breenlab/

Understanding how chromosome 21 gene dosage contributes to neurodevelopmental and systemic phenotypes in trisomy 21 (T21) remains a fundamental challenge. We performed transcriptome-wide RNA sequencing on fetal cortical and hippocampal tissues from 20 T21 cases and 27 euploid controls collected between 13–22 weeks post-conception, a critical period for human brain development. Differential expression analysis revealed 572 dysregulated genes in the prefrontal cortex and 519 in the hippocampus (FDR < 5%), with significant enrichment for chromosome 21 genes. Functional enrichment analyses highlighted disruptions in neurodevelopmental, synaptic, and immune-related pathways. Among the most strongly dysregulated genes was ADARB1, a chromosome 21-encoded RNA editing enzyme, whose overexpression in T21 fetal brain was associated with increased adenosine-to-inosine (A-to-I) editing, including recoding sites in GRIA2 (p.R764G), GRIA3 (p.R775G), and GRIK2 (p.Y571C, p.Q621R). A meta-analysis incorporating nine independent transcriptomic datasets spanning early embryonic and progenitor cell types validated robust chromosome 21 dosage effects, including consistent ADARB1 overexpression. Extending these findings, a meta-analysis of A-to-I editing across datasets revealed widespread over-editing at 3′UTRs and at GRIA3 (p.R775G), a site critical for AMPA receptor desensitization. Together, these results implicate dysregulated RNA editing driven by ADARB1 overexpression as a post-transcriptional mechanism contributing to fetal neuropathology in T21 and provide a framework for understanding the broader molecular consequences of chromosome 21 dosage sensitivity during brain development.

![Abstract (1)](https://github.com/BreenMS/Living-Brain/blob/main/Figure_1.png)

Here we described the main computational code used to generate all results and figures in this body of work.  

All original RNA-sequencing data are publicly at the National Center for Biotechnology Information Gene Expression Omnibus under the following accession number: [GSE301886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301886).  

# This work entails four main levels of analysis:
1. Compute an Alu Editing Index (AEI) from a STAR mapped bam file  [(RNAEditingIndexer v1.0)](https://github.com/a2iEditing/RNAEditingIndexer)<br /> 
2. Quantifying RNA editing sites from STAR mapped bam files using de novo methods [(reditools v2.0)](https://github.com/tizianoflati/reditools2.0) and [(JACUSA2)](https://github.com/dieterich-lab/JACUSA2)<br /> 
3. Quantifying RNA editing from STAR mapped bam files using a list of predefined list of sites (code provided below)<br /> 

An overview of our analytical approach:<br /> 
<img width="500" alt="Screen Shot 2022-07-07 at 10 16 40 AM" src="https://github.com/BreenMS/Living-Brain/blob/main/WorkFlow.png">

# 1. Compute AEI from a STAR mapped bam file:
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


# 2. Quantify RNA editing sites from STAR mapped bam files using de novo methods:
We used already available software from the [reditools v2.0 GitHub account](https://github.com/tizianoflati/reditools2.0) and [JACUSA2 GitHub account](https://github.com/dieterich-lab/JACUSA2) to quantify de novo RNA editing sites based on a STAR mapped bam file. The methods are describe in the original publications: [BMC Bioinformatics (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03562-x) and [Genome Biology (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0). Here, we provide examples of how we executed reditools 2.0 and JACUSA2 on the BrainVar data set. <br /> 

REDITOOLS 2.0
```ruby
mpirun parallel_reditools.py -f Aligned.sortedByCoord.out.bam -r GRCh38.chrom.fa -S -s 2 -ss 5 -mrl 50 -q 10 -bq 20 -C -T 2 -m homopolymeric_sites_hg38.txt -os 5 -Z GRCh38.chrom.fa.fai -G Aligned.sortedByCoord.out.bam.cov -D Aligned.sortedByCoord.out.bam_out
```
<br />  

JACUSA2
```ruby
java -jar $JACUSA2_JAR call-1 -r Aligned.sortedByCoord.out.bam -p 10 -a D,M,Y,E:file=hg38-blacklist.v2_sort.bed:type=BED -s -m 20 -R GRCh38.chrom.fa -P RF-FIRSTSTRAND
```
<br />  



# 3. Quantify RNA editing from STAR mapped bam files using a list of predefined list of sites (based on a predefined list of sites):
It is often of interest to quantify RNA editing sites based on a user defined list of sites. Samtools mpileup has the functionality to execute this task. Here we provide two perl scripts that will achieve this task. The only requirement is installing a recent version of samtools. In the current study, we leveraged lists of known sites from these three resources: [REDIportal](https://academic.oup.com/nar/article/49/D1/D1012/5940507), [cellular and genetic drivers of RNA editing variability in the human brain](https://www.nature.com/articles/s41467-022-30531-0), and [an atlas of human recoding sites](https://www.nature.com/articles/s41467-022-28841-4). <br /> 


query_known_sites.pl= excute mpileup (samtools) to query a list of known editing sites.<br />
parse_pileup_query.pl = a requirement for query_known_sites.pl<br />  
Usage: perl query_known_sites.pl [A predefined list of known editing sites] [STAR mapped bam file] [Output file name]
```ruby
perl query_known_sites.pl CNS_A2G_events.txt SampleName.bam OutputFileName.txt
```
<br />  

# Helpful data files:
CNS_A2G_events.txt = A predefined list of 166,215 A-to-I RNA editing sites detected within each cell population can be found [here](https://github.com/BreenMS/RNA-editing-in-CNS-cell-types).<br /> 

<br />  

# Supplemental Data Tables 1-5:
Supplemental Data 1. Summary of Living Brain Project demographics and RNA editing metrics, including external postmortem datasets. <br />  
Supplemental Data 2. RNA editing summary statistics across the Living Brain Project and secondary in vitro and in vivo experiments. <br />  
Supplemental Data 3. RNA editing summary metrics across secondary postmortem cohorts. <br />  
Supplemental Data 4. Cell-specific cataloging of A-to-I sites in postmortem human cortex. <br />  
Supplemental Data 5. RNA editing quantitative trait loci in the Living Brain Project. <br />  

# Supplemental Note 1
Supplemental Note 1. Description of cell-specific A-to-I sites derived from fluorescence activated nuclei sorted (FANS) neuronal and non-neuronal cell types from the postmortem human cortex. References are included within the supplemental note. <br />



