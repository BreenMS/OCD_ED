# Divergent landscapes of A-to-I editing in postmortem and living human brain

Miguel Rodriguez de los Santos1, Brian H. Kopell1, Ariela Buxbaum Grice1, Gauri Ganesh1, Andy Yang1, Pardis Amini1, Lora E. Liharska1, Eric Vornholt1, John F. Fullard1, Pengfei Dong1, Eric Park1, Sarah Zipkowitz1, Deepak A. Kaji1, Ryan C. Thompson1, Donjing Liu1, You Jeong Park1, Esther Cheng1, Kimia Ziafat1, Emily Moya1, Brian Fennessy1, Lillian Wilkins1, Hannah Silk1, Lisa M. Linares1, Brendan Sullivan1, Vanessa Cohen1, Prashant Kota1, Claudia Feng1, Jessica S. Johnson1, Marysia-Kolbe Rieder1, Joseph Scarpa1, Girish N. Nadkarni1, Minghui Wang1, Bin Zhang1, Pamela Sklar1, Noam D. Beckmann1, Eric E. Schadt1, Panos Roussos1, Alexander W. Charney1, Michael S. Breen1

Affiliations: 1, Icahn School of Medicine at Mount Sinai<br /> 
**Correspondence to: michael[dot] breen [at] mssm [dot] edu<br /> 
[@breenPsychgene](https://twitter.com/breenPsychGene)<br /> 

Adenosine-to-inosine (A-to-I) editing is a prevalent post-transcriptional RNA modification within the brain. Yet, most research has relied on postmortem samples, assuming it is an accurate representation of RNA biology in the living brain. We challenge this assumption by comparing A-to-I editing between postmortem and living prefrontal cortical tissues. Major differences were found, with over 70,000 A-to-I sites showing higher editing levels in postmortem tissues. Increased A-to-I editing in postmortem tissues is linked to higher ADAR and ADARB1 expression, is more pronounced in non-neuronal cells, and indicative of postmortem activation of inflammation and hypoxia. Higher A-to-I editing in living tissues marks sites that are evolutionarily preserved, synaptic, developmentally timed, and disrupted in neurological conditions. Common genetic variants were also found to differentially affect A-to-I editing levels in living versus postmortem tissues. Collectively, these discoveries illuminate the nuanced functions and intricate regulatory mechanisms of RNA editing within the human brain. <br /> <br /> 

![Abstract (1)](https://github.com/BreenMS/Living-Brain/blob/main/Figure_1.png)

Here we described the main computational code used to generate all results and figures in this body of work.  

All RNA editing matrices and sample level RNA editing sites derived from the Living Brain Project will be available shortly on [Synapse.org](https://www.synapse.org/#!Synapse:syn26434508/files/).  

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



