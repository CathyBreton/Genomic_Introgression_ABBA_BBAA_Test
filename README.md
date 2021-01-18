# *Musa_ABBA_BBAA_Introgession*


![Banniere](images/ABBA_Test.png)


## Introduction

<div align="justify">
Hybridization between species represents a major force of evolution, influenced by the external element, as glacial period, migration species or climatic changes, The genome can evolved by providing material for adaptation by natural and or sexual selection. Hybridization can decrease the difference between two species by sharing alleles across the genome, but can also act as a source of variation, impacting adaptation, aiding in evolutionary rescue, promoting range expansion, leading to species divergence, and finally fueling adaptative radiation.  

Since the beginning of  next generation sequencing, and because of the decreasing cost,  the number of whole-genome sequencing increase. Associated to new statistical methods to detect the signature of hybridization at the whole genome or chromosome level, genome sequencing technic provide an information patterns (SNP) across a tree as markers of hybridization. 
<div>


Purpose of Musa_ABBA_BBAA_Introgession
--------------------------------------




<div align="justify">
Due to the explosive expansion in genomic resources, scientist have developed several statistical tests to detect introgression. Patterson’s D-statistic also known as the ABBA-BABA test was developed to quantify the amount of genetic exchange. It considers of an ancestral “A” allele  and derived “B” alleles via mutation across the genome of four taxa. Under the hypothesis  “without introgression” the two allelic patterns “ABBA” or “BABA” occurred with equal frequency (((A,B))B)A) = (((B,A))B)A). An excess of “ABBA” or “BABA” shown by a D-statistic significantly different from zero indicate a gene flow between two taxa. A D-Statistic > 0 means an excess of ABBA indicates an introgression between population P2 to population P3, provided that P1 and P3 are not exchanging gene flow. Whereas D-Statistic < 0 which is an excess of BABA indicate an introgression between P1 and P3. To detect potential past hybridization, we used the ABBA-BABA test ( Martin …..), Patterson’s D test is D = [sum(ABBA) – sum(BABA)] / [sum(ABBA) + sum(BABA)] with ABBA = (1- p1 ) x p2 x p3 x (1- pO ), and BABA = p1 x (1- p2 ) x p3 x (1- pO ). To compute the standard error (Green et al 2010), we used the block jackknife approach. The number of ABBA, BBAA BABA sites was calculate with the workflow suite of Martin described in http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics. 
<div>
  
  
**Workflow - Calculate ABBA_BBA Python test Simon Martin Ref**

**Workflow - Obtain a vcf file **

<div align="justify">
The vcf file previously obtained with GATK version 4 was transformed in a geno file with the script parseVCF.py, some filters --minQual=20 and flag=DP min=5 were applied. On the genotype file, freq.py was used to calculate de frequency of each allele of each population. All the population P1, P2, P3 and outgroup were defined to test the hypothesis on banksii introgression in Papouasie New Guinea (PNG). To calculate the ABBA-BBAA number per windows the script ABBABABAwindows.py was applied with a windows of 10Mb. To compute the variance of D despite non-independence among site we used the jackknife. The block size needs to exceed the distance at which autocorrelation occurs, then we choose 10Mb.  and then to calculate the error and Z score, we used the script calculate_abba_baba_Musa.r adapted to our data. 
<div>



<details>
<summary>Table of content</summary>

## Table of contents

- [**How to cite**](#How-to-cite)
- [**Introduction**](#Introduction)
  - ABBA_BBAA Test
  - Python Test
  - R Test
- [**Workflow - Calculate ABBA_BBA Python test**](#workflow---molecular-karyotype-analysis)
  - Input raw data
  - Read Quality check
  - **Step a : Mapping reads on the reference**
  - DNA Data
  - RNA Data
  - **Step b : Variant discovery**
- [**Workflow - Calculate ABBA_BBA R test**](#workflow---molecular-karyotype-analysis)
  - Merge datasets
  - Filter SNP dataset
  - Split VCF by chromosome
  - Generate molecular karyotype
- [**Authors and acknowledgments**](#authors-and-acknowledgment) 
- [**Contact**](#contact) 

</details>




Dependencies
------------
The tools are developed in Perl, bash, Python3, Java and work on the Linux system and require:

| Tools  | Website | Version |
| ------ | ------- | ------- |
| Bamtools      | https://github.com/pezmaster31/bamtools                         | bamtools/2.4.0 |
| BWA           | http://bio-bwa.sourceforge.net                                  | bwa/0.7.12 |
| Cutadapt      | https://cutadapt.readthedocs.io/en/stable/                      | cutadapt/2.10  |
| FastQC        | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/      | FastQC/0.11.7 |
| GATK V4       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/4.0.5.2 |
| GATK V3       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/3.7-0   |
| Picard Tools  | https://broadinstitute.github.io/picard/                        | picard-tools/2.7.0   |
| sambamba      | https://lomereiter.github.io/sambamba/                          | sambamba/0.6.6 |
| Samtools      | https://github.com/samtools/samtools                            | samtools/1.2  |
| STAR          | https://github.com/alexdobin/STAR                               | STAR/2.5.0b |
| VCFHunter     | https://github.com/SouthGreenPlatform/VcfHunter                 |  |
| Vcftools      | https://vcftools.github.io/index.html                           | vcftools/0.1.14  |



How to cite
-----------
<div align="justify">
Filling the gaps in gene banks: Collecting, characterizing and phenotyping wild banana relatives of Papua new guinea David Eyland, Catherine Breton, Julie Sardos, Simon Kallow, Bart Panis, Rony Swennen, Janet Paofa, François Tardieu, Claude Welcker Steven B. Janssens, Sebastien C. Carpentier. *Crop Science* https://doi.org/10.1002/csc2.20320  
</div>

  
If you use the second workflow.
<div align="justify">
A protocol for detection of large chromosome variations in banana using Next Generation Sequencing. Breton Catherine, Cenci Alberto, Sardos Julie, Chase Rachel, Ruas Max, Rouard Mathieu & Roux Nicolas Published in October. Under review.
</div>

## Authors and acknowledgments

This work is a collaborative work between Catherine Breton, Yann Huber, Mathieu Rouard with the participation, and the use of scripts from Guillaume Martin (CIRAD) who develops and maintains VCFHunter.

## Contact

**Catherine Breton**, Alliance of Bioversity International and CIAT Europe (c.breton@cgiar.org)

The Alliance of Bioversity International and the International Center for Tropical Agriculture (CIAT)
delivers research-based solutions that harness agricultural biodiversity and sustainably transform
food systems to improve people’s lives in a climate crisis.
The Alliance is part of CGIAR, a global research partnership for a food-secure future.
https://www.bioversityinternational.org/       
https://www.ciat.cgiar.org



![Alliance](images/Alliance_logo_wide2.jpg)
