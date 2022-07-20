## Supplementary resources

> Thiago Britto-Borges, Annekathrin Ludt, Etienne Boileau, Enio Gjerga, Federico Marini, Christoph Dieterich
> Magnetique: An interactive web application to explore transcriptome signatures of heart failure

## How to use this repository

This repository contains source codes for the analyses of the MAGNet RNA-seq dataset that are available
via the [Magnetique application](https://shiny.dieterichlab.org/app/magnetique).


Under analyses:

### data

Data preparation for the application is based on the gene abundance table generated with StringTie2 (see manuscript for details).

#### Sample-level filtering

We used RNA-seq data from the Myocardial Applied Genomics Network (MAGNet) (GEO: GSE141910), consisting of whole-transcriptomes of left ventricle (LV) tissues from end-stage heart failure (HF) patients due to dilated cardiomyopathy (DCM, n = 165) or hypertrophic cardiomyopathy (HCM, n = 27), and from unmatched non-failing hearts from organ donors (NFD, n = 162).

#### Gene-level filtering

We applied a threshold on expression (raw counts), on variance (remove 1% of genes with lowest variance), and normalised expression (remove 1% of genes with lowest average expression).

### pvca

Surrogate variable and principal variance component analysis.

### gene_expression

MAGNet_DGE.R: Perform differential gene expression analysis using DESeq (details in the manuscript).
MAGNet_GeneTonic.R: Data wrangling for the application.
figures_supp.R: Supplementary figures (gene enrichment analyses).

### dtu

drimseq_age_sex_race_sva.R: does the transcript isoform pre-filtering for DRIMSeq.  
DRIMSeq_DTU.ipynb: runs DRIMSeq.

### rev_global_test

01_load_oRNAment.R: loads the RBP:RNA interaction file.  
02_rev_globaltest.R: runs the reverse global test and processes the output.
