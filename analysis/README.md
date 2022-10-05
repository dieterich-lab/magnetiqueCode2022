
## Workflow

Reproducing the analyses involves several intermediate steps, mostly requiring working **R **, and **Jupyter** installations with
their respective dependencies (packages are listed at the beginning of each script/notebook). 


### 01data

The first step is [data preparation](01data/data_prep.ipynb). We start from the Stringtie2 gene abundance table (see manuscript for details) available from [Zenodo](https://zenodo.org/record/6854308). Additional metadata is also available. *Depending on where you downloaded the data, you may need to adjust the paths.*

#### Sample-level filtering

We used RNA-seq data from the Myocardial Applied Genomics Network (MAGNet) (GEO: GSE141910), consisting of whole-transcriptomes of left ventricle (LV) tissues from end-stage heart failure (HF) patients due to dilated cardiomyopathy (DCM, n = 165) or hypertrophic cardiomyopathy (HCM, n = 27), and from unmatched non-failing hearts from organ donors (NFD, n = 162). All PPCM patients are moved from the final files.

#### Gene-level filtering

We applied a threshold on expression (raw counts), on variance (remove 1% of genes with lowest variance), and normalised expression (remove 1% of genes with lowest average expression).

### 02gene_expression

The second step is to perform [differential expression](02gene_expression/MAGNet_DGE.R). To reproduce the analyses, you can either use the intermediate results available via

```
wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
unzip -o data.zip
rm data.zip
```

or use the newly generate data, *e.g.* in magnetiqueCode2022/analysis/data. *Depending on where you downloaded the data, or if you are using the newly generated data, you may need to adjust the paths.*


* MAGNet_DGE.R: Perform differential gene expression analysis using DESeq (details in the manuscript).
* MAGNet_GeneTonic.R: Data wrangling for the application.
* figures_supp.R: Supplementary figures (gene enrichment analyses).

### 03dtu

One can then perform differential transcript usage (DTU). The first step is to [prepare the input](03dtu/process_mage_dtu_input.R), [fit the model](03dtu/drimseq_age_sex_race_sva.R), and finally [run DTU analyses](03dtu/DRIMSeq_DTU.ipynb). The input transcript count table and GTF annotations are available from [Zenodo](https://zenodo.org/record/6854308). The file *samples.txt* was generated above. *Depending on where you downloaded the data, you may need to adjust the paths.*


### 04rev_global_test

To generate the results corresponding to the **RBP:RNA View**. To reproduce the analyses, you can either use the intermediate results available via

```
wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
unzip -o data.zip
rm data.zip
```

or use the newly generate data (01data, 02gene_expression, 03dtu). *Depending on where you downloaded the data, or if you are using the newly generated data, you may need to adjust the paths.* See [utils.R](04rev_global_test/utils.R).


* utils.R: load results (DGE, DTU, *etc.* )
* 01_load_oRNAment.R: loads the RBP:RNA interaction file.
* 02_rev_globaltest.R: runs the reverse global test and processes the output.


### 05data_analysis

**Generating CARNIVAL inputs**

Here are provided the scripts used to generate the inputs needed to run the [CARNIVAL analysis](06network_analysis). 
All the results generated here have been provided on the [/output](05data_analysis/output) directory and they consist of:

  + Estimated Transcription Factor (TF) activity scores from the differential gene expression (DGE) data as well their significance.
  + A prior knowledge of signed and directed protein interactions as retreived from the [OmniPath](https://github.com/saezlab/OmnipathR) resource and filtered to contain only expressed genes.

For the analysis, the users must simply run the `prepare_carnival_inputs.R` script which additionally calls the support functions in `estimate_significance.R` script.

The outputs generated here will later be used as inputs for the network analysis with CARNIVAL as presented [here](06network_analysis).


### 06network_analysis

**CARNIVAL analysis**

Here are provided the scripts used for generating the contextualized protein interaction networks by using [CARNIVAL](https://www.nature.com/articles/s41540-019-0118-z) for the three comparisons: _DCMvsHCM_, _DCMvsNFD_ and _HCMvsNFD_.

The analyses take as inputs the data generated in the [data_analysis](05data_analysis) directory and each directory contains the analysis scripts as well as the results for each comparison.

In order to generate the results, the users can simply run the `analysisCARNIVAL_*.R` scripts present on each directory and which additionally uses supporting functions which are present in the `process_results.R` script.

Network results as well as node and edge attributes have been provided in a tabulated _.txt_ format, while visualizations have been provided in the [Cytoscape](https://cytoscape.org/) (_.cys_) as well as graphical (_.png_) format.

**NOTE:** In order to run the _runCARNIVAL()_ functionality, please set the [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) solver path to your correct local path. *You may need to adjust the paths.*


### pvca

Jupyter notebook with additional surrogate variable and principal variance component analyses. To reproduce the analyses, you can use the intermediate results available via

```
wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
unzip -o data.zip
rm data.zip
```



