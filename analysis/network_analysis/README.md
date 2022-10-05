# CARNIVAL Analysis
Here are provided the scripts used for generating the contextualized protein interaction networks by using [CARNIVAL](https://www.nature.com/articles/s41540-019-0118-z) for the three comparisons: _DCMvsHCM_, _DCMvsNFD_ and _HCMvsNFD_.

The analyses take as inputs the data generated in the [data_analysis](https://github.com/dieterich-lab/magnetiqueCode2022/tree/main/analysis/data_analysis) directory and each directory contains the analysis scripts as well as the results for each comparison.

In order to generate the results, the users can simply run the ```analysisCARNIVAL_*.R``` scripts present on each directory and which additionally uses supporting functions which are present in the ```process_results.R``` script.

Network results as well as node and edge attributes have been provided in a tabulated _.txt_ format, while visualizations have been provided in the [Cytoscape](https://cytoscape.org/) (_.cys_) as well as graphical (_.png_) format.

**NOTE:** In order to run the _runCARNIVAL()_ functionality, please set the [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) solver path to your correct local path.
