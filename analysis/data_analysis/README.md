# Generating CARNIVAL inputs
Here are provided the scripts used to generate the inputs needed to run the [CARNIVAL analysis](https://github.com/dieterich-lab/magnetiqueCode2022/tree/main/analysis/network_analysis). 
All the results generated here have been provided on the [/output](https://github.com/dieterich-lab/magnetiqueCode2022/tree/main/analysis/data_analysis/output) directory and they consist of:

  + Estimated Transcription Factor (TF) activity scores from the differential gene expression (DGE) data as well their significance.
  + A prior knowledge of signed and directed protein interactions as retreived from the [OmniPath](https://github.com/saezlab/OmnipathR) resource and filtered to contain only expressed genes.

For the analysis, the users must simply run the ```prepare_carnival_inputs.R``` script which additionally calls the support functions in ```estimate_significance.R``` script.

The outputs generated here will later be used as inputs for the netwrok analysis with CARNIVAL as presented [here](https://github.com/dieterich-lab/magnetiqueCode2022/tree/main/analysis/network_analysis).
