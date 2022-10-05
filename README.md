## Supplementary resources

> Thiago Britto-Borges, Annekathrin Ludt, Etienne Boileau, Enio Gjerga, Federico Marini, Christoph Dieterich  
> Magnetique: An interactive web application to explore transcriptome signatures of heart failure

## How to use this repository

This repository contains source codes for the analyses of the MAGNet RNA-seq dataset that are available
via the [Magnetique application](https://shiny.dieterichlab.org/app/magnetique) (Source code available [here](https://github.com/AnnekathrinSilvia/magnetique/)).

Read counts, associated data, and the final DB dump can be downloaded from [Zenodo](https://zenodo.org/record/6854308).
Intermediate results are also available via

```
wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
unzip -o data.zip
rm data.zip
```

The different analysis scripts are available under [analysis](https://github.com/dieterich-lab/magnetiqueCode2022/tree/main/analysis).
