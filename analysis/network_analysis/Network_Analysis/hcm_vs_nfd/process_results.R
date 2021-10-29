library("readxl")
library("readr")
library("XML")
library("igraph")
library("foreach")
library("doParallel")
library("ranger")
library("palmerpenguins")
library("tidyverse")
library("kableExtra")

# 100
load(file = "res_hcm_vs_nfd_signif_100_gap_01.RData")
load(file = "../../Data_Analysis/output/tfListAll.RData")

network <- res$weightedSIF
attributes <- res$nodesAttributes
TFs <- tfListAll$HCM_vs_NFD$id[which(tfListAll$HCM_vs_NFD$pval<=0.05)]

ctrl <- 1
while(ctrl!=0){
  
  dd <- setdiff(x = setdiff(x = network[, 3], y = network[, 1]), y = TFs)
  idx2rem <- which(network[, 3]%in%dd)
  
  if(length(idx2rem)>0){
    network <- network[-idx2rem, ]
  } else {
    ctrl <- 0
  }
  
}

write.table(x = network, file = "network_hcm_vs_nfd_100.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = attributes, file = "attributes_hcm_vs_nfd_100.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
