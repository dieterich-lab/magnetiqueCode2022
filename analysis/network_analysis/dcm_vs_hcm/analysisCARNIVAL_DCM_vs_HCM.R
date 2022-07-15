#! /usr/bin/env Rscript

set.seed(123)

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
library("CARNIVAL")

load(file = "../../data_analysis/output/tfListAll.RData")
tf_activities_stat <- tfListAll$DCM_vs_HCM
topTF <- tf_activities_stat[which(tf_activities_stat$pval<=0.1), ]
measObj <- matrix(data = , nrow = 1, ncol = nrow(topTF))
measObj[1, ] <- topTF$nes
colnames(measObj) <- topTF$id
measObj <- as.data.frame(measObj)

load(file = "../../data_analysis/keep_genes.RData")
load(file = "../ppi.RData")

idx1 <- which(ppi$source%in%keep_genes)
idx2 <- which(ppi$target%in%keep_genes)
idx <- intersect(x = idx1, y = idx2)
ppi <- ppi[idx, ]

res <- runCARNIVAL(inputObj = NULL,
                   measObj = measObj,
                   netObj = ppi,
                   weightObj = NULL,
                   solverPath = "/beegfs/homes/egjerga/cplex",
                   solver = "cplex",
                   timelimit = 36000,
                   limitPop = 100,
                   poolCap = 100,
                   mipGAP = 0.1,
                   poolrelGAP = 0.1,
                   poolIntensity = 4)

save(res, file = "res_dcm_vs_hcm_signif_100_gap_01.RData")
