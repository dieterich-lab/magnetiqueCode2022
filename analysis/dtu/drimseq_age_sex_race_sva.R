#! /usr/bin/env Rscript
## ---------------------------
## Script name: DTU_mage.R
## Purpose of script: DTU for MAGE project
## Author: Thiago Britto-Borges
## Date Created: 2021-04-06
## Copyright (c) Thiago Britto-Borges, DieterichLab 2021
## Email: thiago.brittoborges@uni-heidelberg.de
## ---------------------------
## Notes:
## ---------------------------
message(commandArgs())
suppressPackageStartupMessages({
  library(DRIMSeq)
  library(BiocParallel)
  library(DEXSeq)
})
message("DRIMSeq age_sex_race_sva start")

register(MulticoreParam(workers=30))

message("DRIMSeq load")
load("mage_dtu_input.RData")

d <- dmFilter(
  d,
  min_samps_gene_expr = 50,
  min_samps_feature_prop = 50,
  min_feature_prop = 0.2,
  min_gene_expr = 50
)

message("Total number of tx after filter: ", sum(elementNROWS(d@counts)))
message("Total number of genes after filter: ", length(d@counts))

sampleData$SV1 <- svseq$sv[, 1]
sampleData$SV2 <- svseq$sv[, 2]
sampleData$Age <- scale(sampleData$Age, center = TRUE, scale = TRUE)

design_full <- model.matrix(~ Etiology + Race + Age + Sex + SV1 + SV2, data = sampleData)

message("DRIMSeq computing precision")
d <- dmPrecision(d, design = design_full, add_uniform = TRUE, BPPARAM = pparam())
saveRDS(d, "drimseq_age_sex_race_sva.RDS")
message("DRIMSeq: fit model")
d <- dmFit(d, design = design_full, verbose = 1, BPPARAM = pparam())
saveRDS(d, "drimseq_age_sex_race_sva.RDS")
message("DRIMSeq: test")
sessionInfo()
