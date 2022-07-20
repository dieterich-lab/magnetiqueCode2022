#! /usr/bin/env Rscript
## ---------------------------
## Script name: load_oRNAment.R
## Purpose of script: load oRNAment data for 
## Author: Thiago Britto-Borges
## Date Created: 2022-05-18
## Copyright (c) Thiago Britto-Borges
## Email: thiago.brittoborges@uni-heidelberg.de
## ---------------------------
## Notes: hs stands for homo sapines
## ---------------------------

message(commandArgs())

suppressPackageStartupMessages({
  library(dplyr)
})

ornament_hs_int2id <- readr::read_csv(
  "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_string_to_int_ID_conversion.csv.gz",
  col_names = c(
    "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name",
    "ensembl_transcript_id_INT", "ensembl_gene_id_INT"))

ornament_hs <- readr::read_csv(
  "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_cDNA_oRNAment.csv.gz",
  col_names = c(
    "ensembl_gene_id", "ensembl_transcript_id", "gene_biotype",
    "transcript_biotype", "transcript_position", "RBP", "score",
    "unpaired_probability", "chromosome", "region", "exon_start", "exon_end"
  )
)
ornament_hs <- ornament_hs %>%
  filter(score > 0.95) # canonical motif 1

ornament_hs <- left_join(
  ornament_hs, 
  ornament_hs_int2id, 
  by=c("ensembl_transcript_id"='ensembl_transcript_id_INT'))
rbp_encoding <- readr::read_csv(
  "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/RBP_id_encoding.csv.gz",
  col_names = c("code", "RBP_name"))

ornament_hs <- ornament_hs %>%
  mutate(RBP = as.numeric(RBP)) %>%
  left_join(rbp_encoding, by=c("RBP"='code'))

ornament_hs$RBP_name <- gsub(x=ornament_hs$RBP_name, '(.*) (.*)', '\\1')

save(rbp.targets, file = "rbp.targets.rda")
