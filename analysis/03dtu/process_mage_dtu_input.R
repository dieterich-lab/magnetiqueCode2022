#! /usr/bin/env Rscript
#' ---------------------------
#' Script name: DTU_mage.R
#' Purpose of script: DTU for MAGE project
#' Author: Thiago Britto-Borges
#' Date Created: 2021-04-06
#' Copyright (c) Thiago Britto-Borges, DieterichLab 2021
#' Email: thiago.brittoborges@uni-heidelberg.de
#' ---------------------------
#' :
#' 
#' ---------------------------
suppressPackageStartupMessages({
  library(tximport)
  library(tidyverse)
  library(rtracklayer)
  library(DRIMSeq)
  library(DEXSeq)
  library(sva) 
})

message("Loading and processing GTF")
gtf <- rtracklayer::import("/prj/MAGE/analysis/baltica/stringtie/merged_guided.gtf")
gtf$gene_id[!is.na(gtf$ref_gene_id)] <- gtf$ref_gene_id[!is.na(gtf$ref_gene_id)]
tx2gene <- gtf %>%
  as.data.frame() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select("gene_name", "transcript_id", "gene_id") %>%
  mutate(gene_name = coalesce(gene_name, gene_id))

  
message("Loading and procesing input.")

# samples.txt generated in 01data
metadata <- read.csv("magnetiqueCode2022/analysis/data/samples.txt")
metadata <- filter(metadata, Etiology != "PPCM")
metadata$Etiology <- factor(metadata$Etiology, ordered = FALSE )
metadata$Sex <- factor(metadata$Sex, ordered = FALSE )
metadata$Race <- factor(metadata$Race, ordered = FALSE )
metadata$Etiology <- relevel(metadata$Etiology, ref="NFD")
metadata$Etiology <- droplevels(metadata$Etiology)

# from Zenodo
counts <- read.csv("magnetiqueCode2022/transcript_count_matrix.csv")
metadata <- plyr::rename(metadata, c("Run" = "sample_id"))

d <- dmDSdata(
  counts = counts,
  samples = metadata[c("sample_id", "Age", "Sex", "Race", "Etiology")]
)
message("Filtering dataset data")
message("Total number of tx: ", sum(elementNROWS(d@counts)))
message("Total number of genes: ", length(d@counts))

# Transcript pre-filtering
d <- dmFilter(
  d,
  min_samps_gene_expr = 8,
  min_samps_feature_prop = 8,
  min_feature_prop = 0.1,
  min_gene_expr = 30
)

message("Total number of tx after filter: ", sum(elementNROWS(d@counts)))
message("Total number of genes after filter: ", length(d@counts))


# Preparing objects DRIMSeq and DEXSeq
sampleData <-  DRIMSeq::samples(d)
countData <- round(as.matrix(d@counts@unlistData), 0)
featureID <- rownames(d@counts)
groupID <- tx2gene[ match(rownames(d@counts), tx2gene$transcript_id), 'gene_id']

dxd <- DEXSeqDataSet(
  countData = countData,
  sampleData = sampleData,
  design = ~sample + exon + Age:exon + Sex:exon + Race:exon + Etiology:exon,
  featureID = rownames(d@counts),
  groupID = tx2gene[ match(rownames(d@counts), tx2gene$transcript_id), 'gene_id']
)
dxd <- estimateSizeFactors(dxd)

# SVA analisys
# SVA requires normalized input data
normCounts <- counts(dxd, normalize=TRUE)
# the first half of the columns contains the features counts
# see for this : other exons in DEXseq manual
normCounts <- normCounts[, 1: (ncol(normCounts)/2)]

n.sv <- 2
mod <- model.matrix(~Etiology+Race+Age+Sex, data=sampleData)
mod0 <- model.matrix(~Race+Age+Sex, data=sampleData) 

message('svaseq')
svseq <- svaseq(
  dat = normCounts,
  mod = mod, 
  mod0 = mod0, 
  n.sv=n.sv)

save.image("mage_dtu_input.RData")
