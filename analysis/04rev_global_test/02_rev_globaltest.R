#! /usr/bin/env Rscript
## ---------------------------
## Script name: rev_globaltest.R
## Purpose of script: use reverve global test to
## tackle the task of identifiy RBP:mRNA targets
## Author: Thiago Britto-Borges
## Date Created: 2022-05-20
## Copyright (c) Thiago Britto-Borges
## Email: thiago.brittoborges@uni-heidelberg.de
## ---------------------------
## Notes: the reverse global test was first described here:
## https://academic.oup.com/nar/article/41/15/e146/2411309
## ---------------------------
library(globaltest)
source("utils.R")

library(SummarizedExperiment)

base_X <- SummarizedExperiment::colData(dge_list$dds)
rownames(base_X) <- base_X$Run
base_X <- base_X[c("Etiology", "Age", "Sex", "Race", "SV1", "SV2")]
cov <- c("Etiology", "Age", "Sex", "Race", "SV1", "SV2")

library(tidyverse)
library(furrr)
library(future)
plan(multicore, workers = 40)

message("Start")
# incidence matrix relating regulators to mRNAs
# target x regulators
W <- load_W()
# matrix of regulators expression profiles
# regulators x samples
Y <- load_Y(dtu)
# matrix of targets expression profiles
# samples x regulators
X <- load_X(dtu)
Y <- t(Y)
common <- intersect(colnames(Y), colnames(W))
Y <- Y[, common]
W <- W[, common]
stopifnot(all(rownames(Y) == rownames(X)))
stopifnot(all(rownames(Y) == rownames(base_X)))
Y <- aggregate(Y, by = list(base_X$Etiology), mean)
rownames(Y) <- Y$Group.1
Y$Group.1 <- NA
X <- aggregate(X, by = list(base_X$Etiology), mean)
rownames(X) <- X$Group.1
X$Group.1 <- NULL

common <- intersect(rownames(W), colnames(X))
X <- X[, common]
W <- W[common, ]
# iterate over targets
message("Start loop")
revgt_result <- furrr::future_map(1:nrow(W), function(i) {
  target <- rownames(W)[i]

  w <- names(W[target, ] == 1)
  x <- Y[, w]
  y <- X[, target]

  gt_obj <- globaltest::gt(
    as.vector(y),
    as.matrix(x),
    directional = 1e-06
  )

  gt_table <- gt_obj %>%
    gtTable(.) %>%
    as.matrix(.) %>%
    reshape2::melt(.) %>%
    as.data.frame(.) %>%
    mutate(transcript_id = target)
  list(gt_obj = gt_obj, gt_table = gt_table)
})

save.image("rev_globaltest.Rdata")

df <- purrr::map(revgt_result, ~ .$gt_table) %>%
  bind_rows() %>%
  rename(c("Var1" = "gene_id")) %>%
  pivot_wider(
    id_cols = c(gene_id, transcript_id),
    values_from = value,
    names_from = Var2
  )

message("Finish sucessifuly.")

gtf <- rtracklayer::import("GRCh38.96.gtf")
tx <- subset(gtf, type == "transcript")
gene <- subset(gtf, type == "gene") %>%
  as_tibble() %>%
  dplyr::select(gene_id, gene_name)

anno <- mcols(tx)[
  c("gene_name", "transcript_name", "transcript_id", "transcript_biotype")
]

df <- df %>% left_join(as_tibble(anno), by = "transcript_id")

df <- df %>%
  left_join(as_tibble(gene), by = "gene_id", suffix = c("", "_regulator")) %>%
  rename(c("gene_id" = "gene_id_regulator")) %>%
  arrange(Pvalue, dplyr::desc(abs(zscores)))

df <- df %>% filter(Pvalue < 0.05)
write_excel_csv(df, "mirna_mrna_revgt_interactions.csv")

# Output CPEB1 targets
library(gt)

dtu_sig <- rowData(dtu)[["DRIMSeq_DCM_vs_HCM"]] %>%
  filter(adj_pvalue <= 0.05) %>%
  pull(feature_id)

x %>%
  filter(gene_name_regulator == "CPEB1") %>%
  mutate(fdr = p.adjust(Pvalue, method = "fdr")) %>%
  filter(fdr <= 0.05, transcript_id %in% dtu_sig) %>%
  select(
    gene_id_regulator,
    gene_name_regulator,
    transcript_name,
    transcript_id,
    transcript_biotype,
    fdr,
    Association
  ) %>%
  gt() %>%
  fmt_scientific(
    columns = fdr,
    decimals = 1
  ) %>%
  fmt_number(
    columns = zscores
  ) %>%
  gtsave("CPEB1_rev_global_test_and_dtu.html")
