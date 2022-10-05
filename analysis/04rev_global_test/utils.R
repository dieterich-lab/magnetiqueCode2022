#' this function is part of the miRNAmRNA package
#' avaiable at 
#' https://git.lumc.nl/mvaniterson/mirnamrna/-/blob/master/R/mirnamrna.R#L10
gtTable <- function(object)
{
  test <- function(set) object@functions$test(set, calculateP=TRUE)

  leaves <- t(sapply(1:size(object), function(i) test(i)))

  ##calculate zscores
  zsc <- (leaves[,"S"]  - leaves[,"ES"]) / leaves[,"sdS"]

  ##association
  positive <- object@functions$positive()

  tbl <- data.frame(
    Pvalue=leaves[,"p"], 
    Association=positive, 
    Weights = weights(object), 
    zscores=zsc)

  rownames(tbl) <- object@functions$cov.names(1:size(object))
  tbl
}
dtu <- readRDS("MAGNetApp/data/DTU/summarized_experiment.RDS")

load_Y <- function(dtu) {
  counts <- SummarizedExperiment::assays(dtu)[["counts"]]

  tx_match <- match(
    SummarizedExperiment::rowData(dtu)$transcript_id, 
    rownames(counts))
  transcript_id <- SummarizedExperiment::rowData(dtu)$transcript_id[tx_match]
  gene_id <- SummarizedExperiment::rowData(dtu)$gene_id[tx_match]
  gene_counts <- aggregate(counts, list(gene_id), sum)
  rownames(gene_counts) <- gene_counts$Group.1
  gene_counts$Group.1 <- NULL
  gene_counts
}
dge_list <- readRDS("MAGNetApp/data/DGE/MAGNet_DCMvsHCM_GeneTonic.rds")
dge <- dge_list$res_de

load('rbp.targets.rda')

load_W <- function() {
  target.mat <- ComplexHeatmap::list_to_matrix(rbp.targets)
  x <- readr::read_csv(
  "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_string_to_int_ID_conversion.csv.gz",
  col_names = c(
    "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name",
    "ensembl_transcript_id_INT", "ensembl_gene_id_INT"))

  new_colnames <- x$ensembl_gene_id[
    match(
      colnames(target.mat), 
      x$external_gene_name)]
  colnames(target.mat) <- new_colnames
  target.mat[, !is.na(colnames(target.mat))]
}


#' X is the mRNA expression
#'
load_X <- function(dtu, type="counts") {
  .x <- SummarizedExperiment::assays(dtu)[[type]]
  .x <- .x[rownames(base_X)]

  .x <- cbind(t(.x), base_X)
  .x
}
