#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./run_scdc.R [config]

library(dplyr)
library(stringr)
library(purrr)

library(Biobase)
library(SCDC)


## I/O, settings

# get params and options
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) { stop("Usage: ./run_scdc.R [config]\n", call.=FALSE) }

params.file <- args[1]
# values are NOT "checked" for consistency, type, etc. before they are passed to SCDC
params <- yaml::read_yaml(params.file)
opts <- params$opts


# check if output path exist, create if not
dir.create(params$dirOutput, recursive=TRUE, showWarnings=FALSE)


## Functions


save_pheatmap_pdf <- function(x, filename, width=8, height=8) {
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

# redefine weighted basis matrix
# see https://github.com/meichendong/SCDC/issues/13
mySCDC_basis <- function(x, ct.sub = NULL, ct.varname, sample, ct.cell.size = NULL){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[,sample])
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  sigma <- sapply(unique(mean.id[, 1]), function(id) {
    y = mean.mat[, mean.id[, 1] %in% id]
    if (is.null(dim(y))){
      res = rep(0, length(y))
      message("Warning: the cell type [", id,"] is only available in at most 1 subject!")
    } else {
      res = apply(y, 1, var, na.rm = TRUE)
    }
    return(res)
  })

  sum.mat2 <- sapply(unique(sample.id), function(sid){
    sapply(unique(ct.id), function(id){
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  # library size factor calculated from the samples:
  if (is.null(ct.cell.size)){
    sum.mat <- rowMeans(sum.mat2, na.rm = T)
  } else {
    if (is.null(names(ct.cell.size))){
      message("Cell size factor vector requires cell type names...")
      break
    } else {
      sum.mat <- ct.cell.size
    }
  }

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x,...){
    y <- apply(x,1,max, na.rm = TRUE)
    #y / median(y, na.rm = T)
  }

  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      apply(y,1,var, na.rm=T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  # q15 <- apply(var.adj,2,quantile, probs = 0.15, na.rm =T)
  q15 <- apply(var.adj, 2, function(zz){
    z1 = min(zz[zz>0])
    z2 = quantile(zz, 0.15, na.rm = T)
    return(max(z1,z2))
  })
  q85 <- apply(var.adj,2,quantile, probs = 0.85, na.rm =T)

  var.adj.q <- t(apply(var.adj, 1, function(y){
                        y[y<q15] <- q15[y<q15]
                         y[y>q85] <- q85[y>q85]
                         return(y)}
                      )
                 ) #+ 1e-4

  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    sid = unlist(strsplit(id,'%'))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[,sid]), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}


## Call

# redefine
assignInNamespace("SCDC_basis",  mySCDC_basis, ns="SCDC")

# sc data
sc.eset.list <- map(params$scData, ~readRDS(.x))
# use all cell types present in the data unless specified
if (all(opts$cellTypes == 'all')) {
    combo <- lapply(sc.eset.list, function(l) {
                if (opts$varName %in% rownames(l@phenoData@varMetadata)) {
                        unlist(unique(l@phenoData@data[opts$varName]), use.names=FALSE)
                    }
                })
    ct.sub <- Reduce(intersect, combo)          
} else {
    ct.sub <- opts$cellTypes
}

if (opts$clustQC$use) {
    sc.eset.list <- map2(sc.eset.list, names(sc.eset.list), function(l, n) {
        qc <- SCDC_qc(l,
                      ct.varname=opts$varName, 
                      sample=opts$sampleName, 
                      scsetname=n,
                      ct.sub=ct.sub, 
                      qcthreshold=opts$clustQC$threshold)
        save_pheatmap_pdf(qc$heatfig, file.path(params$dirOutput, paste(n, '-pheatmap.pdf', sep='')))
        qc$sc.eset.qc
    })
}

# bulk
ex <- strsplit(basename(params$bulkData), split="\\.")[[1]][-1]
if (!(ex == 'loom')) {
    # expressionSet
    bulk.eset <- readRDS(params$bulkData)
    if (!(inherits(bulk.eset, "ExpressionSet"))) {
        stop(paste("Bulk ", params$bulkData, " is not an ExpressionSet", sep=""), call.=FALSE)
    }
} else {
    library(loomR)
    # see https://github.com/mojaveazure/loomR/issues/55
    lfile <- connect(filename=params$bulkData, mode="r+", skip.validate=TRUE)
    exprs <- t(lfile[['matrix']][,])
    rownames(exprs) <- lfile[["row_attrs/Gene"]][] # default loompy column and row attributes 
    colnames(exprs) <- lfile[["col_attrs/CellID"]][]
    # default custom attributes GeneFlag and SampleFlag
    exprs <- exprs[as.logical(as.numeric(lfile[["row_attrs/GeneFlag"]][])),as.logical(as.numeric(lfile[["col_attrs/SampleFlag"]][]))]
    
    pData <- lfile$get.attribute.df(MARGIN = 2, col.names = "CellID")
    pData <- pData[as.logical(as.numeric(lfile[["col_attrs/SampleFlag"]][])),]
    stopifnot(all(rownames(pData)==colnames(exprs)))
    pData$CellID <- NULL
    pData$SampleFlag <- NULL
    for (col in colnames(pData)) {
    if (is(pData[[col]], "character")) {
        pData[[col]] <- as.factor(pData[[col]])
        }
    }
    phenoData <- new("AnnotatedDataFrame", data=pData)
    experimentData <- new("MIAME", title="Bulk RNA-seq")
    
    bulk.eset <- ExpressionSet(assayData=exprs,
                               phenoData=phenoData,
                               experimentData=experimentData)
    lfile$close_all()
}

# tree-based deconvolution
# apply to each sc dataset
if (opts$tree$use) {
    print("Applying tree-based deconvolution...")
    # first assign "meta clusters"
    subcl <- map(sc.eset.list, function(e) {
        for (i in seq_along(params$opts$tree$scdcSubcl)) {
            e$scdcSubcl[e[[opts$varName]] %in% params$opts$tree$scdcSubcl[[i]]] <- names(params$opts$tree$scdcSubcl)[i]
        }
        e
    })
    # then estimate proportions
    subcl <- map(subcl, ~SCDC_prop_subcl_marker(bulk.eset=bulk.eset, 
                                                sc.eset=.x, 
                                                ct.varname=opts$varName, 
                                                fl.varname="scdcSubcl", 
                                                sample=opts$sampleName, 
                                                ct.sub=ct.sub, 
                                                ct.fl.sub=unique(.x$scdcSubcl), 
                                                select.marker=T, 
                                                LFC.lim=5))
    if (params$split$use) {
        print("Applying deconvolution to control and condition separately...")
        deconvoluted <- map(params$split$splitAttrs, ~SCDC_ENSEMBLE(bulk.eset=bulk.eset[,grepl(paste(.x, collapse='|'), bulk.eset[[params$split$nameAttr]])],
                                                                    prop.input=subcl,
                                                                    ct.sub=ct.sub,
                                                                    grid.search=opts$ensembl$gridSearch,
                                                                    search.length=opts$ensembl$searchLen))
    } else {
        deconvoluted <- SCDC_ENSEMBLE(bulk.eset=bulk.eset,
                                      prop.input=subcl,
                                      ct.sub=ct.sub,
                                      grid.search=opts$ensembl$gridSearch,
                                      search.length=opts$ensembl$searchLen)
    }

} else { # ENSEMBLE only

    if (params$split$use) {
        print("Applying deconvolution to control and condition separately...")
        deconvoluted <- map(params$split$splitAttrs, ~SCDC_ENSEMBLE(bulk.eset=bulk.eset[,grepl(paste(.x, collapse='|'), bulk.eset[[params$split$nameAttr]])],
                                                                    sc.eset.list=sc.eset.list,
                                                                    ct.varname=opts$varName,
                                                                    sample=opts$sampleName,
                                                                    ct.sub=ct.sub,
                                                                    grid.search=opts$ensembl$gridSearch,
                                                                    search.length=opts$ensembl$searchLen))
    } else {
        deconvoluted <- SCDC_ENSEMBLE(bulk.eset=bulk.eset,
                                      sc.eset.list=sc.eset.list,
                                      ct.varname=opts$varName,
                                      sample=opts$sampleName,
                                      ct.sub=ct.sub,
                                      grid.search=opts$ensembl$gridSearch,
                                      search.length=opts$ensembl$searchLen)
    }
}

filen <- paste(params$outputName, "RData", sep=".")
filen <- file.path(params$dirOutput, filen, fsep=.Platform$file.sep)
save(params, deconvoluted, file=filen)
print("Results written to disk. Terminating.")


