#!/biosw/R/3.6.2/bin/Rscript

# Usage: ./run_scdc.R [config]

library(dplyr)
library(stringr)
library(purrr)

library(Biobase)
library(SCDC)

# Use SCDC local install with changes to SCDC_basis

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


## Call

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
        print("Applying deconvolution to control and condition(s) separately...")
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


