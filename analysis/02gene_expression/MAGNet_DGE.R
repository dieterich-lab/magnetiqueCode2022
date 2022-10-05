#! /usr/bin/env Rscript

# Data needs to be downloaded first

# ```
# #!/bin/sh
# wget --output-document data.zip 'https://data.dieterichlab.org/s/3gCTLGT4DaAqaqb/download' --header 'Accept: text/html' --header 'Connection: keep-alive'
# unzip -o data.zip
# rm data.zip
# rm MAGNetApp/magnetique.sqlite MAGNetApp/data/colData.txt
# rm -r MAGNetApp/data/DGE MAGNetApp/data/DTU MAGNetApp/data/networks MAGNetApp/data/RBP
# 
# ```

# or use the newly generate data, e.g. in magnetiqueCode2022/analysis/data/

library(DESeq2)
library(IHW)

library(sva)

library("org.Hs.eg.db")

library(dplyr)
library(tibble)

library(openxlsx)


## Params

lfcThreshold.set <- 0 # log2(1.5)
altHypothesis.set <- "greaterAbs"
alpha.set <- 0.05

## Functions

# Find DESEq normalization factors
findDeseqFactorsSingle <- function(count_data)
{
  loggeomeans <- rowMeans(log(count_data))
  deseqFactors <- apply(count_data, 2, deseq, loggeomeans = loggeomeans) 
  deseqFactors
}

deseq <- function(x, loggeomeans) {
  finitePositive <- is.finite(loggeomeans) & x > 0
  if (any(finitePositive)) {
    res <- exp(median((log(x) - loggeomeans)[finitePositive], na.rm = TRUE))
  } else {
    print(utils::head(x))
    stop("Can't normalise accross a condition.
         Too many zero expressed genes. ")
  }
  res
}

# Get DESEq results
# strict filtering to increase power is automatically applied via independent
# filtering on the mean of normalized counts within the results function
deseq_results <- function (dds, num, denom) { 

    res <- results(dds,
                   contrast=c("Etiology", num, denom),
                   lfcThreshold=lfcThreshold.set, 
                   altHypothesis=altHypothesis.set,
                   alpha=alpha.set,
                   filterFun=ihw)
          
    res$SYMBOL <- mapIds(org.Hs.eg.db,
                         keys=rownames(res),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
                         
    # we use ashr, otherwise with apeglm we have to provide
    # coef, and for DCM vs HCM, would have to relevel...
    
    # in fact, for ashr, if res is provided, then coef and contrast are ignored...
    # but res here only contains the contrast of interest.
    
    # lfcThreshold=lfcThreshold.set unused with ashr
    
    res.shrunken <- lfcShrink(dds, 
                              contrast=c("Etiology", num, denom), 
                              res=res,
                              type="ashr")
            
    res.shrunken$SYMBOL <- mapIds(org.Hs.eg.db,
                                  keys=rownames(res.shrunken),
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")
    
    res.tib <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
    
    res.shrunken.tib <- res.shrunken %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
        
    # write to disk, add size factors for reference
    wb <- createWorkbook()

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, sep=""))
    writeDataTable(wb, sheet=1, x=res.tib)

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, "_shrunken", sep=""))
    writeDataTable(wb, sheet=2, x=res.shrunken.tib)
    
    addWorksheet(wb, sheetName="sizeFactors")
    
    sf <- dds$sizeFactor %>%
        data.frame() %>%    
        rownames_to_column(var="sample") %>% 
        dplyr::rename(sizeFactor = ".") %>%
        as_tibble()
            
    writeDataTable(wb, sheet=3, x= sf)
    
    filen <- paste0("condition_", num, "_vs_", denom, ".xlsx", sep="")
    filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
    saveWorkbook(wb, filen, overwrite=TRUE)
    
    res
}

## I/O

# colData, countData
dirloc.in <- 'magnetiqueCode2022'
dirloc.out <- 'magnetiqueCode2022'

load(file.path(dirloc.in, 'MAGNetApp', 'data', 'MAGNet.RData'))

stopifnot(all(rownames(colData) == colnames(countData)))

## SVA

# categorical as factor
colData$AFib <- factor(colData$AFib, levels=c("No", "Yes"))
colData$VTVF <- factor(colData$VTVF, levels=c("No", "Yes"))
colData$Diabetes <- factor(colData$Diabetes, levels=c("No", "Yes"))
colData$Hypertension <- factor(colData$VTVF, levels=c("No", "Yes"))
colData$LibraryPool <- factor(colData$LibraryPool)
colData$TissueSource <- factor(colData$TissueSource)

# use normalized data
normData <- t(t(countData) / findDeseqFactorsSingle(countData))

# based on our analysis, we assume we're not removing any source of biological variability, since
# our observations suggest that surrogate variables mostly correlate with technical artefacts (mostly DuplicationRate, TIN, libraryPool)
# 2 surrogate variables seem enough to capture a lot of the variance that was in LibraryPool and TIN.
n.sv <- 2

# including both the adjustment variables and the variable of interest
mod <- model.matrix(~Etiology+Race+Age+Sex, data=colData)
# including only the adjustment variables
mod0 <- model.matrix(~Race+Age+Sex, data=colData) 
svaObj <- svaseq(normData, mod, mod0, n.sv=n.sv) 

# add to colData
colData$SV1 <- svaObj$sv[,1]
colData$SV2 <- svaObj$sv[,2]

filen <- file.path(dirloc.out, 'MAGNet_SVA.rds')
saveRDS(svaObj, file=filen)

write.table(colData, file=file.path(dirloc.out, 'colData.txt'), sep=",", row.names=TRUE, col.names=TRUE, quote=FALSE)

## DESEq

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~Etiology+Race+Sex+Age+SV1+SV2)
dds <- DESeq(dds)

filen <- file.path(dirloc.out, 'MAGNet_DESeqDataSet.rds')
saveRDS(dds, file=filen)

# loop over all possible contrasts
contrasts <- data.frame(num=c('DCM', 'HCM', 'DCM'), denom=c('NFD', 'NFD', 'HCM'))
results <- apply(contrasts, 1, function(x) {
    deseq_results(dds, as.character(x[1]), as.character(x[2]))
} )
names(results) <- apply(contrasts, 1, paste, collapse="vs")

filen <- file.path(dirloc.out, 'MAGNet_DESeqResults.rds')
saveRDS(results, file=filen)

