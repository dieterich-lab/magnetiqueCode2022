#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./bootstrap_wgcna.R [config] <1/2>
# 1: (Default) Construct network using the original data 
# 2: Perform bootstrapping, soft threshold (beta) is given and fixed

library(WGCNA)
options(stringsAsFactors=FALSE)


#######

# even with sys BLAS... this seems
# to take too much time...?

# Allow multi-threading within WGCNA
# deprecated?
#allowWGCNAThreads()

# + set env?
enableWGCNAThreads(nThreads = 24)

#cpucores <- 16
#require(parallel)
#options("mc.cores"=cpucores)
#
##Set CPU cores for doParallel-related functions
#require(doParallel)
#cores <- makeCluster(detectCores(), type='PSOCK')
#registerDoParallel(cores)

###########


library(flashClust)
library(dplyr)

library("RColorBrewer")

## I/O, settings

# get params and options
args <- commandArgs(trailingOnly=TRUE)
pick <- 1
if (length(args)<1) {
  stop("At least configuration file must be given.\n", call.=FALSE)
} else {
  if (length(args)>1 & as.integer(args[2])<3) {
    pick <- as.integer(args[2])
  }
}

select <- switch(
   pick,
   "constructNetwork",
   "bootstrap")

params.file <- args[1]
# values are NOT checked for consistency, type, etc. before they are passed to wgcna, hclust, etc.
# for correlation, Spearman and bicor are recognised, anything else is Pearson
# some options are hard coded (maxPOutliers, etc. see below, many are WGCNA defaults)
params <- yaml::read_yaml(params.file)
opts <- params$opts
opts$select <- select

verbose <- 5
opts$verbose <- verbose

# used in cutreeDynamic, see below
deepSplit <- 1
opts$pruning$deepSplit <- deepSplit
method <- "hybrid" # between hclust and pam
opts$pruning$method <- method
# Maximum dissimilarity that qualifies modules for merging.
# MEDissThres <- 0.2 # correspond to a correlation of 1-MEDissThres
MEDissThres <- 0.3 # <<< run again, we should include this in config?
opts$pruning$MEDissThres <- MEDissThres

# correlation, defaults to bicor
opts$network$cor <- tolower(opts$network$cor)
if (!opts$network$cor %in% c("pearson", "spearman", "bicor")) {
    opts$network$cor <- "bicor"
}

# final soft threshold
opts$network$beta <- opts$network$betaFinal

# output base name - network options are written here, but
# in case of mismatch, are not changed in the network base name...
# in case of warnings, make sure to adjust config/options!
basename <- paste(params$basename, ".",
                  opts$network$type, "-", 
                  opts$network$cor, ".", 
                  opts$network$clust, ".", sep="")

                  
# check if output path exist, create if not
dir.create(params$networkDir, recursive=TRUE, showWarnings=FALSE)


## Functions

# define functions based on chosen correlation
allAdj <- list()
allAdj[["spearman"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                            power=opts$network$beta, 
                                                            type=opts$network$type, 
                                                            corOptions=list(use='p', method='spearman')) 
                                                }
allAdj[["bicor"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                         power=opts$network$beta, 
                                                         type=opts$network$type, 
                                                         corFnc="bicor", 
                                                         corOptions="use='p', maxPOutliers=0.1")
                                             }
allAdj[["pearson"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                           power=opts$network$beta, 
                                                           type=opts$network$type, 
                                                           corOptions=list(use='p')) 
                                               }                                                
get_adjacency <- allAdj[[tolower(opts$network$cor)]]

allMerge <- list()
allMerge[["spearman"]] <- function(datExpr, colors, eigengenes, opts) { 
                                        mergeCloseModules(datExpr, 
                                                          colors, 
                                                          MEs=eigengenes, 
                                                          cutHeight=opts$pruning$MEDissThres, 
                                                          corOptions=list(use='p', method='spearman'), 
                                                          verbose=opts$verbose) }
allMerge[["bicor"]] <- function(datExpr, colors, eigengenes, opts) { 
                                        mergeCloseModules(datExpr, 
                                                          colors, 
                                                          MEs=eigengenes, 
                                                          cutHeight=opts$pruning$MEDissThres, 
                                                          corFnc="bicor",
                                                          corOptions=list(use='p', maxPOutliers=0.1), 
                                                          verbose=opts$verbose) }
allMerge[["pearson"]] <- function(datExpr, colors, eigengenes, opts) { 
                                        mergeCloseModules(datExpr, 
                                                          colors, 
                                                          MEs=eigengenes, 
                                                          cutHeight=opts$pruning$MEDissThres, 
                                                          corOptions=list(use='p'),
                                                          verbose=opts$verbose) }                                                          
get_merged <- allMerge[[tolower(opts$network$cor)]]                                        

allSubTOM <- list()
# determine TOM from expression data (and not from adjacency), for a subset of genes
allSubTOM[["spearman"]] <- function(datExpr, subset, opts) { 
                                        subsetTOM(datExpr, 
                                                  subset, 
                                                  power=opts$network$beta, 
                                                  networkType=opts$network$type, 
                                                  corOptions=list(use='p', method='spearman'),
                                                  verbose=opts$verbose) }
allSubTOM[["bicor"]] <- function(datExpr, subset, opts) { 
                                        subsetTOM(datExpr, 
                                                  subset, 
                                                  power=opts$network$beta, 
                                                  networkType=opts$network$type, 
                                                  corFnc="bicor",
                                                  corOptions="use='p', maxPOutliers=0.1",
                                                  verbose=opts$verbose) }                                                  
allSubTOM[["pearson"]] <- function(datExpr, subset, opts) { 
                                        subsetTOM(datExpr, 
                                                  subset, 
                                                  power=opts$network$beta, 
                                                  networkType=opts$network$type, 
                                                  ccorOptions="use='p'",
                                                  verbose=opts$verbose) }
get_subTOM <- allSubTOM[[tolower(opts$network$cor)]]

get_tom <- function(get_adjacency, datExpr, opts) {
    # Calculate the similarity of columns (genes) in datExpr by calling the function given in corFnc (for correlation 
    # networks) or distFnc (for distance networks), transforms the similarity according to type and raises it 
    # to power, resulting in a weighted network adjacency matrix. 
    adjacency <- get_adjacency(datExpr, opts)
    # Turn adjacency into topological overlap matrix
    # default TOMDenom = "min" Zhang and Horvath (2005)
    TOM <- TOMsimilarity(adjacency, TOMType=opts$network$type)
    TOM
}

get_eigengenes <- function(datExpr, dTOM, opts) {
    # hierarchical clustering
    geneTree <- flashClust(as.dist(dTOM), method=opts$network$clust)
    # assignments to modules
    # this function has many options, pamRespectsDendro = TRUE, default
    # see https://www.rdocumentation.org/packages/dynamicTreeCut/versions/1.63-1/topics/cutreeDynamic
    # use hybrid between hclust and pam
    # cutHeight ("hybrid") defaults to 99% of the range between the 5th percentile and the maximum 
    # of the joining heights on the dendrogram.
    dynamicMods <- cutreeDynamic(dendro=geneTree, 
                                 distM=dTOM, 
                                 method=opts$pruning$method, 
                                 deepSplit=opts$pruning$deepSplit, 
                                 minClusterSize=opts$network$minClusterSize,
                                 verbose=opts$verbose)
    # convert "branches to colours"
    dynamicColors <- labels2colors(dynamicMods)
    # calculate eigengenes
    # https://www.rdocumentation.org/packages/WGCNA/versions/1.61/topics/moduleEigengenes
    MEList <- moduleEigengenes(datExpr, 
                               colors=dynamicColors, 
                               softPower=opts$network$beta, 
                               scale=TRUE) 
    if (!MEList$allOK) { warning("Not all eigengenes have been calculated correctly!") }
    # eigengene expression for each sample
    MEs <- MEList$eigengenes
    list(geneTree=geneTree, dynamicColors=dynamicColors, MEs=MEs)
}

merge_modules <- function(get_merged, datExpr, colors, eigengenes, opts) {
    # merge automatically modules that are too close as measured by the correlation of their eigengenes
    # https://www.rdocumentation.org/packages/WGCNA/versions/1.61/topics/mergeCloseModules
    # default: useAbs = FALSE (correlation), equalizeQuantiles = FALSE (should be fine under ~50 eigengenes), 
    # quantileSummary = "mean" (consensus), iterate = TRUE (merge), relabel = FALSE (output)
    # re-use MEs (need to call moduleEigengenes beforehand, in principle they could be calculated
    # by mergeCloseModules
    MEMerge <- get_merged(datExpr, colors, eigengenes, opts)
    if (!MEMerge$allOK) { warning("Merge close modules has failed!") }
    mergedColors <- MEMerge$colors
    # eigengenes of the new merged modules
    mergedMEs <- MEMerge$newMEs
    list(mergedColors=mergedColors, mergedMEs=mergedMEs)
}

## Call

print(paste("Calling: ", select, sep=""))

if (select=="constructNetwork") { ### 1. Construct network using the original data set
    # check selected options for consistency - but soft threshold is overwritten by config 
    input <- attach(file.path(params$dataDir, params$inputFile, fsep=.Platform$file.sep))
    inputOpts <- input$opts
    if (!(opts$network$type == inputOpts$network$type)) {
        warning(paste("Setting network type to ", inputOpts$network$type, sep=""))
        opts$network$type <- inputOpts$network$type
    }
    if (!(opts$network$cor == tolower(inputOpts$network$cor))) {
        warning(paste("Setting correlation to ", inputOpts$network$cor, sep=""))
        opts$network$cor <- tolower(inputOpts$network$cor)
    }
    if (!(opts$network$clust == inputOpts$network$clust)) {
        warning(paste("Setting clustering method to ", inputOpts$network$clust, sep=""))
        opts$network$clust <- inputOpts$network$clust
    }
    # run
    refTOM <- get_tom(get_adjacency, datExpr, opts)
    dTOM <- 1-refTOM
    out <- get_eigengenes(datExpr, dTOM, opts)
    refGeneTree <- out$geneTree
    dynamicColors <- out$dynamicColors
    MEs <- out$MEs

    # below only for graphical representation 
    # dissimilarity of module eigengenes
    if (tolower(opts$network$cor)=="spearman"){
        dMEs <- 1-cor(MEs, use="p", method="spearman")
    } else if (tolower(opts$network$cor)=="bicor") { 
        dMEs <- 1-bicor(MEs, use="p", maxPOutliers=0.1)
    } else {# default Pearson
        dMEs <- 1-cor(MEs, use="p")
    }
    # cluster module eigengenes
    METree <- flashClust(as.dist(dMEs), method=opts$network$clust)
    filen <- paste(basename, "clustRefdMEs.pdf", sep="")
    filen <- file.path(params$networkDir, filen, fsep=.Platform$file.sep)
    pdf(filen, width=8, height=6, paper='special')
    sizeGrWindow(8, 6)
    plot(METree, main="Clustering of Module Eigengenes", xlab="", sub="")
    # >>> we already know that we want to cut at MEDissThres (hard coded)
    abline(h=MEDissThres, col="red")
    dev.off()
    
    out <- merge_modules(get_merged, datExpr, dynamicColors, MEs, opts)
    mergedColors <- out$mergedColors
    mergedMEs <- out$mergedMEs
    
    # plot results of detection methods
    filen <- paste(basename, "clustColorsRefMEs.pdf", sep="")
    filen <- file.path(params$networkDir, filen, fsep=.Platform$file.sep)
    pdf(filen, width=12, height=9, paper='special')
    sizeGrWindow(12, 9)
    plotDendroAndColors(refGeneTree, cbind(dynamicColors, mergedColors), 
        c("Dynamic Tree Cut", "Merged Modules"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
    dev.off()

    # finally we save everything with numerical labels corresponding to the colors
    # the index corresponds to the gene (order as in datExpr), the value is the 
    # color (or the index in standardColors(), with "grey"==0)
    refModuleColors <- mergedColors
    colorOrder <- c("grey", standardColors())
    refModuleLabels <- match(refModuleColors, colorOrder)-1
    refMEs <- mergedMEs 
    refOpts <- opts
    refOpts$bootstrap <- NULL
    refOpts$modules <- NULL
    filen <- paste(basename, "refNet.RData", sep="")
    filen <- file.path(params$networkDir, filen, fsep=.Platform$file.sep)
    save(refMEs, refModuleLabels, refModuleColors, refTOM, refGeneTree, refOpts, file=filen)
} else if (select=="bootstrap") { ### 2. Perform bootstrapping 
    nBoot <- opts$bootstrap$nBoot
    refPer <- opts$bootstrap$refPer
    # check if output path exist, create if not
    dirBoot <- file.path(params$networkDir, "bootstrap", fsep=.Platform$file.sep)
    opts$networkDir <- params$networkDir
    opts$dirBoot <- dirBoot
    opts$basename <- basename
    dir.create(dirBoot, recursive=TRUE, showWarnings=FALSE)
    # expression data and traits
    filen <- file.path(params$dataDir, params$inputFile, fsep=.Platform$file.sep)
    # attach - ignore opts, this time we consider refNet options for consistency check...
    input <- attach(filen)
    # load reference TOM, etc.
    load(opts$bootstrap$refNet)
    if (!(ncol(datExpr)==ncol(refTOM))) {
        stop("Number of genes in datExpr does not match that of the reference TOM!.\n", call.=FALSE)
    }
    # check options for consistency
    if (!(opts$network$beta == refOpts$network$beta)) {
        warning(paste("Setting soft threshold power to ", refOpts$network$beta, sep=""))
        opts$network$beta <- refOpts$network$beta
    }
    if (!(opts$network$type == refOpts$network$type)) {
        warning(paste("Setting network type to ", refOpts$network$type, sep=""))
        opts$network$type <- refOpts$network$type
    }
    if (!(opts$network$cor == refOpts$network$cor)) {
        warning(paste("Setting correlation to ", refOpts$network$cor, sep=""))
        opts$network$cor <- refOpts$network$cor
        get_adjacency <- allAdj[[tolower(opts$network$cor)]]
        get_merged <- allMerge[[tolower(opts$network$cor)]]
        get_subTOM <- allSubTOM[[tolower(opts$network$cor)]]
    }
    if (!(opts$network$clust == refOpts$network$clust)) {
        warning(paste("Setting clustering method to ", refOpts$network$clust, sep=""))
        opts$network$clust <- refOpts$network$clust
    }
    # Sample sufficiently large number of TOM entries
    nSamples <- 1/(1-refPer)*1000
    subgSamples <- sample(nrow(refTOM)*(nrow(refTOM)-1)/2, size=nSamples)
    # select reference and define quantile
    # use upper triangle (symmetric) c.f. as.dist() which uses lower triangle...
    refTOVals <- vectorizeMatrix(refTOM)[subgSamples] 
    refQuant <- quantile(refTOVals, probs=refPer, type=8) #approx. median-unbiased
    consTOMVec <- vector("list", nBoot)
    # define subgroups for resampling, based on condition or group
    # <<< datTraits must have rownames matching colnames of datExpr <<<<
    subgroups <- lapply(unique(datTraits[[opts$bootstrap$resamplSet]]), function(x){
        return(datTraits[which(datTraits[[opts$bootstrap$resamplSet]] == x), ]) 
        })
    for (iBoot in 1:nBoot) {
        print(paste("Bootstrap: ", iBoot, "th resampling", sep=""))
        datExprSubList <- list()
        for (idx in seq_along(subgroups)) {
            datTraitsSub <- subgroups[[idx]]
            matchedRows <- rownames(datTraitsSub)
            datExprSub <- datExpr[match(matchedRows, rownames(datExpr)),]
            resample <- sample(1:nrow(datExprSub), size=nrow(datExprSub), replace=TRUE)
            datExprSubList[[idx]] <- datExprSub[resample,]
            # ignore rownames when resampling
            # "reorder" samples (rownames) only to make sure the final expression set
            # matches the "original" datExpr (for eigengene characterisation)
            rownames(datExprSubList[[idx]]) <- matchedRows
        }
        datExprResampled <- do.call(rbind, datExprSubList)
        # << we need to reorder actually... *** all this is not quite efficient...
        # and is this actually correct??? i.e. do we achieve the resampling we aim to do?
        datExprResampled <- datExprResampled[match(rownames(datExpr), rownames(datExprResampled)),]
        stopifnot(all(rownames(datExprResampled) == rownames(datExpr)))
        TOMResampled <- get_tom(get_adjacency, datExprResampled, opts)
        TOVals <- vectorizeMatrix(TOMResampled)[subgSamples] 
        # TOM quantiles
        scaleQuant <- quantile(TOVals, probs=refPer, type=8)
        # scaling powers to match reference TOM values
        scalePowers <- log(refQuant)/log(scaleQuant)
        print(paste("Scaling power: ", scalePowers, sep=""))
        TOMScaled <- TOMResampled^scalePowers
        consTOMVec[[iBoot]] <- TOMScaled
    }
    # save all TOMS (for module quantification)
    filen <- paste(basename, "bootScaledTOMS.RData", sep="")
    filen <- file.path(dirBoot, filen, fsep=.Platform$file.sep)
    save(consTOMVec, opts, file=filen)
    # final consensus
    if (tolower(opts$bootstrap$consensus) == "minimum") {
        print("Taking component-wise minimum of all scaled/resampled TOMs...")
        consTOMMat <- do.call(pmin, consTOMVec)
    } else if (tolower(opts$bootstrap$consensus) == "median") {
        print("Calculating median quantile across all scaled/resampled TOMs...")
        consTOMMat <- matrix(NA, nrow=nrow(refTOM), ncol=ncol(refTOM))
        for (iRow in 1:nrow(refTOM) ) {        
            for (jCol in 1:ncol(refTOM)) {
                TOVals <- vector(mode="numeric", length=nBoot)
                for (iBoot in 1:nBoot){
                    TOVals[iBoot] <- consTOMVec[[iBoot]][iRow, jCol]
                }
            consTOMMat[iRow, jCol] <- median(TOVals)
            }
        }
    }
    dConsTOMMat <- 1-consTOMMat
    out <- get_eigengenes(datExpr, dConsTOMMat, opts)
    geneTree <- out$geneTree
    consDynamicColors <- out$dynamicColors
    consMEs <- out$MEs
    out <- merge_modules(get_merged, datExpr, consDynamicColors, consMEs, opts)
    moduleColors <- out$mergedColors
    # match consensus module colors with "reference" module colors
    moduleColors <- matchLabels(moduleColors, refModuleColors, pThreshold=5e-2)[,1]
    # on the consensus, validate module membership
    # note however that TO is calculated anew from the expression data for each color
    # it would make more sense to subset the consensus TOM...
    # each module is validated independently: if average TOM is greater
    # within module than for random genes, we call this a success.
    # We then use a one-proportion Z-test, i.e. test the population proportion
    # use default (alpha=5%, i.e. |z|>=zHalfAlpha)
    # z = (p_obs-p_ref)/sqrt(p_obs(1-p_obs)/n) with zHalfAlpha = qnorm(1-alpha/2)
    if (opts$modules$validate) {
        print("Validating module membership...")
        alpha <- 0.05 # hard coded
        nGenes <- ncol(datExpr)
        uniqueColors <- unique(moduleColors)
        for (thisColor in uniqueColors) {
            thisModGenes <- which(moduleColors==thisColor)
            modSize <- length(thisModGenes) # use matching size sets
            # use median?
            thisAverageTOM <- mean(get_subTOM(datExpr, thisModGenes, opts))
            nSucces <- 0
            for (iter in 1:opts$modules$nIter) {
                randGenes <- sample(1:nGenes, size=modSize, replace=FALSE)
                randAverageTOM <- mean(get_subTOM(datExpr, randGenes, opts))
                if (thisAverageTOM > randAverageTOM) { nSucces <- nSucces+1 }
            }
            pVal <- prop.test(nSucces, 
                              opts$modules$nIter, 
                              opts$modules$refProp, 
                              correct=FALSE, 
                              alternative="greater")$p.value
            if (pVal > alpha) { 
                # reject H1, module is not significantly different than that generated from random genes
                # only flag moduleColors, NS modules have to be dealt with afterwards
                moduleColors[thisModGenes] <- paste("NS", thisColor, sep="")
            }
        }
    }
    # rename and save to disk
    print("Writing to disk....")
    TOM <- consTOMMat
    MEs <- moduleEigengenes(datExpr, 
                            colors=moduleColors, 
                            softPower=opts$network$beta, 
                            scale=TRUE)$eigengenes
    filen <- paste(basename, "consensusNet.RData", sep="")
    filen <- file.path(dirBoot, filen, fsep=.Platform$file.sep)
    save(MEs, moduleColors, TOM, geneTree, opts, file=filen)
    print("Completed successfully.")
}
