#!/biosw/R/3.5.1/bin/Rscript

# Filter based on connectivity - requires soft threshold
# Recompute soft threshold and write final output

# Usage: ./pick_soft_threshold.R config
# config: configuration file


library(WGCNA)
options(stringsAsFactors=FALSE)

# Allow multi-threading within WGCNA
allowWGCNAThreads()

library(dplyr)
library(stringr)


## Data 

# get params and options
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) { stop("At least configuration file must be given.\n", call.=FALSE) } 

# values are NOT checked for consistency, type, etc. 
params <- yaml::read_yaml(args[1])
opts <- params$opts
verbose <- 5
opts$verbose <- verbose

# Choose a set of potential soft-thresholding powers (beta)
betas <- c(c(1:20), seq(from=22, to=34, by=2))

# ---------------------------------------------------------

## Functions

get_gsg <- function(datExpr) {
    gsg <- goodSamplesGenes(datExpr, verbose=3)
    if (!gsg$allOK) {
        if (sum(!gsg$goodGenes)>0) 
            print(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse=", ")))
        if (sum(!gsg$goodSamples)>0) 
            # TO DO: here we need to also remove corresponding samples from the metadata (phenotypic/trait data).
            # This is not the case now!
            print(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
        # Remove
        datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    }
    datExpr
}

## Call

## restrict to highly connected genes

# correlation, defaults to bicor
opts$network$cor <- tolower(opts$network$cor)
if (!opts$network$cor %in% c("pearson", "spearman", "bicor")) {
    opts$network$cor <- "bicor"
}

# define functions based on chosen correlation
allAdj <- list()
allAdj[["spearman"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                            power=opts$network$betaPrep, 
                                                            type=opts$network$type, 
                                                            corOptions=list(use='p', method='spearman')) 
                                                }
allAdj[["bicor"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                         power=opts$network$betaPrep, 
                                                         type=opts$network$type, 
                                                         corFnc="bicor", 
                                                         corOptions="use='p', maxPOutliers=0.1")
                                            }
allAdj[["pearson"]] <- function(datExpr, opts) { adjacency(datExpr, 
                                                           power=opts$network$betaPrep, 
                                                           type=opts$network$type, 
                                                           corOptions=list(use='p')) 
                                            }                                                


# prepped data
prepped <- attach(file.path(params$dataDir, params$preppedFile, fsep=.Platform$file.sep))
preppedOpts <- prepped$opts

# check selected options for consistency
if (!(opts$network$type == preppedOpts$type)) {
    warning(paste("Setting network type to ", preppedOpts$type, sep=""))
    opts$network$type <- preppedOpts$type
}
if (!(opts$network$cor == tolower(preppedOpts$cor))) {
    warning(paste("Setting correlation to ", preppedOpts$cor, sep=""))
    opts$network$cor <- tolower(preppedOpts$cor)
}
if (!(opts$network$RsquaredCut == preppedOpts$RsquaredCut)) {
    warning(paste("Setting RsquaredCut to ", preppedOpts$RsquaredCut, sep=""))
    opts$network$RsquaredCut <- preppedOpts$RsquaredCut
}

# threshold based on connectivity
get_adjacency <- allAdj[[opts$network$cor]]
adjacency <- get_adjacency(datExpr, opts)
diag(adjacency) <- 0
connectivity <- apply(adjacency, 1, sum)
# scaling to [0,1]
connectivity <- connectivity/max(connectivity)
inNetwork <- connectivity >= quantile(connectivity, probs=opts$connectQuant, type=8)
datExpr <- datExpr[,inNetwork]
datExpr <- get_gsg(datExpr) # see function for datTraits, only datExpr is returned

## estimate soft threshold using final dataset
# default: 
# dataIsExpr = TRUE
if (opts$network$cor=="spearman"){
    sft <- pickSoftThreshold(datExpr, 
                             RsquaredCut=opts$network$RsquaredCut,
                             powerVector=betas, 
                             corOptions=list(use="p", method="spearman"), 
                             networkType=opts$network$type, 
                             blockSize=opts$blockSize,
                             verbose=verbose)
} else if (opts$network$cor=="pearson") { 
    sft <- pickSoftThreshold(datExpr, 
                             RsquaredCut=opts$network$RsquaredCut,
                             powerVector=betas,
                             corOptions=list(use='p'),
                             networkType=opts$network$type,
                             blockSize=opts$blockSize,
                             verbose=verbose)
} else { # default bicor
    sft <- pickSoftThreshold(datExpr, 
                             RsquaredCut=opts$network$RsquaredCut,
                             powerVector=betas, 
                             corFnc="bicor",
                             corOptions=list(use="p", maxPOutliers=0.1), 
                             networkType=opts$network$type, 
                             blockSize=opts$blockSize,
                             verbose=verbose)
}

# plot the results
# scale-free topology fit index should reach values above 0.8 for reasonable powers (<15 for unsigned 
# or hybrid, <30 for signed networks); the mean connectivity remains relatively high
filen <- paste(params$basename, ".input.betas.pdf", sep="")
filen <- file.path(params$dataDir, filen, fsep=.Platform$file.sep)
pdf(filen, width=9, height=5, paper='special')
sizeGrWindow(9, 5)
par(mfrow=c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n",
    main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=betas,cex=cex1,col="blue")
# this line corresponds to using an R^2 cut-off of RsquaredCut
abline(h=opts$network$RsquaredCut,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=betas, cex=cex1, col="blue")
dev.off()

# save data
filen <- paste(params$basename, ".input.RData", sep="")
filen <- file.path(params$dataDir, filen, fsep=.Platform$file.sep)
save(opts, datExpr, datTraits, file=filen)
