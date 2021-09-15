#!/biosw/R/3.5.1/bin/Rscript

# Usage: ./quantify_modules.R [config] <0/1>
# 0: (Default) Consensus network
# 1: Reference network 

library(WGCNA)
options(stringsAsFactors=FALSE)

# Allow multi-threading within WGCNA...? see comments main script!
# allowWGCNAThreads()
enableWGCNAThreads(nThreads = 24)

library(flashClust)
library(dplyr)
library(tibble)

library("RColorBrewer")


## I/O, settings


# get params and options
# ** some options are hard coded (see below)
args <- commandArgs(trailingOnly=TRUE)
useRef <- 0
if (length(args)<1) {
  stop("At least configuration file must be given.\n", call.=FALSE)
} else {
  if (length(args)>1 & as.integer(args[2])<=1) {
    useRef <- as.integer(args[2])
  }
}

params.file <- args[1]
# Read the path/input network, but use options from the results
params <- yaml::read_yaml(params.file)

# output base name
outputBasename <- params$basename
    
# we need to check whether we use the ref or the consensus network
if (useRef == 1) {
    qopts <- params$quantify$refnet
} else {
    qopts <- params$quantify$bootstrap
}
print(paste("Quantify: ", qopts$network, sep=""))

# check if output path exist, create if not
outputPath <- qopts$dirOutput
dir.create(outputPath, recursive=TRUE, showWarnings=FALSE)

# load network, and all data required for analysis
load(file=qopts$network) # MEs, moduleColors, TOM, geneTree, opts

# if we have a reference network...
if (useRef == 1) {
    MEs <- refMEs
    moduleColors <- refModuleColors
    TOM <- refTOM
    geneTree <- refGeneTree
    opts <- refOpts
    rm(refMEs, refModuleLabels, refModuleColors, refTOM, refGeneTree, refOpts)
}

beta <- opts$network$beta
networkType <- opts$network$type
corTyp <- opts$network$cor

# here, we ignore the initial options, and use the final options - refNet or bootstrap
# # datExpr, datTraits
input <- attach(file.path(params$dataDir, params$inputFile, fsep=.Platform$file.sep))

if (!(ncol(datExpr)==length(moduleColors))) {
        stop("Number of genes in datExpr does not match that of the assigned modules!.\n", call.=FALSE)
}

        
## hard coded, params

# use standard pearson here, to correlate MM and GS, the relationship is linear
corFnc <- "cor"
corOptions <- "use = 'p'"

# "hub" genes MM and chosen GS
topQ <- .75
    
# cytoscape
adjThreshold <- 0.1 # adjacency threshold for including edges in the output

# pretty display module-trait correlation heatmap
traitPvalThreshold <- 0.05

# hub gene selection (for module characterisation)
# based on intramodular connectivity z-score threshold 
hubZthreshold <- 2
# after, try to relate these to the highest module 
# membership MM (equivalent to module eigengene-based 
# intramodular connectivity kME).

# Define the correlation (both categorical and numerical).
# The order of factors is as defined in datTraits

# TODO: in this case, we have the first column left, which we need to remove...
datTraits$Run <- NULL

# Design matrix. Since we use also numerical traits, we want to use
# robust correlations (bicor), but this will raise a warning for binary/discrete variables.
# The current behaviour is to revert to standard correlation (Pearson) for these, while
# using bicor for the numerical traits.
robustY <- TRUE # default anyway
# For binary traits, the direction of the correlation is not of immediate interest, we want 
# to see if there is anything significant.

# TODO: hard coded...
# We could have used something like: design <- model.matrix(~0+datTraits$condition+...)
design.condition <- binarizeCategoricalVariable(datTraits$condition,
                                                includePairwise=T,
                                                includeLevelVsAll=T)
design.race <- binarizeCategoricalVariable(datTraits$race,
                                           includePairwise=T,
                                           includeLevelVsAll=F)
design.sex <- binarizeCategoricalVariable(datTraits$sex,
                                          includePairwise=T,
                                          includeLevelVsAll=F)
design <- cbind(datTraits, design.condition, design.race, design.sex)
design$condition <- NULL
design$race <- NULL
design$sex <- NULL


## Functions

prepComma = function(s)
{
  ifelse (s=="", s, paste(",", s));
}

corPvalueStudent = function(cor, nSamples) 
{
  T=sqrt(nSamples-2) * cor/sqrt(1-cor^2)
  2*pt(abs(T),nSamples-2, lower.tail = FALSE)
}

# from WGCNA package, reference
hubGeneSignificance=function(datKME, GS) 
{
  nMEs=dim(as.matrix(datKME))[[2]]
  nGenes= dim(as.matrix(datKME))[[1]]
  if ( dim(GS)[1] !=  nGenes ) 
    stop("Numbers of genes in 'datKME' and 'GS' are not compatible. ")
  Kmax=as.numeric(apply(as.matrix(abs(datKME)),2,max, na.rm = TRUE))
  Kmax[Kmax==0]=1
  datKME=scale(datKME, center=FALSE, scale=Kmax)
  sumKsq=as.numeric(apply(as.matrix(datKME^2) , 2, sum, na.rm = TRUE))
  sumKsq[sumKsq==0]=1
  HGS=as.numeric(apply(I(GS)*datKME, 2, sum,na.rm = TRUE))/ sumKsq
  as.numeric(HGS)
}

## Call

# get module eigengenes and average expression
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# recalculate modules eigengenes with the saved (merged) colors labels, 
# "flag" grey and NS modules, but do not drop
moduleColors[grepl("NS", moduleColors)] <- "grey"
# softPower: only used when the hubgene approximation is necessary because the principal component calculation failed.
calcMEs <- moduleEigengenes(datExpr, colors=moduleColors, softPower=beta, scale=TRUE, excludeGrey=TRUE)
averageExpr <- orderMEs(calcMEs$averageExpr) # average normalized module expression
# reorder eigenvectors such that similar ones (as measured by correlation) are next to each other
# table with rows=samples and columns eigengenes (colors)
orderedMEs <- orderMEs(calcMEs$eigengenes)
# identical(MEs, orderedMEs) TRUE if we ignore grey and NS modules
module_info <- cbind(orderedMEs, averageExpr) %>% rownames_to_column(var="replicate")
filen <- paste(outputBasename, "module-info.csv.gz", sep="")
z <- gzfile(file.path(outputPath, filen, fsep=.Platform$file.sep))
write.csv(module_info, file=z, row.names=FALSE, quote=FALSE)
# variance, not reordered
varExplained <- calcMEs$varExplained
colnames(varExplained) <- colnames(calcMEs$eigengenes)
varExplained <- varExplained[,colnames(orderedMEs)]
filen <- paste(outputBasename, "variance-explained.csv.gz", sep="")
z <- gzfile(file.path(outputPath, filen, fsep=.Platform$file.sep))
write.csv(varExplained, file=z, row.names=FALSE, quote=FALSE)

# calculates Student asymptotic p-values for given correlations
if (corTyp=="spearman"){
    moduleTraitCor <- cor(orderedMEs, 
                          design, 
                          use="p",
                          method="spearman")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
} else if (corTyp=="bicor") {
    bicorAndP <- bicorAndPvalue(orderedMEs, 
                                design, 
                                use="p", 
                                maxPOutliers=0.1, 
                                robustY=robustY) # discrete if FALSE, anyway, binary or discrete, reverts to Pearson...
    moduleTraitCor <- bicorAndP$bicor
    moduleTraitPvalue <- bicorAndP$p
} else {
    moduleTraitCor <- cor(orderedMEs, 
                          design, 
                          use="p")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
}
# no pretty display: correlation heatmap, xLabels <- sapply(strsplit(colnames(moduleTraitCor), designChoice), function(x) x[2])
xLabels <- colnames(moduleTraitCor)
# adjust cell content
# textMatrix = ifelse(moduleTraitPvalue < 0.05, '*', '')
# textMatrix <-  paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
textMatrix <- paste(ifelse(moduleTraitPvalue < traitPvalThreshold, signif(moduleTraitCor, 2), ""), " (", ifelse(moduleTraitPvalue < traitPvalThreshold, signif(moduleTraitPvalue, 1), ""), ")", sep = "")
textMatrix <- replace(textMatrix, textMatrix == " ()", "")

filen <- paste(outputBasename, "heatmap.pdf", sep="")
filen <- file.path(outputPath, filen, fsep=.Platform$file.sep)
pdf(filen, width=10, height=6, paper='special')
sizeGrWindow(10,6)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor,
            xLabels=xLabels,
            yLabels=names(orderedMEs),
            ySymbols=substring(names(orderedMEs), 3),
            yColorLabels=TRUE,
            colors=blueWhiteRed(50, gamma=0.5),
            textMatrix=textMatrix,
            cex.text=0.5,
            zlim=c(-1,1),
            main=paste("Module-Trait Relationships"))
dev.off()

# save the correlation and p-values for reference
colnames(moduleTraitCor) <- paste("cor", colnames(moduleTraitCor), sep=".") 
colnames(moduleTraitPvalue) <- paste("p",colnames(moduleTraitPvalue), sep=".") 
x <- order(c(1:ncol(moduleTraitCor), 1:ncol(moduleTraitPvalue)))
moduleTraitDf <- cbind(moduleTraitCor, moduleTraitPvalue)[,x]
filen <- paste(outputBasename, "module-trait-cor.csv.gz", sep="")
z <- gzfile(file.path(outputPath, filen, fsep=.Platform$file.sep))
write.csv(moduleTraitDf, file=z, row.names=TRUE, quote=FALSE)

# Intramodular analysis: identifying genes with high GS and MM (kME). We quantify associations of individual genes 
# with traits of interest by defining Gene Significance GS as (the absolute value of) the correlation between 
# the gene and the trait. For each module, we also define a quantitative measure of module membership MM as 
# the correlation of the module eigengene and the gene expression profile.

# check signedKME eigengene-based connectivity (signedKME ~ cor(datExpr, orderedMEs))
# cor between expression and module, merging dimension is samples
# same for orderedMEs, which shows the eigengenes for each sample (rows)

# else use corAndPvalue(datExpr, orderedMEs, alternative="greater"), and extract $p and $cor
# default alternative="two.sided", greater ~ positive association, and vice versa

if (corTyp=="spearman"){
    geneModuleMembership <- as.data.frame(cor(datExpr, orderedMEs, use="p", method="spearman"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    geneTraitSignificance <- as.data.frame(cor(datExpr, design, use = "p", method="spearman"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
} else if (corTyp=="bicor") { 
    bicorAndP <- bicorAndPvalue(datExpr, orderedMEs, use="p", maxPOutliers=0.1)
    geneModuleMembership <- bicorAndP$bicor
    MMPvalue <- bicorAndP$p
    bicorAndP <- bicorAndPvalue(datExpr, design, use="p", maxPOutliers=0.1, robustY=robustY)
    geneTraitSignificance <- bicorAndP$bicor
    GSPvalue <- bicorAndP$p
} else {# default Pearson
    geneModuleMembership <- as.data.frame(cor(datExpr, orderedMEs, use="p"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    geneTraitSignificance <- as.data.frame(cor(datExpr, design, use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
}
modNames <- substring(names(orderedMEs), 3)
colnames(geneModuleMembership) <- paste("MM", modNames, sep="")
colnames(MMPvalue) <- paste("p.MM", modNames, sep="")
# no pretty display
#     colnames(geneTraitSignificance) <- paste("GS.", 
#         sapply(strsplit(colnames(geneTraitSignificance), designChoice), function(x) x[2]), sep="")
#     colnames(GSPvalue) <- paste("p.GS.", sapply(strsplit(colnames(GSPvalue), designChoice), function(x) x[2]), sep="")
colnames(geneTraitSignificance) <- paste("GS.", colnames(geneTraitSignificance), sep="")
colnames(GSPvalue) <- paste("p.GS.", colnames(GSPvalue), sep="")

# TODO: hard coded!
# add gene symbol and save to disk    
library(org.Hs.eg.db)

# the unique colors...
selectMod <- unique(moduleColors) # subset using c("selected_color")
#  ... except grey (incl. NS modules)
selectMod <- selectMod[selectMod != "grey"]
# the colors, length(moduleColors)=nGenes, after removing grey we have less so the order is not the same 
colorProbes <- moduleColors[moduleColors != "grey"]
# mask color in selection: assign each color in moduleColors to one found in selectMod
# i.e. grey is NA, so applying is.finite convert to TRUE or FALSE.
# length(inModule)=nGenes
inModule <- is.finite(match(moduleColors, selectMod))
# the genes we want to keep, so now we have colorProbes matching modProbes
modProbes <- colnames(datExpr)[inModule]
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=modProbes,
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
# but this removes duplicate assignments arbitrarily...
# unless we assume the first assignment to be more relevant?!
ens2symbol <- ens2symbol[!duplicated(ens2symbol$ENSEMBL),]
geneInfo <- data.frame(symbol=ens2symbol[match(modProbes, ens2symbol$ENSEMBL), 'SYMBOL'],
                       color=colorProbes)
for (mod in 1:ncol(geneTraitSignificance)){
    oldNames <- colnames(geneInfo)
    geneInfo <- data.frame(geneInfo, 
                           geneTraitSignificance[inModule, mod],
                           GSPvalue[inModule, mod])
    colnames(geneInfo) <- c(oldNames, colnames(geneTraitSignificance)[mod], colnames(GSPvalue)[mod])
}  
rownames(geneInfo) <- rownames(geneTraitSignificance)[inModule]
for (mod in 1:ncol(geneModuleMembership)){
    oldNames <- colnames(geneInfo)
    geneInfo <- data.frame(geneInfo, 
                           geneModuleMembership[inModule, mod],
                           MMPvalue[inModule, mod])
    colnames(geneInfo) <- c(oldNames, colnames(geneModuleMembership)[mod], colnames(MMPvalue)[mod])
}                       
geneInfo <- geneInfo %>% rownames_to_column(var="gene_id")


# connectivity of nodes to other nodes within the same module
# anyway, default is: ignoreColors = if (is.numeric(colors)) 0 else "grey"
# see also softConnectivity (determine connectivity overall, i.e. not intramodular)

# Intramodular connectivity measures how coexpressed a given gene is with respect to genes of the same module. 
# This is a measure equivalent to module membership (kME/MM). 
# default: determines whole network connectivity
# ** Do this for all modules (including grey module)! i.e. dim(iConnect)[1]==nGenes
iConnect <- intramodularConnectivity(TOM,
                                     moduleColors, 
                                     scaleByMax=TRUE)

iConnect_matched <- iConnect[inModule,]         
# add info to gene-info...
geneInfo <- cbind(geneInfo, iConnect_matched[,c(1,2)])
dim.geneInfo <- dim(geneInfo)
geneInfo <- geneInfo[,c(1,2,3,(dim.geneInfo[2]-1):dim.geneInfo[2],4:(dim.geneInfo[2]-2))]
stopifnot(all.equal(dim.geneInfo, dim(geneInfo)))

# Select "hub genes" for each module. It is possible that assigned module and module of highest MM differ (module labels 
# assigned using dynamic tree cut, based on TO measure, module merging, etc.). Genes will have a high MM/kME to 
# their assigned module, and where there is "conflict", such genes can be considered "intermediate" between two modules.
# We rather use intramodular connectivity (k_in), and assess the relationship with MM (and GS) later

# determine intramodular connectivity z-scores per modules 
geneInfo$used <- 0
geneInfo$rank <- NA
for (mod in colnames(geneInfo)[grep("^MM", colnames(geneInfo))]){ # == selectMod
    m.inModule <- geneInfo$color==gsub("^MM", "", mod)
    k <- geneInfo$kWithin[m.inModule]
    k <- scale(k) # (k-mean(k))/sd(k)
    geneInfo$rank[m.inModule] <- rank(-k)
    m.inModule[m.inModule==TRUE] <- k > hubZthreshold
    geneInfo$used[m.inModule] <- 1
}
dim.geneInfo <- dim(geneInfo)
geneInfo <- geneInfo[,c(1,2,3,4,5,(dim.geneInfo[2]-1):dim.geneInfo[2],6:(dim.geneInfo[2]-2))]
stopifnot(all.equal(dim.geneInfo, dim(geneInfo)))

# hub gene selection
# also provide correlation/hub gene significance for every module (MM) against all traits (also use kWithin)

MMs <- colnames(geneInfo)[grepl("^MM",colnames(geneInfo))]
GSs <- colnames(geneInfo)[grepl("^GS",colnames(geneInfo))]

MMList <- list()
cols <- c("cor", "p", "sigMM", "sigKin")
offsets <- seq(1,4*length(GSs),4) # c(1, 5, 9, 13)
nMM <- length(GSs)
nCols <- length(cols)
for (mm in MMs) {
    module <- gsub("MM", "", mm)
    inModule <- is.finite(match(geneInfo$color, module))
    # get MM for genes in given module
    inMM <- abs(geneInfo[inModule,][mm])
    # get k_in (scaled)
    inKin <- geneInfo[inModule,]["kWithin"]
    vec <- numeric(nMM*nCols)
    names_vec <- c()
    for (idx in seq_along(GSs)) {
        gs <- GSs[idx]
        inGS <- abs(geneInfo[inModule,][gs])
        # correlation, use the absolute value, not the sign, i.e. the strengh of association
        corExpr <- parse(text=paste(corFnc, "(inMM, inGS ", prepComma(corOptions), ")"))
        corVal <- signif(eval(corExpr), 2) # truncate brefore calculating significance?
        # Calculates Student asymptotic p-value for given correlations.
        corP <- signif(corPvalueStudent(corVal, sum(is.finite(as.matrix(inMM)) & is.finite(as.matrix(inGS)))), 2)
        # hub gene significance - here not absolute values
        # CHANGED? passed second as vector...
        sigMM <- hubGeneSignificance(geneInfo[inModule,][mm], as.matrix(geneInfo[inModule,][gs])) # redo, but not used in the paper
        sigK <- hubGeneSignificance(inKin, as.matrix(geneInfo[inModule,][gs]))
        jdx <- offsets[idx] + 4 - 1
        vec[offsets[idx]:jdx] <- c(corVal, corP, sigMM, sigK)
        names_vec <- c(names_vec, paste(gs, cols, sep="."))
    }
    names(vec) <- names_vec    
    MMList[[module]] <- vec
}


MMdf_ref <- do.call("rbind", MMList) %>% as.data.frame()

# for each module, get the highest correlation (used for characterising the module)
MMdf <- MMdf_ref %>% dplyr::select(-contains("sig")) %>%
                    dplyr::select(-ends_with(".p"))
# do not bother to filter on the pvalue, the strongest correlations are significant anyway...?
#MMdf_ref$assoc <- colnames(MMdf)[apply(abs(MMdf), 1, which.max)]
MMdf_ref$assoc <- colnames(MMdf)[apply(MMdf, 1, which.max)]
MMdf_ref$assoc <- gsub(".cor", "", MMdf_ref$assoc)
MMdf_ref$sign <- apply(MMdf_ref, 1, function(x) sign(as.numeric(x[paste(x["assoc"], ".sigMM", sep="")])))

# write to disk for reference
filen <- paste(outputBasename, ".kMEs.csv.gz", sep="")
z <- gzfile(file.path(outputPath, filen, fsep=.Platform$file.sep))
write.csv(MMdf_ref, file=z, row.names=TRUE, quote=FALSE)

# now get the top right corner of genes for enrichment
geneInfo$top <- 0
for (module in unique(geneInfo$color)){
    inModule <- is.finite(match(geneInfo$color, module))
    # GS
    col <- MMdf_ref[rownames(MMdf_ref) %in% module,]$assoc
    vals <- as.matrix(abs(geneInfo[inModule,][col])) # for correlation, we used absolute values
    modQuant <- quantile(vals, probs=topQ, type=8) #approx. median-unbiased
    gstop <- abs(geneInfo[col])>modQuant & inModule
    # MM
    col <- paste("MM", module, sep="")
    vals <- as.matrix(abs(geneInfo[inModule,][col])) 
    modQuant <- quantile(vals, probs=topQ, type=8)
    mmtop <- abs(geneInfo[col])>modQuant & inModule
    geneInfo$top[gstop & mmtop] <- 1
}

# write to disk for reference

geneInfo$isHub <- 0
m <- (geneInfo$used == 1) | (geneInfo$top == 1)
geneInfo$isHub[m] <- 1

geneInfo$used <- geneInfo$isHub
geneInfo$top <- NULL
geneInfo$isHub <- NULL
geneInfo <- geneInfo %>% dplyr::rename(isHub=used)


for (mod in colnames(geneInfo)[grep("^MM", colnames(geneInfo))]){ # == selectMod
    m.inModule <- geneInfo$color==gsub("^MM", "", mod)
    k <- geneInfo$isHub[m.inModule]
    msg <- paste("Module: ", mod, ", ", sum(k), " hub genes", sep="")
    print(msg)
}

# ...and write to disk
filen <- paste(outputBasename, "gene-info.csv.gz", sep="")
z <- gzfile(file.path(outputPath, filen, fsep=.Platform$file.sep))
write.csv(geneInfo, file=z, row.names=FALSE, quote=FALSE)


# export to Cytoscape, use TOM 
# DO NOT recalculate topological overlap using intramodularConnectivity.fromExpr!
# i.e. concepts of connectivity/network requires the right TOM, but if we use modules/eigengenes, these are already assigned
# based on the TOM, so their expression, etc. as well as MM (correlation between modules and genes), 
# GS, etc. can be determined from assigned modules and datExpr. 

# check if output path exist, create if not
# all files are overwritten!
print("Writing cytoscape files. Existing files are overwritten!")
filen <- paste(outputBasename, "cytoscape", sep="")
cytoDir <- file.path(outputPath, filen, fsep=.Platform$file.sep)
dir.create(cytoDir, recursive=TRUE, showWarnings=FALSE)

probes <- colnames(datExpr)
for (modCol in selectMod) {  # this excludes grey
    # redefine, using only this module, length(inModule)==nGenes==length(probes)
    inModule <- is.finite(match(moduleColors, modCol))
    # subset genes in this module only (all, not only hub genes)
    modProbes <- probes[inModule]
    # annotation
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key=modProbes,
                                        columns="SYMBOL",
                                        keytype="ENSEMBL")
    # but this removes duplicate assignments arbitrarily...
    ens2symbol <- ens2symbol[!duplicated(ens2symbol$ENSEMBL),]
    # select the corresponding Topological Overlap
    modTOM <- TOM[inModule, inModule]
    dimnames(modTOM) <- list(modProbes, modProbes)
    
    # define node attributes: intramodular connectivity, isHub
    isHub <- geneInfo$isHub[geneInfo$gene_id %in% modProbes]
    nodeAttr <- data.frame(cbind(iConnect$kWithin[inModule], isHub)) %>% dplyr::rename(kIN=V1)
    
    # export the network into edge and node list files for Cytoscape
    eFilen <- paste(outputBasename, "cytoEdges-", paste(modCol, collapse="-"), ".txt", sep="")
    eFilen <- file.path(cytoDir, eFilen, fsep=.Platform$file.sep)
    nFilen <- paste(outputBasename, "cytoNodes-", paste(modCol, collapse="-"), ".txt", sep="")
    nFilen <- file.path(cytoDir, nFilen, fsep=.Platform$file.sep)
    cyt <- exportNetworkToCytoscape(modTOM,
                                    edgeFile=eFilen,
                                    nodeFile=nFilen,
                                    weighted=TRUE,
                                    threshold=adjThreshold, # adjacency threshold for including edges in the output
                                    nodeNames=modProbes,
                                    altNodeNames=ens2symbol$SYMBOL,
                                    nodeAttr=nodeAttr) 
}


# -----------------------------------------------------------
print("Done!")
