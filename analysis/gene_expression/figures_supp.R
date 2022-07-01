#! /usr/bin/env Rscript

# this uses a slightly different set-up than the 
# MAGNet_GeneTonics.R script...

library(topGO)
library(CellPlot)
library(ggplot2)

mapping <- "org.Hs.eg.db"
ID <- "Ensembl"
go.class <- "weight01Count"
# use by default Fisher Test
nodesize <- 10
topNodes <- 250

alpha.set <- 0.05
pvalCutOff <- 0.05 # everywhere 

## 

contrast <- 'DCMvsHCM'
# ontology <- 'BP'
ontology <- 'CC'
# ontology <- 'MF'

## I/O

dirloc <- '/prj/MAGE/analysis/genetonic/DGE/results'
results <- readRDS(file.path(dirloc, 'MAGNet_DESeqResults.rds'))

.x <- results[[contrast]]
.y <- contrast

de <- rownames(GeneTonic::deseqresult2df(.x, FDR=alpha.set))
bg <- rownames(GeneTonic::deseqresult2df(.x))

universe <- bg %in% bg
selection <- bg %in% de
relevant.genes <- factor(as.integer(selection[universe]))
names(relevant.genes) <- bg

topGO.data <- new("topGOdata", 
                  ontology=ontology, 
                  allGenes=relevant.genes, 
                  mapping=mapping, 
                  annotationFun=annFUN.org, 
                  nodeSize=nodesize,
                  ID=ID)
# test statistic
test.stat <- new(go.class, 
                 testStatistic=GOFisherTest, 
                 name="Fisher")
# run Fisher test
sig.groups <- getSigGroups(topGO.data, test.stat)
# output results
fisher.results <- GenTable(topGO.data,
                            pvalCutOff=sig.groups,
                            topNodes=topNodes) #length(topGO.data@graph@nodes))

fisher.results$pvalCutOff <- as.numeric(stringr::str_replace_all(fisher.results$pvalCutOff, "[^0-9e\\-\\.]*", ""))
fisher.results$LogEnriched <- log2(fisher.results$Significant / fisher.results$Expected)

ga <- genesInTerm(topGO.data) # GenesAnnotated | list of genes per go-terms
ga <- ga[fisher.results$GO.ID] # eliminate missing terms
names(ga) <- NULL
fisher.results$GenesAnnotated <- ga
xs <- .x[,c("padj", "log2FoldChange")] # significant stats subset
xs <- subset(xs, padj < pvalCutOff)
fisher.results$GenesSignificant <- lapply(fisher.results$GenesAnnotated, intersect, rownames(xs)) # extract genes
ei.rows <- mclapply(fisher.results$GenesSignificant, function (y) {
if (length(y)) as.list(xs[y,,drop=FALSE])
else as.list(rep(NA_real_, length(xs)))
}, mc.cores = 10)
ei <- mclapply(names(xs), function(z) {
lapply(ei.rows, "[[", z)
}, mc.cores = 10)
ei <- structure(ei, names = names(xs), row.names = seq(nrow(fisher.results)), class = "data.frame")
row.names(ei) <- NULL
fisher.results <- data.frame(fisher.results, ei, stringsAsFactors = FALSE, check.names = FALSE)

x <- head(fisher.results, 15)
x$pvalCutOff <- -log10(x$pvalCutOff)
    
    
filen <- file.path(dirloc, paste(contrast, '_cellplot_pval', ontology, '.pdf', sep=''))  
pdf(filen, width=12, height=8)

cell.plot(x=setNames(x$pvalCutOff, x$Term), 
            cells=x$log2FoldChange, 
            main=contrast, 
            x.mar=c(.5, 0), 
            key.n=7, 
            y.mar=c(.1, 0), 
            cex=1.6, 
            cell.outer=2, 
            bar.scale=.7, 
            space=.2,
            xlab = "-log10(Pval) enrichment")
            
sym.plot(x=setNames(x$pvalCutOff, x$Term), 
            cells=x$log2FoldChange, 
            x.annotated = x$Annotated, 
            main=contrast,
            x.mar=c(.45, 0), 
            key.n=7, 
            cex=1.6, 
            axis.cex=.8, 
            group.cex=.7,
            xlab = "-log10(Pval) enrichment") 

dev.off()
