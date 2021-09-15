#!/biosw/R/4.0.5_deb10/bin/Rscript

library(edgeR)
library(topGO)
library(CellPlot)
library(Glimma)

library(gplots)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)

library(dplyr)
library(purrr)
library(tibble)

library(openxlsx)


# ---------------------------------------------------------

# topGO
mapping <- "org.Hs.eg.db"
ID <- "Ensembl"
# "classicCount" or "weight01Count" extension class algorithms dealing with the GO graph structure
go.class <- "weight01Count"
# use by default Fisher Test
nodesize <- 10
topNodes <- 250

pvalCutOff <- 0.05 # everywhere 
FCCutOff <- 1.2

nTopGenes <- 100 # heatmap

# figures
color_palette <- rev(brewer.pal(11,"RdBu"))
ann_colors = list(group=c("NFD"="gray100", "DCM"="gray10", "HCM"="gray65"))

# ---------------------------------------------------------

## Functions

call_figures <- function(fit, contrast) {

    tt <- topTags(fit, n=nrow(fit))$table
       
    # background - final DGE set
    background <- rownames(tt)
    universe <- background %in% background
    selection <- background %in% rownames(tt[tt$FDR<pvalCutOff,])
    relevant.genes <- factor(as.integer(selection[universe]))
    names(relevant.genes) <- background

    topGO.data <- new("topGOdata", 
                      ontology='BP', # HARD CODED - ONLY BP
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
    xs <- tt[,c("FDR", "logFC")] # significant stats subset
    xs <- subset(xs, FDR < pvalCutOff)
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
    
    # save table
    wb <- createWorkbook()
    addWorksheet(wb, sheetName=paste("TopGO ", contrast, sep=""))
    writeDataTable(wb, sheet=1, x=fisher.results, rowNames=FALSE)
    saveWorkbook(wb, file.path(loc, paste(contrast, "_topGOBP.xlsx", sep=""), fsep=.Platform$file.sep), overwrite=TRUE)
    
    # CellPlot
    x <- head(fisher.results, 15)
    filen <- file.path(loc, paste(contrast, "_cellplot.pdf", sep=""), fsep=.Platform$file.sep)
    pdf(filen, width=12, height=8)
    cell.plot(x=setNames(x$LogEnrich, x$Term), 
              cells=x$logFC, 
              main=contrast, 
              x.mar=c(.5, 0), 
              key.n=7, 
              y.mar=c(.1, 0), 
              cex=1.6, 
              cell.outer=2, 
              bar.scale=.7, 
              space=.2)
    sym.plot(x=setNames(x$LogEnrich, x$Term), 
             cells=x$logFC, 
             x.annotated = x$Annotated, 
             main=contrast,
             x.mar=c(.45, 0), 
             key.n=7, 
             cex=1.6, 
             axis.cex=.8, 
             group.cex=.7) 
    dev.off()
    
    # heatmap
    topGenes <- tt[1:nTopGenes,]
    idx <- match(rownames(topGenes), rownames(lcpm))
    pheatmap(lcpm[idx,],
             scale='row',
             annotation_col=fit$samples[,"group",drop=FALSE],
             show_colnames=FALSE,
             annotation_colors=ann_colors,
             labels_row=annotated[idx,]$hgnc_symbol,
             color=color_palette,
             main=paste("Top", nTopGenes, " DE genes ", contrast, sep=""),
             fontsize_row=3, margin=c(8,6), lhei=c(2,10),
             filename=file.path(loc, paste(contrast, "_pheatmap.pdf", sep=""), fsep=.Platform$file.sep))
             
    
    # volcano - adapted from Tobias' script
    volcano <- tt
    volcano$FDR_log10  <- -log10(volcano$FDR)
    volcano$FDR_log10[which(!is.finite(volcano$FDR_log10))] <- 0
    volcano$FDR_log10[which(volcano$FDR_log10==0)] <- max(volcano$FDR_log10)
    idx <- match(rownames(volcano), annotated$ID)
    volcano$Gene <- annotated[idx,]$hgnc_symbo
    volcano <- volcano %>%
    mutate(color = ifelse(volcano$logFC > FCCutOff & volcano$FDR_log10 > -log10(pvalCutOff),
                          yes = "up_gene",
                          no = ifelse(volcano$logFC < -FCCutOff & volcano$FDR_log10 > -log10(pvalCutOff),
                                      yes = "down_gene",
                                      no = "none")))

    num_genes_down <- nrow(volcano[volcano$color == "down_gene",])
    num_genes_up <- nrow(volcano[volcano$color == "up_gene",])
    if (num_genes_down < 10) {
        num_genes_down <- nrow(volcano[volcano$color == "down_gene",])
        top_labelled_down <- volcano[volcano$color == "down_gene",]
    } else {
        num_genes_down <- 10
        top_labelled_down <- top_n(volcano, n = -num_genes_down, wt = logFC)
    }
    if (num_genes_up < 10) {
        num_genes_up <- nrow(volcano[volcano$color == "up_gene",])
        top_labelled <- volcano[volcano$color == "up_gene",]
    } else {
        num_genes_up <- 10
        top_labelled <- top_n(volcano, n = num_genes_up, wt = logFC)
    }

    colored <- ggplot(volcano, aes(x = logFC, y = FDR_log10)) +
        geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
        theme_bw(base_size = 16) + 
        theme(legend.position = "none") + 
        ggtitle(label = contrast) + 
        xlab(expression(log[2]("FC"))) + 
        ylab(expression(-log[10]("FDR"))) + 
        geom_vline(xintercept = 0, colour = "black") + 
        geom_hline(yintercept = 1.3, colour = "#000000", linetype = "dashed") + 
        scale_color_manual(values = c("up_gene" = rev(color_palette)[2],
                                    "down_gene" = color_palette[2],
                                    "none" = "#636363")) + 
        scale_y_continuous() +
        geom_label_repel(data = top_labelled,
                        mapping = aes(label = Gene),
                        size = 3,
                        fontface = 'bold',
                        color = 'black',
                        nudge_x = 3,
                        force = 1,
                        max.iter = 10000,
                        box.padding = unit(0.25, "lines"),
                        point.padding = unit(0.25, "lines")) +
        geom_label_repel(data = top_labelled_down,
                        mapping = aes(label = Gene),
                        size = 3,
                        fontface = 'bold',
                        color = 'black',
                        force = 1,
                        max.iter = 10000,
                        nudge_x = -3,
                        box.padding = unit(0.25, "lines"),
                        point.padding = unit(0.25, "lines")
        )
    ggsave(file.path(loc, paste(contrast, "_volcano.pdf", sep=""), fsep=.Platform$file.sep))
}


call_glimma <- function(fit, contrast) {

    dt <- decideTests(fit, adjust.method="BH", p.value=pvalCutOff, lfc=FCCutOff)
    glMDPlot(fit, 
             counts=counts,
             groups=fit$samples$group,
             status=dt, 
             transform=TRUE,
             anno=annotated,
             xlab="logMeanExpr",
             ylab="log2FoldChange",
             side.ylab="NormalizedCount",
             path=loc, 
             folder=paste("GlimmaMD_", contrast, sep=""), 
             launch=FALSE)
}


## Call

loc <- "/prj/MAGE/analysis/dge"

objFile <- file.path(loc, "MAGEbaseline_DGEobj.rds", fsep=.Platform$file.sep)
counts <- readRDS(objFile)$counts
lcpm <- cpm(counts, log=TRUE)

fitFile <- file.path(loc, "MAGEbaseline_DGEfit.rds", fsep=.Platform$file.sep)
fits <- readRDS(fitFile)

annotatedFile <- file.path(loc, "MAGEbaseline_DGE.xlsx", fsep=.Platform$file.sep)
annotated <- read.xlsx(annotatedFile, sheet="DCMvsNFD", cols=c(1,2,8,9))
annotated <- annotated[!duplicated(annotated$ID),]
annotated <- annotated[match(rownames(lcpm), annotated$ID),]


map2(fits, names(fits), function(x, y) {
    call_figures(x,y) 
})


map2(fits, names(fits), function(x, y) {
    call_glimma(x,y)
})
