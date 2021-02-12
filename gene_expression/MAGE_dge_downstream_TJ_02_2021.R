#!/usr/bin/env Rscript
### imports

library(edgeR)
library(pheatmap)
library(ca)
library(VennDiagram)
library(lattice)
library(biomaRt)
library(gplots)
library(pathview)
library(gridExtra)
library(topGO)
library(GO.db)
library(grid)
library(gridExtra)
library(parallel)
library(openxlsx)
library(RColorBrewer)
library(viridis)
library(CellPlot)
library(tools)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(Glimma)
library(tm)
library(SnowballC)
library(wordcloud)
library(stringr)

get_data <- function(x, go_obj) {
  scores <- as.numeric(unlist(scoresInTerm(go_obj, x)))
  genes <- genesInTerm(go_obj, whichGO = x)
  sgenes <- sigGenes(go_obj)
  test <- genes[x][[1]]

  if (length(test) > 0) {
    gene_df <- data.frame(test, scores)
    colnames(gene_df) <- c("GeneID", "log2FC")
    found_genes <- (genes[x][[1]][genes[x][[1]] %in% sgenes])
    gene_df <- gene_df[gene_df$log2FC != 0,]

    if (length(found_genes) > 0) {
      geneAnnotation <- id_conversion
      colnames(geneAnnotation) <- c("SYMBOL", "GENEID", "DESCRIPTION")

      df <- data.frame(merge(
        gene_df,
        geneAnnotation,
        by.x = "GeneID",
        by.y = "GENEID",
        all.x = TRUE
      ))

      return(df)
    }
  }
}

run_go_heatmap <- function(go_data, go_table, RPKM, category, direction, focus) {
  head(go_data)

  terms <- go_data$Term
  id <- go_data$GO.ID
  log2 <- go_data$log2
  annotated <- go_data$Annotated

  go_df <- data.frame(id, log2, terms, annotated, stringsAsFactors = FALSE)

  go_df <- go_df[order(-log2),]

  top_terms <- head(go_df$id, 30)
  desc <- head(go_df$terms, 30)

  mylist <- top_terms

  #print(unlist(mylist))

  genes <- lapply(unlist(mylist), function(x) {
    get_data(x, go_table)
  })

  counter <- 1
  for (set in genes) {
    #print(head(set))

    if (!is.null(set)) {
      name <- mylist[counter]
      description <- desc[counter]

      description <- gsub("/", "__", description)

      cpm$GeneID <- rownames(cpm)
      tmp <- merge(cpm, set, by = ("GeneID"))

      tmp[which(is.na(tmp$SYMBOL)), c("SYMBOL")] <- tmp[which(is.na(tmp$SYMBOL)), c("GeneID")]
      tmp$GeneID <- NULL
      tmp$ENTREZID <- NULL
      tmp$LogFC <- NULL

      rownames(tmp) <- NULL
      tmp <- unique(tmp)
      rownames(tmp) <- make.names(tmp$SYMBOL, unique = TRUE)
      gene_names <- data.frame(tmp$SYMBOL)
      colnames(gene_names) <- c("gene")
      #print(head(tmp))
      tmp <- tmp[, labels]
      m <- as.matrix(tmp)


      m <- log(m + 1)
      m[!is.finite(m)] <- NA
      m <- na.omit(m)
      #print(head(m))


      desc_file <- gsub("\\.", "_", description)
      desc_file <- gsub(" ", "_", desc_file)
      desc_file <- gsub(",", "_", desc_file)
      desc_file <- gsub(":", "_", desc_file)

      desc_file <- strtrim(desc_file, 50)



      number <- formatC(counter, width = 3, format = "d", flag = "0")

      if (nrow(m) > 5) {
        pheatmap(
          mat = m,
          main = paste("GO ", category, " ", direction, " in ", focus, "\n ",
                       description, " (", name, ")\nbased on log(RPKM+1) values", sep = ""),
          color = inferno(2000),
          show_colnames = TRUE,
          show_rownames = TRUE,
          filename = paste(baseDir, "GO_heatmaps/", contrast, "_GO_", category, "_", direction,
                           "_", number, "_", desc_file, "_", name, ".pdf", sep = ""),
          fontsize = 8,
          treeheight_row = 0,
          treeheight_col = 40,
          #cellwidth         = 40,
          #cellheight        = 50,
          fontsize_row = 8,
          scale = "row",
          gaps_col = c(3),
          border_color = "black",
          cluster_rows = T,
          cluster_cols = F,
          display_numbers = F,
          #labels_row = gene_names$gene,
          legend.cex = 1,
          #width= 20,
          #height= 15
        )
      }
      counter <- counter + 1
    }
  }
}

wordcloud_function <- function(data = data,
                               removewords = c("process", "activity", "positive", "negative", "response", "regulation"),
                               min.freq = 4,
                               max.words=Inf,
                               random.order=TRUE){
  input <- Corpus(VectorSource(data))

  input <- tm_map(input, content_transformer(tolower))
  input <- tm_map(input, content_transformer(removePunctuation))
  input <- tm_map(input, removeNumbers)
  input <- tm_map(input, stripWhitespace)

  toSpace <- content_transformer(function(x , pattern ) gsub(pattern, " ", x))
  input <- tm_map(input, toSpace, "/")
  input <- tm_map(input, toSpace, "@")
  input <- tm_map(input, toSpace, "\\|")

  input <- tm_map(input, function(x) removeWords(x, stopwords("english")))

  # specify your stopwords as a character vector
  input <- tm_map(input, removeWords, removewords)

  tdm <- TermDocumentMatrix(input)
  m <- as.matrix(tdm)
  v <- sort(rowSums(m),decreasing = TRUE)
  d <- data.frame(word = names(v),freq = v)

  set.seed(1234)
  wordcloud(words = d$word, freq = d$freq, min.freq = min.freq, scale = c(8,.2),
            max.words = max.words, random.order = random.order, rot.per = 0.15,
            colors = brewer.pal(8, "Dark2"))
}


# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

contrast <- args[1]
global_pval <- as.numeric(args[2])
global_fdr <- as.numeric(args[3])
log_threshold <- as.numeric(args[4])
threshold <- as.numeric(args[5])
focus <- args[6]
## fix "not find zip" error from openxlsx
Sys.setenv(R_ZIPCMD = "/usr/bin/zip")
require(openxlsx)
# visible R code
options(echo = TRUE)

gene_counts <- read.csv("new_data_02_2021/gene_count_matrix.csv",row.names="gene_id",as.is=T)
colnames(gene_counts)<-gsub("_stringtieRef","",colnames(gene_counts))

library(stringr)
library(edgeR)

meta<-readRDS("new_data_02_2021/meta.rds")

dedup=read.table("new_data_02_2021/mage_dedup_metrics_final.txt",as.is=T,header=F)
dedup=dedup[,c(1,3)]
dedup[,1]=gsub(".metrics.txt-Unknown","",dedup[,1])
colnames(dedup)=c("Run","DuplicationRate")
meta=merge(meta,dedup,by.x=1,by.y=1)
rownames(meta)<-meta[,1]

ins=intersect(rownames(meta),colnames(gene_counts))

meta <- meta[ins, ]
gene_counts <- gene_counts[, ins]

#library(edgeR)
DGEObj <- DGEList(counts=gene_counts)

design <- model.matrix(~ meta$etiology + meta$race + meta$sex + meta$Age + meta$DuplicationRate)
design = design[,-c(5,7,9)]

keep <- filterByExpr(DGEObj, design=design[, 2:4])
DGEObj <- DGEObj[keep, , keep.lib.sizes=FALSE]


DGEObj <- calcNormFactors(DGEObj)


DGEObj <- estimateGLMCommonDisp(DGEObj, design, verbose=TRUE)
DGEObj <- estimateGLMTagwiseDisp(DGEObj, design)

#Todo save edgeR y object

fit <- glmQLFit(DGEObj, design, robust=TRUE)

                                        #Collect different contrast comparisons.

if (contrast == "DCM"){

  qlf <- glmQLFTest(fit,contrast=c(0,1,0,-1,0,0,0,0)) #DCM vs Normal

} else if (contrast == "HCM"){

  qlf <- glmQLFTest(fit,contrast=c(0,0,1,-1,0,0,0,0)) #HCM vs Normal

} else if (contrast == "DCM_HCM"){

  qlf <- glmQLFTest(fit,contrast=c(0,1,-1,0,0,0,0,0)) #DCM vs HCM

}

counts <- DGEObj$counts

annotation <- "hsapiens_gene_ensembl"

ensembl <- useMart("ensembl", dataset = annotation, host = "useast.ensembl.org")
# set base output directory
baseDir <- paste("./",contrast,"_edgeR__p", global_pval, "__FDR_", global_fdr, "__thresh_", threshold, "/", sep = "")
# baseDir <- paste("./edgeR_p", global_pval, "__FDR_", global_fdr, "__thresh_", threshold, "/", sep = "")
dir.create(baseDir, showWarnings = TRUE, recursive = TRUE, mode = "0755")
dir.create(paste(baseDir, "/GO_heatmaps/", sep = ""), showWarnings = TRUE, recursive = TRUE, mode = "0755")
dir.create(paste(baseDir, "/GO_gene_lists/", sep = ""), showWarnings = TRUE, recursive = TRUE, mode = "0755")


# RPKM <- rpkm(DGEObj, DGEObj$genes$Length, prior.count = 0)
cpm <- cpm(DGEObj)
# get list of gene names with rpkm > 1 in all samples

go_background_num <- nrow(cpm[apply(cpm[, -1], MARGIN = 1, function(x) all(x > 1)),])
go_background_blank <- as.vector(rep(0, go_background_num))

names(go_background_blank) <- row.names(cpm[apply(cpm[, -1], MARGIN = 1, function(x) all(x > 1)),])

message(paste("background gene list contains: ", go_background_num, " genes", sep = ""))

id_conversion <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), values = names(go_background_blank), mart = ensembl)

head(go_background_blank)

write.table(go_background_blank, file = paste(output = baseDir, contrast, "_background_genes.csv", sep = ""), row.names = TRUE, na = "", col.names = TRUE, sep = "\t", quote = FALSE)

write.table(cpm, file = paste(output = baseDir, contrast, "_cpm.csv", sep = ""), row.names = TRUE, na = "", col.names = TRUE, sep = "\t", quote = FALSE)

cpm <- data.frame(cpm)
head(DGEObj$genes)

dim(DGEObj)

write.table(cpm(DGEObj), file = paste(output = baseDir, contrast, "_cpm.csv", sep = ""), row.names = TRUE, na = "", col.names = TRUE, sep = "\t", quote = FALSE)

# MDS plot
pdf(paste(baseDir, "plot_mds.pdf", sep = ""), title = "MA plot")
col.status <- c("blue", "red", "green", "cyan", "black")[DGEObj$samples$group]
plotMDS(DGEObj, col = col.status)
dev.off()


isDE <- as.logical(decideTestsDGE(qlf, p.value = global_pval))


DEnames <- rownames(DGEObj)[isDE]
summary(isDE)
# main data table of edgeR results

save.image(file = paste(contrast, "_", global_pval, ".RData", sep = ""))
q()

genemap3 <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "wikigene_description",  "gene_biotype"), values = rownames(qlf$table), mart = ensembl)
edgerTable <- topTags(qlf, n = nrow(DGEObj))$table

edgerTable.idx <- match(rownames(qlf$table), genemap3$ensembl_gene_id)
qlf$table$gene <- genemap3$external_gene_name[edgerTable.idx]
qlf$genes$desc <- genemap3$wikigene_description[edgerTable.idx]
qlf$genes$go <- genemap3$gene_biotype[edgerTable.idx]

print(head(qlf))
print(head(qlf$genes))

colnames(qlf$genes) <- c("GeneID", "Description", "Type")

head(qlf$table)

gene_names <- data.frame(qlf$table$gene)
colnames(gene_names) <- c("Symbol")

#glMDPlot(   qlf,
#            anno = gene_names,
#            status = decideTestsDGE(qlf, p.value = global_pval),
#            counts = DGEObj,
#            groups = groups,
#            transform = TRUE,
#            launch = FALSE,
#            folder= paste(contrast,"_glimma",sep=""),
#            path = baseDir,
#            html=paste(contrast,"_MD",sep=""),
#            side.main = "Symbol",
#            main = "Gene expression"
#)




DGEObj <- cpm(DGEObj, log = TRUE, prior.count = 1)


selected_cpm_rows <- DGEObj[edgerTable$GeneID[edgerTable$FDR < global_fdr & abs(edgerTable$logFC) > log_threshold],]

selected_logFC_rows <- edgerTable$logFC[edgerTable$FDR < global_fdr & abs(edgerTable$logFC) > log_threshold]
selected_logFC_rows <- sort(selected_logFC_rows)

nrow(selected_cpm_rows)
# colnames(selected_cpm_rows) <- c("ATF6 1","ATF6 2","ATF6 3","ATF6 4","Con 1","Con 2","Con 3")
# plots distribution of log fold changes from logFC vector
shade <- function(x, y, col, n = 500, xlab = 'x', ylab = 'y', ...) {
  # x, y: the x and y coordinates
  # col: a vector of colours (hex, numeric, character), or a colorRampPalette
  # n: the vertical resolution of the gradient
  # ...: further args to plot()
  plot(x, y, type = 'n', las = 1, xlab = xlab, ylab = ylab, ...)
  e <- par('usr')
  height <- diff(e[3:4]) / (n - 1)
  y_up <- seq(0, e[4], height)
  y_down <- seq(0, e[3], -height)
  ncolor <- max(length(y_up), length(y_down))
  pal <- if (!is.function(col)) colorRampPalette(col)(ncolor) else col(ncolor)
  # plot rectangles to simulate colour gradient
  sapply(seq_len(n),
         function(i) {
           rect(min(x), y_up[i], max(x), y_up[i] + height, col = 'red', border = NA)
           rect(min(x), y_down[i], max(x), y_down[i] - height, col = 'blue', border = NA)
         })
  # plot white polygons representing the inverse of the area of interest
  polygon(c(min(x), x, max(x), rev(x)),
          c(e[4], ifelse(y > 0, y, 0),
            rep(e[4], length(y) + 1)), col = 'white', border = NA)
  polygon(c(min(x), x, max(x), rev(x)),
          c(e[3], ifelse(y < 0, y, 0),
            rep(e[3], length(y) + 1)), col = 'white', border = NA)
  lines(x, y)
  abline(h = 0)
  box()
}
# start code for heatmap
pdf(paste(output = baseDir, contrast, "_heatmap.pdf", sep = ""), title = paste(contrast, " - heatmap", sep = ""), width = 100, height = 14)
# get log10 transformed RPKMs from edger
m <- log(vpm(DGEObj)[rownames(selected_cpm_rows),] + 1)
# only keep the selected ones
m <- m[rownames(selected_cpm_rows),]

print(head(m))

m <- m[, labels]

print(m)
# remove infinite numbers before plotting
# i.e. set them NA so that pheatmap ignores them anyway
m[!is.finite(m)] <- NA

pheatmap(
  mat = m,
  color = inferno(2000),
  show_colnames = FALSE,
  show_rownames = FALSE,
  fontsize = 14,
  treeheight_row = 0,
  treeheight_col = 40,
  cellwidth = 20,
  cellheigth = 50,
  scale = "row",
  border_color = "black",
  cluster_rows = T,
  cluster_cols = T,
  display_numbers = F,
  main = paste("Global heatmap\nlog10 CPM values, pval and FDR <", global_fdr)
)

dev.off()
####################################################################################
pdf(paste(baseDir, contrast, ".pdf", sep = ""), title = paste(contrast, " results"))
plotSmear(qlf, xlab = "Log Concentration", ylab = "Log Fold-Change", smooth.scatter = FALSE, lowess = FALSE, de.tags = DEnames, pch = 19, cex = 0.4, main = paste(contrast, " expression", sep = ""))
abline(h = c((-1) * threshold, threshold), col = "blue")

dev.off()

go_list_up <- go_background_blank
go_list_down <- go_background_blank
go_list_both <- go_background_blank

edgerResults <- edgerTable[edgerTable$PValue <= global_pval,]


tmp <- edgerResults[edgerResults$FDR <= global_fdr,]
go_list_both[tmp$GeneID] <- tmp$logFC

tmp <- edgerResults[edgerResults$FDR <= global_fdr & edgerResults$logFC > log_threshold,]
go_list_up[tmp$GeneID] <- tmp$logFC

tmp <- edgerResults[edgerResults$FDR <= global_fdr & edgerResults$logFC < -1 * (log_threshold),]
go_list_down[tmp$GeneID] <- tmp$logFC

length(go_list_both[go_list_both != 0])
length(go_list_up[go_list_up != 0])
length(go_list_down[go_list_down != 0])

allGenes <- edgerResults$logFC
names(allGenes) <- edgerResults$GeneID

topDiffGenesUp <- function(allScore) {
  return(allScore > 0)
}
topDiffGenesDown <- function(allScore) {
  return(allScore < 0)
}
topDiffGenesBoth <- function(allScore) {
  return(allScore > 0 | allScore < 0)
}

pdf(paste(output = baseDir, contrast, "_topgo.pdf", sep = ""), title = paste(contrast, " - topGO", sep = ""), width = 16, height = 10)


paintChart <- function(input, title, color) {
  input$Term <- factor(input$Term, levels = input$Term[order(input$log2)])
  gfx <- ggplot(data = input, aes(x = Term, y = log2, fill = Term)) +
    geom_bar(colour = "black", fill = color, width = .8, stat = "identity") +
    guides(fill = FALSE) +
    theme(text = element_text(size = 14)) +
    xlab("GO Term") +
    ylab(expression('-log'[2] * '(p-value)')) +
    ggtitle(paste(" ", contrast, " | ", title, " in ", focus, sep = "")) +
    coord_flip()
  write.table(input$GO.ID, file = paste(output = baseDir, gsub(" ", "_", title), ".csv", sep = ""), row.names = F, na = "", col.names = TRUE, sep = "\t", quote = FALSE)
  return(gfx)
}


go_annotations <- getBM(attributes = c("ensembl_gene_id", "go_id"), values = names(go_background_blank), mart = ensembl)


go_annotations <- go_annotations %>%
  group_by(ensembl_gene_id) %>%
  mutate(go_id = paste0(go_id, collapse = ", ")) %>%
  distinct()
write.table(go_annotations, file = "go_mapping_tmp", row.names = F, na = "", col.names = F, sep = "\t", quote = FALSE)
geneID2GO <- readMappings(file = "go_mapping_tmp")

test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = global_pval)

go_up_bp <- new("topGOdata", nodeSize = 10, ontology = "BP", description = contrast, geneSel = topDiffGenesUp, allGenes = go_list_up, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_down_bp <- new("topGOdata", nodeSize = 10, ontology = "BP", description = contrast, geneSel = topDiffGenesDown, allGenes = go_list_down, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_up_mf <- new("topGOdata", nodeSize = 10, ontology = "MF", description = contrast, geneSel = topDiffGenesUp, allGenes = go_list_up, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_down_mf <- new("topGOdata", nodeSize = 10, ontology = "MF", description = contrast, geneSel = topDiffGenesDown, allGenes = go_list_down, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_up_cc <- new("topGOdata", nodeSize = 10, ontology = "CC", description = contrast, geneSel = topDiffGenesUp, allGenes = go_list_up, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_down_cc <- new("topGOdata", nodeSize = 10, ontology = "CC", description = contrast, geneSel = topDiffGenesDown, allGenes = go_list_down, annot = annFUN.gene2GO, gene2GO = geneID2GO)

go_bp <- new("topGOdata", nodeSize = 10, ontology = "BP", description = contrast, geneSel = topDiffGenesBoth, allGenes = go_list_both, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_mf <- new("topGOdata", nodeSize = 10, ontology = "MF", description = contrast, geneSel = topDiffGenesBoth, allGenes = go_list_both, annot = annFUN.gene2GO, gene2GO = geneID2GO)
go_cc <- new("topGOdata", nodeSize = 10, ontology = "CC", description = contrast, geneSel = topDiffGenesBoth, allGenes = go_list_both, annot = annFUN.gene2GO, gene2GO = geneID2GO)


resultElimFisBPup <- getSigGroups(go_up_bp, test.stat)
resultElimFisBPdown <- getSigGroups(go_down_bp, test.stat)
resultElimFisBP <- getSigGroups(go_bp, test.stat)

resultElimFisMFup <- getSigGroups(go_up_mf, test.stat)
resultElimFisMFdown <- getSigGroups(go_down_mf, test.stat)
resultElimFisMF <- getSigGroups(go_mf, test.stat)

resultElimFisCCup <- getSigGroups(go_up_cc, test.stat)
resultElimFisCCdown <- getSigGroups(go_down_cc, test.stat)
resultElimFisCC <- getSigGroups(go_cc, test.stat)


allRes_bp <- GenTable(go_bp, elimFisher = resultElimFisBP, topNodes = 200, numChar = 100)
allRes_mf <- GenTable(go_mf, elimFisher = resultElimFisMF, topNodes = 200, numChar = 100)
allRes_cc <- GenTable(go_cc, elimFisher = resultElimFisCC, topNodes = 200, numChar = 100)

allRes_mf[, 6] <- gsub("<", "", allRes_mf[, 6])
allRes_bp[, 6] <- gsub("<", "", allRes_bp[, 6])
allRes_cc[, 6] <- gsub("<", "", allRes_cc[, 6])

print(head(allRes_mf))
allRes_mf[, 6] <- as.numeric(allRes_mf[, 6]);
allRes_mf$log2 <- (-1) * log(allRes_mf[, 6], base = 2)

print(head(allRes_bp))
allRes_bp[, 6] <- as.numeric(allRes_bp[, 6]);
allRes_bp$log2 <- (-1) * log(allRes_bp[, 6], base = 2)

print(head(allRes_cc))
allRes_cc[, 6] <- as.numeric(allRes_cc[, 6]);
allRes_cc$log2 <- (-1) * log(allRes_cc[, 6], base = 2)

allRes_down_bp <- GenTable(go_down_bp, elimFisher = resultElimFisBPdown, topNodes = 200, numChar = 100)
allRes_down_bp[, 6] <- gsub("<", "", allRes_down_bp[, 6])
allRes_down_bp[, 6] <- as.numeric(allRes_down_bp[, 6]);
allRes_down_bp$log2 <- (-1) * log(allRes_down_bp[, 6], base = 2)
paintChart(head(allRes_down_bp, threshold), "GO enrichment: Biological process - down", "#DD8888")

allRes_up_bp <- GenTable(go_up_bp, elimFisher = resultElimFisBPup, topNodes = 200, numChar = 100)
allRes_up_bp[, 6] <- gsub("<", "", allRes_up_bp[, 6])
allRes_up_bp[, 6] <- as.numeric(allRes_up_bp[, 6]);
allRes_up_bp$log2 <- (-1) * log(allRes_up_bp[, 6], base = 2)
paintChart(head(allRes_up_bp, threshold), "GO enrichment: Biological process - up", "#379b20")

allRes_down_mf <- GenTable(go_down_mf, elimFisher = resultElimFisMFdown, topNodes = 200, numChar = 100)
allRes_down_mf[, 6] <- gsub("<", "", allRes_down_mf[, 6])
allRes_down_mf[, 6] <- as.numeric(allRes_down_mf[, 6]);
allRes_down_mf$log2 <- (-1) * log(allRes_down_mf[, 6], base = 2)
paintChart(head(allRes_down_mf, threshold), "GO enrichment: Molecular function - down", "#DD8888")

allRes_up_mf <- GenTable(go_up_mf, elimFisher = resultElimFisMFup, topNodes = 200, numChar = 100)
allRes_up_mf[, 6] <- gsub("<", "", allRes_up_mf[, 6])
allRes_up_mf[, 6] <- as.numeric(allRes_up_mf[, 6]);
allRes_up_mf$log2 <- (-1) * log(allRes_up_mf[, 6], base = 2)
paintChart(head(allRes_up_mf, threshold), "GO enrichment: Molecular function - up", "#379b20")

allRes_down_cc <- GenTable(go_down_cc, elimFisher = resultElimFisCCdown, topNodes = 200, numChar = 100)
allRes_down_cc[, 6] <- gsub("<", "", allRes_down_cc[, 6])
allRes_down_cc[, 6] <- as.numeric(allRes_down_cc[, 6])
allRes_down_cc$log2 <- (-1) * log(allRes_down_cc[, 6], base = 2)
paintChart(head(allRes_down_cc, threshold), "GO enrichment: Cellular component - down", "#DD8888")


allRes_up_cc <- GenTable(go_up_cc, elimFisher = resultElimFisCCup, topNodes = 200, numChar = 100)
allRes_up_cc[, 6] <- gsub("<", "", allRes_up_cc[, 6])
allRes_up_cc[, 6] <- as.numeric(allRes_up_cc[, 6]);
allRes_up_cc$log2 <- (-1) * log(allRes_up_cc[, 6], base = 2)
paintChart(head(allRes_up_cc, threshold), "GO enrichment: Cellular component - up", "#379b20")

dev.off()

run_go_heatmap(allRes_bp, go_down_bp, cpm, "BP", "down", focus)
run_go_heatmap(allRes_cc, go_up_cc, cpm, "CC", "up", focus)
run_go_heatmap(allRes_cc, go_down_cc, cpm, "CC", "down", focus)
run_go_heatmap(allRes_mf, go_up_mf, cpm, "MF", "up", focus)
run_go_heatmap(allRes_mf, go_down_mf, cpm, "MF", "down", focus)
run_go_heatmap(allRes_bp, go_up_bp, cpm, "BP", "up", focus)
## done excel export

genemap <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "wikigene_description", "uniprot_gn_id","gene_biotype"), values = edgerTable$GeneID, mart = ensembl)
genemap2 <- getBM(attributes = c("ensembl_gene_id", "go_id"), filter = "ensembl_gene_id", values = edgerTable$GeneID, mart = ensembl)


edgerTable.idx <- match(edgerTable$GeneID, genemap$ensembl_gene_id)
edgerTable.idx2 <- match(edgerTable$GeneID, genemap2$ensembl_gene_id)
# add biomart data to results tablesTRUE
edgerTable$external_gene_name <- genemap$external_gene_name[edgerTable.idx]
edgerTable$Uniprot <- genemap$uniprot_gn_id[edgerTable.idx]
edgerTable$description <- genemap$wikigene_description[edgerTable.idx]
edgerTable$GO <- genemap2$go_id[edgerTable.idx2]
edgerTable$ensembl_gene_id <- genemap$ensembl_gene_id[edgerTable.idx]
edgerTable$Type <- genemap$gene_biotype[edgerTable.idx]


edgerTable <- edgerTable[c("external_gene_name", "GeneID", "logFC", "logCPM", "PValue", "FDR", "GO", "description", "Uniprot","Type")]

write.table(edgerTable, file = paste(output = baseDir, contrast, "_list.csv", sep = ""), row.names = TRUE, na = "", col.names = TRUE, sep = "\t", quote = FALSE)

edgerResults <- edgerTable[edgerTable$PValue <= global_pval,]

edgerFDR <- edgerResults[edgerResults$FDR <= global_fdr,]

edgerFiltered <- edgerResults[edgerResults$FDR <= global_fdr & abs(edgerResults$logFC) > log_threshold,]

edgerSortedByFC <- edgerFDR[order(edgerFiltered$logFC),]

edgerSortedByFCselection <- edgerSortedByFC[edgerSortedByFC$logFC <= (-1) * log_threshold | edgerSortedByFC$logFC >= log_threshold,]
edgerSortedByFCselection <- edgerSortedByFCselection[order(edgerFiltered$logFC),]

## create excel workbook and file
wb <- createWorkbook()

addWorksheet(wb, sheetName = paste("p<", global_pval, sep = ""));
writeDataTable(wb, sheet = 1, x = edgerResults, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR<", global_fdr, sep = ""));
writeDataTable(wb, sheet = 2, x = edgerFDR, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR<", global_fdr, ", by FC", sep = ""));
writeDataTable(wb, sheet = 3, x = edgerSortedByFC, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("p<", global_pval, ", FDR< ", global_fdr, ", logFC>", log_threshold, sep = ""));
writeDataTable(wb, sheet = 4, x = edgerSortedByFCselection, rowNames = FALSE);

addWorksheet(wb, sheetName = "CPM");
writeDataTable(wb, sheet = 5, x = cpm, rowNames = TRUE);

addWorksheet(wb, sheetName = paste("TopGO up MF, top ", threshold));
writeDataTable(wb, sheet = 6, x = allRes_up_mf, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO down MF, top ", threshold));
writeDataTable(wb, sheet = 7, x = allRes_down_mf, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO up CC, top ", threshold));
writeDataTable(wb, sheet = 8, x = allRes_up_cc, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO down CC, top ", threshold));
writeDataTable(wb, sheet = 9, x = allRes_down_cc, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO up BP, top ", threshold));
writeDataTable(wb, sheet = 10, x = allRes_up_bp, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO down BP, top ", threshold));
writeDataTable(wb, sheet = 11, x = allRes_down_bp, rowNames = FALSE);

addWorksheet(wb, sheetName = paste("TopGO background list"));
writeDataTable(wb, sheet = 12, x = data.frame(names(go_background_blank)), rowNames = FALSE);

saveWorkbook(wb, paste(output = baseDir, contrast, "_edger.xlsx", sep = ""), overwrite = TRUE)

save.image(file = paste(contrast, "_", global_pval, ".RData", sep = ""))

pdf(paste(baseDir, contrast, "_volcano.pdf", sep = ""), title = paste("Volcano plot - ", contrast, sep = ""))

volcano <- edgerTable

# colnames(volcano) <- c("logFC","FDR","status")

volcano$status = 0
volcano$status[ volcano$logFC < 0 &  volcano$FDR < global_fdr ] <-  -1
volcano$status[ volcano$logFC > 0 &  volcano$FDR < global_fdr ] <-  1

volcano$FDR_log10  <- -log10(volcano$FDR)
volcano$FDR_log10[which(!is.finite(volcano$FDR_log10))] <- 0
volcano$FDR_log10[which(volcano$FDR_log10==0)] <- max(volcano$FDR_log10)


volcano$Uniprot <- paste("<a href=\"https://www.uniprot.org/uniprot/",volcano$Uniprot,"\" target=\"_blank\">",volcano$Uniprot,"</a>", sep="")
# volcano$RGD <- paste("<a href=\"https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=",volcano$RGD,"\" target=\"_blank\">",volcano$RGD,"</a>", sep="")

# volcano$uniprot <- c("test")


names(volcano)[names(volcano) == "external_gene_name"] <- "Gene"
names(volcano)[names(volcano) == "description"] <- "Description"

# volcano <- volcano[, c("Gene","GeneID","Uniprot","Description","FDR_log10","logFC","status")]
volcano_plot <- volcano[, c("FDR_log10","logFC","status")]

volcano$status = "not significant"
volcano$status[ volcano$logFC < 0 &  volcano$FDR < global_fdr ] <-  "downregulated"
volcano$status[ volcano$logFC > 0 &  volcano$FDR < global_fdr ] <-  "upregulated"


volcano <- volcano[, c("Gene","GeneID","Uniprot","Type","Description","FDR","status","logFC","FDR_log10")]

names(volcano)[names(volcano) == "status"] <- "Status"


#glXYPlot(   x=volcano_plot$logFC,
#            y=volcano_plot$FDR_log10,
#            xlab="logFC",
#            ylab="log10(FDR)",
#            status = volcano_plot$status,
#            anno = volcano,
#            counts = DGEObj,
#            groups = groups,
#            main = paste(contrast, ": Volcano plot", sep = ""),
#            display.columns = c("Gene","GeneID","Uniprot","Type","Description","FDR","Status"),
#            #sample.cols = inferno(length(unique(groups))),
#            #side.gridstep=TRUE,
#            side.main = c("Gene"),
#            launch=FALSE,
#            folder= paste(contrast,"_glimma",sep=""),
#            path = baseDir,
#            html=paste(contrast,"_volcano",sep="")
#)

print(head(volcano))

volcano <- volcano %>%
  # Modify dataset to add new coloumn of colors
  mutate(color = ifelse(volcano$logFC > 1 & volcano$FDR_log10 > 1.3,
                        yes = "up_gene",
                        no = ifelse(volcano$logFC < -1 & volcano$FDR_log10 > 1.3,
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

print(top_labelled)
print(top_labelled_down)
# Subset table to only show certain gene labels
# Color corresponds to fold change directionality
colored <- ggplot(volcano, aes(x = logFC, y = FDR_log10)) +
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
  theme_bw(base_size = 16) + # clean up theme
  theme(legend.position = "none") + # remove legend
  ggtitle(label = "Volcano Plot", subtitle = contrast) + # add title
  xlab(expression(log[2]("fold change"))) + # x-axis label
  #xlab("log[2] fold change") + # x-axis label
  ylab(expression(-log[10]("FDR"))) + # y-axis label
  geom_vline(xintercept = 0, colour = "black") + # add line at 0
  geom_hline(yintercept = 1.3, colour = "#000000", linetype = "dashed") + # p(0.05) = 1.3
  # annotate(geom = "text",
  #          label = "Down-regulated\ngenes",
  #          x = -2, y = 85,
  #          size = 7, colour = "#3182bd") + # add down_gene text
  # annotate(geom = "text",
  #          label = "Up-regulated\ngenes",
  #          x = 2, y = 85,
  #          size = 7, colour = "#E64B35") + # add up_gene text
  scale_color_manual(values = c("up_gene" = "#E64B35",
                                "down_gene" = "#3182bd",
                                "none" = "#636363")) + # change colors
  #scale_x_continuous(limits=c(-12,12), breaks=seq(-12,12,2)) +
  scale_y_continuous() +
  # Add layer of text annotation to volcano plot.
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
# Plot figure
colored



dev.off()




test <- function(x, go_obj) {
  scores <- as.numeric(unlist(scoresInTerm(go_obj, x)))
  return(scores)
}

plot_cell <- function(go_cat) {
  go_data <- NULL

  if (go_cat == "MF") {
    go_data <- allRes_mf
    goobj <- go_mf
  }

  if (go_cat == "CC") {
    go_data <- allRes_cc
    goobj <- go_cc
  }

  if (go_cat == "BP") {
    go_data <- allRes_bp
    goobj <- go_bp
  }


  terms <- go_data$Term
  id <- go_data$GO.ID
  log2 <- go_data$log2
  annotated <- go_data$Annotated



  df <- data.frame(id, log2, terms, annotated)
  colnames(df) <- c("id", "log2", "term", "Annotated")

  field <- lapply(unlist(id), function(x) { test(x, goobj) })
  df$scores <- field

  df <- head(df, threshold)

  cell.plot(x = setNames(df$log2, df$term),
            cells = df$scores,
            main = paste(contrast, "- GO enrichment -", go_cat, "- top", threshold, "GO terms"),
            x.mar = c(.7, 0),
            key.n = 7,
            y.mar = c(.1, 0),
            cex = 1.6,
            cell.outer = 2,
            cell.limit = 8,
            cell.lwd = 0.05,
            bar.scale = .7,
            #space = .2
  )
  sym.plot(x = setNames(df$log2, df$term),
           cells = df$scores,
           x.annotated = df$Annotated,
           main = paste(contrast, "- GO enrichment -", go_cat, "- top", threshold, "GO terms"),
           x.mar = c(.7, 0),
           # key.n = 7,
           # cex = 1.6,
           #space= 0.2,
           group.cex = .7)
}
pdf(paste(baseDir, contrast, "_cellplot.pdf", sep = ""), title = paste("Cellplot - ", contrast, sep = ""), width = 12, height = 8)

plot_cell("MF")
plot_cell("BP")
plot_cell("CC")

#### Genes per Go category tables come here

allRes_combined_bp <- allRes_bp
allRes_combined_mf <- allRes_mf
allRes_combined_cc <- allRes_cc

go_combined_bp <- go_bp
go_combined_mf <- go_mf
go_combined_cc <- go_cc

pdf(paste(baseDir, contrast, "_GO_BP_word_cloud.pdf", sep = ""), title = paste("GO BP term word cloud - ", contrast, sep = ""))
    wordcloud_function(data = allRes_bp$Term)
dev.off()

pdf(paste(baseDir, contrast, "_GO_CC_word_cloud.pdf", sep = ""), title = paste("GO CC term word cloud - ", contrast, sep = ""))
    wordcloud_function(data = allRes_cc$Term)
dev.off()

pdf(paste(baseDir, contrast, "_GO_MF_word_cloud.pdf", sep = ""), title = paste("GO MF term word cloud - ", contrast, sep = ""))
    wordcloud_function(data = allRes_mf$Term)
dev.off()

categories <- c("bp", "cc", "mf")

test_list <- lapply(unlist(categories), function(go_category) {
  none <- lapply(unlist(c("up", "down", "combined")), function(z) {
    go_data <- get(paste("allRes_", z, "_", go_category, sep = ""))

    terms <- go_data$Term
    id <- go_data$GO.ID
    log2 <- go_data$log2
    annotated <- go_data$Annotated

    main_df <- data.frame(id, log2, terms, annotated, stringsAsFactors = FALSE)
    main_df <- main_df[order(-log2),]

    top_terms <- head(main_df$id, threshold)

    mylist <- top_terms

    goterms <- Term(GOTERM)
    df <- data.frame(c("gene", "go"))

    mylist <- top_terms

    df_list <- lapply(unlist(mylist), function(x) {
      get_data(x, get(paste("go_", z, "_", go_category, sep = "")))
    })
    ## create excel workbook and file
    wb <- createWorkbook()
    counter <- 1
    for (i in seq(1:length(df_list))) {
      if (length(df_list[i][[1]]) > 0) {
        current_term <- gsub("\n", "", strtrim(paste(gsub(":", " ", main_df$id[i]), " ", main_df$terms[i], sep = ""), 31))
        current_term <- gsub(":", " ", current_term)
        current_term <- gsub("/", "_", current_term)

        file_term <- gsub("\n", "", strtrim(paste(gsub(":", "_", main_df$id[i]), "__", main_df$terms[i], sep = ""), 90))
        file_term <- gsub(":", "_", file_term)
        file_term <- gsub(" ", "_", file_term)
        file_term <- gsub(",", "_", file_term)
        file_term <- gsub("/", "_", file_term)

        print_pos <- formatC(counter, width = 3, format = "d", flag = "0")


        addWorksheet(wb, sheetName = current_term)

        writeDataTable(wb, sheet = counter, x = df_list[i][[1]], rowNames = FALSE)
        write.table(df_list[i][[1]], file = paste(baseDir, "GO_gene_lists/", "GO_", toupper(go_category), "_genes_", z, "_pos_", print_pos, "_", file_term, ".csv", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
        counter <- counter + 1
      }
    }

    saveWorkbook(wb, paste(baseDir, "GO_gene_lists/", contrast, "_", go_category, "_", z, ".xlsx", sep = ""), overwrite = TRUE)
  })
})
