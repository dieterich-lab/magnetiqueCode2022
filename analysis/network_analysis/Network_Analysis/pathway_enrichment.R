set.seed(1234)

library("readr")
library("piano")
library("dplyr")
library("ggplot2")
library("tibble")
library("tidyr")
library("dplyr")
library("scales")
library("plyr")
library("GSEABase")
library("network")
library("reshape2")
library("cowplot")
library("pheatmap")
library("ggraph")
library("tidygraph")
library("fgsea")
library("readxl")
library("openxlsx")
library("VennDiagram")

dir.create("output")

gsc = loadGSC("c2.cp.v7.4.symbols.gmt")
geneList <- gsc$gsc
# geneList <- geneList[which(grepl(pattern = "REACTOME_", x = names(geneList), fixed = TRUE))]
# names(geneList) <- gsub(pattern = "REACTOME_", replacement = "", x = names(geneList), fixed = TRUE)

gsc = loadGSC("h.all.v7.4.symbols.gmt")
hallmarks <- gsc$gsc

oraListPx <- list()
oraListHlm <- list()

## DCM vs NFD
# CARNIVAL ORA Px
sif <- read.delim(file = "dcm_vs_nfd/network_dcm_vs_nfd_100.txt")
attr <- read.delim(file = "dcm_vs_nfd/attributes_dcm_vs_nfd_100.txt")

genes <- unique(c(sif[, 1], sif[, 3]))
universe <- attr[, 1]

fgseaOver <- fora(pathways = geneList, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.05))
pdf(file = "output/dcm_vs_nfd_carnival_ora_px.pdf", width = 12, height = 16)
fgseaOver %>% 
  top_n(cnt, wt=-padj) %>% 
  mutate(hitsPerc=overlap*100/size) %>% 
  ggplot(aes(x=hitsPerc, 
             y=pathway, 
             colour=padj, 
             size=overlap)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
dev.off()

fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)

oraListPx[[length(oraListPx)+1]] <- fgseaOver

# CARNIVAL ORA Px
fgseaOver <- fora(pathways = hallmarks, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.1))
if(cnt>0){
  
  pdf(file = "output/dcm_vs_nfd_carnival_ora_hlm.pdf", width = 8, height = 8)
  fgseaOver %>% 
    top_n(cnt, wt=-padj) %>% 
    mutate(hitsPerc=overlap*100/size) %>% 
    ggplot(aes(x=hitsPerc, 
               y=pathway, 
               colour=padj, 
               size=overlap)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
  dev.off()
  
  fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)
  
  oraListHlm[[length(oraListHlm)+1]] <- fgseaOver
  
}


## HCM vs NFD
# CARNIVAL ORA Px
sif <- read.delim(file = "hcm_vs_nfd/network_hcm_vs_nfd_100.txt")
attr <- read.delim(file = "hcm_vs_nfd/attributes_hcm_vs_nfd_100.txt")

genes <- unique(c(sif[, 1], sif[, 3]))
universe <- attr[, 1]

fgseaOver <- fora(pathways = geneList, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.05))
pdf(file = "output/hcm_vs_nfd_carnival_ora_px.pdf", width = 12, height = 18)
fgseaOver %>% 
  top_n(cnt, wt=-padj) %>% 
  mutate(hitsPerc=overlap*100/size) %>% 
  ggplot(aes(x=hitsPerc, 
             y=pathway, 
             colour=padj, 
             size=overlap)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
dev.off()

fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)

oraListPx[[length(oraListPx)+1]] <- fgseaOver

# CARNIVAL ORA Hlm
fgseaOver <- fora(pathways = hallmarks, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.1))
if(cnt>0){
  
  pdf(file = "output/hcm_vs_nfd_carnival_ora_hlm.pdf", width = 8, height = 8)
  fgseaOver %>% 
    top_n(cnt, wt=-padj) %>% 
    mutate(hitsPerc=overlap*100/size) %>% 
    ggplot(aes(x=hitsPerc, 
               y=pathway, 
               colour=padj, 
               size=overlap)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
  dev.off()
  
  fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)
  
  oraListHlm[[length(oraListHlm)+1]] <- fgseaOver
  
}

## DCM vs HCM
# CARNIVAL ORA Pathways
sif <- read.delim(file = "dcm_vs_hcm/network_dcm_vs_hcm_100.txt")
attr <- read.delim(file = "dcm_vs_hcm/attributes_dcm_vs_hcm_100.txt")

genes <- unique(c(sif[, 1], sif[, 3]))
universe <- attr[, 1]

fgseaOver <- fora(pathways = geneList, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.05))
pdf(file = "output/dcm_vs_hcm_carnival_ora_px.pdf", width = 12, height = 16)
fgseaOver %>% 
  top_n(cnt, wt=-padj) %>% 
  mutate(hitsPerc=overlap*100/size) %>% 
  ggplot(aes(x=hitsPerc, 
             y=pathway, 
             colour=padj, 
             size=overlap)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
dev.off()

fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)

oraListPx[[length(oraListPx)+1]] <- fgseaOver

# CARNIVAL ORA Hallmarks
fgseaOver <- fora(pathways = hallmarks, genes = genes, universe = universe, minSize = 5, maxSize = Inf)

cnt <- length(which(fgseaOver$padj<=0.1))
if(cnt>0){
  
  pdf(file = "output/dcm_vs_hcm_carnival_ora_hlm.pdf", width = 14, height = 12)
  fgseaOver %>% 
    top_n(length(idx2keep), wt=-padj) %>% 
    mutate(hitsPerc=overlap*100/size) %>% 
    ggplot(aes(x=hitsPerc, 
               y=pathway, 
               colour=padj, 
               size=overlap)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="Pathway", colour="adjusted p-value", size="Count")
  dev.off()
  
  fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)
  
  oraListHlm[[length(oraListHlm)+1]] <- fgseaOver
  
} else {
  
  fgseaOver$overlapGenes <- as.character(fgseaOver$overlapGenes)
  
  oraListHlm[[length(oraListHlm)+1]] <- fgseaOver
  
}

## Save workbooks
# CARNIVAL ORA Px
names(oraListPx) <- c("DCM_vs_NFD", "HCM_vs_NFD", "DCM_vs_HCM")

wb <- createWorkbook()
addWorksheet(wb, "DCM_vs_NFD")
addWorksheet(wb, "HCM_vs_NFD")
addWorksheet(wb, "DCM_vs_HCM")

writeDataTable(wb, "DCM_vs_NFD", x = oraListPx$DCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "HCM_vs_NFD", x = oraListPx$HCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "DCM_vs_HCM", x = oraListPx$DCM_vs_HCM, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "output/carnival_nodes_over_representation_analysis_gene_sets.xlsx", overwrite = TRUE)
save(oraListPx, file = "output/carnival_nodes_over_representation_analysis_gene_sets.RData")

# CARNIVAL ORA Hlm
names(oraListHlm) <- c("DCM_vs_NFD", "HCM_vs_NFD", "DCM_vs_HCM")

wb <- createWorkbook()
addWorksheet(wb, "DCM_vs_NFD")
addWorksheet(wb, "HCM_vs_NFD")
addWorksheet(wb, "DCM_vs_HCM")

writeDataTable(wb, "DCM_vs_NFD", x = oraListHlm$DCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "HCM_vs_NFD", x = oraListHlm$HCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "DCM_vs_HCM", x = oraListHlm$DCM_vs_HCM, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "output/carnival_nodes_over_representation_analysis_hallmark_sets.xlsx", overwrite = TRUE)

## Check Overlap and Different Sets
setDCM <- oraListPx$DCM_vs_NFD$pathway[which(oraListPx$DCM_vs_NFD$padj<=0.05)]
setHCM <- oraListPx$HCM_vs_NFD$pathway[which(oraListPx$HCM_vs_NFD$padj<=0.05)]
setCommon <- intersect(x = setDCM, y = setHCM)

venn.diagram(list(DCM_vs_NFD = setDCM, HCM_vs_NFD = setHCM), fill = c("purple", "yellow"), 
             alpha = c(0.5, 0.5), lwd =0, "output/venn_diagram.tiff", height = 5000, width = 5000, 
             main = "Overlap of Significant REACTOME Sets")
