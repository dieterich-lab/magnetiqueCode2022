#
set.seed(1234)

#
library(readr)
library(Rsubread)
library(vsn)
library(fgsea)
library(GSA)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(dorothea)
library(progeny)
library(viper)
library(stringr)
library(RColorBrewer)
library(gplots)
library(EnhancedVolcano)
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(CARNIVAL)
library(circlize)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)
library(readr)
library(piano)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(scales)
library(plyr)
library(GSEABase)
library(network)
library(reshape2)
library(cowplot)
library(pheatmap)
library(ggraph)
library(tidygraph)
library(piano)
library(OmnipathR)
library(pals)

load(file = "output/tfList50.RData")
dcmTF <- as.character(tfList$DCM_vs_Healthy$GeneID)
hcmTF <- as.character(tfList$HCM_vs_Healthy$GeneID)

source("~/Documents/src/plotDistributions.R")
source("~/Documents/src/createRegulonList.R")
source("~/Documents/src/support_functions.R")

load(file = "output/output_carnival/dcm_vs_healthy/res.RData")
dcmNetwork <- res$weightedSIF
dcmAttrib <- res$nodesAttributes
idx <- 1
while(length(idx)>0){
  
  diffdiff <- setdiff(x = dcmNetwork[, 3], y = dcmNetwork[, 1])
  diffdiff <- setdiff(x = diffdiff, dcmTF)
  idx2rem <- which(dcmNetwork[, 3]%in%diffdiff)
  if(length(idx2rem)>0){
    dcmNetwork <- dcmNetwork[-idx2rem, ]
  }
  idx <- diffdiff
  
}
write.table(x = dcmNetwork, file = "output/output_carnival/dcm_vs_healthy/res.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = dcmAttrib, file = "output/output_carnival/dcm_vs_healthy/act.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

load(file = "output/output_carnival/hcm_vs_healthy/res.RData")
hcmNetwork <- res$weightedSIF
hcmAttrib <- res$nodesAttributes
idx <- 1
while(length(idx)>0){
  
  diffdiff <- setdiff(x = hcmNetwork[, 3], y = hcmNetwork[, 1])
  diffdiff <- setdiff(x = diffdiff, dcmTF)
  idx2rem <- which(hcmNetwork[, 3]%in%diffdiff)
  if(length(idx2rem)>0){
    hcmNetwork <- hcmNetwork[-idx2rem, ]
  }
  idx <- diffdiff
  
}
write.table(x = hcmNetwork, file = "output/output_carnival/hcm_vs_healthy/res.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = hcmAttrib, file = "output/output_carnival/hcm_vs_healthy/act.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## Identify alternated protein activity levels
res1 <- dcmNetwork
res2 <- hcmNetwork

act1 <- dcmAttrib
act2 <- hcmAttrib

allProt1 <- unique(c(dcmNetwork[, 1], dcmNetwork[, 3]))
allProt1 <- allProt1[-which(allProt1=="Perturbation")]

allProt2 <- unique(c(hcmNetwork[, 1], hcmNetwork[, 3]))
allProt2 <- allProt2[-which(allProt2=="Perturbation")]

act1 <- act1[which(act1[, 1]%in%allProt1), ]
act2 <- act2[which(act2[, 1]%in%allProt2), ]

commonProt <- intersect(x = act1[, 1], y = act2[, 1])
protSgn <- matrix(data = , nrow = length(commonProt), ncol = 2)
rownames(protSgn) <- commonProt
colnames(protSgn) <- c("DCM_vs_Healthy", "HCM_vs_Healthy")
for(ii in 1:nrow(protSgn)){
  
  protSgn[ii, 1] <- sign(as.numeric(act1[, 5])[which(act1[, 1]==commonProt[ii])])
  protSgn[ii, 2] <- sign(as.numeric(act2[, 5])[which(act2[, 1]==commonProt[ii])])
  
}
altProt <- rownames(protSgn)[which(rowSums(protSgn)==0)]

save(altProt, file = "output/output_carnival/altProt.RData")

## Get common interactions
res1 <- dcmNetwork
res2 <- hcmNetwork

idxRes1 <- c()
for(ii in 1:nrow(res1)){
  
  idx1 <- which(res2[, 1]==res1[ii, 1])
  idx2 <- which(res2[, 2]==res1[ii, 2])
  idx3 <- which(res2[, 3]==res1[ii, 3])
  
  idx <- intersect(x = idx1, y = intersect(x = idx2, y = idx3))
  
  if(length(idx)>0){
    idxRes1 <- c(idxRes1, idx)
  }
  
}

idxRes2 <- c()
for(ii in 1:nrow(res2)){
  
  idx1 <- which(res1[, 1]==res2[ii, 1])
  idx2 <- which(res1[, 2]==res2[ii, 2])
  idx3 <- which(res1[, 3]==res2[ii, 3])
  
  idx <- intersect(x = idx1, y = intersect(x = idx2, y = idx3))
  
  if(length(idx)>0){
    idxRes2 <- c(idxRes2, idx)
  }
  
}

commonInt <- unique(rbind(res2[idxRes1, 1:3], res1[idxRes2, 1:3]))
commonInt <- as.data.frame(commonInt)

colnames(commonInt) = c("from", "sign", "to")
labels_edge = c("-1" = "inhibition", "1" = "activation")
nodes = data.frame(union(commonInt$from, commonInt$to))
colnames(nodes) = "nodes"
nodes$label = nodes$nodes

tidygraph::tbl_graph(nodes = nodes, edges = commonInt) %>%
  ggraph::ggraph(layout = "nicely") + 
  geom_node_point(color = "#C0C0C0", size = 8) +
  geom_edge_link(arrow = arrow(), aes(edge_colour=as.factor(sign))) +
  theme_graph() +
  geom_node_text(aes(label = label), vjust = 0.4)

write.table(x = commonInt, file = "output/output_carnival/carnival_common_interactions.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

## ORA of the common network module
exP <- loadGSC(file = "~/Downloads/human_pathways.gmt")
exP <- exP$gsc

universeDCM <- dcmAttrib[, 1]
allNodesDCM <- unique(c(dcmNetwork[, 1], dcmNetwork[, 3]))
allNodesDCM <- allNodesDCM[-which(allNodesDCM=="Perturbation")]
stats <- rep(0, length(allNodesDCM))
names(stats) <- allNodesDCM
for(ii in 1:length(stats)){
  stats[ii] <- sign(as.numeric(dcmAttrib[which(dcmAttrib[, 1]==allNodesDCM[ii]), 5]))
}
dcmFGSEA <- fgseaSimple(pathways = exP, stats = stats, nperm = 100000, minSize = 5, maxSize = Inf)
topDCM <- dcmFGSEA$pathway[which(dcmFGSEA$pval<=0.1)]

universeHCM <- hcmAttrib[, 1]
allNodesHCM <- unique(c(hcmNetwork[, 1], hcmNetwork[, 3]))
allNodesHCM <- allNodesDCM[-which(allNodesHCM=="Perturbation")]
stats <- rep(0, length(allNodesHCM))
names(stats) <- allNodesHCM
for(ii in 1:length(stats)){
  stats[ii] <- sign(as.numeric(hcmAttrib[which(hcmAttrib[, 1]==allNodesHCM[ii]), 5]))
}
hcmFGSEA <- fgseaSimple(pathways = exP, stats = stats, nperm = 100000, minSize = 5, maxSize = Inf)
topHCM <- hcmFGSEA$pathway[which(hcmFGSEA$pval<=0.1)]

significantPathways <- unique(topDCM, topHCM)

mm <- matrix(data = , nrow = length(significantPathways), ncol = 5)
# rownames(mm) <- significantPathways
colnames(mm) <- c("pathway", "dcm_nes", "dcm_pval", "hcm_nes", "hcm_pval")
for(ii in 1:nrow(mm)){
  
  # mm[ii, 1] <- oraDCM[which(rownames(oraDCM)==significantPathways[ii]), 2]
  # mm[ii, 2] <- oraHCM[which(rownames(oraHCM)==significantPathways[ii]), 2]
  mm[ii, 1] <- significantPathways[ii]
  mm[ii, 2] <- dcmFGSEA$NES[which(dcmFGSEA$pathway==significantPathways[ii])]
  mm[ii, 3] <- dcmFGSEA$pval[which(dcmFGSEA$pathway==significantPathways[ii])]
  mm[ii, 4] <- hcmFGSEA$NES[which(hcmFGSEA$pathway==significantPathways[ii])]
  mm[ii, 5] <- hcmFGSEA$pval[which(hcmFGSEA$pathway==significantPathways[ii])]
  
}

idx <- c()
for(ii in 1:nrow(mm)){
  
  if(sign(as.numeric(mm[ii, 2]))==((-1)*sign(as.numeric(mm[ii, 4])))){
    idx <- c(idx, ii)
  }
  
}
altPx <- mm[idx, 1]

df1 <- matrix(data = , nrow = length(altPx), ncol = 4)
colnames(df1) <- c("pathway", "nes", "pval", "comparison")
df1[, 1] <- altPx
for(ii in 1:nrow(df1)){
  df1[ii, 2] <- dcmFGSEA$NES[which(dcmFGSEA$pathway==altPx[ii])]
  df1[ii, 3] <- dcmFGSEA$pval[which(dcmFGSEA$pathway==altPx[ii])]
  df1[ii, 4] <- "dcm_vs_healthy"
}
df1 <- as.data.frame(df1)
df1$nes <- as.numeric(df1$nes)
df1$pval <- as.numeric(df1$pval)

df2 <- matrix(data = , nrow = length(altPx), ncol = 4)
colnames(df2) <- c("pathway", "nes", "pval", "comparison")
df2[, 1] <- altPx
for(ii in 1:nrow(df2)){
  df2[ii, 2] <- hcmFGSEA$NES[which(hcmFGSEA$pathway==altPx[ii])]
  df2[ii, 3] <- hcmFGSEA$pval[which(hcmFGSEA$pathway==altPx[ii])]
  df2[ii, 4] <- "hcm_vs_healthy"
}
df2 <- as.data.frame(df2)
df2$nes <- as.numeric(df2$nes)
df2$pval <- as.numeric(df2$pval)

df <- rbind(df1, df2)
colnames(df) <- c("Pathway", "NES", "pval", "Comparison")

pdf(file = "output/output_carnival/altered_pathways.pdf", 
    width = 12, height = 7)
pp <- ggplot(df,aes(x=Comparison,y=NES,fill=-log10(pval))) + 
  geom_col(position="dodge",width=0.8) +
  coord_flip() + #scale_fill_viridis(trans='log10',option="B")+
  scale_fill_gradientn(colours=coolwarm(100), guide = "colourbar")+
  facet_grid(Pathway~.)+
  theme(strip.text.y = element_text(angle = 0, size = 12))
plot(pp)
dev.off()
