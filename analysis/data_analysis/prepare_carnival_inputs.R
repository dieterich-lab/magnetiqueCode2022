#! /usr/bin/env Rscript

#
set.seed(66901)

#
dir.create("output")

#
library("readr")
library("Rsubread")
library("vsn")
library("fgsea")
library("GSA")
library("pheatmap")
library("edgeR")
library("limma")
library("dplyr")
library("tidyr")
library("dorothea")
library("progeny")
library("viper")
library("stringr")
library("RColorBrewer")
library("gplots")
library("EnhancedVolcano")
library("readxl")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("CARNIVAL")
library("circlize")
library("ComplexHeatmap")
library("VennDiagram")
library("tibble")
library("viridis")
library("BiRewire")
library("ggpubr")
library("foreach")
library("doParallel")
library("openxlsx")
library("biomaRt")
library("OmnipathR")

gene_counts <- read.csv("../Data/gene_count_matrix.csv",row.names="gene_id",as.is=T)
colnames(gene_counts)<-gsub("_stringtieRef","",colnames(gene_counts))

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org");
resMArt=getBM(attributes=c("ensembl_gene_id","external_gene_name","pfam"),mart=mart)

MAGNet_data <- load('../Data/MAGNet_data.RData')
meta <- as.data.frame(cell.attrs)
gene_counts <- cts[as.logical(as.numeric(levels(gene.attrs$GeneFlag))[gene.attrs$GeneFlag]),
                   as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag])]
meta <- meta[as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag]),]
meta$Etiology <- droplevels(meta$Etiology)

ins=intersect(rownames(meta),colnames(gene_counts))

meta <- meta[ins, ]
gene_counts <- gene_counts[, ins]

y <- DGEList(counts=gene_counts)

design <- model.matrix(~ meta$Etiology + meta$Race + meta$Sex + meta$Age + meta$DuplicationRate)

keep <- filterByExpr(y, design=design[,2:4])
ensg_keep <- names(keep)[which(keep)]
keep_genes <- unique(resMArt$external_gene_name[which(resMArt$ensembl_gene_id%in%ensg_keep)])
save(keep_genes, file = "keep_genes.RData")

source("../src/support_functions.R")
source("estimate_significance.R")

dirloc <- "/prj/MAGE/analysis/genetonic/results" 

results <- readRDS(file.path(dirloc, "MAGNet_DCMvsNFD_GeneTonic.rds"))
DCM_vs_NFD <- matrix(data = , nrow = length(results$res_de$baseMean), ncol = length(results$res_de))
DCM_vs_NFD[, 1] <- results$res_de$baseMean
DCM_vs_NFD[, 2] <- results$res_de$log2FoldChange
DCM_vs_NFD[, 3] <- results$res_de$lfcSE
DCM_vs_NFD[, 4] <- results$res_de$stat
DCM_vs_NFD[, 5] <- results$res_de$pvalue
DCM_vs_NFD[, 6] <- results$res_de$padj
DCM_vs_NFD[, 7] <- results$res_de$weight
DCM_vs_NFD[, 8] <- results$res_de$SYMBOL
colnames(DCM_vs_NFD) <- names(results$res_de)
rownames(DCM_vs_NFD) <- names(results$res_de$SYMBOL)
DCM_vs_NFD <- as.data.frame(DCM_vs_NFD)
DCM_vs_NFD$baseMean <- as.numeric(DCM_vs_NFD$baseMean)
DCM_vs_NFD$log2FoldChange <- as.numeric(DCM_vs_NFD$log2FoldChange)
DCM_vs_NFD$lfcSE <- as.numeric(DCM_vs_NFD$lfcSE)
DCM_vs_NFD$stat <- as.numeric(DCM_vs_NFD$stat)
DCM_vs_NFD$pvalue <- as.numeric(DCM_vs_NFD$pvalue)
DCM_vs_NFD$padj <- as.numeric(DCM_vs_NFD$padj)
DCM_vs_NFD$weight <- as.numeric(DCM_vs_NFD$weight)
DCM_vs_NFD$SYMBOL <- as.character(DCM_vs_NFD$SYMBOL)

results <- readRDS(file.path(dirloc, "MAGNet_HCMvsNFD_GeneTonic.rds"))
HCM_vs_NFD <- matrix(data = , nrow = length(results$res_de$baseMean), ncol = length(results$res_de))
HCM_vs_NFD[, 1] <- results$res_de$baseMean
HCM_vs_NFD[, 2] <- results$res_de$log2FoldChange
HCM_vs_NFD[, 3] <- results$res_de$lfcSE
HCM_vs_NFD[, 4] <- results$res_de$stat
HCM_vs_NFD[, 5] <- results$res_de$pvalue
HCM_vs_NFD[, 6] <- results$res_de$padj
HCM_vs_NFD[, 7] <- results$res_de$weight
HCM_vs_NFD[, 8] <- results$res_de$SYMBOL
colnames(HCM_vs_NFD) <- names(results$res_de)
rownames(HCM_vs_NFD) <- names(results$res_de$SYMBOL)
HCM_vs_NFD <- as.data.frame(HCM_vs_NFD)
HCM_vs_NFD$baseMean <- as.numeric(HCM_vs_NFD$baseMean)
HCM_vs_NFD$log2FoldChange <- as.numeric(HCM_vs_NFD$log2FoldChange)
HCM_vs_NFD$lfcSE <- as.numeric(HCM_vs_NFD$lfcSE)
HCM_vs_NFD$stat <- as.numeric(HCM_vs_NFD$stat)
HCM_vs_NFD$pvalue <- as.numeric(HCM_vs_NFD$pvalue)
HCM_vs_NFD$padj <- as.numeric(HCM_vs_NFD$padj)
HCM_vs_NFD$weight <- as.numeric(HCM_vs_NFD$weight)
HCM_vs_NFD$SYMBOL <- as.character(HCM_vs_NFD$SYMBOL)

results <- readRDS(file.path(dirloc, "MAGNet_DCMvsHCM_GeneTonic.rds"))
DCM_vs_HCM <- matrix(data = , nrow = length(results$res_de$baseMean), ncol = length(results$res_de))
DCM_vs_HCM[, 1] <- results$res_de$baseMean
DCM_vs_HCM[, 2] <- results$res_de$log2FoldChange
DCM_vs_HCM[, 3] <- results$res_de$lfcSE
DCM_vs_HCM[, 4] <- results$res_de$stat
DCM_vs_HCM[, 5] <- results$res_de$pvalue
DCM_vs_HCM[, 6] <- results$res_de$padj
DCM_vs_HCM[, 7] <- results$res_de$weight
DCM_vs_HCM[, 8] <- results$res_de$SYMBOL
colnames(DCM_vs_HCM) <- names(results$res_de)
rownames(DCM_vs_HCM) <- names(results$res_de$SYMBOL)
DCM_vs_HCM <- as.data.frame(DCM_vs_HCM)
DCM_vs_HCM$baseMean <- as.numeric(DCM_vs_HCM$baseMean)
DCM_vs_HCM$log2FoldChange <- as.numeric(DCM_vs_HCM$log2FoldChange)
DCM_vs_HCM$lfcSE <- as.numeric(DCM_vs_HCM$lfcSE)
DCM_vs_HCM$stat <- as.numeric(DCM_vs_HCM$stat)
DCM_vs_HCM$pvalue <- as.numeric(DCM_vs_HCM$pvalue)
DCM_vs_HCM$padj <- as.numeric(DCM_vs_HCM$padj)
DCM_vs_HCM$weight <- as.numeric(DCM_vs_HCM$weight)
DCM_vs_HCM$SYMBOL <- as.character(DCM_vs_HCM$SYMBOL)

#
DCM_vs_NFD <- DCM_vs_NFD[complete.cases(DCM_vs_NFD), ]
HCM_vs_NFD <- HCM_vs_NFD[complete.cases(HCM_vs_NFD), ]
DCM_vs_HCM <- DCM_vs_HCM[complete.cases(DCM_vs_HCM), ]

ttopList <- list()
ttopList[[length(ttopList)+1]] <- DCM_vs_NFD
ttopList[[length(ttopList)+1]] <- HCM_vs_NFD
ttopList[[length(ttopList)+1]] <- DCM_vs_HCM

names(ttopList) <- c("DCM_vs_NFD", "HCM_vs_NFD", "DCM_vs_HCM")

save(ttopList, file = "output/ttopList.RData")

wb <- createWorkbook()
addWorksheet(wb, "DCM_vs_NFD")
addWorksheet(wb, "HCM_vs_NFD")
addWorksheet(wb, "DCM_vs_HCM")

writeDataTable(wb, "DCM_vs_NFD", x = ttopList$DCM_vs_NFD, colNames = TRUE, rowNames = TRUE)
writeDataTable(wb, "HCM_vs_NFD", x = ttopList$HCM_vs_NFD, colNames = TRUE, rowNames = TRUE)
writeDataTable(wb, "DCM_vs_HCM", x = ttopList$DCM_vs_HCM, colNames = TRUE, rowNames = TRUE)

saveWorkbook(wb, "output/DGE_List.xlsx", overwrite = TRUE)

# Most regulated genes
for(ii in 1:length(ttopList)){

  temp <- ttopList[[ii]]
  idx <- intersect(x = which(abs(temp$log2FoldChange)>1), y = which(temp$padj<0.05))
  idx1 <- intersect(x = which(temp$log2FoldChange>1), y = which(temp$padj<0.05))
  idx2 <- intersect(x = which(temp$log2FoldChange < -1), y = which(temp$padj<0.05))
  pdf(file = paste0("output/volcano_", names(ttopList)[ii], ".pdf"), width = 12, height = 8)
  EnhancedVolcano(toptable = temp, lab = temp$SYMBOL, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1,
                  pointSize = 3, title = "Volcano plot of differentially expressed genes",
                  subtitle = paste0(length(idx), " significantly enriched genes (pCutoff = 0.05, Log-FCcutoff = 1): ", length(idx1), " up & ", length(idx2), " down"))
  dev.off()

  # idx1 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][1])
  # idx2 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][2])
  # # top50 <- temp$ID[order(abs(temp$log2FoldChange), decreasing = TRUE)[1:50]]
  # top50 <- temp$geneID[order(abs(temp$log2FoldChange), decreasing = TRUE)[1:50]]
  # select <- which(rownames(data)%in%top50)
  # mycol <- colorRamp2(c(-2,0,2), c("dodgerblue", "black", "yellow"))
  # mm <- as.matrix(data)[select, c(idx1, idx2)]
  # pdf(file = paste0("output/top50_expr_genes_scaled_", names(ttopList)[ii], ".pdf"),
  #     width = 60, height = 20)
  # Heatmap(scale(mm), col=mycol, cluster_columns = FALSE, cluster_rows = TRUE,
  #         column_title = "Expression heatmap of top 50 regulated genes")
  # dev.off()

}

sets <- list()
for(ii in 1:length(ttopList)){
  temp <- ttopList[[ii]]
  sets[[length(sets)+1]] <- as.character(temp$SYMBOL[order(abs(temp$log2FoldChange), decreasing = TRUE)[1:100]])
}
DCM_vs_NFD <- sets[[1]]
HCM_vs_NFD <- sets[[2]]
DCM_vs_HCM <- sets[[3]]

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(DCM_vs_NFD, HCM_vs_NFD, DCM_vs_HCM),
  category.names = c("DCM vs NFD" , "HCM vs NFD" , "DCM vs HCM"),
  filename = 'output/venn_diagramm_top100_genes.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

## DoRothEA Analysis
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C", "D"))

idx1 <- which(regulons$tf%in%keep_genes)
idx2 <- which(regulons$target%in%keep_genes)

regulons <- regulons[intersect(x = idx1, y = idx2), ]

tfList <- list()
tfListAll <- list()

for(ii in 1:3){
  
  print(paste0("Step ---- ", ii, "/", 3))
  
  stats <- ttopList[[ii]]$log2FoldChange[which(ttopList[[ii]]$SYMBOL%in%keep_genes)]
  names(stats) <- ttopList[[ii]]$SYMBOL[which(ttopList[[ii]]$SYMBOL%in%keep_genes)]
  uID <- unique(names(stats))
  ss <- matrix(data = , nrow = length(uID), ncol = 1)
  rownames(ss) <- uID
  colnames(ss) <- names(ttopList)[ii]
  for(jj in 1:nrow(ss)){
    idx <- which(uID==rownames(ss)[jj])
    ss[jj, 1] <- median(as.numeric(stats[idx]))
  }
  
  tf_activities_stat <- estimate_significance(expr = ss, regulons = regulons, nperm = 1000)
  
  tfListAll[[length(tfListAll)+1]] <- tf_activities_stat
  
  tmp <- as.matrix(tf_activities_stat[order(abs(tf_activities_stat$nes), decreasing = TRUE), ])
  tmp <- as.data.frame(tmp)
  write_excel_csv2(x = tmp, 
                   file = paste0("output/tf_activities_all_", names(ttopList)[ii], ".xls"),
                   col_names = TRUE)
  
  # save(tf_activities_stat, file = paste0("output/tf_activities_all_", names(ttopList)[ii], ".RData"))
  
  tmp <- tf_activities_stat[order(tf_activities_stat$nes, decreasing = TRUE), ]
  tmp$significance <- "unsignif"
  tmp$significance[intersect(x = which(tmp$pval<=0.1), y = which(tmp$nes>0))] <- "signif up"
  tmp$significance[intersect(x = which(tmp$pval<=0.1), y = which(tmp$nes<0))] <- "signif dn"
  
  pdf(file = paste0("output/", tolower(names(ttopList)[ii]), ".pdf"), width = 40, height = 14)
  ggdotchart(tmp, x = "id", y = "nes",
             color = "significance",                                # Color by groups
             palette = c("#3333FF", "#FF0000", "#808080"), # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             # group = "significance",                                # Order by groups
             dot.size = 6,                                 # Large dot size
             label = round(tmp$nes, 2),                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr()                        # ggplot2 theme
  )+
    geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
  dev.off()
  
  # tf_activities_stat_signif <- tmp[-which(tmp$significance=="unsignif"), ] 
  # pdf(file = paste0("output/", tolower(names(ttopList)[ii]), "_signif.pdf"), width = 18, height = 12)
  # ggplot(tf_activities_stat_signif,aes(x = reorder(id, nes), y = nes)) +
  #   geom_bar(aes(fill = nes), stat = "identity") +
  #   scale_fill_gradient2(low = "midnightblue", high = "firebrick4",
  #                        mid = "whitesmoke", midpoint = 0) +
  #   theme_minimal() +
  #   theme(axis.title = element_text(face = "bold", size = 12),
  #         axis.text.x =
  #           element_text(angle = 45, hjust = 1, size =10, face= "bold"),
  #         axis.text.y = element_text(size =10, face= "bold"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   xlab("Transcription Factors")
  # dev.off()
  tmp2 <- tmp[-which(tmp$significance=="unsignif"), ]
  rownames(tmp2) <- tmp2$id
  require(dplyr)
  tf_activities_stat_top25 <- tmp2 %>%
    as.data.frame() %>%
    dplyr::top_n(nrow(tmp2), wt = abs(nes)) %>%
    dplyr::arrange(nes) %>%
    dplyr::mutate(GeneID = factor(id))
  
  pdf(file = paste0("output/", tolower(names(ttopList)[ii]), "_signif.pdf"), width = 14, height = 7)
  ggplot(tf_activities_stat_top25,aes(x = reorder(id, nes), y = nes)) +
    geom_bar(aes(fill = nes), stat = "identity") +
    # scale_fill_gradient2(low = "darkblue", high = "indianred",
    #                      mid = "whitesmoke", midpoint = 0) +
    scale_fill_viridis_c(option = 'virdis', direction = -1) +
    xlab("Transcription Factor") + ylab("NES") +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x =
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  dev.off()
  
  tfList[[length(tfList)+1]] <- tf_activities_stat_signif
  
  # save(tf_activities_stat_signif, file = paste0("tf_activities_stat_signif_", tolower(names(ttopList)[ii]), ".RData"))
  
}

names(tfList) <- names(ttopList)
names(tfListAll) <- names(ttopList)

save(tfList, file = "output/tfList.RData")
save(tfListAll, file = "output/tfListAll.RData")

wb <- createWorkbook()
addWorksheet(wb, "DCM_vs_NFD")
addWorksheet(wb, "HCM_vs_NFD")
addWorksheet(wb, "DCM_vs_HCM")

writeDataTable(wb, "DCM_vs_NFD", x = tfListAll$DCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "HCM_vs_NFD", x = tfListAll$HCM_vs_NFD, colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, "DCM_vs_HCM", x = tfListAll$DCM_vs_HCM, colNames = TRUE, rowNames = FALSE)

saveWorkbook(wb, "output/TF_Activity_List.xlsx", overwrite = TRUE)

## Prepare the background network
interactions <- import_omnipath_interactions()
interactions <- interactions[which(interactions$is_directed==1), ]
sums <- interactions$is_stimulation+interactions$is_inhibition
interactions <- interactions[which(sums==1), ]
ppi <- matrix(data = "1", nrow = nrow(interactions), ncol = 3)
ppi[, 1] <- interactions$source_genesymbol
ppi[, 3] <- interactions$target_genesymbol
ppi[which(interactions$is_inhibition==1), 2] <- "-1"
colnames(ppi) <- c("source", "sign", "target")
ppi <- as.data.frame(ppi)

save(ppi, "output/ppi.RData")
