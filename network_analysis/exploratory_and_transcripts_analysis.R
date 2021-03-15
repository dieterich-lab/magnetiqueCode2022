#
set.seed(1234)

#
dir.create("output_analysis1")

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
# library(xlsx)
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(CARNIVAL)
library(circlize)
library(ComplexHeatmap)
library(VennDiagram)
library(tibble)

# source("~/Documents/src/limmaWrapper.R")
# source("~/Documents/src/pianoWrapper.R")
source("~/Documents/src/plotDistributions.R")
source("~/Documents/src/createRegulonList.R")
source("~/Documents/src/support_functions.R")

MAGE_metadata <- read.csv("~/Desktop/MAGE-Project/Data/MAGE_metadata.txt")

gene_count_matrix <- read_csv("~/Desktop/MAGE-Project/Data/gene_count_matrix.csv")
rNames <- gene_count_matrix$gene_id
gene_count_matrix <- gene_count_matrix[, 2:ncol(gene_count_matrix)]

counts <- gene_count_matrix
rownames(counts) <- rNames
colnames(counts) <- gsub(pattern = "_stringtieRef", replacement = "", x = colnames(counts), fixed = TRUE)

ensgID <- rownames(counts)
MappedID_MAGE <- read.delim("~/Desktop/MAGE-Project/Data/MappedID_MAGE_Stringtie.txt", header=FALSE)

counts <- as.data.frame(counts)

varQuant <- 0.01
countQuant <- 0.01

chk <- counts[rowSums(counts) > 10, ]
var.all <- apply(chk, 1, var)
chk <- chk[var.all > quantile(var.all, probs=varQuant, type=8), ]
sf <- DESeq2::estimateSizeFactorsForMatrix(chk)
cts.norm <- t(t(chk)/sf)
keep <- rowMeans(cts.norm) > quantile(rowMeans(cts.norm), probs=countQuant, type=8)
chk <- chk[keep,]

counts <- chk

##
targets <- matrix(data = , nrow = ncol(counts), ncol = 2)
for(ii in 1:ncol(counts)){
  
  targets[ii, 1] <- colnames(counts)[ii]
  idx <- which(MAGE_metadata$Run==colnames(counts)[ii])
  if(length(idx) > 0){
    
    targets[ii, 2] <- MAGE_metadata$etiology[idx]
    
  }
  
}

targets <- as.data.frame(targets)
targets$V1 <- as.character(targets$V1)
targets$V2 <- as.character(targets$V2)
targets <- targets[-which(targets$V2==""), ]
idx2rem <- which(colnames(counts)%in%setdiff(x = colnames(counts), y = targets$V1))
counts <- counts[, -idx2rem]
colnames(targets) <- c("sample", "condition")

targets$condition[which(grepl(pattern = "Donor", x = targets$condition))] <- "Healthy"
targets$condition[which(grepl(pattern = "DCM", x = targets$condition))] <- "DCM"
targets$condition[which(grepl(pattern = "HCM", x = targets$condition))] <- "HCM"
targets$condition[which(grepl(pattern = "PPCM", x = targets$condition))] <- "PPCM"

targets$sample[which(targets$condition=="Healthy")] <- 
  paste0("Healthy_", 1:length(which(targets$condition=="Healthy")))

targets$sample[which(targets$condition=="DCM")] <- 
  paste0("DCM_", 1:length(which(targets$condition=="DCM")))

targets$sample[which(targets$condition=="HCM")] <- 
  paste0("HCM_", 1:length(which(targets$condition=="HCM")))

targets$sample[which(targets$condition=="PPCM")] <- 
  paste0("PPCM_", 1:length(which(targets$condition=="PPCM")))

idx2rem <- which(targets$condition=="PPCM")
targets <- targets[-idx2rem, ]
counts <- counts[, -idx2rem]

##
DGE <- DGEList(counts) %>%  calcNormFactors() 
data <- voom(DGE)
df <- data$E

df <- as.data.frame(df)
data <- df

plotDistributions(data = df, path = "output_analysis1/", width = 40, height = 15)

## All groups
data.pca <- t(df)
data.pca <- cbind(data.pca, as.matrix(targets))
data.pca <- as.data.frame(data.pca)
data.pca[, 1:(ncol(data.pca)-2)] <- lapply(data.pca[, 1:(ncol(data.pca)-2)], function(x) as.numeric(as.character(x)))

res.pca <- prcomp(data.pca[, -c(ncol(data.pca)-1, ncol(data.pca))], scale. = TRUE)
res.plot <-  as.data.frame(cbind(res.pca$x[, 1], res.pca$x[, 2], as.character(data.pca$condition), rownames(data.pca)))
res.plot[, 1:2] <- lapply(res.plot[, 1:2], function(x) as.numeric(as.character(x)))
res.plot[, 3:4] <- lapply(res.plot[, 3:4], function(x) as.character(x))
colnames(res.plot) <- c("pc1", "pc2", "Group", "sample")
percentages <- ((res.pca$sdev)^2 / sum(res.pca$sdev^2)*100)[1:2]

pdf(file = "output_analysis1/pca_of_samples.pdf", width = 30, height = 30)
pp <- ggplot(res.plot, aes(x=pc1, y=pc2, color=Group)) +
  geom_point(size=10, alpha = 0.5) +
  #geom_text(size=6) +
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=40,face="bold")) +
  #geom_path(arrow=arrow()) +
  #theme_minimal() +
  xlab(paste0("PC1 (", round(x = percentages[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(x = percentages[2], digits = 2), "%)")) +
  # xlim(c(-max(abs(res.pca$x[, 1])),max(abs(res.pca$x[, 1])))) +
  # ylim(c(-max(abs(res.pca$x[, 2])),max(abs(res.pca$x[, 2])))) + 
  xlim(c(-220, 220)) +
  ylim(c(-220, 220)) +
  theme(legend.position = "none") +
  geom_text_repel(data = res.plot, aes(label=sample), box.padding = 1.5, size=5)
plot(pp)
dev.off()

## Correlation Heatmap
pheatmap(mat = cor(x = df), cluster_cols = TRUE, cluster_rows = TRUE, width = 32, 
         height = 32, filename = "output_analysis1/correlation_of_samples_plot.pdf", main = "Correlation of samples")

##
ttopList <- list()

limmaRes <- runLimma(measurements = df, targets = targets, comparisons = list(c(3, -1)))
ttop <- topTable(fit = limmaRes[[1]], coef = 1, number = nrow(data), adjust.method = "fdr")
ttop$ID <- rownames(ttop)
ttopList[[length(ttopList)+1]] <- ttop

limmaRes <- runLimma(measurements = df, targets = targets, comparisons = list(c(2, -1)))
ttop <- topTable(fit = limmaRes[[1]], coef = 1, number = nrow(data), adjust.method = "fdr")
ttop$ID <- rownames(ttop)
ttopList[[length(ttopList)+1]] <- ttop

names(ttopList) <- c("DCM_vs_Healthy", "HCM_vs_Healthy")

MappedID_MAGE_Stringtie <- read.delim("~/Desktop/MAGE-Project/Data/MappedID_MAGE_Stringtie.txt", header=FALSE)

ttopTemp <- list()
for(ii in 1:length(ttopList)){
  tmp <- ttopList[[ii]]
  for(jj in 1:nrow(tmp)){
    idx <- which(MappedID_MAGE_Stringtie$V1==tmp$ID[jj])
    if(length(idx)>0){
      vv <- MappedID_MAGE_Stringtie$V2[idx]
      if(vv==""){
        tmp$ID[jj] <- rownames(tmp)[jj]
      } else {
        tmp$ID[jj] <- MappedID_MAGE_Stringtie$V2[idx]
      }
    }
  }
  ttopTemp[[length(ttopTemp)+1]] <- tmp
}
ttopList <- ttopTemp
names(ttopList) <- c("DCM_vs_Healthy", "HCM_vs_Healthy")

for(ii in 1:length(ttopList)){
  
  temp <- ttopList[[ii]]
  write_excel_csv(x = temp, file = paste0("output_analysis1/diff_gene_expr_", names(ttopList)[ii], ".xls"), 
                  col_names = TRUE, delim = "\t")
  # write.csv(x = temp, file = paste0("output_analysis1/diff_gene_expr_", names(ttopList)[ii], ".csv"), 
  #           quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}

for(ii in 1:length(ttopList)){
  
  temp <- ttopList[[ii]]
  idx <- intersect(x = which(abs(temp$logFC)>2), y = which(temp$adj.P.Val<0.05))
  idx1 <- intersect(x = which(temp$logFC>2), y = which(temp$adj.P.Val<0.05))
  idx2 <- intersect(x = which(temp$logFC < -2), y = which(temp$adj.P.Val<0.05))
  pdf(file = paste0("output_analysis1/volcano_", names(ttopList)[ii], ".pdf"), width = 12, height = 8)
  EnhancedVolcano(toptable = temp, lab = temp$ID, x = "logFC", y = "adj.P.Val", pCutoff = 0.05, FCcutoff = 1,
                  pointSize = 3, title = "Volcano plot of differentially expressed genes",
                  subtitle = paste0(length(idx), " significantly enriched genes (pCutoff = 0.05, Log-FCcutoff = 2): ", length(idx1), " up & ", length(idx2), " down"))
  dev.off()
  
  # idx1 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][1])
  # idx2 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][2])
  # top50 <- temp$ID[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
  # hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  # select <- which(rownames(data)%in%top50)
  # pdf(file = paste0("output/top50_expr_genes_", names(ttopList)[ii], ".pdf"), width = 40, height = 30)
  # heatmap.2(t(as.matrix(data)[select, c(idx1, idx2)]), col = hmcol,
  #           Rowv = FALSE, Colv = FALSE, scale="none",
  #           dendrogram="none", trace="none", margin=c(12,8), 
  #           main = "Expression heatmap of top 50 regulated genes", cexRow = .9)
  # dev.off()
  
  idx1 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][1])
  idx2 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][2])
  # top50 <- temp$ID[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
  top50 <- rownames(temp)[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
  select <- which(rownames(data)%in%top50)
  mycol <- colorRamp2(c(-3,0,3), c("dodgerblue", "black", "yellow"))
  mm <- as.matrix(data)[select, c(idx1, idx2)]
  pdf(file = paste0("output_analysis1/top50_expr_genes_scaled_", names(ttopList)[ii], ".pdf"), 
      width = 60, height = 20)
  Heatmap(scale(mm), col=mycol, cluster_columns = FALSE, cluster_rows = TRUE,
          column_title = "Expression heatmap of top 50 regulated genes")
  dev.off()
  
}

sets <- list()
for(ii in 1:length(ttopList)){
  temp <- ttopList[[ii]]
  sets[[length(sets)+1]] <- temp$ID[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
}
DCM_vs_Healthy <- sets[[1]]
HCM_vs_Healthy <- sets[[2]]

## Gene set enrichment analysis
gmt <- GSA.read.gmt(filename = "~/Downloads/human_ontology.gmt")
exP <- gmt$genesets
names(exP) <- gmt$geneset.names

sets <- list()
for(ii in 1:length(ttopList)){
  
  temp <- ttopList[[ii]]
  stats <- temp$logFC
  names(stats) <- temp$ID
  
  fgseaRes <- fgseaSimple(pathways = exP, stats = stats, minSize = 1, maxSize = Inf, nperm = 10000)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  fgseaRes <- fgseaRes[order(fgseaRes$padj, decreasing = FALSE), ]
  write_excel_csv(x = fgseaRes[, 1:7], file = paste0("output_analysis1/all_ontology_sets_", names(ttopList)[ii],  ".xls"), 
                  col_names = TRUE, delim = "\t")
  
  # create a theme for dot plots, which can be reused
  theme_dotplot <- theme_bw(14) +
    theme(axis.text.y = element_text(size = rel(.7)),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = rel(.7)),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.5),
          panel.grid.minor.x = element_blank())
  
  ff <- fgseaRes[order(fgseaRes$padj, decreasing = FALSE), ]
  ff <- ff[intersect(x = which(ff$padj<=0.05), y = which(abs(ff$NES)>=2)), ]
  sets[[length(sets)+1]] <- ff$pathway
  level_order <- ff$pathway
  pdf(file = paste0("output_analysis1/top_ontology_sets_", names(ttopList)[ii],  ".pdf"), width = 20, height = 20)
  pp <- ggplot(ff, aes(x = NES, y = factor(pathway, level_order))) + theme_dotplot
  pp <- pp + geom_point(aes(colour = padj, size = NES)) + xlab(label = "Normalized Enrichment Scores") + ylab(label = "Genesets")
  pp <- pp + ggtitle(paste0("Top Enriched Ontology sets - ", names(ttopList)[ii], " (padj <= 0.05 && abs-NES>=2)"))
  plot(pp)
  dev.off()
  
}

DCM_vs_Healthy <- sets[[1]]
HCM_vs_Healthy <- sets[[2]]

## PROGENy Analysis
data <- as.matrix(df)

genesMapped <- rep("", nrow(data))
for(ii in 1:length(genesMapped)){
  genesMapped[ii] <- MappedID_MAGE_Stringtie$V2[which(MappedID_MAGE_Stringtie$V1==rownames(data)[ii])]
}

idx2rem <- which(genesMapped=="")
data <- data[-idx2rem, ]
genesMapped <- genesMapped[-idx2rem]
uGenes <- unique(genesMapped)

dd <- matrix(data = , nrow = length(uGenes), ncol = ncol(data))
rownames(dd) <- uGenes
colnames(dd) <- colnames(data)

for(ii in 1:nrow(dd)){
  
  idx <- which(genesMapped==uGenes[ii])
  if(length(idx)==1){
    dd[ii, ] <- data[idx, ]
  } else {
    dd[ii, ] <- colMedians(x = data[idx, ])
  }
  
}

PathwayActivity_counts <- progeny(dd, scale=TRUE, organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

idxOrder <- c(which(targets$condition=="Healthy"), 
              which(targets$condition=="DCM", 
              which(targets$condition=="HCM")))
pheatmap(t(PathwayActivity_counts[idxOrder, ]),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy Samplewise Activities", angle_col = NULL,
                         treeheight_col = 0,  border_color = NA)

for(ii in 1:length(ttopList)){
  
  if(ii==1){
    currSample <- "DCM"
  } else {
    currSample <- "HCM"
  }
  
  progeny_sample <- PathwayActivity_counts
  
  idx1 <- which(targets$condition=="Healthy")
  rownames(progeny_sample)[idx1] <- paste0("Healthy_", 1:length(idx1))
  
  idx2 <- which(targets$condition==currSample)
  rownames(progeny_sample)[idx2] <- paste0(currSample, 1:length(idx2))
  
  idx2rem <- setdiff(x = 1:nrow(progeny_sample), y = c(idx1, idx2))
  
  progeny_sample <- progeny_sample[-idx2rem, ]  
    
    bb <- matrix(data = , nrow = nrow(progeny_sample)*ncol(progeny_sample), ncol = 3)
    colnames(bb) <- c("Pathway", "Condition", "Activity")
    cnt <- 1
    for(ll in 1:nrow(progeny_sample)){

      for(jj in 1:ncol(progeny_sample)){

        bb[cnt, 1] <- colnames(progeny_sample)[jj]

        if(grepl(pattern = "Healthy", x = rownames(progeny_sample)[ll], fixed = TRUE)){
          bb[cnt, 2] <- "Healthy"
        } else{
          bb[cnt, 2] <- strsplit(x = names(ttopList)[ii], split = "_", fixed = TRUE)[[1]][1]
        }

        bb[cnt, 3] <- as.numeric(progeny_sample[ll, jj])

        cnt <- cnt + 1
      }

    }
    bb <- as.data.frame(bb)
    bb$Pathway <- as.character(bb$Pathway)
    bb$Condition <- as.character(bb$Condition)
    bb$Activity <- as.numeric(as.character(bb$Activity))

    pdf(file = paste0("output_analysis1/progeny_samplewise_activity_boxplot_", names(ttopList)[ii], ".pdf"), width = 18, height = 6)
    pp <- ggplot(bb, aes(x=Pathway, y=Activity, fill=Condition)) +
      geom_boxplot(position=position_dodge(1))
    plot(pp)
    dev.off()
  
}


for(ii in 1:length(ttopList)){
  
  stats <- ttopList[[ii]]$logFC
  names(stats) <- ttopList[[ii]]$ID
  uID <- unique(names(stats))
  ss <- matrix(data = , nrow = length(uID), ncol = 1)
  rownames(ss) <- uID
  colnames(ss) <- names(ttopList)[ii]
  for(jj in 1:nrow(ss)){
    idx <- which(uID==rownames(ss)[jj])
    ss[jj, 1] <- median(as.numeric(stats[idx]))
  }
  PathwayActivity_zscore <- progeny(ss, scale=FALSE, organism="Human", 
                                    top = 100, perm = 10000, z_scores = TRUE) %>%
    t()
  colnames(PathwayActivity_zscore) <- "NES"
  
  PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))
  
  pdf(file = paste0("output_analysis1/progeny_differential_activities_", names(ttopList)[ii], ".pdf"), 
      width = 12, height = 5)
  ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("Pathways")
  dev.off()
  
  prog_matrix <- getModel("Human", top=100) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")
  
  ss <- ss %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")
  
  scat_plots <- progeny::progenyScatter(df = ss, 
                                        weight_matrix = prog_matrix, 
                                        statName = "t_values", verbose = FALSE)
  
  pdf(file = paste0("output_analysis1/jak_stat_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`JAK-STAT`)
  dev.off()
  
  pdf(file = paste0("output_analysis1/androgen_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`Androgen`)
  dev.off()
  
}

## DoRothEA analysis
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tfList <- list()
for(ii in 1:length(ttopList)){
  
  stats <- ttopList[[ii]]$logFC
  names(stats) <- ttopList[[ii]]$ID
  uID <- unique(names(stats))
  ss <- matrix(data = , nrow = length(uID), ncol = 1)
  rownames(ss) <- uID
  colnames(ss) <- names(ttopList)[ii]
  for(jj in 1:nrow(ss)){
    idx <- which(uID==rownames(ss)[jj])
    ss[jj, 1] <- median(as.numeric(stats[idx]))
  }
  
  tf_activities_stat <- dorothea::run_viper(ss, regulons,
                                            options =  list(minsize = 5, eset.filter = FALSE, 
                                                            cores = 1, verbose = FALSE, nes = TRUE))
  colnames(tf_activities_stat) <- "t"
  
  tmp <- as.matrix(tf_activities_stat[order(abs(tf_activities_stat[, 1]), decreasing = TRUE), ])
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- "t"
  tmp$id <- rownames(tmp)
  tmp <- tmp[, c(2, 1)]
  tmp$t <- as.character(tmp$t)
  write_excel_csv2(x = tmp, 
                   file = paste0("output_analysis1/tf_activities_all_", names(ttopList)[ii], ".xls"),
                   col_names = TRUE)
  
  tf_activities_stat_top25 <- tf_activities_stat %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t") %>%
    dplyr::top_n(25, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))
  
  pdf(file = paste0("output_analysis1/top25_diff_tf_act_", names(ttopList)[ii], ".pdf"), 
      width = 10, height = 7)
  ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
    geom_bar(aes(fill = NES), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "indianred", 
                         mid = "whitesmoke", midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 12),
          axis.text.x = 
            element_text(angle = 45, hjust = 1, size =10, face= "bold"),
          axis.text.y = element_text(size =10, face= "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    xlab("Transcription Factors")
  dev.off()
  
  
  tf_activities_stat_top50 <- tf_activities_stat %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t") %>%
    dplyr::top_n(50, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))
  
  tfList[[length(tfList)+1]] <- tf_activities_stat_top50
  
}
names(tfList) <- names(ttopList)
save(tfList, file = "output_analysis1/tfList50.RData")


# ## Differential TF activity analysis
# regulonsTab <- dorothea_hs
# regulonsTab <- regulonsTab[which(regulonsTab$confidence%in%c("A", "B", "C")), ]
# 
# xx <- as.matrix(regulonsTab[, c(1, 4, 3)])
# xx[, 2] <- gsub(pattern = " ", replacement = "", x = xx[, 2])
# xx[, 2] <- as.numeric(xx[, 2])
# xx <- as.data.frame(xx)
# xx$tf <- as.character(xx$tf)
# xx$mor <- as.numeric(as.character(xx$mor))
# xx$target <- as.character(xx$target)
# regulons <- createRegulonList(regulon_table = xx)
# 
# top50TF <- list()
# for(ii in 1:length(ttopList)){
#   
#   ttop <- ttopList[[ii]]
#   stats <- ttop$logFC
#   names(stats) <- ttop$ID
#   
#   tf_activities <- viper(eset = stats, regulon = regulons, nes = T, method = "none", minsize = 4, eset.filter = F)
#   tf_activities <- as.data.frame(tf_activities)
#   tf_activities$V1 <- as.numeric(as.character(tf_activities$V1))
#   
#   save(tf_activities, file = paste0("output_analysis1/diff_tf_act_", names(ttopList)[ii], ".RData"))
#   
#   temp <- data.frame(y = rownames(tf_activities), x = round(x = tf_activities$V1, digits = 3))
#   pdf(file = paste0("output_analysis1/diff_tf_act_", names(ttopList)[ii], ".pdf"), width = 10, height = 30)
#   p <- ggplot(temp, aes(x, y)) +
#     geom_bar(stat = "identity", aes(fill = x), show.legend = FALSE) +
#     geom_text(aes(label = x))
#   p <- p + ggtitle("TF activitiy scores based on gene regulation levels - All (DoRothEA A, B & C levels)") + labs(x = "Activity scores", y = "Transcription Factor")
#   plot(p)
#   dev.off()
#   
#   dtf <- temp[order(abs(temp$x), decreasing = TRUE)[1:50], ]
#   top50TF[[length(top50TF)+1]] <- dtf$y
#   pdf(file = paste0("output_analysis1/diff_tf_act_top50_", names(ttopList)[ii], ".pdf"), width = 10, height = 15)
#   p <- ggplot(dtf, aes(x, y)) +
#     geom_bar(stat = "identity", aes(fill = x), show.legend = FALSE) +
#     geom_text(aes(label = x))
#   p <- p + ggtitle("TF activitiy scores based on gene regulation levels - Top 50") + labs(x = "Activity scores", y = "Transcription Factor")
#   plot(p)
#   dev.off()
#   
# }
# 
# # universe <- unique(unlist(exP))
# # for(ii in 1:length(ttopList)){
# #   genes <- top50TF[[ii]]
# #   ff <- fora(pathways = exP, genes = genes, minSize = 1, maxSize = Inf, universe = universe)
# #   ff$overlapGenes <- as.character(ff$overlapGenes)
# #   ff$Ratio <- ff$overlap/50
# #   ff <- ff[which(ff$padj<=0.05), ]
# #   
# #   level_order <- ff$pathway
# #   pdf(file = paste0("output_analysis1/top_ora_sets_", names(ttopList)[ii],  ".pdf"), width = 15, height = 30)
# #   pp <- ggplot(ff, aes(x = Ratio, y = factor(pathway, level_order))) + theme_dotplot
# #   pp <- pp + geom_point(aes(colour = padj, size = Ratio)) + xlab(label = "Ratio Members (Overlap/Size)") + ylab(label = "Genesets")
# #   pp <- pp + ggtitle(paste0("Ontology enrichment results - ", names(ttopList)[ii]))
# #   plot(pp)
# #   dev.off()
# #   
# # }
# 
# ## Prepare normalized matrix
# zz <- scale(x = df, center = TRUE, scale = TRUE)
# 
# rr <- c()
# for(ii in 1:nrow(zz)){
#   idx <- which(MappedID_MAGE_Stringtie$V1==rownames(zz)[ii])
#   if(length(idx)>0){
#     rr <- c(rr, MappedID_MAGE_Stringtie$V2[idx])
#   } else {
#     rr <- c(rr, "")
#   }
# }
# 
# idx2rem <- unique(c(which(rr==""), which(duplicated(rr))))
# rr <- rr[-idx2rem]
# zz <- zz[-idx2rem, ]
# 
# rownames(zz) <- rr
# colnames(zz)[which(targets$condition=="Healthy")] <- paste0("Healthy_", 1:length(which(targets$condition=="Healthy")))
# colnames(zz)[which(targets$condition=="DCM")] <- paste0("DCM_", 1:length(which(targets$condition=="DCM")))
# colnames(zz)[which(targets$condition=="HCM")] <- paste0("HCM_", 1:length(which(targets$condition=="HCM")))
# 
# ## PROGENy analysis
# for(ii in 1:length(ttopList)){
#   
#   ttop <- ttopList[[ii]]
#   stats <- ttop$logFC
#   names(stats) <- ttop$ID
#   
#   uNames <- unique(names(stats))
#   ss <- rep(0, length(uNames))
#   
#   for(jj in 1:length(uNames)){
#     ss[jj] <- median(x = stats[which(names(stats)==uNames[jj])], na.rm = TRUE)
#   }
#   
#   stats <- ss
#   names(stats) <- uNames
#   
#   progeny_diff <- progeny(expr = as.matrix(stats), scale = TRUE, organism = "Human", perm = 10000)
#   rownames(progeny_diff) <- c("diff_expr")
#   pa <- as.data.frame(progeny_diff)
#   dtf <- data.frame(x = colnames(progeny_diff), y = progeny_diff[1, ])
#   pdf(file = paste0("output_analysis1/progeny_differential_activities_", names(ttopList)[ii], ".pdf"), width = 12, height = 8)
#   p <- ggplot(dtf, aes(x, y)) +
#     geom_bar(stat = "identity", aes(fill = x), show.legend = FALSE) +
#     geom_text(aes(label = y))
#   plot(p)
#   dev.off()
#   
#   # ss <- strsplit(x = dirs[ii], split = "/", fixed = TRUE)[[1]][3]
#   # ss <- strsplit(x = ss, split = "_", fixed = TRUE)[[1]]
#   # idx1 <- which(tolower(targets$condition)==paste0(ss[1], "_", ss[2]))
#   # idx2 <- which(tolower(targets$condition)==paste0(ss[1], "_", ss[4]))
#   # idx <- c(idx1, idx2)
#   idx1 <- which(targets$condition=="Healthy")
#   idx2 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_", fixed = TRUE)[[1]][1])
#   idx <- c(idx1, idx2)
#   
#   tt <- zz[, idx]
#   progeny_sample <- progeny(expr = tt, scale = TRUE, organism = "Human", perm = 10000)
#   rownames(progeny_sample) <- colnames(zz)[idx]
#   pheatmap(mat = t(progeny_sample), cluster_cols = FALSE, cluster_rows = TRUE, 
#            filename = paste0("output_analysis1/progeny_samplewise_activities_", names(ttopList)[ii], ".pdf"), 
#            main = paste0("Sample-wise pathway activities - ", names(ttopList)[ii]), 
#            width = 52, height = 10)
#   
#   bb <- matrix(data = , nrow = nrow(progeny_sample)*ncol(progeny_sample), ncol = 3)
#   colnames(bb) <- c("Pathway", "Condition", "Activity")
#   cnt <- 1
#   for(ll in 1:nrow(progeny_sample)){
#     
#     for(jj in 1:ncol(progeny_sample)){
#       
#       bb[cnt, 1] <- colnames(progeny_sample)[jj]
#       
#       if(grepl(pattern = "Healthy", x = rownames(progeny_sample)[ll], fixed = TRUE)){
#         bb[cnt, 2] <- "Healthy"
#       } else{
#         bb[cnt, 2] <- strsplit(x = names(ttopList)[ii], split = "_", fixed = TRUE)[[1]][1]
#       }
#       
#       bb[cnt, 3] <- as.numeric(progeny_sample[ll, jj])
#       
#       cnt <- cnt + 1
#     }
#     
#   }
#   bb <- as.data.frame(bb)
#   bb$Pathway <- as.character(bb$Pathway)
#   bb$Condition <- as.character(bb$Condition)
#   bb$Activity <- as.numeric(as.character(bb$Activity))
#   
#   pdf(file = paste0("output_analysis1/progeny_samplewise_activity_boxplot_", names(ttopList)[ii], ".pdf"), width = 18, height = 6)
#   pp <- ggplot(bb, aes(x=Pathway, y=Activity, fill=Condition)) +
#     geom_boxplot(position=position_dodge(1))
#   plot(pp)
#   dev.off()
#   
#   
# }
