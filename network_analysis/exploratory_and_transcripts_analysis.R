#
set.seed(1234)

#
dir.create("output")

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
MappedID_MAGE_Stringtie <- read.delim("~/Desktop/MAGE-Project/Data/MappedID_MAGE_Stringtie.txt", header=FALSE)

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

DGE <- DGEList(counts) %>%  calcNormFactors() 
data <- voom(DGE)
df <- data$E

df <- as.data.frame(df)
data <- df

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

#
DCM_vs_Healthy <- read_excel("output/DCM_vs_Healthy.xlsx")
HCM_vs_Healthy <- read_excel("output/HCM_vs_Healthy.xlsx")
DCM_vs_HCM <- read_excel("output/DCM_vs_HCM.xlsx")

for(ii in 1:nrow(DCM_vs_Healthy)){
  
  idx <- which(MappedID_MAGE_Stringtie$V1==DCM_vs_Healthy$ID[ii])
  if(length(idx)==1){
    DCM_vs_Healthy$hgnc_symbol[ii] <- MappedID_MAGE_Stringtie$V2[idx]
  }
  
}

for(ii in 1:nrow(DCM_vs_HCM)){
  
  idx <- which(MappedID_MAGE_Stringtie$V1==DCM_vs_HCM$ID[ii])
  if(length(idx)==1){
    DCM_vs_HCM$hgnc_symbol[ii] <- MappedID_MAGE_Stringtie$V2[idx]
  }
  
}

for(ii in 1:nrow(HCM_vs_Healthy)){
  
  idx <- which(MappedID_MAGE_Stringtie$V1==HCM_vs_Healthy$ID[ii])
  if(length(idx)==1){
    HCM_vs_Healthy$hgnc_symbol[ii] <- MappedID_MAGE_Stringtie$V2[idx]
  }
  
}

ttopList <- list()
ttopList[[length(ttopList)+1]] <- DCM_vs_Healthy
ttopList[[length(ttopList)+1]] <- HCM_vs_Healthy
ttopList[[length(ttopList)+1]] <- DCM_vs_HCM

names(ttopList) <- c("DCM_vs_Healthy", "HCM_vs_Healthy", "DCM_vs_HCM")

for(ii in 1:length(ttopList)){
  
  temp <- ttopList[[ii]]
  temp$hgnc_symbol <- as.character(temp$hgnc_symbol)
  idx <- intersect(x = which(abs(temp$logFC)>2), y = which(temp$FDR<0.05))
  idx1 <- intersect(x = which(temp$logFC>2), y = which(temp$FDR<0.05))
  idx2 <- intersect(x = which(temp$logFC < -2), y = which(temp$FDR<0.05))
  pdf(file = paste0("output/volcano_", names(ttopList)[ii], ".pdf"), width = 12, height = 8)
  EnhancedVolcano(toptable = temp, lab = temp$hgnc_symbol, x = "logFC", y = "FDR", pCutoff = 0.05, FCcutoff = 1,
                  pointSize = 3, title = "Volcano plot of differentially expressed genes",
                  subtitle = paste0(length(idx), " significantly enriched genes (pCutoff = 0.05, Log-FCcutoff = 2): ", length(idx1), " up & ", length(idx2), " down"))
  dev.off()
  
  idx1 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][1])
  idx2 <- which(targets$condition==strsplit(x = names(ttopList)[ii], split = "_vs_", fixed = TRUE)[[1]][2])
  # top50 <- temp$ID[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
  top50 <- temp$hgnc_symbol[order(abs(temp$logFC), decreasing = TRUE)[1:50]]
  select <- which(rownames(data)%in%top50)
  mycol <- colorRamp2(c(-3,0,3), c("dodgerblue", "black", "yellow"))
  mm <- as.matrix(data)[select, c(idx1, idx2)]
  pdf(file = paste0("output/top50_expr_genes_scaled_", names(ttopList)[ii], ".pdf"), 
      width = 60, height = 20)
  Heatmap(scale(mm), col=mycol, cluster_columns = FALSE, cluster_rows = TRUE,
          column_title = "Expression heatmap of top 50 regulated genes")
  dev.off()
  
}

sets <- list()
for(ii in 1:length(ttopList)){
  temp <- ttopList[[ii]]
  sets[[length(sets)+1]] <- as.character(temp$hgnc_symbol[order(abs(temp$logFC), decreasing = TRUE)[1:50]])
}
DCM_vs_Healthy <- sets[[1]]
HCM_vs_Healthy <- sets[[2]]
DCM_vs_HCM <- sets[[3]]

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(DCM_vs_Healthy, HCM_vs_Healthy, DCM_vs_HCM),
  category.names = c("DCM vs Healthy" , "HCM vs Healthy" , "DCM vs HCM"),
  filename = 'output/venn_diagramm_top50_genes.png',
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


## Gene set enrichment analysis
gmt <- GSA.read.gmt(filename = "~/Downloads/human_ontology.gmt")
exP <- gmt$genesets
names(exP) <- gmt$geneset.names

sets <- list()
for(ii in 1:length(ttopList)){
  
  temp <- ttopList[[ii]]
  stats <- temp$logFC
  names(stats) <- temp$hgnc_symbol
  
  fgseaRes <- fgseaSimple(pathways = exP, stats = stats, minSize = 1, maxSize = Inf, nperm = 10000)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  fgseaRes <- fgseaRes[order(fgseaRes$padj, decreasing = FALSE), ]
  write_excel_csv(x = fgseaRes[, 1:7], file = paste0("output/all_ontology_sets_", names(ttopList)[ii],  ".xls"), 
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
  pdf(file = paste0("output/top_ontology_sets_", names(ttopList)[ii],  ".pdf"), width = 20, height = 20)
  pp <- ggplot(ff, aes(x = NES, y = factor(pathway, level_order))) + theme_dotplot
  pp <- pp + geom_point(aes(colour = padj, size = NES)) + xlab(label = "Normalized Enrichment Scores") + ylab(label = "Genesets")
  pp <- pp + ggtitle(paste0("Top Enriched Ontology sets - ", names(ttopList)[ii], " (padj <= 0.05 && abs-NES>=2)"))
  plot(pp)
  dev.off()
  
}

## PROGENy Analysis
df <- data
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

PathwayActivity_counts <- progeny(dd, scale=TRUE, organism="Human", top = 200)
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

ttList <- list()
for(ii in 1:length(ttopList)){
  
  tt <- c()
  
  ss <- strsplit(x = names(ttopList)[ii], split = "_vs_")[[1]]
  
  progeny_sample <- PathwayActivity_counts
  
  idx1 <- which(targets$condition==ss[2])
  rownames(progeny_sample)[idx1] <- paste0(ss[2], 1:length(idx1))
  
  idx2 <- which(targets$condition==ss[1])
  rownames(progeny_sample)[idx2] <- paste0(ss[1], 1:length(idx2))
  
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
          if(grepl(pattern = "DCM", x = rownames(progeny_sample)[ll], fixed = TRUE)){
            bb[cnt, 2] <- "DCM"
          } else {
            bb[cnt, 2] <- "HCM"
          }
        }

        bb[cnt, 3] <- as.numeric(progeny_sample[ll, jj])

        cnt <- cnt + 1
      }

    }
    bb <- as.data.frame(bb)
    bb$Pathway <- as.character(bb$Pathway)
    bb$Condition <- as.character(bb$Condition)
    bb$Activity <- as.numeric(as.character(bb$Activity))

    pdf(file = paste0("output/progeny_samplewise_activity_boxplot_", names(ttopList)[ii], ".pdf"), width = 18, height = 6)
    pp <- ggplot(bb, aes(x=Pathway, y=Activity, fill=Condition)) +
      geom_boxplot(position=position_dodge(1))
    plot(pp)
    dev.off()
    
    px <- c("Androgen", "EGFR", "Estrogen", "Hypoxia",  "JAK-STAT", "MAPK", "NFkB", 
            "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
    
    for(jj in 1:length(px)){
      
      act1 <- bb$Activity[intersect(x = which(bb$Pathway==px[jj]), y = which(bb$Condition==ss[1]))]
      act2 <- bb$Activity[intersect(x = which(bb$Pathway==px[jj]), y = which(bb$Condition==ss[2]))]
      
      p <- t.test(x = act1, y = act2)
      tt <- c(tt, p$p.value)
      
    }
    
    names(tt) <- px
    ttList[[length(ttList)+1]] <- tt
  
}
names(ttList) <- names(ttopList)

progenyDFList <- list()
for(ii in 1:length(ttopList)){
  
  stats <- ttopList[[ii]]$logFC
  names(stats) <- ttopList[[ii]]$hgnc_symbol
  uID <- unique(names(stats))
  ss <- matrix(data = , nrow = length(uID), ncol = 1)
  rownames(ss) <- uID
  colnames(ss) <- names(ttopList)[ii]
  for(jj in 1:nrow(ss)){
    idx <- which(uID==rownames(ss)[jj])
    ss[jj, 1] <- median(as.numeric(stats[idx]))
  }
  PathwayActivity_zscore <- progeny(ss, scale=FALSE, organism="Human", 
                                    top = 200, perm = 10000, z_scores = TRUE) %>%
    t()
  colnames(PathwayActivity_zscore) <- "NES"
  
  PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
    rownames_to_column(var = "Pathway") %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Pathway = factor(Pathway))
  
  progenyDFList[[length(progenyDFList)+1]] <- PathwayActivity_zscore_df
  
  pdf(file = paste0("output/progeny_differential_activities_", names(ttopList)[ii], ".pdf"), 
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
  
  prog_matrix <- getModel("Human", top=200) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")
  
  ss <- ss %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")
  
  scat_plots <- progeny::progenyScatter(df = ss, 
                                        weight_matrix = prog_matrix, 
                                        statName = "t_values", verbose = FALSE)
  
  pdf(file = paste0("output/jak_stat_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`JAK-STAT`)
  dev.off()
  
  pdf(file = paste0("output/hypoxia_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`Hypoxia`)
  dev.off()
  
  pdf(file = paste0("output/egfr_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`EGFR`)
  dev.off()
  
  pdf(file = paste0("output/androgen_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`Androgen`)
  dev.off()
  
  pdf(file = paste0("output/estrogen_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`Estrogen`)
  dev.off()
  
  pdf(file = paste0("output/tnfa_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`TNFa`)
  dev.off()
  
  pdf(file = paste0("output/wnt_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`WNT`)
  dev.off()
  
  pdf(file = paste0("output/trail_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`Trail`)
  dev.off()
  
  pdf(file = paste0("output/tgfb_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`TGFb`)
  dev.off()
  
  pdf(file = paste0("output/mapk_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`MAPK`)
  dev.off()
  
  pdf(file = paste0("output/vegf_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`VEGF`)
  dev.off()
  
  pdf(file = paste0("output/nfkb_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`NFkB`)
  dev.off()
  
  pdf(file = paste0("output/p53_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`p53`)
  dev.off()
  
  pdf(file = paste0("output/pi3k_most_responsive_genes_", tolower(names(ttopList)[ii]), ".pdf"), width = 10, height = 10)
  plot(scat_plots[[1]]$`PI3K`)
  dev.off()
  
}
names(progenyDFList) <- names(ttopList)

## volcano progeny
temp <- matrix(data = , nrow = length(progenyDFList)*nrow(progenyDFList$DCM_vs_Healthy), ncol = 3)
cnt <- 1
for(ii in 1:length(progenyDFList)){
  
  currDF <- progenyDFList[[ii]]
  currDF$Pathway <- as.character(currDF$Pathway)
  currDF$NES <- as.numeric(as.character(currDF$NES))
  currTT <- ttList[[ii]]
  
  currDF <- currDF[order(as.character(currDF$Pathway), decreasing = TRUE), ]
  currTT <- currTT[order(names(currTT), decreasing = TRUE)]
  
  for(jj in 1:length(currTT)){
    
    temp[cnt, 1] <- paste0(names(currTT)[jj], "_", names(progenyDFList)[ii])
    temp[cnt, 2] <- currDF$NES[jj]
    temp[cnt, 3] <- as.numeric(currTT)[jj]
    
    cnt <- cnt + 1
    
  }
  
}
temp <- as.data.frame(temp)
colnames(temp) <- c("case", "nes", "pval")
temp$nes <- as.numeric(temp$nes)
temp$pval <- as.numeric(temp$pval)

idx <- intersect(x = which(abs(temp$nes)>1), y = which(temp$pval<0.1))
idx1 <- intersect(x = which(temp$nes>1), y = which(temp$pval<0.1))
idx2 <- intersect(x = which(temp$nes < -1), y = which(temp$pval<0.1))
pdf(file = paste0("output/progeny_deregulated_pathways.pdf"), width = 12, height = 8)
EnhancedVolcano(toptable = temp, lab = temp$case, x = "nes", y = "pval", pCutoff = 0.1, FCcutoff = 1,
                pointSize = 3, title = "Volcano plot of PROGENy scores",
                subtitle = paste0(length(idx), " significantly regulated pathways (pCutoff = 0.1, NES = 1): ", length(idx1), " up & ", length(idx2), " down"), 
                xlab = "NES", ylim = c(0, 5), labSize = 1.5)
dev.off()

## DoRothEA analysis
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tfList <- list()
for(ii in 1:length(ttopList)){
  
  stats <- ttopList[[ii]]$logFC
  names(stats) <- ttopList[[ii]]$hgnc_symbol
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
                   file = paste0("output/tf_activities_all_", names(ttopList)[ii], ".xls"),
                   col_names = TRUE)
  
  tf_activities_stat_top25 <- tf_activities_stat %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "t") %>%
    dplyr::top_n(50, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))
  
  pdf(file = paste0("output/top50_diff_tf_act_", names(ttopList)[ii], ".pdf"), 
      width = 12, height = 7)
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
save(tfList, file = "output/tfList50.RData")
