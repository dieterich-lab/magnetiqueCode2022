library("stringr")
library("edgeR")
library("proBatch")
library("biomaRt")

dir.create("output")

gene_counts <- read.csv("../Data/gene_count_matrix.csv",row.names="gene_id",as.is=T)
colnames(gene_counts)<-gsub("_stringtieRef","",colnames(gene_counts))

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org");
resMArt=getBM(attributes=c("ensembl_gene_id","external_gene_name","pfam"),mart=mart)

# meta<-readRDS("../../../Data/meta.rds")
# 
# dedup=read.table("../../../Data/mage_dedup_metrics_final.txt",as.is=T,header=F)
# dedup=dedup[,c(1,3)]
# dedup[,1]=gsub(".metrics.txt-Unknown","",dedup[,1])
# colnames(dedup)=c("Run","DuplicationRate")
# meta=merge(meta,dedup,by.x=1,by.y=1)
# rownames(meta)<-meta[,1]
MAGNet_data <- load('../Data/MAGNet_data.RData')
meta <- as.data.frame(cell.attrs)
gene_counts <- cts[as.logical(as.numeric(levels(gene.attrs$GeneFlag))[gene.attrs$GeneFlag]),
                   as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag])]
meta <- meta[as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag]),]
meta$Etiology <- droplevels(meta$Etiology)

ins=intersect(rownames(meta),colnames(gene_counts))

meta <- meta[ins, ]
gene_counts <- gene_counts[, ins]

#library(edgeR)
y <- DGEList(counts=gene_counts)

design <- model.matrix(~ meta$Etiology + meta$Race + meta$Sex + meta$Age + meta$DuplicationRate)

keep <- filterByExpr(y, design=design[,2:4])
ensg_keep <- names(keep)[which(keep)]
keep_genes <- unique(resMArt$external_gene_name[which(resMArt$ensembl_gene_id%in%ensg_keep)])

save(keep_genes, file = "keep_genes.RData")

y <- y[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(y)


y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

#Todo save edgeR y object

fit <- glmQLFit(y,design,robust=TRUE)

data <- log_transform_dm(data_matrix = fit$fitted.values, log_base = 2)

tt1 <- c()
tt2 <- c()
tt3 <- c()
for(ii in 1:nrow(data)){
  
  ss1 <- rownames(meta)[which(meta$Etiology=="DCM")]
  ss2 <- rownames(meta)[which(meta$Etiology=="NFD")]
  tt <- t.test(x = as.numeric(data[ii, which(colnames(data)%in%ss1)]),
               y = as.numeric(data[ii, which(colnames(data)%in%ss2)]))
  tt1 <- c(tt1, tt$statistic)
  
  ss1 <- rownames(meta)[which(meta$Etiology=="HCM")]
  ss2 <- rownames(meta)[which(meta$Etiology=="NFD")]
  tt <- t.test(x = as.numeric(data[ii, which(colnames(data)%in%ss1)]),
               y = as.numeric(data[ii, which(colnames(data)%in%ss2)]))
  tt2 <- c(tt2, tt$statistic)
  
  ss1 <- rownames(meta)[which(meta$Etiology=="DCM")]
  ss2 <- rownames(meta)[which(meta$Etiology=="HCM")]
  tt <- t.test(x = as.numeric(data[ii, which(colnames(data)%in%ss1)]),
               y = as.numeric(data[ii, which(colnames(data)%in%ss2)]))
  tt3 <- c(tt3, tt$statistic)
  
}

names(tt1) <- rownames(data)
names(tt2) <- rownames(data)
names(tt3) <- rownames(data)

dgeFit <- readRDS(file = "../Data/MAGE_DGEFit_baseline.rds")
DCM_vs_Healthy <- dgeFit$DCMvsNFD$table
HCM_vs_Healthy <- dgeFit$HCMvsNFD$table
DCM_vs_HCM <- dgeFit$DCMvsHCM$table

library("biomaRt")
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="apr2019.archive.ensembl.org")
resMArt=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description"),mart=mart)

## Fix DCM_vs_Healthy
geneID <- rownames(DCM_vs_Healthy)
mm <- matrix(data = , nrow = length(geneID), ncol = 7)
colnames(mm) <- c(colnames(DCM_vs_Healthy), "geneID", "hgnc_symbol", "t_val")
for(ii in 1:length(geneID)){
  
  mm[ii, 1] <- DCM_vs_Healthy$logFC[ii]
  mm[ii, 2] <- DCM_vs_Healthy$logCPM[ii]
  mm[ii, 3] <- DCM_vs_Healthy$F[ii]
  mm[ii, 4] <- DCM_vs_Healthy$PValue[ii]
  mm[ii, 5] <- geneID[ii]
  mm[ii, 6] <- resMArt$hgnc_symbol[which(resMArt$ensembl_gene_id==geneID[ii])[1]]
  mm[ii, 7] <- as.numeric(tt1[which(names(tt1)==geneID[ii])])
  
}
DCM_vs_Healthy <- as.data.frame(mm)
DCM_vs_Healthy$logFC <- as.numeric(DCM_vs_Healthy$logFC)
DCM_vs_Healthy$logCPM <- as.numeric(DCM_vs_Healthy$logCPM)
DCM_vs_Healthy$F <- as.numeric(DCM_vs_Healthy$F)
DCM_vs_Healthy$PValue <- as.numeric(DCM_vs_Healthy$PValue)
DCM_vs_Healthy$t_val <- as.numeric(DCM_vs_Healthy$t_val)

rownames(DCM_vs_Healthy) <- DCM_vs_Healthy$geneID

DCM_vs_Healthy <- DCM_vs_Healthy[, 1:6]

DCM_vs_Healthy$logFC <- as.numeric(DCM_vs_Healthy$logFC)
DCM_vs_Healthy$logCPM <- as.numeric(DCM_vs_Healthy$logCPM)
DCM_vs_Healthy$F <- as.numeric(DCM_vs_Healthy$F)
DCM_vs_Healthy$PValue <- as.numeric(DCM_vs_Healthy$PValue)

save(DCM_vs_Healthy, file = "output/DCM_vs_Healthy.RData")

## Fix HCM_vs_Healthy
geneID <- rownames(HCM_vs_Healthy)
mm <- matrix(data = , nrow = length(geneID), ncol = 7)
colnames(mm) <- c(colnames(HCM_vs_Healthy), "geneID", "hgnc_symbol", "t_val")
for(ii in 1:length(geneID)){
  
  mm[ii, 1] <- HCM_vs_Healthy$logFC[ii]
  mm[ii, 2] <- HCM_vs_Healthy$logCPM[ii]
  mm[ii, 3] <- HCM_vs_Healthy$F[ii]
  mm[ii, 4] <- HCM_vs_Healthy$PValue[ii]
  mm[ii, 5] <- geneID[ii]
  mm[ii, 6] <- resMArt$hgnc_symbol[which(resMArt$ensembl_gene_id==geneID[ii])[1]]
  mm[ii, 7] <- as.numeric(tt2[which(names(tt2)==geneID[ii])])
  
}
HCM_vs_Healthy <- as.data.frame(mm)
HCM_vs_Healthy$logFC <- as.numeric(HCM_vs_Healthy$logFC)
HCM_vs_Healthy$logCPM <- as.numeric(HCM_vs_Healthy$logCPM)
HCM_vs_Healthy$F <- as.numeric(HCM_vs_Healthy$F)
HCM_vs_Healthy$PValue <- as.numeric(HCM_vs_Healthy$PValue)
HCM_vs_Healthy$t_val <- as.numeric(HCM_vs_Healthy$t_val)

rownames(HCM_vs_Healthy) <- HCM_vs_Healthy$geneID

HCM_vs_Healthy <- HCM_vs_Healthy[, 1:6]

HCM_vs_Healthy$logFC <- as.numeric(HCM_vs_Healthy$logFC)
HCM_vs_Healthy$logCPM <- as.numeric(HCM_vs_Healthy$logCPM)
HCM_vs_Healthy$F <- as.numeric(HCM_vs_Healthy$F)
HCM_vs_Healthy$PValue <- as.numeric(HCM_vs_Healthy$PValue)

save(HCM_vs_Healthy, file = "output/HCM_vs_Healthy.RData")

## Fix DCM_vs_HCM
geneID <- rownames(DCM_vs_HCM)
mm <- matrix(data = , nrow = length(geneID), ncol = 7)
colnames(mm) <- c(colnames(DCM_vs_HCM), "geneID", "hgnc_symbol", "t_val")
for(ii in 1:length(geneID)){
  
  mm[ii, 1] <- DCM_vs_HCM$logFC[ii]
  mm[ii, 2] <- DCM_vs_HCM$logCPM[ii]
  mm[ii, 3] <- DCM_vs_HCM$F[ii]
  mm[ii, 4] <- DCM_vs_HCM$PValue[ii]
  mm[ii, 5] <- geneID[ii]
  mm[ii, 6] <- resMArt$hgnc_symbol[which(resMArt$ensembl_gene_id==geneID[ii])[1]]
  mm[ii, 7] <- as.numeric(tt3[which(names(tt3)==geneID[ii])])
  
}
DCM_vs_HCM <- as.data.frame(mm)
DCM_vs_HCM$logFC <- as.numeric(DCM_vs_HCM$logFC)
DCM_vs_HCM$logCPM <- as.numeric(DCM_vs_HCM$logCPM)
DCM_vs_HCM$F <- as.numeric(DCM_vs_HCM$F)
DCM_vs_HCM$PValue <- as.numeric(DCM_vs_HCM$PValue)
DCM_vs_HCM$t_val <- as.numeric(DCM_vs_HCM$t_val)

rownames(DCM_vs_HCM) <- DCM_vs_HCM$geneID

DCM_vs_HCM <- DCM_vs_HCM[, 1:6]

DCM_vs_HCM$logFC <- as.numeric(DCM_vs_HCM$logFC)
DCM_vs_HCM$logCPM <- as.numeric(DCM_vs_HCM$logCPM)
DCM_vs_HCM$F <- as.numeric(DCM_vs_HCM$F)
DCM_vs_HCM$PValue <- as.numeric(DCM_vs_HCM$PValue)

save(DCM_vs_HCM, file = "output/DCM_vs_HCM.RData")
