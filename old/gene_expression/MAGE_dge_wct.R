#!/usr/bin/Rscript

## DGE using cell type proportion as covariate

library(stringr)
library(dplyr)
library(edgeR)

# use reference data
dataDir <- '/mnt/beegfs/prj/MAGE/analysis/data/stringtie'
MAGNet_data <- load(file.path(dataDir, 'MAGNet_data.RData'))

meta <- as.data.frame(cell.attrs)
# already filtered 
gene_counts <- cts[as.logical(as.numeric(levels(gene.attrs$GeneFlag))[gene.attrs$GeneFlag]),
                   as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag])]
# filter meta after filtering counts!
meta <- meta[as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag]),]
# then drop unused levels
meta$Etiology <- droplevels(meta$Etiology)

# add cell-type ratios
filen <- 'scaden_deconvolution_MAGNet_counts_trained_3set_by_sample_1000cells_2000samples.txt'
filen <- file.path('/mnt/beegfs/prj/MAGE/analysis/deconvolution/scaden/results', filen)
ctr <- read.table(filen)
meta <- cbind(meta, ctr[match(rownames(meta), rownames(ctr)),])

form <- as.formula(paste("~ ", paste(names(meta)[c(3:ncol(meta))], collapse="+")))
design <- model.matrix(form, data=meta)

# since cell type proportions sum to 1, design matrix may not be full column rank
# drop smallest fraction -> Lymphocyte
idx <- ncol(design) - ncol(ctr) + 1
m <- colMeans(design[,c(idx:ncol(design))])
idx <- idx + match(m[which.min(m)], m) - 1
design <- design[,-idx]

print(colnames(design))


y <- DGEList(counts=gene_counts, group=meta$Etiology)
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

# save DGE object
saveRDS(y, file="MAGE_DGEObj.rds")

fit <- glmQLFit(y,design,robust=TRUE)

QlfList<-list();
QlfList[["DCMvsNFD"]] <- glmQLFTest(fit, coef=2) # DCM vs Normal
QlfList[["HCMvsNFD"]] <- glmQLFTest(fit, coef=3)  # HCM vs Normal

contrast <- rep(0, ncol(design))
contrast[2] <- 1
contrast[3] <- -1
QlfList[["DCMvsHCM"]]<- glmQLFTest(fit,contrast=contrast) # DCM vs HCM

# save fit
saveRDS(QlfList, file="MAGE_DGEFit.rds")

# add annotations
library(biomaRt)
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="apr2019.archive.ensembl.org")
resMArt=getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype"), mart=mart)

#also use code for diseases..
#Open Targets database genetic_association_mendelian_somatic_mutation_combined.bed Nov 25 2020

Disease=read.delim("/mnt/beegfs/prj/MAGE/database/gene2disease.txt",header=F,as.is=T)
colnames(Disease)=c("ID","GeneSymbol","Disease")
tmp=tapply(Disease$Disease,Disease$ID,paste,collapse=",")
Disease=data.frame(ID=names(tmp),Disease=tmp)

#RBP - from Hentze paper
library(openxlsx)

rbps=readWorkbook("/mnt/beegfs/prj/MAGE/database/TableS1.xlsx",2)
rbps=as.character(subset(rbps$ID,rbps$found_in_at_least_2_Hs_RIC_studies=="YES" | rbps$knownRBPs=="YES" | rbps$"3470_human_RBPs"=="YES"))
rbps=data.frame(ID=rbps,RBP=rep(1,length(rbps)))

#Some TF paper    
tfs=read.delim("/mnt/beegfs/prj/MAGE/database/annotated_human_TFs_slim.txt",header=F,as.is=T)
tfs=unique(tfs[,2])
tfs=data.frame(ID=tfs,TF=rep(1,length(tfs)))


wb <- createWorkbook()
n<-1;

for(k in names(QlfList))
{
    DfQlf=as.data.frame(topTags(QlfList[[k]],n=nrow(QlfList[[k]])))

    DfQlf=merge(data.frame(ID=rownames(DfQlf),DfQlf),resMArt,by.x=1,by.y=1,all.x=T)
    DfQlf=merge(DfQlf,Disease,by.x=1,by.y=1,all.x=T)
    
    DfQlf=DfQlf[order(DfQlf$PValue),]
    DfQlf=DfQlf[DfQlf$FDR<0.05,]

    DfQlf=merge(DfQlf,rbps,by.x=1,by.y=1,all.x=T)
    DfQlf=merge(DfQlf,tfs,by.x=7,by.y=1,all.x=T)
    
    addWorksheet(wb, sheetName = k);
    writeDataTable(wb, sheet = n, x=DfQlf);
    n=n+1;
}

# adjust output path
saveWorkbook(wb, "MAGE_DGE.xlsx", overwrite = TRUE)
