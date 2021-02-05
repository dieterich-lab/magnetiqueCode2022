setwd("/Volumes/prj/MAGE/Christoph/stringtie")
#read stringtie gene expression data

gene_counts <- read.csv("gene_count_matrix.csv",row.names="gene_id",as.is=T)
colnames(gene_counts)<-gsub("_stringtieRef","",colnames(gene_counts))

library(stringr)
library(edgeR)

meta<-readRDS("/Volumes/prj/MAGE/Christoph/meta.rds")

dedup=read.table("/Volumes/prj/MAGE/Christoph/mage_dedup_metrics_final.txt",as.is=T,header=F)
dedup=dedup[,c(1,3)]
dedup[,1]=gsub(".metrics.txt-Unknown","",dedup[,1])
colnames(dedup)=c("Run","DuplicationRate")
meta=merge(meta,dedup,by.x=1,by.y=1)
rownames(meta)<-meta[,1]

ins=intersect(rownames(meta),colnames(gene_counts))

meta <- meta[ins, ]
gene_counts <- gene_counts[, ins]

#library(edgeR)
y <- DGEList(counts=gene_counts)

#time=c("4","4","4","4","4","4","12","12","12","12","12","12","12","12")
#genotype=c("wt","wt","wt","het_ko","het_ko","het_ko","wt","wt","wt","het_ko","het_ko","het_ko","het_ko","het_ko")

design <- model.matrix(~ meta$etiology + meta$race + meta$sex + meta$Age + meta$DuplicationRate)
design = design[,-c(5,7,9)]

keep <- filterByExpr(y, design=design[,2:4])
y <- y[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(y)

#y <- estimateDisp(y)
                                        #Zwichenspeichern

y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
#y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

                                        #TODO: save y object

fit <- glmQLFit(y,design,robust=TRUE)
                                        #plotQLDisp(fit)
QlfList<-list();

QlfList[["DCMvsNormal"]] <- glmQLFTest(fit,contrast=c(0,1,0,-1,0,0,0,0)) #DCM vs Normal
QlfList[["DCMvsHCM"]] <- glmQLFTest(fit,contrast=c(0,1,-1,0,0,0,0,0)) #DCM vs HCM
QlfList[["HCMvsNormal"]]<- glmQLFTest(fit,contrast=c(0,0,1,-1,0,0,0,0)) #HCM vs Normal

                                        #add annotation
library(biomaRt)
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="apr2019.archive.ensembl.org")
resMArt=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","description"),mart=mart)

#also use code for diseases..

Disease=read.delim("/Volumes//prj/MAGE/database/gene2disease.txt",header=F,as.is=T)
colnames(Disease)=c("ID","GeneSymbol","Disease")
tmp=tapply(Disease$Disease,Disease$ID,paste,collapse=",")
Disease=data.frame(ID=names(tmp),Disease=tmp)


#RBP - redo with Netze stuiff
library(openxlsx)

rbps=readWorkbook("/Volumes//prj/MAGE/database/TableS1.xlsx",2)
rbps=as.character(subset(rbps$ID,rbps$found_in_at_least_2_Hs_RIC_studies=="YES" | rbps$knownRBPs=="YES" | rbps$"3470_human_RBPs"=="YES"))
rbps=data.frame(ID=rbps,RBP=rep(1,length(rbps)))
#Mmmh ???
    
tfs=read.delim("~/Downloads/annotated_human_TFs_slim.txt",header=F,as.is=T)
tfs=unique(tfs[,2])
tfs=data.frame(ID=tfs,TF=rep(1,length(tfs)))

#exit(0);
library(openxlsx)
wb <- createWorkbook()
n<-1;

for(k in names(QlfList))
{
    DfQlf=as.data.frame(topTags(QlfList[[k]],n=22360))

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

saveWorkbook(wb, "MAGE_DGE_Nov2020.xlsx", overwrite = TRUE)

exit(0);

                                        #RBPs

rbpsDGE=subset(DfQlf,DfQlf$hgnc_symbol %in% rbps);

#TFs

todo=as.character(DfQlf[grep("MEF2C-AS1",DfQlf$hgnc_symbol),1])

for(goI in todo)
{
idx=grep(goI,as.character(DfQlf[,1]))
countTab=cpm(y)[goI,]
                                        #cool.. but should remove outliers
#pdf(paste0("/Volumes//prj/MAGE/Christoph/plots/",DfQlf[idx,"hgnc_symbol"],".pdf"))
boxplot(list(DCM=countTab[which(design[,2]==1)],HCM=countTab[which(design[,3]==1)],Normals=countTab[which(design[,4]==1)]),ylab="cpm",xlab="Etiology",main=paste0("Gene: ",DfQlf[idx,"hgnc_symbol"]," FDR: ",DfQlf[idx,"FDR"]))
#dev.off()
}

