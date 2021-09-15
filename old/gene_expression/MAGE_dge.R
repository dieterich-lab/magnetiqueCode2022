#!/biosw/R/4.0.5_deb10/bin/Rscript

library(stringr)
library(edgeR)

# use reference data
dataDir <- '/prj/MAGE/analysis/data/stringtie'
MAGNet_data <- load(file.path(dataDir, 'MAGNet_data.RData'))

meta <- as.data.frame(cell.attrs)
# already filtered 
gene_counts <- cts[as.logical(as.numeric(levels(gene.attrs$GeneFlag))[gene.attrs$GeneFlag]),
                   as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag])]
# filter meta after filtering counts!
meta <- meta[as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag]),]
# then drop unused levels
meta$Etiology <- droplevels(meta$Etiology)


design <- model.matrix(~ meta$Etiology + meta$Race + meta$Sex + meta$Age + meta$DuplicationRate)

y <- DGEList(counts=gene_counts, group=meta$Etiology)
keep <- filterByExpr(y, group=y$samples$group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTagwiseDisp(y, design)

# save DGE object
saveRDS(y, file="MAGEbaseline_DGEobj.rds")

fit <- glmQLFit(y,design,robust=TRUE)

QlfList<-list();
QlfList[["DCMvsNFD"]] <- glmQLFTest(fit, coef=2) # DCM vs Normal
QlfList[["HCMvsNFD"]] <- glmQLFTest(fit, coef=3)  # HCM vs Normal
QlfList[["DCMvsHCM"]]<- glmQLFTest(fit,contrast=c(0,1,-1,0,0,0,0)) # DCM vs HCM

# save fit
saveRDS(QlfList, file="MAGEbaseline_DGEfit.rds")

# add annotations
library(biomaRt)
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="apr2019.archive.ensembl.org")
resMArt=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","gene_biotype","description"),mart=mart)

#also use code for diseases..
#Open Targets database genetic_association_mendelian_somatic_mutation_combined.bed Nov 25 2020

Disease=read.delim("/prj/MAGE/database/gene2disease.txt",header=F,as.is=T)
colnames(Disease)=c("ID","GeneSymbol","Disease")
tmp=tapply(Disease$Disease,Disease$ID,paste,collapse=",")
Disease=data.frame(ID=names(tmp),Disease=tmp)

#RBP - from Hentze paper
library(openxlsx)

rbps=readWorkbook("/prj/MAGE/database/TableS1.xlsx",2)
rbps=as.character(subset(rbps$ID,rbps$found_in_at_least_2_Hs_RIC_studies=="YES" | rbps$knownRBPs=="YES" | rbps$"3470_human_RBPs"=="YES"))
rbps=data.frame(ID=rbps,RBP=rep(1,length(rbps)))

#Some TF paper    
tfs=read.delim("/prj/MAGE/database/annotated_human_TFs_slim.txt",header=F,as.is=T)
tfs=unique(tfs[,2])
tfs=data.frame(ID=tfs,TF=rep(1,length(tfs)))


wb <- createWorkbook()
n<-1;

for(k in names(QlfList))
{
    DfQlf=as.data.frame(topTags(QlfList[[k]],n=dim(QlfList[[k]])[1]))

    DfQlf=merge(data.frame(ID=rownames(DfQlf),DfQlf),resMArt,by.x=1,by.y=1,all.x=T)
    DfQlf=merge(DfQlf,Disease,by.x=1,by.y=1,all.x=T)
    
    DfQlf=DfQlf[order(DfQlf$PValue),]
    # DfQlf=DfQlf[DfQlf$FDR<0.05,] 

    DfQlf=merge(DfQlf,rbps,by.x=1,by.y=1,all.x=T)
    DfQlf=merge(DfQlf,tfs,by.x=7,by.y=1,all.x=T)
    
    addWorksheet(wb, sheetName = k);
    writeDataTable(wb, sheet = n, x=DfQlf);
    n=n+1;
}

# adjust output path
saveWorkbook(wb, "/prj/MAGE/analysis/dge/MAGEbaseline_DGE.xlsx", overwrite = TRUE)
