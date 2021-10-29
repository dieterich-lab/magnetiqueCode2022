#! /usr/bin/env Rscript

#
set.seed(66901)

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
library("topGO")
library("CellPlot")
library("openxlsx")

source("cellplot_functions.R")

load(file = "output/ttopList.RData")

cases <- names(ttopList)

for(ii in 1:length(ttopList)){
  
  ttop <- ttopList[[ii]]
  
  ttop$ID <- rownames(ttop)
  
  fgSet<-subset(ttop$ID,ttop$padj<=0.05)
  bgSet<-unique(ttop$ID);
  
  selec<-numeric(length=length(bgSet))
  names(selec)<-bgSet
  
  selec[fgSet]<-1
  
  testGoBP<-new("topGOdata",description=cases[ii],ontology="BP",
                allGenes=factor(selec),nodeSize=5, annot=annFUN.org, 
                mapping = "org.Hs.eg.db", ID = "ensembl")
  
  test.stat1 <- new("elimCount", testStatistic = GOFisherTest,cutOff=0.05)
  resultElimFis <- getSigGroups(testGoBP, test.stat1)
  
  allResBP<-GenTable(testGoBP, elimFisher = resultElimFis, topNodes=100)
  allResBP[,6]<-as.numeric(allResBP[,6])
  
  write.csv(subset(allResBP,allResBP[,6]<=0.05),file=paste0("output/diffGenes_", tolower(cases[ii]), ".csv"))
  ga<-getGenesFromGOandExpression(paste0("output/diffGenes_", tolower(cases[ii]), ".csv"),testGoBP,selec);
  all.ga<-getGenesFromGOandCuff(paste0("output/diffGenes_", tolower(cases[ii]), ".csv"),testGoBP,selec);
  
  UpDownCellPlot=lapply(sapply(ga,function(x){subset(as.numeric(ttop[,"log2FoldChange"]),ttop[,"ID"] %in% x)}),
                        function(x){table(ifelse(x>0,"Upregulated","Downregulated"))})
  tmp=data.frame(ID=names(UpDownCellPlot),Downregulated=0,Upregulated=0)
  rownames(tmp)<-names(UpDownCellPlot)
  for(k in names(UpDownCellPlot))
  {
    tmp2=UpDownCellPlot[[k]]
    for(t in names(tmp2))
    {
      tmp[k,t]<-tmp2[t]
    }
  }
  UpDownCellPlot=tmp;
  log2foldChangeCellPlot=sapply(ga,function(x){subset(ttop[,"log2FoldChange"],ttop[,"ID"] %in% x)})
  log2foldChangeCellPlot=data.frame(ID=names(log2foldChangeCellPlot),fc=I(as.list(log2foldChangeCellPlot)))
  padjChangeCellPlot=sapply(ga,function(x){subset(ttop[,"padj"],ttop[,"ID"] %in% x)})
  padjChangeCellPlot=data.frame(ID=names(padjChangeCellPlot),fc=I(as.list(padjChangeCellPlot)))
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = paste0("GO terms BP - ", paste0(strsplit(x = cases[ii], split = "_", fixed = TRUE)[[1]], collapse = " ")));
  writeDataTable(wb, sheet = 1, x=merge(allResBP,data.frame(ID=names(ga),Genes=sapply(ga,paste,collapse=", ")),by.x=1,by.y=1));
  saveWorkbook(wb, file.path(paste0('output/DGE_', cases[ii], '.xlsx')), overwrite = TRUE)
  
  require(CellPlot)
  forCellPlot=readWorkbook(file.path(paste0('output/DGE_', cases[ii], '.xlsx')),1)
  forCellPlot=merge(forCellPlot,log2foldChangeCellPlot,by.x=1,by.y=1)
  forCellPlot=forCellPlot[order(forCellPlot$elimFisher),]
  
  backup=forCellPlot
  forCellPlot=forCellPlot[1:10,]
  #pdf(paste0("CellPlot_",labelCond,".pdf"))
  labelCond=cases[ii]
  pdf(file = file.path(paste0('output/DGE_', cases[ii], '.pdf')), width = 8, height = 8)
  cell.plot(x = setNames( -log10(forCellPlot[,"elimFisher"]), forCellPlot[,"Term"]), 
            cells = sapply(sapply(sapply(forCellPlot$fc,as.character),strsplit,","),as.numeric), 
            main =paste0("GO BP enrichment (CellPlot_",labelCond,")"), 
            x.mar = c(.4, 0), 
            key.n = 7, 
            y.mar = c(.1, 0), 
            cex = 1.6, 
            cell.outer = 3, 
            bar.scale = .7, 
            space = .2,
            xlab = "GO Term Enrichment (-log10 pvalue)",
            key.lab = "Differential Expression (log2 fold change)")
  dev.off();
  
}
  