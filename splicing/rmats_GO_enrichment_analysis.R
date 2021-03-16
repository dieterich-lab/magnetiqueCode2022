suppressPackageStartupMessages({
  library(CellPlot)
  library(topGO)
  library(dplyr)
  library(stringr)
  library(org.Hs.eg.db )
  library(readr)
})

setwd('/prj/MAGE/analysis/baltica/')
files  <- Sys.glob('rmats/*/*.MATS.JC.txt')
x <- lapply(files, readr::read_delim, '\t')
names(x) <- gsub(files, pattern = "rmats/([HD]CM-vs-CTRL)/(.*).MATS.JC.txt", replacement = '\\1_\\2')

x2 <- bind_rows(x, .id='type')
x2 <- x2 %>% tidyr::separate(type, into=c('comparison','as_type'), sep ='_')


run.topgo <- function(df, comparison, mapping='org.Hs.eg.db', desc = '', ont='MF') {
  df <- x2 %>% 
    filter(comparison == comparison) %>% 
    group_by(GeneID) %>% 
    dplyr::select(GeneID, FDR, IncLevelDifference) %>% 
    slice(which.min(FDR))
  
  fg <- unique(subset(df, FDR < 0.05)[['GeneID']])
  bg <- unique(df[['GeneID']])
  in_universe <- bg %in% bg
  in_selection <- bg %in% fg
  
  all_genes <- factor(as.integer(in_selection[in_universe]))
  names(all_genes) <- bg
  
  g <- new(
    "topGOdata",
    ontology = as.character(ont),
    description = '',
    allGenes = all_genes,
    mapping = mapping, # TODO generalise 
    annotationFun = annFUN.org,
    ID = "Ensembl"
  )
  t <- new("elimCount",
           testStatistic = GOFisherTest,
           name = "Fisher test") # test definition
  s <- getSigGroups(g, t) # run F-test
  r <- GenTable(g,
                pvalCutOff = s,
                topNodes = length(g@graph@nodes)) # return data.frame
  
  r$pvalCutOff <- as.numeric(str_replace_all(r$pvalCutOff, "[^0-9e\\-\\.]*", ""))
  r$LogEnrich <- log2(r$Significant / r$Expected)

  
  anno <- stack(genesInTerm(g))
  colnames(anno) <- c('GeneID', 'term')
  anno$term <- as.character(anno$term)
  anno <- subset(anno, term %in% r$GO.ID)
  anno$gene_significant <- anno$GeneID %in% fg
  anno <- left_join(anno, df, ) %>%
    distinct()

  

  write.csv(r, str_glue('topgo_{ont}.csv'))
  write.csv(anno, str_glue('GO_annotation_{ont}.csv'))
  
  r <- subset(r, pvalCutOff <= 0.05 & Significant > 5)
  r <- head(r[order(-r$LogEnrich),], n=20)
  
  try({
    cell.plot.x <- setNames(r$LogEnrich, r$Term)
    cell.plot.cells <- split(anno$IncLevelDifference, anno$term)
    png(filename = paste(ont, 'cellplot.png', sep = '_'), width = 923, height = 508)
    cell.plot(x = cell.plot.x,
              #lapply(na.omit(r$padj), function(x) -log10(x)),
              cells = cell.plot.cells[r$GO.ID],
              key.lab = "IncLevelDifference",
              main = str_glue('{ont} enrichment'),
              x.mar = c(.3, .01),
              key.n = 7,
              y.mar = c(.1, -.1),
              cex = 1.6,
              cell.outer = 1,
              bar.scale = .2,
              space = .1)
    dev.off()
  })
  
}
run.topgo(x2, ont='MF', comparison = 'HCM-vs-CTRL')
run.topgo(x2, ont='CC', comparison = 'HCM-vs-CTRL')
run.topgo(x2, ont='BP', comparison = 'HCM-vs-CTRL')
