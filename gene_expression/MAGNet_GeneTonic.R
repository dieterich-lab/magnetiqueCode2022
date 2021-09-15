#! /usr/bin/env Rscript

# prep results for loading into GeneTonic

library("GeneTonic")

library("org.Hs.eg.db")

library("topGO")

## Params

alpha.set <- 0.05
mapping <- "org.Hs.eg.db"
topTablerows <- 500
topGO_method2 <- 'weight01'
p_value_column <- 'p.value_weight01'

## I/O

dirloc <- '/prj/MAGE/analysis/genetonic/results'

dds <- readRDS(file.path(dirloc, 'MAGNet_DESeqDataSet.rds'))
results <- readRDS(file.path(dirloc, 'MAGNet_DESeqResults.rds'))

## Run

purrr::map2(results, names(results), function(.x,.y) {

    anno_df <- data.frame(gene_id=rownames(.x),
                          gene_name=.x$SYMBOL,
                          stringsAsFactors=FALSE,
                          row.names=rownames(.x))
                          
    de <- deseqresult2df(.x, FDR=alpha.set)$SYMBOL
    bg <- deseqresult2df(.x)$SYMBOL
    ont <- lapply(c("BP", "MF", "CC"), function(o) {
        topgo <- pcaExplorer::topGOtable(de,
                                         bg,
                                         ontology=o,
                                         mapping=mapping,
                                         geneID="symbol",
                                         topGO_method2=topGO_method2,
                                         topTablerows=topTablerows)
        res_enrich <- shake_topGOtableResult(topgo, 
                                             p_value_column=p_value_column)
        get_aggrscores(res_enrich=res_enrich,
                       res_de=.x,
                       annotation_obj=anno_df,
                       aggrfun=mean)
    })
    names(ont) <- c("BP", "MF", "CC")
    
    l <- list(dds=dds, res_de=.x, res_enrich=ont, annotation_obj=anno_df, project_id=.y)
    filen <- file.path(dirloc, paste('MAGNet_', .y, '_GeneTonic.rds', sep=''))
    saveRDS(l, file=filen)
    
})
