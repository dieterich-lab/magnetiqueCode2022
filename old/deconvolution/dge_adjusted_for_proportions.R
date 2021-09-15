#!/biosw/R/4.0.5_deb10/bin/Rscript

## Using estimates for cell type composition, we derive gene sets 
## that are likely to be differentially regulated at the cell type level.

## We fit 2 models of gene expression differences, a standard model and an extended model that includes the cell type proportions.
## Genes that are likely to be differentially regulated at the cell type level are those whose p-values decrease
## when adjusting for cell type composition


library(stringr)
library(dplyr)
library(edgeR)

library(biomaRt)

library(openxlsx)

library(ggplot2)
library(RColorBrewer)


# output
output_dir <- '/prj/MAGE/analysis/deconvolution/scaden/analysis'

all_contrasts <- c("DCMvsNFD", "HCMvsNFD", "DCMvsHCM")
basename_standard <- 'MAGNet_reduced'
basename_extended <- 'MAGNet_full'

# visualization
FDRThreshold <- 0.05

# cell type composition
filen <- 'MAGNet_cell_composition.txt'
filen <- file.path('/prj/MAGE/analysis/deconvolution/scaden/results', filen)
ctr <- read.table(filen)

# use reference data
data_dir <- '/prj/MAGE/analysis/data/stringtie'
MAGNet_data <- load(file.path(data_dir, 'MAGNet_data.RData'))

meta <- as.data.frame(cell.attrs)
# already filtered 
gene_counts <- cts[as.logical(as.numeric(levels(gene.attrs$GeneFlag))[gene.attrs$GeneFlag]),
                   as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag])]
# filter meta after filtering counts!
meta <- meta[as.logical(as.numeric(levels(meta$SampleFlag))[meta$SampleFlag]),]
# then drop unused levels
meta$Etiology <- droplevels(meta$Etiology)

# gene annotation
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                dataset="hsapiens_gene_ensembl",
                host="apr2019.archive.ensembl.org")
resMArt <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype"), 
                 mart=mart)
                 
# redefine the order of the matrix
genes <- merge(data.frame(ensembl_gene_id=rownames(gene_counts)), resMArt, by.x=1, by.y=1, all.x=T)
gene_counts <- gene_counts[match(genes$ensembl_gene_id, rownames(gene_counts)),]


## functions

wrapResults <- function(l, basename) {
    
    # output all predefined fits to spreadsheet
    wb <- createWorkbook()
    purrr::imap(unname(l), function(.x, .y) {
        d <- as.data.frame(topTags(.x, n=nrow(.x)))
        d <- d[order(d$PValue),]
        name <- names(l)[[.y]]
        addWorksheet(wb, sheetName=name)
        writeDataTable(wb, sheet=.y, x=d)
    })
    filen <- file.path(output_dir, paste(basename, ".xlsx", sep=""))
    saveWorkbook(wb, filen, overwrite=FALSE)
}

                 
fitAll <- function(x, group, genes, design, basename){
    
    # fit predefined contrats 
    
    y <- DGEList(counts=x, group=group, genes=genes)
    keep <- filterByExpr(y, group=y$samples$group)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
    y <- estimateGLMTagwiseDisp(y, design)
    
    filen <- file.path(output_dir, paste(basename, "_DGEobj.rds", sep=""))
    saveRDS(y, file=filen)

    fit <- glmQLFit(y, design, robust=TRUE)
    QlfList <- list()
    QlfList[["DCMvsNFD"]] <- glmQLFTest(fit, coef=2) # DCM vs Normal
    QlfList[["HCMvsNFD"]] <- glmQLFTest(fit, coef=3)  # HCM vs Normal

    contrast <- rep(0, ncol(design))
    contrast[2] <- 1
    contrast[3] <- -1
    QlfList[["DCMvsHCM"]]<- glmQLFTest(fit, contrast=contrast) # DCM vs HCM

    filen <- file.path(output_dir, paste(basename, "_DGEfit.rds", sep=""))
    saveRDS(QlfList, file=filen)

    wrapResults(QlfList, basename)
}

# adapted from https://github.com/shenorrLab/bseqsc
pvalueScatter_ <- function(formula, data, pval.th = 0.05, pval.th.y = pval.th, pval.th.x = pval.th, label.th = Inf, title){
  
  # parse formula
  t <- terms(formula, data = data)
  vars <- as.character(attr(t, 'variables'))[-1L]
  data$x <- data[[vars[[2L]]]]
  data$y <- data[[vars[[1L]]]]
  data <- data[data$x <= pval.th.x | data$y <= pval.th.y, , drop = FALSE]
  
  # add layer variables
  data$doLabel <- -log10(data$x) >= label.th & data$x <= data$y
  data$label <- ifelse(data$doLabel, data$hgnc_symbol, '') #<<<<<<<<<
  data$Improved <- data$Improved %||% data$x <= data$y
  
  # plot
  ggplot(data, aes_string(y = 'y', x = 'x', color = 'Improved')) + geom_point(size = .3) + geom_abline(slope = 1) + 
      scale_color_manual(values = c('TRUE' = '#aa2222', 'FALSE' = 'black'), guide = FALSE) + 
      scale_x_continuous(trans = reverselog_trans(10)) + scale_y_continuous(trans = reverselog_trans(10)) +
      geom_text(aes_string(label = 'label'), size = 2, hjust = 0.3, vjust = 1.5) + 
      xlab(paste(vars[[2L]], " FDR", sep="")) + ylab(paste(vars[[1L]], " FDR", sep="")) +
      ggtitle(title) +
      theme_classic() +
      theme(text=element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12))

}


`%||%` <- function(a, b) if( is.null(a) ) b else a


reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
      scales::log_breaks(base = base), 
      domain = c(1e-100, Inf))
}


basisfit <- function(x) x$basisfit


plot_fdr_csfit <- function(x, subset = NULL, alternative = 'two.sided', by = c('cell-type', 'alternative')
                          , xlab='# called', ylab= 'FDR', ylim=c(0,1), ...){
  
  if( length(extra <- list(...)) )
    stop("Unused arguments: ", str_out(names(extra), Inf))
  
  # extract fitted model
  
  fit <- basisfit(x)
  
  if( !is.list(fit) ) stop("Invalid input data: must be a list.")
  if( is.null(rhat <- fit$stat) ){
    warning("Nothing to plot: data does not contain any group comparison data [was csSAMfit run with a group variable in argument `data`?]")
    return(invisible())
  }
  
  ncell <- nrow(rhat)
  cellnames <- rownames(rhat)
  if( is.null(cellnames) ){
    if( is.character(subset) ){
      if( length(subset) > ncell ){
        stop("Invalid number of cell types (", length(subset), "):"
            , " should have at most the same number of cell type data in `x` (", ncell, ")")
      }
      cellnames <- subset
    }else cellnames <- as.character(seq(ncell))
  }
  if( is.null(subset) ) subset <- seq(ncell)
  if( is.character(subset) ){ # partial match the types
    subset <- pmatch(subset, cellnames)
    subset <- subset[!is.na(subset)]
    if( !length(subset) )
      stop("None of the cell types ", str_out(subset), " matched [available: ", str_out(cellnames),"]")
  }
  

  if( is.null(FDR <- fit$FDR) )
    stop("Could not find FDR data in fit object: check if csSAM was run with group data and nperms > 0.")
#            stop("Could not find FDR data in fit object: provide reference data, or, alternatively, run csSAM with group data and nperms > 0.")
  
  alt_set <- names(FDR)
  if( !'all' %in% alternative ){
    ialt <- charmatch(alternative, alt_set, nomatch = 0L)
    FDR <- FDR[ialt]
  }
  if( !length(FDR) ){
    stop("No FDR data found for alternative ", str_out(alternative, Inf)
        , "\n  Available hypothesis: ", str_out(alt_set, Inf))   
  }
  
  df <- plyr::ldply(setNames(subset, cellnames[subset]), function(i){                
        df <- plyr::ldply(FDR, function(x){
              n <- ncol(x$ncall.g)
              j <- which(x$ncall.g[i, ]<=0)
              if( length(j) != n ){
                nc <- rev(x$ncall.g[i, -j])
                ifirst <- !duplicated(nc)
                val <- rev(x$fdr.g[i, -j])[ifirst]
                nc <- nc[ifirst]
                val <- as.numeric(rbind(val, c(val[-1], tail(val, 1L))))
                nc <- as.numeric(rbind(nc, nc))
                
                
                # complete plot with data for first called
                if( nc[1L] != 1 ){
                  nc <- c(1, nc)
                  val <- c(val[1], val) 
                }
                
                data.frame(x = nc, y = val)
              }else
                data.frame(x = seq(1, ncol(rhat), length.out = n), y = rev(x$fdr.g[i,]))
              
            }, .id = 'Alternative')
        df
      }, .id = 'Cell.type')
  
  # order cell type levels as input columns
  df[['Cell.type']] <- factor(df[['Cell.type']], levels = cellnames)
  
  
  # compute breaks in log scale
#  breaks <- .log10_break(df$x)
  
  by <- match.arg(by)
  if( by == 'alternative' ){
    by.colour <- 'Cell.type'
    by.facet <- 'Alternative'
  }else if( by == 'cell-type' ){
    by.colour <- 'Alternative'
    by.facet <- 'Cell.type'
  }
  
  aes_specs <- aes_string(colour = by.colour) 
  
  label_pretty <- as_labeller(function(labels){
    x <- gsub('.', ' ', labels, fixed = TRUE)
    if( any(w <- grepl(' ', x)) ) x[w] <- paste0(toupper(substr(x[w], 1, 1)), substring(x[w], 2))
    x
  })
  # create ggplot plot
  p <- ggplot(df, aes_string(x = 'x', y = 'y')) + 
      geom_line(aes_specs, size=1) +
      scale_y_continuous(ylab, limits=ylim) + 
#      scale_x_continuous(xlab, breaks = breaks) +
      scale_x_log10(labels = scales::label_number(accuracy = 1)) +
      scale_colour_manual(values = brewer.pal(ncell, "Set1")) +
      xlab('Number of regulated genes') +
      theme_minimal() +
      theme(
            legend.position="bottom",
            legend.title = element_blank()) +
      facet_grid(paste0('~ ', by.facet), labeller = label_pretty)  
  
  ltitle <- label_pretty(by.colour)
  p <- p + labs(colour = ltitle, linetype = ltitle)
  p
}


plotEffectSize_ <- function(x){
    # reorder feature levels
    new_feat <- ~ stats::reorder(factor(Feature), rank(as.numeric(factor(Cell.type)) * max(abs(t)) + abs(t)))
    cst <- dplyr::mutate_(x, .dots = setNames(list(new_feat), 'Feature'))
    
    # match color from previous plot
    
    # use csfit...
    # cellnames <- rownames(csfit$basisfit$stat)
    # cst[['Cell.type']] <- factor(cst[['Cell.type']], levels = cellnames)

    cst$Cell.type <- factor(cst$Cell.type, levels=c("Ventricular_cardiomyocyte", "Endothelial", "Pericyte", "Fibroblast", "Smooth_muscle", "Immune"))
    cols <- brewer.pal(length(unique(cst$Cell.type)), "Set1")
    values <-  c("Ventricular_cardiomyocyte"=cols[1], 
                 "Endothelial"=cols[2],
                 "Pericyte"=cols[3],
                 "Fibroblast"=cols[4],
                 "Smooth_muscle"=cols[5],
                 "Immune"=cols[6])
    
    # plot
    ggplot(cst, aes_string(x = 'Feature', y = 't', fill = 'Cell.type')) + 
         geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey30') + 
         ylab('Effect size') + xlab('') + 
         geom_bar(stat = 'identity', position = 'dodge') + coord_flip() +
         scale_fill_manual(values=values) +
         theme_classic() +
         theme(legend.position="bottom",
               legend.title = element_blank(),
               axis.text.y = element_text(size=6)) +
         guides(fill=guide_legend(nrow=2, byrow=TRUE))
}



## run

# standard or baseline model ~Etiology + Race + Age + Sex + DuplicationRate
form <- as.formula(paste("~ ", paste(names(meta)[c(3:ncol(meta))], collapse="+")))
design <- model.matrix(form, data=meta)
# fitAll(gene_counts, meta$Etiology, genes, design, basename_standard)

# extended model - with cell type composition
meta <- cbind(meta, ctr[match(rownames(meta), rownames(ctr)),])
form <- as.formula(paste("~ ", paste(names(meta)[c(3:ncol(meta))], collapse="+")))
design <- model.matrix(form, data=meta)
# we might run into rank issues... so just drop this Unknown cell type
design <- design[,!colnames(design) %in% 'Unknown']

# fitAll(gene_counts, meta$Etiology, genes, design, basename_extended)

# visualize p-values comparison between the 2 models: base model (y-axis) and extended model (x-axis)
lapply(all_contrasts, function(c) {

    filen <- file.path(output_dir, paste(basename_standard, ".xlsx", sep=""))
    base <- read.xlsx(filen, sheet=c)
    base$hgnc_symbol[base$hgnc_symbol==""] <- NA
    base$hgnc_symbol <- ifelse(is.na(base$hgnc_symbol), base$ensembl_gene_id, base$hgnc_symbol)
    base <- base[order(base$FDR),]
    
    filen <- file.path(output_dir, paste(basename_extended, ".xlsx", sep=""))
    adj <- read.xlsx(filen, sheet=c)
    adj$hgnc_symbol[adj$hgnc_symbol==""] <- NA
    adj$hgnc_symbol <- ifelse(is.na(adj$hgnc_symbol), adj$ensembl_gene_id, adj$hgnc_symbol)
    adj <- adj[order(adj$FDR),]
    
    pvals <- cbind(base[, c('hgnc_symbol', 'ensembl_gene_id'), drop=FALSE], 
                   Base=base$FDR,
                   Adjusted=adj$FDR[match(base$hgnc_symbol, adj$hgnc_symbol)],
                   stringsAsFactors=FALSE)
    pvals <- mutate(pvals, Regulated=Adjusted <= FDRThreshold & Adjusted < Base)
    p <- pvalueScatter_(Base ~ Adjusted, pvals, 
                        pval.th=FDRThreshold, 
                        label.th=3.5, 
                        title=paste(c, " FDR<", FDRThreshold, sep=""))
    filen <- file.path(output_dir, paste(c, "_pvals.csv", sep=""))
    write.csv(pvals, file=filen, row.names=F)
    filen <- file.path(output_dir, paste(c, "_pvals_scatter.png", sep=""))
    ggsave(filen, plot=p)
    
})

# cell type-specific differential expression using subsets of selected genes identified above
# we use csSAM (wrapped by bseqsc): for each gene, in each cell type, calculate the contrast in deconvoluted expression between 
# groups, and repeat the procedure with permuted group-label data. Differences in gene expression between two deconvolved cell 
# types (FDR) are calculated as the ratio of genes whose contrast exceeds a given threshold in the real dataset compared with 
# the average number of genes exceeding the same threshold in the permuted dataset. 
# For up and down separately, perform one-tailed tests.

library(bseqsc)

# version 2 (factor effect with more than 2 levels, using lm interaction model), does not work...
# but in all cases model does not handle variables of class [numeric], and does not seem
# to allow for multiple covariates
# so to simplify the interpretation we subsample and run the model 3 times 
# once for each contrast 

# counts
filen <- file.path(output_dir, paste(basename_extended, "_DGEobj.rds", sep=""))
counts <- readRDS(filen)$counts
# lcpm <- cpm(counts, log=TRUE) # use log?? or raw counts?
cpm <- cpm(counts, normalized.lib.sizes=TRUE, log=FALSE)

genes <- genes[match(rownames(cpm), genes$ensembl_gene_id),]
genes$hgnc_symbol[genes$hgnc_symbol==""] <- NA
genes$hgnc_symbol <- ifelse(is.na(genes$hgnc_symbol), genes$ensembl_gene_id, genes$hgnc_symbol)
rownames(cpm) <- genes$hgnc_symbol

cpm <- cpm[!is.na(rownames(cpm)),]
cpm <- cpm[!duplicated(rownames(cpm)),]

# we could have done that above....
ctr <- ctr[,!colnames(ctr) %in% 'Unknown']

lapply(all_contrasts, function(c) {
    
    # pvals from above
    filen <- file.path(output_dir, paste(c, "_pvals.csv", sep=""))
    pvals <- read.csv(filen)
    selected_genes <- subset(pvals, Regulated)$hgnc_symbol
    
    # requires eset?
    # actually, the fitting is not well-documented, and I have a hard time matching what is described
    # in the paper and the source code...
    # ... even when passing explicitely the cell type proportions the results do not seem to change?!
    if( !all(rownames(meta)==colnames(cpm)) ){ stop("Mismatch!") }
    eset <- ExpressionSet(assayData=cpm,
                          phenoData=new("AnnotatedDataFrame", data=meta))
    selected <- eset[rownames(eset@assayData$exprs) %in% selected_genes, ]
    selected <- selected[,selected$Etiology %in% unlist(str_split(c, 'vs'))] 
    #csfit <- bseqsc_csdiff(selected ~ Etiology  | Ventricular_cardiomyocyte + Endothelial + Pericyte + Fibroblast + Smooth_muscle + Immune , 
    #                       verbose=2, 
    #                       nperms=5000, 
    #                       .rng = 12345)
    
    csfit <- bseqsc_csdiff(selected ~ Etiology,
                           data=t(ctr[rownames(ctr) %in% colnames(selected@assayData$exprs),]),
                           verbose=2, 
                           nperms=5000, 
                           .rng = 12345)
    # or 
    #csfit2 <- csSAMfit(selected@assayData$exprs,
    #                   t(ctr[rownames(ctr) %in% colnames(selected@assayData$exprs),]),
    #                   y=selected@phenoData@data$Etiology,
    #                   verbose=2, 
    #                   nperms=5000, 
    #                   .rng = 12345)

    filen <- file.path(output_dir, paste(c, "_csfit.rds", sep=""))
    saveRDS(csfit, file=filen)
    
    
    p <- plot_fdr_csfit(csfit, alt = c('great', 'less'), by = 'alt', xlab = 'Number of significant genes', ylab = 'FDR cutoff')
    p$data$Alternative <- c('greater' = 'Up-regulated', less = 'Down-regulated')[as.character(p$data$Alternative)]
    p <- p + geom_hline(yintercept=FDRThreshold)
    filen <- file.path(output_dir, paste(c, "_fdr.png", sep=""))
    ggsave(filen, plot=p)
    
    # call csSAM - returns cell type-level statistics (expression, etc.)
    # globally filter for genes with the strongest signal, but omit Unknown cells
    # n disable any filtering/ordering? n maximum number of top features to retain -- in each cell type
    # **after** ordering according to arguments `sort.by` and `decreasing`. If passing Inf may not work... and different values of
    # n may lead to more or less results...
    # Default is to return only the cell types for which a `threshold` value was provided. If `NULL`, as default, then no subsetting is performed.
    # use min = minimum FDR across alternatives, default is sort.by = 'FDR'
    # ** default using threshold/subset does not work... 
    cst <- csTopTable(csfit, alt='min', merge=TRUE, n=50, nodups=TRUE)
    # but we just filter afterwards anyway...
    cst <- cst[cst$FDR<0.05,]
    cst <- cst[order(cst$FDR, cst$se),]
    filen <- file.path(output_dir, paste(c, "_ctstats.csv", sep=""))
    write.csv(cst, file=filen, row.names=F)
    # but keep only top for plotting
    cst <- cst %>% group_by(Cell.type) %>% arrange(FDR, se, .by_group = TRUE) %>% slice_head(n=5)
    p <- plotEffectSize_(cst)
    filen <- file.path(output_dir, paste(c, "_effect_size.png", sep=""))
    ggsave(filen, plot=p)
})
