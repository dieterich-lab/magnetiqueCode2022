create_variables <- function(background.network = background.network){
  
  #
  cc1 <- paste0("xb", 1:nrow(background.network))
  cc1Exp <- paste0("reaction ", background.network$pfam_source, "=", 
                   background.network$pfam_target, " of ", 
                   background.network$gene_source, "=", background.network$gene_target)
  
  #
  ppi <- unique(background.network[, c("gene_source", "gene_target")])
  cnt <- 1 + length(cc1)
  cc2 <- paste0("xb", cnt:(nrow(ppi)+cnt-1))
  cc2Exp <- paste0("interaction ", ppi$gene_source, "=", ppi$gene_target)
  
  #
  nodes <- unique(c(ppi$gene_source, ppi$gene_target))
  cnt <- 1 + length(cc1) + length(cc2)
  cc3 <- paste0("xb", cnt:(length(nodes)+cnt-1))
  cc3Exp <- paste0("node ", nodes)
  
  #
  ppi <- unique(background.network[, c("gene_source", "gene_target", "pfam_source", "pfam_target")])
  cnt <- 1 + length(cc1) + length(cc2) + length(cc3)
  cc4Exp <- unique(c(paste0("domain ", ppi$pfam_source, " of protein ", ppi$gene_source),
                     paste0("domain ", ppi$pfam_target, " of protein ", ppi$gene_target)))
  cc4 <- paste0("xb", cnt:(length(cc4Exp)+cnt-1))
  
  #
  cnt <- 1 + length(cc1) + length(cc2) + length(cc3) + length(cc4)
  cc5 <- paste0("xb", cnt:(length(nodes)+cnt-1))
  cc5Exp <- paste0("dist ", nodes)
  
  variables <- list()
  variables[[length(variables)+1]] <- c(cc1, cc2, cc3, cc4, cc5)
  variables[[length(variables)+1]] <- c(cc1Exp, cc2Exp, cc3Exp, cc4Exp, cc5Exp)
  names(variables) <- c("var", "var_exp")
  return(variables)
  
}