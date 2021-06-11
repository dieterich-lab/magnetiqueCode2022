write_mediation_constraints <- function(variables = variables, 
                                        background.network = background.network){
  
  ppi <- unique(background.network[, c("gene_source", "gene_target", 
                                       "pfam_source", "pfam_target")])
  
  cc1 <- rep("", nrow(ppi))
  cc2 <- rep("", nrow(ppi))
  
  for(ii in 1:nrow(ppi)){
    
    ## print(paste0("Progress ", ii, "/", nrow(ppi)))
    
    idx1 <- which(variables$var_exp==paste0("domain ",
                                            ppi$pfam_source[ii],
                                            " of protein ",
                                            ppi$gene_source[ii]))
    
    idx2 <- which(variables$var_exp==paste0("domain ",
                                            ppi$pfam_target[ii],
                                            " of protein ",
                                            ppi$gene_target[ii]))
    
    idx3 <- which(variables$var_exp==paste0("reaction ",
                                            ppi$pfam_source[ii],
                                            "=",
                                            ppi$pfam_target[ii],
                                            " of ",
                                            ppi$gene_source[ii],
                                            "=",
                                            ppi$gene_target[ii]))
    
    cc1[ii] <- paste0(variables$var[idx1], " + ", variables$var[idx2], 
                          " - 2 ", variables$var[idx3], " >= 0")
    
    idx1 <- which(variables$var_exp==paste0("node ", 
                                            ppi$gene_source[ii]))
    
    idx2 <- which(variables$var_exp==paste0("node ", 
                                            ppi$gene_target[ii]))
    
    idx3 <- which(variables$var_exp==paste0("interaction ", 
                                            ppi$gene_source[ii],
                                            "=",
                                            ppi$gene_target[ii]))
    
    cc2[ii] <- paste0(variables$var[idx1], " + ", variables$var[idx2], 
                      " - 2 ", variables$var[idx3], " >= 0")
    
  }
  
  return(c(cc1, unique(cc2)))
  
}