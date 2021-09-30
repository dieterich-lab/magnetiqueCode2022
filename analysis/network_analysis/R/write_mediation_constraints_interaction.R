write_mediation_constraints_interaction  <- function(variables = variables, 
                                                     background.network = background.network){
  
  ppi <- unique(background.network[, c("gene_source", "gene_target")])
  
  cc <- rep("", nrow(background.network))
  
  for(ii in 1:nrow(ppi)){
    
    ## print(paste0("Progress ", ii, "/", nrow(ppi)))
    
    idx1 <- which(variables$var_exp==paste0("node ", 
                                            ppi$gene_source[ii]))
    
    idx2 <- which(variables$var_exp==paste0("node ", 
                                            ppi$gene_target[ii]))
    
    idx3 <- which(variables$var_exp==paste0("interaction ", 
                                            ppi$gene_source[ii],
                                            "=",
                                            ppi$gene_target[ii]))
    
    cc[ii] <- paste0(variables$var[idx1], " + ", variables$var[idx2], 
                     " - 2 ", variables$var[idx3], " >= 0")
    
  }
  
  return(cc)
  
}