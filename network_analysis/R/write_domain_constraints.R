write_domain_constraints <- function(background.network = background.network){
  
  ppi <- unique(background.network[, c("gene_source", "gene_target")])
  
  c1 <- rep("", nrow(ppi))
  c2 <- rep(list(c()), nrow(ppi))
  for(ii in 1:nrow(ppi)){
    
    ## print(paste0("Progress --- ", ii, "/", nrow(ppi)))
    
    #
    idx <- intersect(x = which(background.network$gene_source==ppi$gene_source[ii]), 
                     y = which(background.network$gene_target==ppi$gene_target[ii]))
    
    c1[ii] <- paste0(paste0(variables$var[which(grepl(pattern = "reaction", x = variables$var_exp))[idx]], collapse = " + "),
                     " - ", variables$var[which(grepl(pattern = "interaction", x = variables$var_exp))[ii]], " >= 0")
    
    #
    c2[[ii]] <- paste0(variables$var[which(grepl(pattern = "reaction", x = variables$var_exp))[idx]], 
                       " - ", variables$var[which(grepl(pattern = "interaction", x = variables$var_exp))[ii]], " <= 0")
    
  }
  
  return(c(c1, unlist(c2)))
  
}