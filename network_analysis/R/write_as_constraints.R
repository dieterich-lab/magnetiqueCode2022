write_as_constraints <- function(background.network = background.network,
                                 variables = variables,
                                 lThresh = NULL){
  
  if(is.null(lThresh)){
    
    return(NULL)
    
  } else {
    
    idx <- which(background.network$min_score<=lThresh)
    
    if(length(idx)==0){
      
      return(NULL)
      
    } else {
      
      cc <- rep("", length(idx))
      
      for(ii in 1:length(idx)){
        
        idx_var <- which(variables$var_exp==paste0("reaction ", 
                                                   background.network$pfam_source[idx[ii]],
                                                   "=",
                                                   background.network$pfam_target[idx[ii]],
                                                   " of ",
                                                   background.network$gene_source[idx[ii]],
                                                   "=",
                                                   background.network$gene_target[idx[ii]]))
        
        cc[ii] <- paste0(variables$var[idx_var], " = 0")
        
      }
      
      return(cc)
      
    }
    
  }
  
}