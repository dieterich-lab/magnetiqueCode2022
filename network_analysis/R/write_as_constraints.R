write_as_constraints <- function(background.network = background.network,
                                 variables = variables,
                                 pValThresh = NULL){
  
  if(is.null(pValThresh)){
    
    return(NULL)
    
  } else {
    
    idx1 <- which(background.network$source_fdr<=pValThresh)
    idx2 <- which(background.network$source_score<=0)
    idx <- intersect(x = idx1, y = idx2)
    
    if(length(idx)==0){
      cc1 <- NULL
      
    } else {
      
      cc1 <- rep("", length(idx))
      
      for(ii in 1:length(idx)){
        
        idx_var <- which(variables$var_exp==paste0("domain ",
                                                   background.network$pfam_source[idx[ii]],
                                                   " of protein ",
                                                   background.network$gene_source[idx[ii]]))
        
        cc1[ii] <- paste0(variables$var[idx_var[1]], " = 0")
        
      }
      
    }
    
    idx1 <- which(background.network$target_fdr<=pValThresh)
    idx2 <- which(background.network$target_score<=0)
    idx <- intersect(x = idx1, y = idx2)
    
    if(length(idx)==0){
      cc2 <- NULL
      
    } else {
      
      cc2 <- rep("", length(idx))
      
      for(ii in 1:length(idx)){
        
        idx_var <- which(variables$var_exp==paste0("domain ",
                                                   background.network$pfam_target[idx[ii]],
                                                   " of protein ",
                                                   background.network$gene_target[idx[ii]]))
        
        cc2[ii] <- paste0(variables$var[idx_var[1]], " = 0")
        
      }
      
    }
    
    return(c(cc1, cc2))
    
  }
  
}