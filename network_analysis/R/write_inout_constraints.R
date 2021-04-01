write_inout_constraints <- function(background.network = background.network,
                                    input.node = input.node,
                                    input.scores = input.scores){
  
  allNodes <- unique(c(background.network$gene_source, 
                       background.network$gene_target))
  pSet <- setdiff(x = allNodes, y = input.node)
  
  ppi <- unique(background.network[, c("gene_source", "gene_target")])
  
  ccin <- rep("", length(pSet))
  ccout <- rep("", length(setdiff(pSet, input.scores$id)))
  for(ii in 1:length(pSet)){
    
    if(pSet[ii]%in%input.scores$id){
      
      # #
      # idx <- which(ppi$gene_source==pSet[ii])
      # varout <- variables$var[which(grepl(pattern = "interaction", 
      #                                     x = variables$var_exp))[idx]]
      # ccout[ii] <- paste0(paste0(varout, collapse = " + "), " - ", 
      #                     variables$var[which(variables$var_exp==paste0("node ", pSet[ii]))],
      #                     " >= 0")
      
      #
      idx <- which(ppi$gene_target==pSet[ii])
      varin <- variables$var[which(grepl(pattern = "interaction", 
                                          x = variables$var_exp))[idx]]
      ccin[ii] <- paste0(paste0(varin, collapse = " + "), " - ", 
                         variables$var[which(variables$var_exp==paste0("node ", pSet[ii]))],
                         " >= 0")
      
    } else {
      
      #
      idx <- which(ppi$gene_source==pSet[ii])
      varout <- variables$var[which(grepl(pattern = "interaction", 
                                          x = variables$var_exp))[idx]]
      ccout[ii] <- paste0(paste0(varout, collapse = " + "), " - ", 
                          variables$var[which(variables$var_exp==paste0("node ", pSet[ii]))],
                          " >= 0")
      
      #
      idx <- which(ppi$gene_target==pSet[ii])
      varout <- variables$var[which(grepl(pattern = "interaction", 
                                          x = variables$var_exp))[idx]]
      ccin[ii] <- paste0(paste0(varout, collapse = " + "), " - ", 
                         variables$var[which(variables$var_exp==paste0("node ", pSet[ii]))],
                         " >= 0")
      
    }
    
  }
  
  return(c(ccin, ccout))
  
}