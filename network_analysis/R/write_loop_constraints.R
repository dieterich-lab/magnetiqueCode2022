write_loop_constraints <- function(variables = variables,
                                   background.network = background.network){
  
  sif <- unique(background.network[, c("gene_source", "gene_target")])
  
  species <- unique(c(sif[, 1], sif[, 2]))
  
  speciesVar <- variables$var[grepl(pattern = "node", x = variables$var_exp)]
  speciesExp <- variables$var_exp[grepl(pattern = "node", x = variables$var_exp)]
  reacVar <- variables$var[grepl(pattern = "interaction", x = variables$var_exp)]
  distVar <- variables$var[grepl(pattern = "dist ", x = variables$var_exp)]
  
  cc1 <- paste0(speciesVar, " - ", distVar, " <= 0")
  
  cc2 <- paste0(distVar, " <= ", 1001)
  
  cc3 <- rep("", nrow(sif))
  for(ii in 1:nrow(sif)){
    
    var1 <- variables$var[which(variables$var_exp==paste0("dist ", sif[ii, 2]))]
    var2 <- variables$var[which(variables$var_exp==paste0("dist ", sif[ii, 1]))]
    var3 <- variables$var[which(variables$var_exp==paste0("interaction ", sif[ii, 1], "=", sif[ii, 2]))]
    
    cc3[ii] <- paste0(var1, " - ", var2, " + 1001 ", var3, " <= ", -(1-1001))
    
  }
  
  
  return(c(cc1, cc2, cc3))
  
}