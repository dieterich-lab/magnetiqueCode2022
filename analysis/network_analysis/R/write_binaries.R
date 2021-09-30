write_binaries <- function(variables = variables){
  
  binaries <- variables$var[setdiff(x = 1:length(variables$var), 
                                    y = which(grepl(pattern = "dist ", x = variables$var_exp)))]
  
  return(binaries)
  
}