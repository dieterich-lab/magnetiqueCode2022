write_bounds <- function(variables = variables){
  
  bounds1 <- variables$var[setdiff(x = 1:length(variables$var), 
                                   y = which(grepl(pattern = "dist ", x = variables$var_exp)))]
  
  bounds2 <- variables$var[which(grepl(pattern = "dist ", x = variables$var_exp))]
  
  bounds <- c(paste0("0 <= ", bounds1, " <= 1"), paste0("0 <= ", bounds2, " <= 1001"))
  
  return(bounds)
  
}