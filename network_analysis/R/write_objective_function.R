write_objective_function <- function(background.network = background.network, 
                                     variables = variables, 
                                     input.scores = input.scores, 
                                     lambda1=10, 
                                     lambda2=0.01){
  
  print("Writing the objective function and constraints. This might take a bit of time..")
  
  # Write main objective - Minimizing the TF scores
  tfNodes <- input.scores$id
  idx <- which(variables$var_exp%in%paste0("node ", tfNodes))
  if(length(idx)==0){
    stop("No input protein mapped to the background network. Check inputs/background network")
  }
  
  objective.function <- ""
  cnt <- 1
  for(ii in 1:nrow(input.scores)){
    
    idx <- which(variables$var_exp==paste0("node ", input.scores$id[ii]))
    if(length(idx)==1){
      
      if(cnt == 1){
        
        if(input.scores$bin[ii]==1){
          objective.function <- paste0(objective.function, lambda1, " ", variables$var[idx])
        } else {
          objective.function <- paste0(objective.function, "- ", lambda1, " ", variables$var[idx])
        }
        cnt <- cnt + 1
        
      } else {
        
        if(input.scores$bin[ii]==1){
          objective.function <- paste0(objective.function, " + ", lambda1, " ", variables$var[idx])
        } else {
          objective.function <- paste0(objective.function, " - ", lambda1, " ", variables$var[idx])
        }
        
      }
      
    }
    
  }
  
  # Write second objective - AS scores
  idx <- which(background.network$min_score!=0)
  obj <- rep("", length(idx))
  for(ii in 1:length(idx)){
    
    currScore <- background.network$min_score[idx[ii]]
    idx_var <- which(variables$var_exp==paste0("reaction ", 
                                               background.network$pfam_source[idx[ii]],
                                               "=",
                                               background.network$pfam_target[idx[ii]],
                                               " of ",
                                               background.network$gene_source[idx[ii]],
                                               "=",
                                               background.network$gene_target[idx[ii]]))
    if(currScore > 0){
      obj[ii] <- paste0(" - ", currScore, " ", variables$var[idx_var])
    } else {
      obj[ii] <- paste0(" + ", abs(currScore), " ", variables$var[idx_var])
    }
  }
  
  obj <- paste0(obj, collapse = "")
  
  objective.function <- paste0(objective.function, obj)
  
  # Write third objective - size penalty factor
  idx <- which(grepl(pattern = "interaction", x = variables$var_exp, fixed = TRUE))
  obj <- paste0(" + ", lambda2, " ", variables$var[idx])
  obj <- paste0(obj, collapse = "")
  objective.function <- paste0(objective.function, obj)
  
  return(objective.function)
  
}