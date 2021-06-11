solve_with_cplex <- function(solverPath = solverPath, 
                             variables = variables,
                             condition = 1){
  
  if (Sys.info()[1]=="Windows") {
    file.copy(from = solverPath,to = getwd())
    system(paste0("cplex.exe -f cplexCommand_", condition, ".txt"))
    file.remove("cplex.exe")
  } else {
    system(paste0(solverPath, " -f cplexCommand_", condition, ".txt"))
  }
  
}