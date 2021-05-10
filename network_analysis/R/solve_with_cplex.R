solve_with_cplex <- function(solverPath = solverPath, 
                             variables = variables){
  
  if (Sys.info()[1]=="Windows") {
    file.copy(from = solverPath,to = getwd())
    system(paste0("cplex.exe -f cplexCommand.txt"))
    file.remove("cplex.exe")
  } else {
    system(paste0(solverPath, " -f cplexCommand.txt"))
  }
  
}