computeILP <- function(variables = variables,
                       objective.function = objective.function,
                       background.network = bg, 
                       solverPath = solverPath,
                       mipgap = 0.05,
                       relgap = 0.05,
                       timelimit = 3600,
                       condition = condition,
                       allC = allC){
  
  allC <- allC[sample(1:length(allC))]
  
  bounds <- write_bounds(variables = variables)
  binaries <- write_binaries(variables = variables)
  
  # write the .lp file
  data = paste0("testFile_", condition, ".lp")
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  write(objective.function, data, append = TRUE)
  write("Subject To", data, append = TRUE)
  write(allC, data, append = TRUE)
  write("Bounds", data, append = TRUE)
  write(bounds, data, append = TRUE)
  write("Binaries", data, append = TRUE)
  write(binaries, data, append = TRUE)
  write("End", data, append = TRUE)
  
  # write cplexCommand file
  data2 = paste0("cplexCommand_", condition, ".txt")
  write(paste0("read testFile_", condition, ".lp"), data2)
  write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
  write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
  write(paste0("set mip pool replace ", replace), data2, append = TRUE)
  write(paste0("set mip limits populate ", 2), data2, append = TRUE)
  write(paste0("set mip pool capacity ", 1), data2, append = TRUE)
  write(paste0("set mip pool intensity ", 0), data2, append = TRUE)
  write(paste0("set timelimit ", timelimit), data2, append = TRUE)
  write("populate", data2, append = TRUE)
  write(paste0("write results_", condition, ".txt sol all"), data2, append = TRUE)
  write("quit", data2, append = TRUE)
  
  solve_with_cplex(solverPath = solverPath, variables = variables, 
                   condition = condition)
  
  sif <- read_solution_cplex(variables = variables, 
                             background.network = bn$background.network, 
                             condition = condition)
  
  cleanupILP(condition = condition)
  
  return(sif)
  
}