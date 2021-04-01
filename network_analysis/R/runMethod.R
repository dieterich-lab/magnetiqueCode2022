res <- runMethod(input.scores = input.scores, 
                 rmats.input = rmats.input, 
                 background.network = bg, 
                 map.table = map.table,
                 solverPath = "/usr/bin/cplex",
                 input.node = NULL,
                 lThresh = NULL,
                 top = 50,
                 lambda1 = 10,
                 lambda2 = 0.01,
                 mipgap = 0.05,
                 relgap = 0.05,
                 populate = 1000,
                 nSolutions = 100,
                 intensity = 0,
                 timelimit = 3600,
                 process_log = FALSE){
  
  options(scipen=999)
  
  input.scores <- bin_measurements(input.scores = input.scores, top = 50)
  bn <- prepare_bn(background.network = background.network, 
                   rmats.input = rmats.input, input.scores = input.scores,
                   map.table = map.table, input.node = input.node)
  variables <- create_variables(background.network = bn$background.network)
  objective.function <- write_objective_function(background.network = bn$background.network, 
                                                 variables = variables, 
                                                 input.scores = input.scores, 
                                                 lambda1 = lambda1, 
                                                 lambda2 = lambda2)
  cc1 <- write_input_constraints(variables = variables, 
                                 input.node = bn$input.node)
  cc2 <- write_as_constraints(background.network = bn$background.network, 
                              variables = variables, 
                              lThresh = lThresh)
  cc3 <- write_mediation_constraints(variables = variables, 
                                     background.network = bn$background.network)
  cc4 <- write_inout_constraints(background.network = bn$background.network, 
                                 input.node = bn$input.node, 
                                 input.scores = input.scores)
  cc5  <- write_domain_constraints(background.network = bn$background.network)
  allC <- c(cc1, cc2, cc3, cc4, cc5)
  idx2rem <- c(which(allC==""), which(is.na(allC)))
  if(length(idx2rem)>0){
    allC <- allC[-idx2rem]
  }
  # allC <- paste0("c", 1:length(allC), ": ", allC)
  # allC <- paste("c", 1:length(allC), ":\t", allC, "\t \t", sep = "")
  
  bounds <- paste0("0 <= ", variables$var, " <= 1")
  binaries <- variables$var
  
  # write the .lp file
  data = "testFile.lp"
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
  data2 = "cplexCommand.txt"
  write("read testFile.lp", data2)
  write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
  write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
  write(paste0("set mip pool replace ", replace), data2, append = TRUE)
  write(paste0("set mip limits populate ", populate), data2, append = TRUE)
  write(paste0("set mip pool capacity ", nSolutions), data2, append = TRUE)
  write(paste0("set mip pool intensity ", intensity), data2, append = TRUE)
  write(paste0("set timelimit ", timelimit), data2, append = TRUE)
  write("populate", data2, append = TRUE)
  write("write results1.txt sol all", data2, append = TRUE)
  write("quit", data2, append = TRUE)
  
  solve_with_cplex(solverPath = solverPath, variables = variables)
  
  
}