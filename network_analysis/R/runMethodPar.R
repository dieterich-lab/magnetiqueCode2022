runMethodPar <- function(input.scores = input.scores, 
                        rmats.input = rmats.input, 
                        background.network = bg, 
                        map.table = map.table,
                        solverPath = "/usr/bin/cplex",
                        input.node = NULL,
                        pValThresh = NULL,
                        top = 50,
                        lambda1 = 10,
                        lambda2 = 0.01,
                        mipgap = 0.05,
                        relgap = 0.05,
                        populate = 1000,
                        nSolutions = 100,
                        intensity = 0,
                        timelimit = 3600,
                        constraints3 = NULL,
                        constraints4 = NULL,
                        constraints5 = NULL,
                        constraints6 = NULL){
  
  options(scipen=999)
  
  input.scores <- bin_measurements(input.scores = input.scores, top = 50)
  bn <- prepare_bn(background.network = background.network, 
                   rmats.input = rmats.input, map.table = map.table, 
                   input.node = input.node)
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
                              pValThresh = pValThresh)
  if(is.null(constraints3)){
    cc3 <- write_mediation_constraints(variables = variables,
                                       background.network = bn$background.network)
    # save(cc3, file = "data_test/cc3.RData")
  } else {
    load(file = constraints3)
  }
  if(is.null(constraints4)){
    cc4 <- write_inout_constraints(background.network = bn$background.network,
                                   input.node = bn$input.node,
                                   input.scores = input.scores)
    # save(cc4, file = "data_test/cc4.RData")
  } else {
    load(file = constraints4)
  }
  if(is.null(constraints5)){
    cc5  <- write_domain_constraints(background.network = bn$background.network)
    # save(cc5, file = "cc5.RData")
  } else {
    load(file = constraints5)
  }
  if(is.null(constraints6)){
    cc6 <- write_loop_constraints(variables = variables, 
                                  background.network = bn$background.network)
    # save(cc6, file = "cc6.RData")
  } else {
    load(file = constraints6)
  }
  allC <- c(cc1, cc2, cc3, cc4, cc5, cc6)
  idx2rem <- c(which(allC==""), which(is.na(allC)))
  if(length(idx2rem)>0){
    allC <- allC[-idx2rem]
  }
  save(allC, file = "data_test/allC.RData")
  
  solutionList <- list()
  allC_ctrl <- allC
  for(ii in 1:nSolutions){
    
    allC <- allC_ctrl[sample(1:length(allC_ctrl))]
    
    bounds <- write_bounds(variables = variables)
    binaries <- write_binaries(variables = variables)
    
    # write the .lp file
    data = paste0("testFile_", ii, ".lp")
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
    write(paste0("read testFile_", ii, ".lp"), data2)
    write(paste0("set mip tolerances mipgap ", mipgap), data2, append = TRUE)
    write(paste0("set mip pool relgap ", relgap), data2, append = TRUE)
    write(paste0("set mip pool replace ", replace), data2, append = TRUE)
    write(paste0("set mip limits populate ", 3), data2, append = TRUE)
    write(paste0("set mip pool capacity ", 1), data2, append = TRUE)
    write(paste0("set mip pool intensity ", intensity), data2, append = TRUE)
    write(paste0("set timelimit ", timelimit), data2, append = TRUE)
    write("populate", data2, append = TRUE)
    write("write results1.txt sol all", data2, append = TRUE)
    write("quit", data2, append = TRUE)
    
    solve_with_cplex(solverPath = solverPath, variables = variables)
    
    sif <- read_solution_cplex(variables = variables, 
                               background.network = bn$background.network)
    
    solutionList[[length(solutionList)+1]] <- sif
    
    save(solutionList, file = "solutionList.RData")
    
    cleanupILP(condition = ii)
    
  }
  
  # combSol <- combineSolutions(solutionList = solutionList, input.scores = input.scores, weightThresh = 0)
  
  return(solutionList)
  
}