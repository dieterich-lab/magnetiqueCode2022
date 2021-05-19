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
                         nSolutions = 100,
                         timelimit = 3600,
                         constraints3 = NULL,
                         constraints4 = NULL,
                         constraints5 = NULL,
                         constraints6 = NULL){
  
  options(scipen=999)
  
  input.scores <- bin_measurements(input.scores = input.scores, top = top)
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
  
  # solutionList <- vector("list", length = nSolutions)
  # allC_ctrl <- allC
  solutionList <- foreach(ii = 1:nSolutions) %dopar% {
    
    computeILP(variables = variables,
               background.network = bn$background.network, 
               solverPath = solverPath, 
               mipgap = mipgap, 
               relgap = relgap, 
               timelimit = timelimit, 
               condition = ii, 
               allC = allC, 
               objective.function = objective.function)
    
  }
  
  # combSol <- combineSolutions(solutionList = solutionList, input.scores = input.scores, weightThresh = 0)
  
  return(solutionList)
  
}