reduce_network_ilp <- function(sif = sif, lambda = 30, solverPath){
  
  #
  create_variables_net <- function(sif){
    variables = list()
    variables[[length(variables)+1]] <- paste0("xb", 1:nrow(sif))
    variables[[length(variables)+1]] <- paste0("reaction ", sif[, 1], "=", sif[, 3])
    names(variables) <- c("var", "var_exp")
    return(variables)
  }
  
  #
  create_bounds_net <- function(variables){
    return(paste0("0 <= ", variables$var, " <= 1"))
  }
  
  #
  write_objective_function_net <- function(variables, sif, lambda){
    
    of1 <- paste0(paste0( " + ", as.numeric(sif[, 2]), " ", variables$var), collapse = "")
    of1 <- substr(x = of1, start = 4, stop = nchar(of1))
    
    of2 <- paste0(paste0( " - ", lambda, " ", variables$var), collapse = "")
    
    return(paste0(of1, of2))
  }
  
  #
  write_constraints_net <- function(sif = sif){
    
    uSpecies <- unique(c(sif[, 1], sif[, 3]))
    
    cc1 <- c()
    for(ii in 1:length(uSpecies)){
      
      idx <- which(sif[, 1]==uSpecies[ii])
      if(length(idx)>0){
        
        for(jj in 1:length(idx)){
          
          if(jj==1){
            
            cc <- variables$var[idx[1]]
            
          } else {
            
            cc <- paste0(cc, " + ", variables$var[idx[jj]])
            
          }
          
        }
        
        cc1 <- c(cc1, cc)
        
      }
      
    }
    
    cc2 <- c()
    for(ii in 1:length(uSpecies)){
      
      idx <- which(sif[, 3]==uSpecies[ii])
      if(length(idx)>0){
        
        for(jj in 1:length(idx)){
          
          if(jj==1){
            
            cc <- variables$var[idx[1]]
            
          } else {
            
            cc <- paste0(cc, " + ", variables$var[idx[jj]])
            
          }
          
        }
        
        cc2 <- c(cc2, cc)
        
      }
      
    }
    
    return(c(paste0(cc1, " >= 1"), paste0(cc2, " >= 1")))
    
  }
  
  #
  write_problem_net <- function(of, variables, bounds, constraints){
    
    # write the .lp file
    data = paste0("testFile.lp")
    write("enter Problem", data)
    write("", data, append = TRUE)
    write("Maximize", data, append = TRUE)
    write(of, data, append = TRUE)
    write("Subject To", data, append = TRUE)
    write(constraints, data, append = TRUE)
    write("Bounds", data, append = TRUE)
    write(bounds, data, append = TRUE)
    write("Binaries", data, append = TRUE)
    write(variables$var, data, append = TRUE)
    write("End", data, append = TRUE)
    
    # write cplexCommand file
    data2 = "cplexCommand.txt"
    write(paste0("read testFile.lp"), data2)
    write(paste0("set mip tolerances mipgap ", 0), data2, append = TRUE)
    write(paste0("set mip pool relgap ", 0), data2, append = TRUE)
    write(paste0("set mip pool replace ", 2), data2, append = TRUE)
    write(paste0("set mip limits populate ", 3), data2, append = TRUE)
    write(paste0("set mip pool capacity ", 1), data2, append = TRUE)
    write(paste0("set mip pool intensity ", 0), data2, append = TRUE)
    write(paste0("set timelimit ", 3600), data2, append = TRUE)
    write("populate", data2, append = TRUE)
    write("write results1.txt sol all", data2, append = TRUE)
    write("quit", data2, append = TRUE)
    
    if (Sys.info()[1]=="Windows") {
      file.copy(from = solverPath,to = getwd())
      system(paste0("cplex.exe -f cplexCommand.txt"))
      file.remove("cplex.exe")
    } else {
      system(paste0(solverPath, " -f cplexCommand.txt"))
    }
    
  }
  
  #
  read_solution_cplex_net <- function(variables, sif){
    
    cplexSolutionData <- xmlParse("results1.txt")
    cplexSolution <- xmlToList(cplexSolutionData)
    cplexSolution <- cplexSolution$CPLEXSolution[[4]]
    functional_reactions <- c()
    for(ii in 1:length(cplexSolution)){
      
      currEdge <- cplexSolution[[ii]]
      if(currEdge[3]=="1"){
        functional_reactions <- c(functional_reactions, 
                                  as.numeric(gsub(pattern = "xb", replacement = "", x = currEdge[1])))
      }
      
    }
    
    return(sif[functional_reactions, ])
    
  }
  
  #
  cleanup_ilp_net <- function(){
    
    if(file.exists(paste0("testFile.lp"))){
      file.remove(paste0("testFile.lp"))
    }
    
    if(file.exists(paste0("results1.txt"))){
      file.remove(paste0("results1.txt"))
    }
    
    if(file.exists("cplex.log")){
      file.remove("cplex.log")
    }
    
    if(file.exists(paste0("cplexCommand.txt"))){
      file.remove(paste0("cplexCommand.txt"))
    }
    
    AllFiles <- list.files()
    CloneFiles <- which(grepl(pattern = "clone",x = AllFiles,fixed = TRUE))
    if (length(CloneFiles)>0) {
      for (counter in seq_len(length(CloneFiles))) {
        file.remove(AllFiles[CloneFiles[counter]])
      }
    }
    
  }
  
  #
  variables <- create_variables_net(sif = sif)
  bounds <- create_bounds_net(variables = variables)
  of <- write_objective_function_net(variables = variables, sif = sif, lambda = lambda)
  constraints <- write_constraints_net(sif = sif)
  write_problem_net(of, variables, bounds, constraints)
  reduced_sif <- read_solution_cplex_net(variables = variables, sif = sif)
  cleanup_ilp_net()
  
}