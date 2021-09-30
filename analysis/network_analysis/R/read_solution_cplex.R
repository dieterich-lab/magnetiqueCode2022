read_solution_cplex <- function(variables = variables,
                                background.network = background.network,
                                condition = 1){
  
  cplexSolutionFileName <- paste0("results_", condition, ".txt")
  
  reacVar <- variables$var[which(grepl(pattern = "reaction",
                                       x = variables$var_exp))]
  intVar <- variables$var[which(grepl(pattern = "interaction",
                                      x = variables$var_exp))]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  print("xmlParse Successful")
  cplexSolution <- xmlToList(cplexSolutionData)
  print("xmlToList Successful")
  
  rm(cplexSolutionData)
  
  cplexSolutionIntAll <- list()
  cplexSolutionReacAll <- list()
  
  sifAll1 <- list()
  sifAll2 <- list()
  
  # 2:(length(cplexSolution)-1)
  for(ii in 1:(length(cplexSolution)-1)){
    
    sif1 <- matrix(data = NA, nrow = 1, ncol = 3)
    sif2 <- matrix(data = NA, nrow = 1, ncol = 3)
    
    currSolution <- cplexSolution[[ii]][[4]]
    
    varvar <- unlist(lapply(currSolution, '[', 1))
    valval <- as.numeric(unlist(lapply(currSolution, '[', 3)))
    
    idxReac <- intersect(x = which(varvar%in%reacVar), y = which(valval>=0.99))
    idxInt <- intersect(x = which(varvar%in%intVar), y = which(valval>=0.99))
    
    if(length(idxReac)>0 && length(idxInt)>0){
      
      for(jj in 1:length(idxReac)){
        
        currVar <- varvar[idxReac[jj]]
        currReac <- variables$var_exp[which(variables$var==currVar)]
        currReac <- strsplit(x = currReac, split = " ", fixed = TRUE)[[1]][2]
        
        sif2bind <-
          as.matrix(c(strsplit(x = currReac, split = "=")[[1]][1],
                      "1", strsplit(x = currReac, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif1 <- rbind(sif1, sif2bind)
        
      }
      sifAll1[[length(sifAll1)+1]] <- unique(sif1[-1, ])
      
      for(jj in 1:length(idxInt)){
        
        currVar <- varvar[idxInt[jj]]
        currReac <- variables$var_exp[which(variables$var==currVar)]
        currReac <- strsplit(x = currReac, split = " ", fixed = TRUE)[[1]][2]
        
        sif2bind <-
          as.matrix(c(strsplit(x = currReac, split = "=")[[1]][1],
                      "1", strsplit(x = currReac, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif2 <- rbind(sif2, sif2bind)
        
      }
      sifAll2[[length(sifAll2)+1]] <- unique(sif2[-1, ])
      
    }
    
  }
  
  res1 <- sifAll1
  res2 <- sifAll2
  
  for(ii in 1:length(res2)){
    
    if(ii == 1){
      
      combSIF <- unique(res2[[ii]])
      
    } else {
      
      currSIF <- unique(res2[[ii]])
      
      for(jj in 1:nrow(currSIF)){
        
        ss <- currSIF[jj, 1]
        tt <- currSIF[jj, 3]
        
        idx1 <- which(combSIF[, 1]==ss)
        idx2 <- which(combSIF[, 3]==tt)
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx)>0){
          
          combSIF[idx, 2] <- as.character(as.numeric(combSIF[idx, 2])+1)
          
        } else {
          
          toBind <- t(as.matrix(c(ss, "1", tt)))
          colnames(toBind) <- colnames(combSIF)
          
          combSIF <- unique(rbind(combSIF, toBind))
          
        }
        
      }
      
    }
    
  }
  
  domains <- rep("", nrow(combSIF))
  for(ii in 1:nrow(combSIF)){
    
    ss <- combSIF[ii, 1]
    tt <- combSIF[ii, 3]
    
    idx1 <- which(bg$gene_source==ss)
    idx2 <- which(bg$gene_target==tt)
    
    idx <- intersect(x = idx1, y = idx2)
    cnt <- 1
    
    if(length(idx)>0){
      
      for(jj in 1:length(idx)){
        
        ss2 <- bg$pfam_source[idx[jj]]
        tt2 <- bg$pfam_target[idx[jj]]
        
        counts <- 0
        for(kk in 1:length(res1)){
          
          idx3 <- which(res1[[kk]][, 1]==ss2)
          idx4 <- which(res1[[kk]][, 3]==tt2)
          
          if(length(intersect(x = idx3, y = idx4))>0){
            
            counts <- counts + 1
            
          }
          
        }
        
        if(counts > 0){
          
          if(cnt == 1){
            
            pp <- paste0(bg$pfam_source[idx[jj]], "=", bg$pfam_target[idx[jj]])
            cnt <- cnt + 1
            
          } else {
            
            pp <- paste0(pp, "; ", bg$pfam_source[idx[jj]], "=", bg$pfam_target[idx[jj]])
            
          }
          
        }
        
      }
      
      domains[ii] <- pp
      
      pp <- ""
      
    }
    
  }
  
  sif <- matrix(data = , nrow = nrow(combSIF), ncol = 4)
  sif[, 1:3] <- combSIF
  sif[, 4] <- domains
  
  idx2rem <- setdiff(x = which(sif[, 4]==""), y = which(sif[, 1]=="Perturbation"))
  if(length(idx2rem)>0){
    sif <- sif[-idx2rem, ]
  }
  
  if(nrow(sif)>0){
    
    colnames(sif) <- c("source", "weight", "target", "reaction")
    
    # return(sif)
    returnList <- list()
    returnList[[length(returnList)+1]] <- sifAll1
    returnList[[length(returnList)+1]] <- sifAll2
    returnList[[length(returnList)+1]] <- sif
    
    names(returnList) <- c("domain_interactions", "protein_interactions", "combined_interactions")
    
    return(returnList)
  } else {
    
    return(NULL)
    
  }
  
}