read_solution_cplex <- function(cplexSolutionFileName="results1.txt", 
                                variables = variables, 
                                background.network = background.network){
  
  reacVar <- variables$var[which(grepl(pattern = "reaction", 
                                       x = variables$var_exp))]
  intVar <- variables$var[which(grepl(pattern = "interaction", 
                                      x = variables$var_exp))]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)
  
  cplexSolutionIntAll <- list()
  cplexSolutionReacAll <- list()
  
  sifAll1 <- list()
  sifAll2 <- list()
  
  for(ii in 1:(length(cplexSolution)-2)){
    
    sif1 <- matrix(data = NA, nrow = 1, ncol = 3)
    sif2 <- matrix(data = NA, nrow = 1, ncol = 3)
    
    currSolution <- cplexSolution[[ii]][[4]]
    
    for(jj in 1:length(currSolution)){
      
      if((currSolution[[jj]][1]%in%reacVar) && 
         round(as.numeric(currSolution[[jj]][3])==1)){
        
        reaction <- 
          strsplit(variables[[2]][which(variables[[1]]==currSolution[[jj]][1])], 
                   split = " ")[[1]][2]
        
        sif2bind <- 
          as.matrix(c(strsplit(x = reaction, split = "=")[[1]][1], 
                      "1", strsplit(x = reaction, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif1 <- rbind(sif1, sif2bind)
        
      }
      
      if((currSolution[[jj]][1]%in%intVar) && 
         round(as.numeric(currSolution[[jj]][3])==1)){
        
        reaction <- 
          strsplit(variables[[2]][which(variables[[1]]==currSolution[[jj]][1])], 
                   split = " ")[[1]][2]
        
        sif2bind <- 
          as.matrix(c(strsplit(x = reaction, split = "=")[[1]][1], 
                      "1", strsplit(x = reaction, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif2 <- rbind(sif2, sif2bind)
        
      }
      
    }
    
    save(sif1, file = paste0("sif1_", ii, ".RData"))
    save(sif2, file = paste0("sif2_", ii, ".RData"))
    
    sifAll1[[length(sifAll1)+1]] <- unique(sif1[-1, ])
    sifAll2[[length(sifAll2)+1]] <- unique(sif2[-1, ])
    
  }
  
  sifAll <- sifAll2
  for(ii in 1:length(sifAll)){
    
    if(ii==1){
      
      sif <- sifAll[[1]]
      
    } else {
      
      # sif <- unique(rbind(sif, sifAll[[ii]]))
      if(nrow(sifAll[[ii]]) > 0){
        
        for(jj in 1:nrow(sifAll[[ii]])){
          
          idx1 <- which(sif[, 1]==sifAll[[ii]][jj, 1])
          idx2 <- which(sif[, 3]==sifAll[[ii]][jj, 3])
          
          idx <- intersect(x = idx1, y = idx2)
          
          if(length(idx) > 0){
            
            sif[idx, 2] <- as.character(as.numeric(sif[idx, 2])+1)
            
          } else {
            
            sif <- rbind(sif, sifAll[[ii]][jj, ])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  sif <- sif[complete.cases(sif), ]
  
  temp <- matrix(data = , nrow = nrow(sif), ncol = 4)
  temp[, 1:3] <- sif
  
  for(ii in 1:nrow(temp)){
    
    idx1 <- which(background.network$gene_source==temp[ii, 1])
    idx2 <- which(background.network$gene_target==temp[ii, 3])
    
    idx <- intersect(x = idx1, y = idx2)
    
    numVec <- c()
    nn <- rep("", length(idx))
    for(jj in 1:length(idx)){
      
      pfam1 <- background.network$pfam_source[idx[jj]]
      pfam2 <- background.network$pfam_target[idx[jj]]
      
      nn[jj] <- paste0(pfam1, "=", pfam2)
      
      cnt <- 0
      for(kk in 1:length(sifAll1)){
        
        index1 <- which(sifAll1[[kk]][, 1]==pfam1)
        index2 <- which(sifAll1[[kk]][, 1]==pfam2)
        index <- intersect(x = index1, index2)
        if(length(index)>0){
          cnt <- cnt + 1
        }
        
      }
      
      numVec <- c(numVec, cnt)
      
    }
    
    
  }
  
  return(sif)
  
}