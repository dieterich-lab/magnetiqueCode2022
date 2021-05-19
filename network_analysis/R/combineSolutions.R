combineSolutions <- function(solutionList = solutionList, 
                             input.scores = input.scores, 
                             weightThresh = 5){
  
  for(ii in 1:length(solutionList)){
    
    if(ii==1){
      
      combSIF <- solutionList[[1]]
      
    } else {
      
      currSIF <- solutionList[[ii]]
      
      for(jj in 1:nrow(currSIF)){
        
        idx1 <- which(combSIF[, 1]==currSIF[jj, 1])
        idx2 <- which(combSIF[, 3]==currSIF[jj, 3])
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx)>0){
          
          combSIF[idx, 2] <- as.character(as.numeric(combSIF[idx, 2])+1)
          
          ss <- strsplit(x = currSIF[jj, 4], split = ";", fixed = TRUE)[[1]]
          
          if(length(ss)>0){
            
            for(kk in 1:length(ss)){
              
              if(!grepl(pattern = ss[kk], x = combSIF[idx, 4])){
                
                combSIF[idx, 4] <- paste0(combSIF[idx, 4], ";", ss[kk])
                
              }
              
            }
            
          }
          
        } else {
          
          combSIF <- rbind(combSIF, currSIF[jj, ])
          
        }
        
      }
      
    }
    
  }
  
  combSIF <- combSIF[which(as.numeric(combSIF[, 2])>=weightThresh), ]
  
  uSpecies <- unique(c(combSIF[, 1], combSIF[, 3]))
  attributes <- matrix(data = , nrow = length(uSpecies), ncol = 2)
  attributes[, 1] <- uSpecies
  attributes[, 2] <- "P"
  attributes[which(uSpecies%in%input.scores$id), 2] <- "T"
  attributes[which(uSpecies=="Perturbation"), 2] <- "I"
  colnames(attributes) <- c("species", "type")
  
  returnList <- list()
  returnList[[length(returnList)+1]] <- combSIF
  returnList[[length(returnList)+1]] <- attributes
  names(returnList) <- c("combine_solution", "node_attributes")
  
  return(returnList)
  
  
}