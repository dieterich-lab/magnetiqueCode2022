reduce_network <- function(sif = sif, 
                           input.nodes = "Perturbation", 
                           input.scores){
  
  df <- as.data.frame(sif[, c(1, 3)])
  gg <- graph_from_data_frame(d = df, directed = TRUE)
  adj <- get.adjacency(graph = gg)
  
  ww <- (as.numeric(sif[, 2]) - min(as.numeric(sif[, 2])))/(max(as.numeric(sif[, 2]))-min(as.numeric(sif[, 2])))
  ww <- 1 - ww
  
  tf <- input.scores$id
  
  toBind <- matrix(data = , nrow = 1, ncol = 4)
  colnames(toBind) <- colnames(sif)
  for(ii in 1:length(input.nodes)){
    
    for(jj in 1:length(tf)){
      
      idx1 <- which(rownames(adj)==input.nodes[ii])
      idx2 <- which(rownames(adj)==tf[jj])
      
      if((length(idx1)>0) && (length(idx2)>0)){
        
        sP <- all_shortest_paths(graph = gg, from = idx1, to = idx2, weights = ww)[[1]]
        for(kk in 1:length(sP)){
          
          currPath <- sP[[kk]]
          
          zz <- matrix(data = , nrow = length(currPath)-1, ncol = 4)
          colnames(zz) <- colnames(sif)
          for(ll in 1:(length(currPath)-1)){
            
            for(mm in (ll+1):length(currPath)){
              
              zz[ll, 1] <- rownames(adj)[currPath[ll]]
              zz[ll, 3] <- rownames(adj)[currPath[mm]]
              
              ii1 <- which(sif[, 1]==zz[ll, 1])
              ii2 <- which(sif[, 3]==zz[ll, 3])
              index <- intersect(x = ii1, y = ii2)
              
              zz[ll, 2] <- sif[idx, 2]
              zz[ll, 4] <- sif[idx, 4]
              
            }
            
          }
          
          toBind <- unique(rbind(toBind, zz))
          
        }
        
      }
      
    }
    
  }
  
}