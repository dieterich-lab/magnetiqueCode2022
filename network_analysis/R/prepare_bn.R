prepare_bn <- function(background.network = background.network,
                       rmats.input = rmats.input,
                       map.table = map.table,
                       input.node = NULL, 
                       input.scores = input.scores){
  
  print("Processing the Background Network..")
  
  integrate_scores_in_bn <- function(rmats.input = rmats.input, 
                                     background.network = background.network, 
                                     map.table = map.table){
    
    # print("Integrating AS scores in the Background Network..")
    source_score <- rep(0, nrow(background.network))
    target_score <- rep(0, nrow(background.network))
    
    for(ii in 1:nrow(rmats.input)){
      
      idx1 <- which(paste0("chr", map.table$Chromosome.scaffold.name)==rmats.input$chr[ii])
      if(rmats.input$strand[ii]=="-"){
        idx2 <- which(map.table$Strand==-1)
      } else {
        idx2 <- which(map.table$Strand==1)
      }
      
      idx3 <- which(map.table$Genomic.coding.start>=rmats.input$exonStart_0base[ii])
      idx4 <- which(map.table$Genomic.coding.end<=rmats.input$exonEnd[ii])
      
      idx <- intersect(x = intersect(x = idx1, y = idx2), 
                       y = intersect(x = idx3, y = idx4))
      
      if(length(idx)>0){
        
        for(jj in 1:length(idx)){
          
          currTranscript <- map.table$Transcript.stable.ID[idx[jj]]
          idx_source <- which(background.network$transcript_source==currTranscript)
          idx_target <- which(background.network$transcript_target==currTranscript)
          
          if(length(idx_source)>0){
            source_score[idx_source] <- rmats.input$IncLevelDifference[ii]
          }
          
          if(length(idx_target)>0){
            target_score[idx_target] <- rmats.input$IncLevelDifference[ii]
          }
          
        }
        
      }
      
    }
    
    min_score <- pmin(source_score, target_score)
    
    background.network$source_score <- source_score
    background.network$target_score <- target_score
    background.network$min_score <- min_score
    
    return(background.network)
    
  }
  
  # remove first self-interactions
  idx2rem <- which(background.network$gene_source==background.network$gene_target)
  if(length(idx2rem)>0){
    background.network <- background.network[-idx2rem, ]
  }
  
  if(is.null(input.node)){
    
    ppi <- unique(background.network[, c("gene_source", "gene_target")])
    
    top.proteins <- unique(setdiff(x = ppi$gene_source, y = ppi$gene_target))
    
    for(ii in 1:length(top.proteins)){
      
      upfam <- unique(background.network$pfam_source[which(background.network$gene_source==top.proteins[ii])])
      utrans <- unique(background.network$transcript_source[which(background.network$gene_source==top.proteins[ii])])
      uentrez <- unique(background.network$entrez_source[which(background.network$gene_source==top.proteins[ii])])
      
      toBind <- matrix(data = , nrow = length(upfam), ncol = ncol(background.network))
      colnames(toBind) <- colnames(background.network)
      for(jj in 1:nrow(toBind)){
        
        toBind[jj, 1] <- "Perturbation"
        toBind[jj, 2] <- utrans[jj]
        toBind[jj, 3] <- "Perturbation"
        toBind[jj, 4] <- upfam[jj]
        toBind[jj, 5] <- "Perturbation"
        toBind[jj, 6] <- uentrez[jj]
        toBind[jj, 7] <- "Perturbation"
        toBind[jj, 8] <- top.proteins[ii]
        
      }
      
      toBind <- as.data.frame(toBind)
      background.network <- rbind(background.network, toBind)
      
    }
    
    background.network <- unique(background.network)
    input.node <- "Perturbation"
    
  }
  
  background.network <- integrate_scores_in_bn(rmats.input = rmats.input, 
                                               background.network = background.network, 
                                               map.table = map.table)
  
  ppi <- unique(background.network[, c("gene_source", "gene_target" )])
  targets <- input.scores$id[which(input.scores$id%in%unique(c(ppi$gene_source, ppi$gene_target)))]
  gg <- graph_from_data_frame(d = ppi, directed = TRUE)
  adj <- get.adjacency(graph = gg)
  ppiRed <- matrix(data = , nrow = 1, ncol = 2)
  for(ii in 1:length(input.node)){
    
    for(jj in 1:length(targets)){
      
      sP <- get.all.shortest.paths(graph = gg, 
                                   from = which(rownames(adj)==input.node[ii]), 
                                   to = which(rownames(adj)==targets[jj]))
      
      if(length(sP[[1]]) > 0){
        
        sP <- sP[[1]]
        
        for(kk in 1:length(sP)){
          
          interactors <- sP[[kk]]
          if(length(interactors)>1){
            
            toBind <- matrix(data = , nrow = length(interactors)-1, ncol = 2)
            for(ll in 1:(length(interactors)-1)){
              
              toBind[ll, 1] <- rownames(adj)[interactors[ll]]
              toBind[ll, 2] <- rownames(adj)[interactors[ll+1]]
              
            }
            
            ppiRed <- unique(rbind(ppiRed, toBind))
            
          }
          
        }
        
      }
      
    }
    
  }
  
  ppiRed <- ppiRed[-c(1), ]
  idx2keep <- c()
  for(ii in 1:nrow(ppiRed)){
    
    idx1 <- which(background.network$gene_source==ppiRed[ii, 1])
    idx2 <- which(background.network$gene_target==ppiRed[ii, 2])
    
    idx2keep <- c(idx2keep, intersect(x = idx1, y = idx2))
    
  }
  
  background.network <- background.network[idx2keep, ]
  
  returnList <- list()
  returnList[[length(returnList)+1]] <- background.network
  returnList[[length(returnList)+1]] <- input.node
  names(returnList) <- c("background.network", "input.node")
  
  return(returnList)
  
}