prepare_bn <- function(background.network = background.network,
                       rmats.input = rmats.input,
                       map.table = map.table,
                       input.node = NULL){
  
  print("Processing the Background Network..")
  
  integrate_scores_in_bn <- function(rmats.input = rmats.input, 
                                     background.network = background.network, 
                                     map.table = map.table){
    
    # print("Integrating AS scores in the Background Network..")
    source_score <- rep(0, nrow(background.network))
    target_score <- rep(0, nrow(background.network))
    source_fdr <- rep(1, nrow(background.network))
    target_fdr <- rep(1, nrow(background.network))
    
    uTranscripts <- unique(rmats$transcript_id)
    for(ii in 1:length(uTranscripts)){
      
      score <- min(rmats$IncLevelDifference[which(rmats$transcript_id==uTranscripts[ii])])
      fdr <- min(rmats$FDR[which(rmats$transcript_id==uTranscripts[ii])])
      
      idx <- which(background.network$exon_source==uTranscripts[ii])
      if(length(idx)>0){
        source_score[idx] <- score
        source_fdr[idx] <- fdr
      }
      
      idx <- which(background.network$exon_target==uTranscripts[ii])
      if(length(idx)>0){
        target_score[idx] <- score
        target_fdr[idx] <- fdr
      }
      
    }
    
    min_score <- pmin(source_score, target_score)
    min_fdr <- pmin(source_fdr, target_fdr)
    
    background.network$source_score <- source_score
    background.network$target_score <- target_score
    background.network$min_score <- min_score
    
    background.network$source_fdr <- source_fdr
    background.network$target_fdr <- target_fdr
    background.network$min_fdr <- min_fdr
    
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
      utrans <- unique(background.network$exon_source[which(background.network$gene_source==top.proteins[ii])])
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
  
  returnList <- list()
  returnList[[length(returnList)+1]] <- background.network
  returnList[[length(returnList)+1]] <- input.node
  names(returnList) <- c("background.network", "input.node")
  
  return(returnList)
  
}