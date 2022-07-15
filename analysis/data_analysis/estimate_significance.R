estimate_significance <- function(expr = expr, regulons = regulons, nperm=1000){
  
  tf_activities_stat <- dorothea::run_viper(expr, regulons,
                                            options =  list(minsize = 3, eset.filter = FALSE, 
                                                            cores = 1, verbose = FALSE, nes = TRUE))
  tf_activities_stat <- as.data.frame(tf_activities_stat)
  input.scores <- matrix(data = , nrow = nrow(tf_activities_stat), ncol = 2)
  input.scores[, 1] <- rownames(tf_activities_stat)
  input.scores[, 2] <- tf_activities_stat[, 1]
  colnames(input.scores) <- c("id", "nes")
  input.scores <- as.data.frame(input.scores)
  input.scores$id <- as.character(input.scores$id)
  input.scores$nes <- as.numeric(input.scores$nes)
  
  sif <- matrix(data = "+", nrow = nrow(regulons), ncol = 3)
  sif[, 1] <- regulons$tf
  sif[which(regulons$mor==-1), 2] <- "-"
  sif[, 3] <- regulons$target
  
  random_scores <- matrix(data = , nrow = nrow(input.scores), ncol = nperm)
  rownames(random_scores) <- input.scores$id[order(input.scores$id)]
  colnames(random_scores) <- paste0("perm_", 1:nperm)
  for(ii in 1:nperm){
    
    model.dsg <- birewire.induced.bipartite(sif, delimitators = list(negative="-", positive="+"))
    tmp <- birewire.rewire.dsg(model.dsg, verbose = FALSE, delimitators = list(negative="-", positive="+"))
    
    idx <- which(x = tmp$positive==1, arr.ind = TRUE)
    randomPos <- matrix(data = , nrow = nrow(idx), ncol = 4)
    randomPos[, 1] <- rownames(tmp$positive)[idx[, 1]]
    randomPos[, 2] <- "A"
    randomPos[, 3] <- colnames(tmp$positive)[idx[, 2]]
    randomPos[, 4] <- "1"
    
    idx <- which(x = tmp$negative==1, arr.ind = TRUE)
    randomNeg <- matrix(data = , nrow = nrow(idx), ncol = 4)
    randomNeg[, 1] <- rownames(tmp$positive)[idx[, 1]]
    randomNeg[, 2] <- "A"
    randomNeg[, 3] <- colnames(tmp$positive)[idx[, 2]]
    randomNeg[, 4] <- "-1"
    
    random <- unique(rbind(randomPos, randomNeg))
    colnames(random) <- colnames(regulons)
    random <- as.data.frame(random)
    random$mor <- as.numeric(random$mor)
    
    tf <- dorothea::run_viper(expr, random,
                              options =  list(minsize = 3, eset.filter = FALSE, 
                                              cores = 1, verbose = FALSE, nes = TRUE))
    
    random_scores[which(rownames(random_scores)%in%intersect(x = rownames(random_scores), y = rownames(tf))), ii] <- 
      as.numeric(tf)[which(rownames(tf)%in%intersect(x = rownames(random_scores), y = rownames(tf)))]
    
  }
  
  ind <- apply(random_scores, 1, function(x) all(is.na(x)))
  random_scores <- random_scores[!ind, ]
  
  input.scores <- input.scores[which(input.scores$id%in%rownames(random_scores)), ]
  
  # ecdf_fun <- function(x,perc) ecdf(x)(perc)
  
  significance <- rep(1, nrow(input.scores))
  for(ii in 1:nrow(input.scores)){
    
    ind <- which(rownames(random_scores)==input.scores$id[ii])
    ss <- as.numeric(random_scores[ind, ])
    idx2rem <- which(is.na(ss))
    if(length(idx2rem)>0){ss <- ss[-idx2rem]}
    
    require(Hmisc)
    k <- sum(abs(as.numeric(ss)) >= abs(input.scores$nes[ii]))
    significance[ii] <- zapsmall(binconf(k, length(as.numeric(ss)), method='exact'))[1]
    
  }
  
  input.scores$pval <- significance
  
  return(input.scores)
  
}