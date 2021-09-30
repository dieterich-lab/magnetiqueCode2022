bin_measurements <- function(input.scores = input.scores, top = 50){
  
  input.scores$bin <- rep(1, nrow(input.scores))
  input.scores$bin[order(abs(input.scores$t), decreasing = TRUE)[1:top]] <- -1
  
  return(input.scores)
  
}