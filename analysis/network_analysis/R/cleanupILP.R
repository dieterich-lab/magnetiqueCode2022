cleanupILP <- function(condition=condition){
  
  if(file.exists(paste0("testFile_",condition,".lp"))){
    file.remove(paste0("testFile_",condition,".lp"))
  }
  
  if(file.exists(paste0("results_", condition, ".txt"))){
    file.remove(paste0("results_", condition, ".txt"))
  }
  
  if(file.exists("cplex.log")){
    file.remove("cplex.log")
  }
  
  if(file.exists(paste0("cplexCommand_", condition, ".txt"))){
    file.remove(paste0("cplexCommand_", condition, ".txt"))
  }
  
  AllFiles <- list.files()
  CloneFiles <- which(grepl(pattern = "clone",x = AllFiles,fixed = TRUE))
  if (length(CloneFiles)>0) {
    for (counter in seq_len(length(CloneFiles))) {
      file.remove(AllFiles[CloneFiles[counter]])
    }
  }
  
}