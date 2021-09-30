cleanup_files <- function(process_log = FALSE){
  
  if(file.exists(paste0("testFile.lp"))){
    file.remove(paste0("testFile.lp"))
  }
  
  if(file.exists(paste0("results1.txt"))){
    file.remove(paste0("results_cplex.txt"))
  }
  
  if(process_log == FALSE){
    if(file.exists("cplex.log")){
      file.remove("cplex.log")
    }
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