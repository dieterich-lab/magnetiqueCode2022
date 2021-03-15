library(CARNIVAL)
library(OmnipathR)

dir.create("output")

ff <- c("output/dcm_vs_healthy/", "output/hcm_vs_healthy/")

interactions <- import_Omnipath_Interactions()
idx1 <- which(interactions$is_directed==1)
idx2 <- which((interactions$is_stimulation+interactions$is_inhibition)==1)
idx2keep <- intersect(x = idx1, y = idx2)
interactions <- interactions[idx2keep, ]

netObj <- matrix(data = , nrow = nrow(interactions), ncol = 3)
netObj[, 1] <- interactions$source_genesymbol
netObj[which(interactions$is_stimulation==1), 2] <- "1"
netObj[which(interactions$is_inhibition==1), 2] <- "-1"
netObj[, 3] <- interactions$target_genesymbol
netObj <- as.data.frame(netObj)
netObj$V1 <- as.character(netObj$V1)
netObj$V2 <- as.character(netObj$V2)
netObj$V3 <- as.character(netObj$V3)
colnames(netObj) <- c("Source", "Sign", "Target")

load(file = "TF/tfList50.RData")

for(ii in 1:length(ff)){
  
  dir.create(ff[ii])
  
  inputObj <- matrix(data = , nrow = 1, ncol = 50)
  inputObj[1, ] <- tfList[[ii]]$NES
  colnames(inputObj) <- tfList[[ii]]$GeneID
  inputObj <- as.matrix(inputObj)
  res <- runCARNIVAL(measObj = inputObj, netObj = netObj, solverPath = "~/Desktop/cplex", 
                     solver = "cplex", timelimit = 10800,
                     dir_name = ff[ii])
  
  save(res, file = paste0(ff[ii], "res.RData"))
  
}
