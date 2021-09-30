set.seed(1234)

library(readxl)
library(readr)
library(XML)
library(igraph)
library(foreach)
library(doParallel)
library(ranger)
library(palmerpenguins)
library(tidyverse)
library(kableExtra)

## Source R functions
files.sources = paste0("R/", list.files("R/"))
sapply(files.sources, source)

## Load data
load(file = "data_set/bg.RData")
load(file = "data_set/map.table.RData")

#
tfActList <- list()

tf_activities_all_DCM_vs_Healthy <- read_delim("data_set/tf_activities_all_DCM_vs_Healthy.csv", 
                                               ";", escape_double = FALSE, trim_ws = TRUE)
tfActList[[length(tfActList)+1]] <- tf_activities_all_DCM_vs_Healthy
tfActList[[length(tfActList)+1]] <- tf_activities_all_DCM_vs_Healthy

tf_activities_all_HCM_vs_Healthy <- read_delim("data_set/tf_activities_all_HCM_vs_Healthy.csv", 
                                               ";", escape_double = FALSE, trim_ws = TRUE)
tfActList[[length(tfActList)+1]] <- tf_activities_all_HCM_vs_Healthy
tfActList[[length(tfActList)+1]] <- tf_activities_all_HCM_vs_Healthy

#
rmatsList <- list()
pValThresh <- 0.05

load(file = "data_set/rmats_dcm.RData")
rmats$FDR <- 1
rmats$IncLevelDifference <- 0
rmatsList[[length(rmatsList)+1]] <- rmats

load(file = "data_set/rmats_dcm.RData")
rmats$IncLevelDifference[which(rmats$FDR>pValThresh)] <- 0
rmatsList[[length(rmatsList)+1]] <- rmats

load(file = "data_set/rmats_hcm.RData")
rmats$FDR <- 1
rmats$IncLevelDifference <- 0
rmatsList[[length(rmatsList)+1]] <- rmats

load(file = "data_set/rmats_hcm.RData")
rmats$IncLevelDifference[which(rmats$FDR>pValThresh)] <- 0
rmatsList[[length(rmatsList)+1]] <- rmats

## Set parameters
top = 50
lambda1 <- 10
lambda2 <- 1
input.node <- NULL
mipgap = 0.05
relgap = 0.05
populate = 1000
nSolutions = 100
intensity = 0
timelimit = 36000
process_log = FALSE
replace = 1
solverPath <- "/beegfs/homes/egjerga/cplex"
# solverPath <- "~/Downloads/cplex"
weightThreshold <- 10
background.network = bg
map.table = map.table

all_constraints = c("data_set/dcm_constraints/allC_ctrl.RData", 
                    "data_set/dcm_constraints/allC_as.RData",
                    "data_set/hcm_constraints/allC_ctrl.RData", 
                    "data_set/hcm_constraints/allC_as.RData")

background_network_path = c("data_set/dcm_constraints/bg_post_ctrl.RData", 
                            "data_set/dcm_constraints/bg_post_as.RData",
                            "data_set/hcm_constraints/bg_post_ctrl.RData",
                            "data_set/hcm_constraints/bg_post_as.RData")

objective_function_path = c("data_set/dcm_constraints/objective_function_ctrl.RData",
                            "data_set/dcm_constraints/objective_function_as.RData",
                            "data_set/hcm_constraints/objective_function_ctrl.RData",
                            "data_set/hcm_constraints/objective_function_as.RData")

registerDoParallel(cores=4)

## Obtain and solve solutions - 1:4
res <- foreach(ii = 1:4) %dopar% {
  
  runMethod(input.scores = tfActList[[ii]], 
            rmats.input = rmatsList[[ii]], 
            background.network = bg, 
            map.table = map.table, 
            solverPath = solverPath, 
            input.node = NULL, 
            pValThresh = pValThresh, 
            top = top, 
            lambda1 = lambda1, 
            lambda2 = lambda2, 
            mipgap = mipgap, 
            relgap = relgap, 
            populate = populate, 
            nSolutions = nSolutions, 
            timelimit = timelimit, 
            all_constraints = all_constraints[ii], 
            background_network_path = background_network_path[ii], 
            objective_function_path = NULL, 
            condition = ii)
  
}

save(res, file = "all_solutions_1_2_3_4.RData")

