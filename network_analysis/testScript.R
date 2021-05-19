library(readxl)
library(readr)
library(XML)
library(igraph)

## Source R functions
files.sources = paste0("R/", list.files("R/"))
sapply(files.sources, source)

## Load data
load(file = "data_test/bg.RData")
load(file = "data_test/rmats.RData")
load(file = "data_test/map.table.RData")
# load(file = "~/Desktop/MAGE-Project/Data_Exploration_Stringtie/github_scripts/v8/output/tfList50.RData")
tf_activities_all_DCM_vs_Healthy <- read_delim("data_test/tf_activities_all_DCM_vs_Healthy.csv", 
                                           ";", escape_double = FALSE, trim_ws = TRUE)

## Set parameters
top = 50
pValThresh <- 0.05
lambda1 <- 10
lambda2 <- 0.1
input.node <- NULL
mipgap = 0.05
relgap = 0.05
populate = 3
nSolutions = 100
intensity = 0
timelimit = 3600
process_log = FALSE
replace = 2
solverPath <- "/beegfs/homes/egjerga/cplex"
# solverPath <- "~/Downloads/cplex"
weightThreshold <- 10
input.scores = tf_activities_all_DCM_vs_Healthy
top = 50 
rmats.input = rmats
background.network = bg
map.table = map.table
input.node = input.node
constraints3 = "data_test/cc3.RData"
constraints4 = "data_test/cc4.RData"
constraints5 = "data_test/cc5.RData"
constraints6 = "data_test/cc6.RData"

## Obtain and slve solutions
solutionList <- runMethodPar(input.scores = input.scores, rmats.input = rmats.input, 
                             background.network = background.network, map.table = map.table, 
                             solverPath = solverPath, input.node = input.node, pValThresh = 0.05, 
                             top = 50, lambda1 = lambda1, lambda2 = lambda2, 
                             nSolutions = nSolutions, constraints3 = constraints3, 
                             constraints4 = constraints4, constraints5 = constraints5, 
                             constraints6 = constraints6)

save(solutionList, file = "solutionList.RData")
