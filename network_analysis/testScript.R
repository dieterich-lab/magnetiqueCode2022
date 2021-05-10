library(readxl)
library(readr)
library(XML)
library(igraph)

load(file = "data_test/bg.RData")
# load(file = "../../DIGGER/DIGGER_Resource/dataInput.RData")
# rmats1 <- read_excel("~/Desktop/Alternative_Splicing/DIGGER/rmats_results/DCM-vs-CTRL_RI.xlsx")
# rmats2 <- read_excel("~/Desktop/Alternative_Splicing/DIGGER/rmats_results/DCM-vs-CTRL_A3SS.xlsx")
# rmats3 <- read_excel("~/Desktop/Alternative_Splicing/DIGGER/rmats_results/DCM-vs-CTRL_A5SS.xlsx")
# rmats4 <- read_excel("~/Desktop/Alternative_Splicing/DIGGER/rmats_results/DCM-vs-CTRL_MXE.xlsx")
# rmats4 <- rmats4[, -c(6,7)]
# rmats5 <- read_excel("~/Desktop/Alternative_Splicing/DIGGER/rmats_results/DCM-vs-CTRL_SE.xlsx")
# colnames(rmats1) <- colnames(rmats5)
# colnames(rmats2) <- colnames(rmats5)
# colnames(rmats3) <- colnames(rmats5)
# colnames(rmats4) <- colnames(rmats5)
# asInput <- unique(rbind(rmats1, rmats2, rmats3, rmats4, rmats5))
# save(asInput, file = "data_test/asInput.RData")
load(file = "data_test/asInput.RData")
load(file = "data_test/map.table.RData")
# load(file = "~/Desktop/MAGE-Project/Data_Exploration_Stringtie/github_scripts/v8/output/tfList50.RData")
tf_activities_all_DCM_vs_Healthy <- read_delim("data_test/tf_activities_all_DCM_vs_Healthy.csv", 
                                           ";", escape_double = FALSE, trim_ws = TRUE)


top = 50
lThresh <- -0.1
lambda1 <- 10
lambda2 <- 0.001
input.node <- NULL
mipgap = 0.05
relgap = 0.05
populate = 100
nSolutions = 20
intensity = 0
timelimit = 3600
process_log = FALSE
replace = 2
solverPath <- "~/Documents/cplex"
weightThreshold <- 10

res <- runMethod(input.scores = tf_activities_all_DCM_vs_Healthy, top = 50, 
                 rmats.input = asInput, background.network = bg, 
                 map.table = map.table, input.node = input.node)

