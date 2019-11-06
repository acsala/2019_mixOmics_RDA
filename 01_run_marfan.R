# initialize packages ####

#load mixOmics
library(mixOmics)

library(parallel)
library(elasticnet)
#source all sRDA code
setwd("./sRDA/R/")
file.sources = list.files(pattern = "*.r", ignore.case = T)
sapply(file.sources,source,.GlobalEnv)
setwd("../../")
getwd()


# analyse marfan data ####
load("/trinity/home/cscala/01attila/06_msPLS/02_RealData_analyses/01_RealDataAnalysis/datasetXYZ.RData")

#run sRDA / sCCA from sRDA pacakge
system.time(
  res_sRDA <- sRDA(X, Z, nonzero = c(25), multiple_LV = T, nr_LVs = 2, penalization = "ust")
)


#run spls from mixOmics
system.time(
  res_spls <- spls(X,Z, keepX = c(25, 25), keepY = c(dim(Z)[2],dim(Z)[2]))  
) 
class(res_spls)


