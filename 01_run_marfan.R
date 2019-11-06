# initialize packages ####

#for benchmarking
library(profvis)

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

# helper functions ####
source("./00_helper_functions.R")

# load data ####
#load("/trinity/home/cscala/01attila/06_msPLS/02_RealData_analyses/01_RealDataAnalysis/datasetXYZ.RData")
load("./Marfan/datasetXYZ.RData")

# Benchmarking ####

profvis({
  res_sRDA <- sRDAccp(predictor = Y, predicted = Z, nonzero = c(25), 
                      multiple_LV = T, nr_LVs = 2, penalization = "ust")
  sum(2)
})

microbenchmark(res_sRDA <- sRDAccp(predictor = Y, predicted = Z,
                                   nonzero = c(25), 
                                   multiple_LV = T, nr_LVs = 2, 
                                   penalization = "ust"), 
               spls(X = Y,Y = Z, keepX = c(25, 25), 
                    keepY = c(dim(Z)[2],dim(Z)[2]), 
                    ncomp = 2), 
               times = 10, unit = "s")

# analyse marfan data ####
#run sRDA / sCCA from sRDA pacakge
system.time(
  res_sRDA <- sRDAccp(predictor = Y, predicted = Z, nonzero = c(25), 
                      multiple_LV = T, nr_LVs = 2, penalization = "ust")
)

#run spls from mixOmics
system.time(
  res_spls <- spls(X = Y,Y = Z, keepX = c(25, 25), keepY = c(dim(Z)[2],dim(Z)[2]), ncomp = 2)  
) 
class(res_spls)
#str(res_spls)

# after obtaining results, put sRDA outputs in mixOmics' "mixo_spls" class 
res_sRDA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                            old_rda_output = res_sRDA)


#plotting ####
X11()
plotIndiv(res_spls)     ## sample plot     
plotVar(res_spls)       ## variable plot
plotLoadings(res_spls, comp = 1, size.name = rel(0.5))
X11()
cim(res_spls, comp = 1)

X11()
plotIndiv(res_sRDA)
plotVar(res_sRDA)
plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))
X11()
cim(res_sRDA, comp = 1)

# look at explained variances ####
res_spls$explained_variance
res_sRDA$explained_variance

# examine correlations of latent variates / scores
cor(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cor(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

cov(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cov(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

# examine sum abs correlation between latent variates and outcome variables
sum(abs(res_sRDA$loadings[["Y"]][,1]))
sum(abs(res_spls$loadings[["Y"]][,1]))

sum(abs(res_sRDA$loadings[["Y"]][,2]))
sum(abs(res_spls$loadings[["Y"]][,2]))

res_sRDA$loadings[["Y"]][,1][1:4]
cor(res_sRDA$variates[["X"]][,1],Z[,1:4])

sum(abs(cor(res_sRDA$variates[["X"]][,1],Z)))
sum(abs(cor(res_spls$variates[["X"]][,1],Z)^2))

