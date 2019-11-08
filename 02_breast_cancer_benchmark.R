# initialize packages ####
#load mixOmics
library(mixOmics)

library(PMA)
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
data(breastdata)

X <- t(breastdata$dna)
Y <- t(breastdata$rna)
dim(X); dim(Y)
