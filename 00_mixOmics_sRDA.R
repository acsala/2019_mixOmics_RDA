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


# helper functions ####
source("./00_helper_functions.R")

# run with breast cancer data ####
data(breast.TCGA)
# extract training data and name each data frame
X <- breast.TCGA$data.train$mrna
Y <- breast.TCGA$data.train$protein


#run spls from mixOmics
system.time(
  res_spls <- spls(X,Y, keepX = c(25, 25), keepY = c(dim(Y)[2],dim(Y)[2]))  
) 
class(res_spls)

#run sRDA / sCCA from sRDA pacakge
system.time(
  res_sRDA <- sRDA(X, Y, nonzero = c(25), multiple_LV = T, nr_LVs = 2, penalization = "ust")
)



# after obtaining results, put sRDA outputs in mixOmics' "mixo_spls" class 
res_sRDA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                            old_rda_output = res_sRDA)

#plotting
plotIndiv(res_spls)     ## sample plot     
plotVar(res_spls)       ## variable plot
plotLoadings(res_spls, comp = 1, size.name = rel(0.5))
plotIndiv(res_sRDA)
plotVar(res_sRDA)
plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))

?cim

# look at explained variances
res_spls$explained_variance
res_sRDA$explained_variance

# examine correlations of latent variates / scores
cor(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cor(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

cov(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cov(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

# examine sum squared correlation between latent variates and outcome variables
sum(cor(res_sRDA$variates[["X"]][,1], Y)^2)
sum(cor(res_spls$variates[["X"]][,1], Y)^2)

sum(cor(res_sRDA$variates[["X"]][,2], Y)^2)
sum(cor(res_spls$variates[["X"]][,2], Y)^2)



str(res_spls)
str(res_sRDA)


# run it on simple data ####
data(nutrimouse)
X <- nutrimouse$gene  
Y <- nutrimouse$lipid
dim(X); dim(Y)

#run spls from mixOmics
system.time(
  res_spls <- spls(X,Y, keepX = c(25, 25), keepY = c(dim(Y)[2],dim(Y)[2]))  
) 
class(res_spls)

#run sRDA / sCCA from sRDA pacakge
system.time(
  res_sRDA <- sRDA(X, Y, nonzero = c(25), multiple_LV = T, nr_LVs = 2, penalization = "ust")
)


# after obtaining results, put sRDA outputs in mixOmics' "mixo_spls" class 
res_sRDA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                            old_rda_output = res_sRDA)

#plotting
plotIndiv(res_spls)     ## sample plot     
plotVar(res_spls)       ## variable plot
plotLoadings(res_spls, comp = 1, size.name = rel(0.5))
plotIndiv(res_sRDA)
plotVar(res_sRDA)
plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))

# look at explained variances
res_spls$explained_variance
res_sRDA$explained_variance

# examine correlations of variates / scores
cor(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cor(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

var(res_spls$variates[["X"]], res_spls$variates[["Y"]])
calc_variance(res_sRDA$variates[["X"]][,1], res_sRDA$variates[["Y"]][,1])/ (dim(X)[1]-1)
var(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])

str(res_spls)
str(res_sRDA)

example(spls)

