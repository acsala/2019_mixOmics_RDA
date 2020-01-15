# ***************************************************************
# Code to compare sPLS from mixOMics and sRDA and sCCA from sRDA on
# breast cancer data from mixOmics'package
#
# 2019-11-18

# initialize packages ####

#load mixOmics
library(mixOmics)

library(microbenchmark)
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

# run sRDA and sPLS with breast cancer data from mixOmics package####
data(breast.TCGA)

# ***************************************************************
# Start data analysis
# ***************************************************************

# extract training data and name each data frame
X <- breast.TCGA$data.train$mrna
Y <- breast.TCGA$data.train$protein
dim(X);dim(Y)
# [1] 150 200
# [1] 150 142

# Run it on other data ####
#info on the data
#https://cran.r-project.org/web/packages/PMA/PMA.pdf

library(PMA)
#prepare data
data(breastdata)

X <- t(breastdata$dna)
Y <- t(breastdata$rna)

dim(X);dim(Y)
# [1] 89 2149
# [1] 89 19672

sum(apply(X,2,sd)>0.3)
sum(apply(Y,2,sd)>0.99)

X_cleaned <- X[,apply(X,2,sd)>0.3]
Y_cleaned <- Y[,apply(Y,2,sd)>0.99]

ncomp <- 5
nr_nonz <- 10

res_all <- get_PLS_CCA_RDA_results(X = X, Y = Y,
                                   nr_nonz = nr_nonz,
                                   nr_comp = ncomp,
                                   pls_mode = "regression",
                                   penalty_mode = "ust",
                                   CCA = F)


res_sCCA <- res_all$res_sCCA
res_sRDA <- res_all$res_sRDA
res_spls <- res_all$res_spls


# Compare the variables selection ####
res_sRDA <- get_nonzero_variables(res_object = res_sRDA)
res_spls <- get_nonzero_variables(res_object = res_spls)
res_sCCA <- get_nonzero_variables(res_object = res_sCCA)

res_sRDA$nz_loading_names[["X"]]
res_spls$nz_loading_names[["X"]]
res_sCCA$nz_loading_names[["X"]]

plot(1:ncomp,get_nr_of_common_components(res_sRDA, res_spls,
                                         ncomp = ncomp),
     ylab = "Nr of common variables", xlab = "Component",
     ylim = c(0,nr_nonz), pch = 19, lwd = 2)


# we can look at explained variances, they are about the same
sum(res_spls$explained_variance$Y)
sum(res_sRDA$explained_variance$Y)
sum(res_sCCA$explained_variance$Y)

sum(res_spls$explained_variance$X)
sum(res_sRDA$explained_variance$X)
sum(res_sCCA$explained_variance$X)

# examine correlations of latent variates / scores
# these are very similar too in both methods
cor(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cor(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])


# covariance is higher in sPLS and it is standardized in RDA
# so RDA's cov is equal to its var
cov(res_spls$variates[["X"]], res_spls$variates[["Y"]])
cov(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])
cov(res_sCCA$variates[["X"]], res_sCCA$variates[["Y"]])

# examine sum squared correlation between latent variates and ####
# outcome variables, these will be the loadings at RDA and they
# are standardized for sPLS (equal to 1)
cor(res_sRDA$variates$X[,1],Y[,1:4])

res_sRDA$loadings$Y[,1][1:4]
sum(res_sRDA$loadings$Y[,1]^2)
sum(res_spls$loadings$Y[,1]^2)
sum(res_sCCA$loadings$Y[,1]^2)


# plot results #####

#X11()
plotIndiv(res_spls)     ## sample plot
plotVar(res_spls)       ## variable plot
plotLoadings(res_spls, comp = 1, size.name = rel(0.5))
cim(res_spls, comp = 1)

#X11()
plotIndiv(res_sRDA)
plotVar(res_sRDA)
plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))
cim(res_sRDA, comp = 1)

#X11()
plotIndiv(res_sCCA)
plotVar(res_sCCA)
plotLoadings(res_sCCA, comp = 1, size.name = rel(0.5))
cim(res_sCCA, comp = 1)



