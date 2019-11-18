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

#run spls from mixOmics
res_spls <- spls(X,Y,
                 keepX = c(25, 25),
                 keepY = c(dim(Y)[2],dim(Y)[2]))

#run sRDA / sCCA from sRDA pacakge
res_sRDA <- sRDAccp(X, Y,
                    nonzero = c(25),
                    multiple_LV = T, nr_LVs = 2,
                    penalization = "ust")


#plot the results from the breast cancer data####
# after obtaining the results on the breast cancer data from mixOmics,
# we can plot with the results of sRDA with mixOmics' plots too
class(res_spls)

# after obtaining results, put sRDA outputs in mixOmics' "mixo_spls" class 
res_sRDA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                            old_rda_output = res_sRDA)

X11()
plotIndiv(res_spls)     ## sample plot     
plotVar(res_spls)       ## variable plot
plotLoadings(res_spls, comp = 1, size.name = rel(0.5))
cim(res_spls, comp = 1)

X11()
plotIndiv(res_sRDA)
plotVar(res_sRDA)
plotLoadings(res_sRDA, comp = 1, size.name = rel(0.5))
cim(res_sRDA, comp = 1)


# we can look at explained variances, they are about the same
res_spls$explained_variance
# $X
# comp 1    comp 2 
# 0.1744936 0.1244934 
# 
# $Y
# comp 1    comp 2 
# 0.1259336 0.1266823 

res_sRDA$explained_variance
# $X
# comp 1    comp 2 
# 0.1792986 0.1242932 
# 
# $Y
# comp 1    comp 2 
# 0.1264332 0.1325913 

# examine correlations of latent variates / scores
# these are very similar too in both methods
cor(res_spls$variates[["X"]], res_spls$variates[["Y"]])
#         comp1         comp2
# comp1  0.90191517 -2.711205e-17
# comp2 -0.07072779  5.618405e-01
cor(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])
#           comp1      comp2
# comp1 0.90898771 0.08819028
# comp2 0.04641389 0.56621434

# covariance is higher in sPLS and it is standardized in RDA
# so RDA's cov is equal to its var
cov(res_spls$variates[["X"]], res_spls$variates[["Y"]])
#         comp1         comp2
# comp1 10.9556301 -2.947167e-16
# comp2 -0.7918976  5.629406e+00

cov(res_sRDA$variates[["X"]], res_sRDA$variates[["Y"]])
#           comp1      comp2
# comp1 0.90898771 0.08819028
# comp2 0.04641389 0.56621434

# examine sum squared correlation between latent variates and ####
# outcome variables, these will be the loadings at RDA and they
# are standardized for sPLS (equal to 1)
cor(res_sRDA$variates$X[,1],Y[,1:4])
#         14-3-3_epsilon    4E-BP1 4E-BP1_pS65 4E-BP1_pT37
# [1,]     -0.1234524 0.2773921    0.424042   0.2000258
res_sRDA$loadings$Y[,1][1:4]
# 14-3-3_epsilon         4E-BP1    4E-BP1_pS65    4E-BP1_pT37 
# -0.1234524      0.2773921      0.4240420      0.2000258

sum(res_sRDA$loadings$Y[,1]^2)
# [1] 14.26879
sum(res_spls$loadings$Y[,1]^2)
# [1] 1

# run CCA ######
res_sCCA <- sCCA(X, Y,
                    nonzero = c(25),
                    multiple_LV = T, nr_LVs = 2,
                    penalization = "ust")

# after obtaining results, put sCCA outputs in mixOmics' "mixo_spls" class 
res_sCCA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                            old_rda_output = res_sCCA)

X11()
plotIndiv(res_sCCA)
plotVar(res_sCCA)
plotLoadings(res_sCCA, comp = 1, size.name = rel(0.5))
cim(res_sCCA, comp = 1)

cov(res_sCCA$variates[["X"]], res_sCCA$variates[["Y"]])
#         comp1      comp2
# comp1 0.90898771 0.08819028
# comp2 0.04641389 0.56621434

