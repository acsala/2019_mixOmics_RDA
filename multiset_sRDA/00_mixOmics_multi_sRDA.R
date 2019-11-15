################################################################################
#
#                                                                                     
#   Filename    :	  analyse_multi_sRDA.R    												  
#                                                                                     
#   Manuscript  :   BiomJ article Multiset sparse redundancy analysis for high 
#                   dimensional omics data" by Attila Csala, Michel H. Hof, and 
#                   Aeilko H. Zwinderman
#
#   Date        :   12-01-2018
#
################################################################################

library(mixOmics)
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna)

Y <- breast.TCGA$data.train$protein

list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5))

set.seed(4)
MyResult.diablo <- block.spls(X, Y, keepX=list.keepX)

X11()
plotIndiv(MyResult.diablo) ## sample plot
plotArrow (MyResult.diablo) ## variable plot
plotVar(MyResult.diablo)
plotLoadings(MyResult.diablo)
network(MyResult.diablo)


#run multi sRDA ####
# load functions
sapply(list.files(pattern="[.]R$", path="./functions/", full.names=TRUE), source)

X1 <- breast.TCGA$data.train$mrna
X2 <- breast.TCGA$data.train$mirna
X3 <- breast.TCGA$data.train$protein
      
Data <- cbind(X1,X2,X3)
EXPL_X = c(0,0,0)
RESP_Y = c(1,0,0)
EXPL_Z = c(1,1,0)
path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z)
      
# blocks of outer model
blocks = list(1:dim(X1)[2], 
                    (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
                    (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2])
                    )
      
      
modes = c("predictive","predictive", "predicted")
      
#if the analysis takes too long, reduce the nonzero grid
time_data <- system.time(
  s_satpls <- multi_sRDA(Data, path_matrix, blocks, modes,
                         scaled=T, penalization = "ust", nonzero =  c(16,17), 
                         lambda = Inf, maxiter = 100, cross_validate = F)
  )

X1_residual <- get_residuals(Dataset = X1, LV = s_satpls$scores[,1])
X2_residual <- get_residuals(Dataset = X2, LV = s_satpls$scores[,2])

time_data <- system.time(
  s_satpls_residual1 <- multi_sRDA(cbind(X1_residual,X2_residual,X3), 
                                   path_matrix, blocks, modes,
                         scaled=T, penalization = "ust", nonzero =  c(16,17), 
                         lambda = Inf, maxiter = 100, cross_validate = F)
)

#latent variabels are orthogonal to each other
cor(s_satpls$scores[,1], s_satpls_residual1$scores[,1])

str(s_satpls)
str(MyResult.diablo)

cor(s_satpls$scores[,1],s_satpls$scores[,2])^2
MyResult.diablo$variates$mRNA

cor(MyResult.diablo$variates$mRNA[,1],
    MyResult.diablo$variates$miRNA[,1])^2



x

      
      #calculate sensitivity and specificity / TPR and TNR
      
      nzero_X_positiong <- which(abs(s_satpls$weights[s_satpls$weights[,2]=="EXPL_X",3])>0)
      nzero_Y_positiong <- which(abs(s_satpls$weights[s_satpls$weights[,2]=="RESP_Y",3])>0)
      nzero_Z_positiong <- which(abs(s_satpls$weights[s_satpls$weights[,2]=="EXPL_Z",3])>0)
      
      zero_X_positiong <- which((s_satpls$weights[s_satpls$weights[,2]=="EXPL_X",3])==0)
      zero_Y_positiong <- which((s_satpls$weights[s_satpls$weights[,2]=="RESP_Y",3])==0)
      zero_Z_positiong <- which((s_satpls$weights[s_satpls$weights[,2]=="EXPL_Z",3])==0)
      
      #sensitivity/ TPR and TAVI
      sensX = sum(nzero_X_positiong %in% (p0[1]+1):(p1+p0)[1])/min(p1[1],length(nzero_X_positiong))
      sensY = sum(nzero_Y_positiong %in% (p0[2]+1):(p1+p0)[2])/min(p1[2],length(nzero_Y_positiong))
      sensZ = sum(nzero_Z_positiong %in% (p0[3]+1):(p1+p0)[3])/min(p1[3],length(nzero_Z_positiong))
      
      sens2X = sum(nzero_X_positiong %in% (p0[1]+1):(p1+p0)[1])/p1[1]
      sens2Y = sum(nzero_Y_positiong %in% (p0[2]+1):(p1+p0)[2])/p1[2]
      sens2Z = sum(nzero_Z_positiong %in% (p0[3]+1):(p1+p0)[3])/p1[3]
      
      
      #specificiy/TNR
      specX = sum(zero_X_positiong %in% 1:p0[1])/p0[1]
      specY = sum(zero_Y_positiong %in% 1:p0[2])/p0[2]
      specZ = sum(zero_Z_positiong %in% 1:p0[3])/p0[3]
      
      
      sens.m[i,]  =c(sensX,sensY,sensZ)
      sens2.m[i,] =c(sens2X,sens2Y,sens2Z)
      spec.m[i,]  =c(specX,specY,specZ)
      
      
      iter    = s_satpls$model$iter
      nonzero = s_satpls$nonzero
      lambda  = s_satpls$lambda
      

get_residuals <- function(Dataset, LV){
        
        # calculate the residuals
        calcres = function(Xcol)
          Xcol - solve(t(LV)%*%LV) %*% t(LV) %*% Xcol %*% t(LV)
        
        Res_data = apply(Dataset, 2, calcres)
        
        return(Res_data)
        
}


str(s_satpls)
str(MyResult.diablo)

s_satpls$weights[s_satpls$weights[,2]==levels(s_satpls$weights$block)[1],3]
s_satpls$model


reshape_multisRDA_output_to_mixOmics <- 
  function(mix_omics_output,
           old_rda_output, 
           components = 1){
  # function takes the mix_omics_output output and overwrites 
  #  loadings, 
  #  variates, 
  #  explained_variance
  # with old_rda_output.
  
  #overwrite loadings
  names(old_rda_output$ALPHA[[1]]) <- colnames(mix_omics_output$X)
  names(old_rda_output$ALPHA[[2]]) <- colnames(mix_omics_output$X)
  
  loadings <- list()
  loadings[["X"]] <- cbind(old_rda_output$ALPHA[[1]], old_rda_output$ALPHA[[2]])
  colnames(loadings[["X"]]) <- cbind("comp1", "comp2")
  rownames(loadings$X) <- rownames(mix_omics_output$loadings$X)
  loadings[["Y"]] <- cbind(old_rda_output$BETA[[1]], old_rda_output$BETA[[2]])
  colnames(loadings[["Y"]]) <- cbind("comp1", "comp2")
  rownames(loadings$Y) <- rownames(mix_omics_output$loadings$Y)
  mix_omics_output$loadings <- loadings
  
  #overwrite variates (scores)
  variates <- list()
  variates[["X"]] <- cbind(old_rda_output$XI[[1]], old_rda_output$XI[[2]])
  colnames(variates[["X"]]) <- cbind("comp1", "comp2")
  rownames(variates$X) <- rownames(mix_omics_output$variates$X)
  variates[["Y"]] <- cbind(old_rda_output$ETA[[1]], old_rda_output$ETA[[2]])
  colnames(variates[["Y"]]) <- cbind("comp1", "comp2")
  rownames(variates$Y) <- rownames(mix_omics_output$variates$Y)
  mix_omics_output$variates <- variates
  
  # variates explained vairance after variates and loadings are replaced
  explained_variance <- list()
  explained_variance[["X"]] <- explained_variance(mix_omics_output$X,
                                                  mix_omics_output$variates$X,
                                                  mix_omics_output$ncomp)
  explained_variance[["Y"]] <- explained_variance(mix_omics_output$Y,
                                                  mix_omics_output$variates$Y,
                                                  mix_omics_output$ncomp)
  
  mix_omics_output$explained_variance <- explained_variance
  
  return(mix_omics_output)
  
}
