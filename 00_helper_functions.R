# helper functions ####
library(matrixStats)
library(Rfast)
library(Rcpp)
sourceCpp("test.cpp")

reshape_sRDA_output_to_mixOmics <- function(mix_omics_output, old_rda_output){
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


get_residuals <- function(Dataset, LV){
  
  # calculate the residuals
  calcres = function(Xcol)
    Xcol - solve(t(LV)%*%LV) %*% t(LV) %*% Xcol %*% t(LV)
  
  Res_data = apply(Dataset, 2, calcres)
  
  return(Res_data)
  
}
calc_variance <- function(vec_1,vec_2){
  
  return(vec_1 %*% vec_2)
  
}

calc_sd <- function(vec_1){
  return(sqrt(calc_variance(vec_1, vec_1)))
}

calc_cor <- function(vec_1,vec_2){
  
  return(calc_variance(vec_1,vec_2) / (calc_sd(vec_1) * calc_sd(vec_2)))
  
}

colScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  ################
  # Get the column means
  ################
  cm = colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = colSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}


deflation = function(X, y){
  # Computation of the residual matrix R
  # Computation of the vector p.
  
  #is.na.tX <- is.na(t(X))
  #save(list=ls(),file="temp3.Rdata")

  p <- crossprod(X,y) / as.vector(crossprod(y))
  
  R <- X - tcrossprod(y,p)
  return(list(p=p,R=R))
}