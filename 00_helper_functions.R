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
  
  comp = old_rda_output$nr_latent_variables
  for (i in 1:comp){
    names(old_rda_output$ALPHA[[i]]) <- colnames(mix_omics_output$X)
  }
  
  
  #overwrite loadings
  loadings <- list()
  for (i in 1:comp){
    loadings[["X"]] <- cbind(loadings[["X"]], old_rda_output$ALPHA[[i]])
    loadings[["Y"]] <- cbind(loadings[["Y"]], old_rda_output$BETA[[i]])
  }
  
  comp_names <- c("comp1")
  if (comp > 1){
    for (i in 2:comp){
      comp_names <- cbind(comp_names,paste0("comp",i))
    }
  }
  
  colnames(loadings[["X"]]) <- comp_names
  colnames(loadings[["Y"]]) <- comp_names
  
  rownames(loadings$X) <- rownames(mix_omics_output$loadings$X)
  rownames(loadings$Y) <- rownames(mix_omics_output$loadings$Y)
  
  mix_omics_output$loadings <- loadings
  
  #overwrite variates (scores)
  variates <- list()
  for (i in 1:comp){
    variates[["X"]] <- cbind(variates[["X"]], old_rda_output$XI[[i]])
    variates[["Y"]] <- cbind(variates[["Y"]], old_rda_output$ETA[[i]])
  }
  

  colnames(variates[["X"]]) <- comp_names
  colnames(variates[["Y"]]) <- comp_names
  
  rownames(variates$X) <- rownames(mix_omics_output$variates$X)
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

get_nonzero_variables <- function(res_object){
  
  
  nonzr_positionsX <- apply(res_object$loadings[["X"]], 
                           2, function(x) which(x!=0))
  
  nonzr_positionsY <- apply(res_object$loadings[["Y"]], 
                            2, function(x) which(x!=0))
  
  res_object$nzX_pos <- nonzr_positionsX
  res_object$nzY_pos <- nonzr_positionsX
  
  comp <- res_object$ncomp

  nz_loadings <- list()
  for (i in 1:comp){
    nz_loadings[["X"]] <- cbind(nz_loadings[["X"]],
                                names(res_object$loadings[["X"]][,i][nonzr_positionsX[,i]]))
    nz_loadings[["Y"]] <- cbind(nz_loadings[["Y"]],
                                names(res_object$loadings[["Y"]][,i][nonzr_positionsY[,i]]))
  }

  res_object$nz_loading_names <- nz_loadings
  
  return(res_object)
}

get_nr_of_common_components <- function(res_object1, res_object2, ncomp){
  nr_of_common_comps <- c()
  for (i in 1:ncomp){
    if (i == 1){
      nr_of_common_comps <- cbind(sum(res_object1$nz_loading_names[["X"]][,i] %in% 
                                        res_object2$nz_loading_names[["X"]][,i]))
    } else {
      nr_of_common_comps <- cbind(nr_of_common_comps, sum(res_object1$nz_loading_names[["X"]][,i] %in% 
                                                            res_object2$nz_loading_names[["X"]][,i]))
    }
    
  }
  
  return(nr_of_common_comps)
  
}

get_PLS_CCA_RDA_results <- function(X, Y, 
                                    nr_nonz = 10, 
                                    nr_comp = 3,
                                    pls_mode = "canonical",
                                    penalty_mode = "ust",
                                    RDA = T, PLS = T, CCA = T){
  
  ncomp = nr_comp
  
  res_spls <- c("not_called")
  res_sRDA <- c("not_called")
  res_sCCA <- c("not_called")
  
  if(PLS){
    #run spls from mixOmics
    res_spls <- spls(X,Y,
                     keepX = rep(nr_nonz, ncomp),
                     keepY = rep(dim(Y)[2], ncomp),
                     ncomp = ncomp, mode = pls_mode)
  }
  
  if(RDA){
    #run sRDA / sCCA from sRDA pacakge
    res_sRDA <- sRDA(X, Y,
                     nonzero = c(nr_nonz),
                     multiple_LV = T, 
                     nr_LVs = ncomp,
                     penalization = penalty_mode)
    # after obtaining results, put sRDA outputs in mixOmics' "mixo_spls" class 
    res_sRDA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                                old_rda_output = res_sRDA)
  }
  
  if(CCA){
    # run CCA ######
    res_sCCA <- sCCA(X, Y,
                     nonzero = c(nr_nonz),
                     multiple_LV = T, 
                     nr_LVs = ncomp,
                     penalization = penalty_mode)
    # after obtaining results, put sCCA outputs in mixOmics' "mixo_spls" class 
    res_sCCA <- reshape_sRDA_output_to_mixOmics(mix_omics_output = res_spls,
                                                old_rda_output = res_sCCA)
  }
  
  
  return <- list(res_sRDA = res_sRDA,
                 res_sCCA = res_sCCA,
                 res_spls = res_spls)
  
}
