#helper functions ####
library(matrixStats)
library(colorspace)
#library(Rfast)
#library(Rcpp)
 #sourceCpp("test.cpp")
library(RColorBrewer)

regression_line_col <-  adjustcolor( "red", alpha.f = 0.5)

reshape_sRDA_output_to_mixOmics <- function(X, Y, mix_omics_output, old_rda_output){
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
    explained_variance[["X"]] <- apply(variates$X,2,function(x) sum(cor(X, x)^2)/dim(X)[2])
    explained_variance[["Y"]] <- explained_variance(Y, variates$X, ncomp)
    #explained_variance[["Y"]] <- apply(variates$X,2,function(x) sum(cor(Y, x)^2)/dim(Y)[2])

    

    
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
    res_sRDA <- reshape_sRDA_output_to_mixOmics(X, Y, mix_omics_output = res_spls,
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
    res_sCCA <- reshape_sRDA_output_to_mixOmics(X, Y, mix_omics_output = res_spls,
                                                old_rda_output = res_sCCA)
  }
  
  
  return <- list(res_sRDA = res_sRDA,
                 res_sCCA = res_sCCA,
                 res_spls = res_spls)
  
}


sRDA_mixOmics = function(X,
                         Y,
                         ncomp = 2,
                         keepX = 10,
                         keepY = dim(Y)[2],
                         scale = TRUE,
                         tol = 1e-06,
                         max.iter = 100,
                         penalty_mode = c("none", "enet", "ust"),
                         ridge_penalty = Inf,
                         cross_validate = FALSE,
                         nr_CVfolds = 10,
                         logratio = "none"
                         )
{
    #tester values
    ## X
    ## Y
    ## ncomp = 5
    ## keepX = c(20)
    ## keepY = dim(Y)[2]
    ## scale = TRUE
    ## tol = 1e-06
    ## max.iter = 100
    ## penalty_mode = c("enet")
    ## ridge_penalty = Inf
    ## cross_validate = FALSE
    ## logratio = "none"
    ## nr_CVfolds = 10

    input.X = X # save the checked X, before logratio/multileve/scale
    input.Y = Y
    
    if(ncomp > 1){
        multiple_LV = TRUE
    }else{
        multiple_LV = FALSE
    }
    
    if(scale == TRUE) {
        X <- scale(X)
        Y <- scale(Y)
    }

    # call to 'sRDA'
    result = sRDA(predictor = X,
                  predicted = Y,
                  penalization = penalty_mode,
                  ridge_penalty = ridge_penalty,
                  nonzero = keepX,
                  multiple_LV = multiple_LV,
                  nr_LVs = ncomp,
                  cross_validate = cross_validate,
                  nr_subsets = nr_CVfolds)

    # create correct names structure
    colnames <- list("X" = colnames(X), "Y" = colnames(Y))  
    names <- list ("sample" = rownames(X), "colnames" = colnames, "blocks" = c("X", "Y"))

    #make lists
    if(typeof(result$ALPHA)!="list"){
        result$ALPHA <- list(result$ALPHA)
        result$BETA <- list(result$BETA)
        }

    # create correct loadings structure
    loadings <- list("X" = do.call(cbind,result$ALPHA), "Y" =  do.call(cbind,result$BETA))
    colnames(loadings$X) <- paste0("comp", seq_len(ncol(loadings$X)))
    colnames(loadings$Y) <- paste0("comp", seq_len(ncol(loadings$Y)))

    loadings.star <- list("X" = apply(loadings$X,2,scale), "Y" = apply(loadings$Y,2,scale))
    rownames(loadings.star$X) <- rownames(loadings$X)
    rownames(loadings.star$Y) <- rownames(loadings$Y)
    
    rownames(loadings$X) = rownames(loadings.star$X) = colnames$X
    rownames(loadings$Y) = rownames(loadings.star$Y) = colnames$Y

    if(typeof(result$XI)!="list"){
        result$XI = list(result$XI)
        result$ETA = list(result$ETA)
    }
   
    # create correct variates structure
    variates <- list("X" = do.call(cbind,result$XI), "Y" =  do.call(cbind,result$ETA))
    colnames(variates$X) <- paste0("comp", seq_len(ncol(variates$X)))
    colnames(variates$Y) <- paste0("comp", seq_len(ncol(variates$Y)))

    rownames(variates$X) = rownames(variates$Y) = names$sample

    # variates explained vairance after variates and loadings are replaced
    explained_variance <- list("X" = apply(variates$X,2,function(x) sum(cor(input.X, x)^2)/dim(X)[2]),
                               "Y" = explained_variance(input.Y,
                                                        variates$X,
                                                        ncomp))

    
    # log transform from mixOmics
    #-----------------------------#
    #-- logratio transformation --#
    
    X = logratio.transfo(X=X, logratio=logratio)
    
    #as X may have changed
    if (ncomp > min(ncol(X), nrow(X)))
    stop("'ncomp' should be smaller than ", min(ncol(X), nrow(X)),
    call. = FALSE)
    
    #-- logratio transformation --#
    #-----------------------------#

    # choose the desired output from 'result'
    out = list(
        call = match.call(), 
        X = X, 
        Y = Y, 
        ncomp = ncomp,
        mode = penalty_mode,
        keepX = result$nr_nonzeros,
        ridgePenaltySelected = result$ridge_penalty,
        keepY = keepY,
        variates = variates,
        loadings = loadings,
        loadings.star = loadings.star,
        names = names,
        tol = tol,
        iter = result$nr_iterations,
        max.iter = max.iter,
        scale = scale,
        logratio = logratio,
        explained_variance = explained_variance,
        input.X = input.X,
        CV_results = result$CV_results
        #nzv = result$nzv, #not implemented
        #mat.c = result$mat.c # not implemented, comes from internal_mint.block,
        #prob for multilevel studies
    )

    class(out) = c("mixo_spls")
    # output if multilevel analysis
    
    return(invisible(out))
}


## plot_CV_results <- function(res_object,
##                             xlim = c(),
##                             ylim = c(),
##                             ylab = "y",
##                             xlab = "x"){

##     plot_cus(res_object$CV_results$mean_abs_cors[,3],
##              xlim = xlim,
##              ylim = ylim,
##              ylab = ylab,
##              xlab = xlab)

## }

## plot_cus <- function(x = c(),y = c(), 
##                      xlim = c(),
##                      ylim = c(),
##                      xlab = "x",
##                      ylab = "y",
##                      bg = "#555555AA",
##                      col = "white",
##                      plot_Yaxis = "T",
##                      plot_Xaxis = "T",
##                      xaxt = "n",
##                      yaxt = "n",
##                      type = "n",
##                      axes = FALSE,
##                      plot_points = TRUE){

##     plot(x = x, y = y, 
##          xlim = xlim,
##          ylim = ylim,
##          xlab = xlab,
##          ylab = ylab,
##          axes = axes,
##          xaxt = xaxt,
##          yaxt = yaxt,
##          type = type
##          )
##     if(plot_Xaxis) axis(side = 1, las = 1, lwd = 2)
##     if(plot_Yaxis) axis(side = 2, las = 1, lwd = 2)
    
##     if(plot_points) points(x, y, pch = 21, cex = 1.5,
##                            col = col, bg = bg,
##                            lwd = 1)
## }


plot_CV_results <- function(res_object,
                            xlim = c(),
                            ylim = c(),
                            ylab = "Sum of absolute correlations",
                            xlab = "Number of nonzeros",
                            labels_Xaxis = c(),
                            labels_Yaxis = c(),
                            labels_Xat = c(),
                            labels_Yat = c(),
                            spline_all_knots = T,
                            spline_df = T){
  
  x = 1:length(res_object$CV_results$mean_abs_cors[,3])
  y = res_object$CV_results$mean_abs_cors[,3]
  plot_cus(y = y,
           x = x,
           xlim = xlim,
           ylim = ylim,
           ylab = ylab,
           xlab = xlab,
           labels_Xaxis = res_object$CV_results$mean_abs_cors[,2],
           labels_Yaxis = labels_Yaxis,
           labels_Xat = x,
           labels_Yat = labels_Yat)

  if(spline_df == T){
    spline_df = length(x)
  }

  spline_res <- smooth.spline(x, y,
                              all.knots = spline_all_knots,
                              df= spline_df)

  lines(predict(spline_res, x)$y, col = "grey", lwd = 2)

}

plot_cus <- function(x = c(),y = c(), 
                     xlim = c(),
                     ylim = c(),
                     xlab = "x",
                     ylab = "y",
                     bg = "#555555AA",
                     col = "white",
                     plot_Yaxis = "T",
                     plot_Xaxis = "T",
                     labels_Xaxis = c(),
                     labels_Yaxis = c(),
                     labels_Xat = c(),
                     labels_Yat = c(),
                     xaxt = "n",
                     yaxt = "n",
                     type = "n",
                     axes = FALSE,
                     plot_points = TRUE){

  plot(x = x, y = y, 
       xlim = xlim,
       ylim = ylim,
       xlab = xlab,
       ylab = ylab,
       axes = axes,
       xaxt = xaxt,
       yaxt = yaxt,
       type = type
       )
    
  if(plot_points) points(x, y, pch = 21, cex = 1.5,
                         col = col, bg = bg,
                         lwd = 1)

  if(plot_Xaxis) axis(side = 1, las = 1,
                      lwd = 2, labels = labels_Xaxis, at = labels_Xat)
  if(plot_Yaxis) axis(side = 2, las = 1,
                      lwd = 2, labels = labels_Yaxis, at = labels_Yat)

  
}

