# helper functions ####
reshape_sRDA_output_to_mixOmics <- function(mix_omics_output, old_rda_output){
  # function takes the mix_omics_output output and overwrites 
  #  loadings, 
  #  variates, 
  #  explained_variance
  # with old_rda_output.
  
  #overwrite loadings
  rownames(old_rda_output$ALPHA[[1]]) <- colnames(mix_omics_output$X)
  rownames(old_rda_output$ALPHA[[2]]) <- colnames(mix_omics_output$X)
  
  loadings <- list()
  loadings[["X"]] <- cbind(old_rda_output$ALPHA[[1]], old_rda_output$ALPHA[[2]])
  colnames(loadings[["X"]]) <- cbind("comp1", "comp2")
  loadings[["Y"]] <- cbind(old_rda_output$BETA[[1]], old_rda_output$BETA[[2]])
  colnames(loadings[["Y"]]) <- cbind("comp1", "comp2")
  mix_omics_output$loadings <- loadings
  
  #overwrite variates (scores)
  variates <- list()
  variates[["X"]] <- cbind(old_rda_output$XI[[1]], old_rda_output$XI[[2]])
  colnames(variates[["X"]]) <- cbind("comp1", "comp2")
  variates[["Y"]] <- cbind(old_rda_output$ETA[[1]], old_rda_output$ETA[[2]])
  colnames(variates[["Y"]]) <- cbind("comp1", "comp2")
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

plot_double <- function(plot_1,plot_2){
  
  par(mfrow = c(2,1))
  plot(plot_1)
  plot(plot_2)
  par(mfrow = c(1,1))
  
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

