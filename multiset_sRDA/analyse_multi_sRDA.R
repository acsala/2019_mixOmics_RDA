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

### 1, Replicate the simulation study as described in Section 3.2 in the manuscript, including Table 1 ####

  # load functions
  sapply(list.files(pattern="[.]R$", path="./functions/", full.names=TRUE), source)
  
  set.seed(4)
  
  #Repeat simulation 20 times, which will take about 40 minutes to run. Change this to 1000 if you'd like to get 
  #the exact results of the manuscript
  nr_of_simulations <- 20
  
  # create table template as in manuscript Table 1
  table_1 <- matrix(data = rep(0,3*3*3),ncol = 3)
  colnames(table_1) <- c("n=100", "n=250", "n=500")
  rownames(table_1) <- c("TPR_X1", "TPR_X2", "TPR_X3",
                         "TAVI_X1", "TAVI_X2", "TAVI_X3",
                         "TNR_X1", "TNR_X2", "TNR_X3")
  
  # 3 studies with different sample sizes
  n = c(100,250,500)
  
  elapsed_time <- system.time(
  for (n_th in 1:3){
    
    ## parameters for data generation
    N = n[n_th]   # number of individuals/ sample size
    k = 4      # number of datasets
    m = 2      # number of latent variables (LV's) per dataset
    # number of irrelevant variables per dataset
    p0=c(1000,500,200,10)
    # number of variables associated with the LV's
    p1=c(10,8,8,2)
    
    sens.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
    sens2.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
    spec.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
    ridge_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations),ncol = 4)
    lasso_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations),ncol = 4)
    
    start_pos <- 1
    
    for (i in start_pos:nr_of_simulations){
      
      print("nr of simulation")
      print(i)    
      
      # for replication
      nr_seed = runif(1, 10, 10^8)
      set.seed(nr_seed)
      
      multiblockdata <- multiset_data_generator(N, k,m,
                                                p0,
                                                p1)
      X1 <- multiblockdata$X1
      X2 <- multiblockdata$X2
      X3 <- multiblockdata$X3
      X4 <- multiblockdata$X4
      
      Data <- cbind(X1,X2,X3,X4)
      EXPL_X = c(0,0,0,0)
      RESP_Y = c(1,0,0,0)
      EXPL_Z = c(0,1,0,0)
      RESP_V = c(0,0,1,0)
      path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)
      
      # blocks of outer model
      blocks = list(1:dim(X1)[2], 
                    (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
                    (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]),
                    (dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+1):
                      (dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+dim(X4)[2]))
      
      
      modes = c("predictive","predictive", "predictive", "predicted")
      
      #if the analysis takes too long, reduce the nonzero grid
      time_data <- system.time(
        s_satpls <- multi_sRDA(Data, path_matrix, blocks, modes,
                               scaled=T, penalization = "ust", nonzero =  c(15,10,8,5), 
                               lambda = Inf, maxiter = 100, cross_validate = T)
      )
      
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
      
    }
    
    print("Mean TPR:")
    print(apply(sens2.m,2,mean))
    table_1[1:3,n_th] <- apply(sens2.m,2,mean)
    
    print("Mean TAVI:")
    print(apply(sens.m,2,mean))
    table_1[4:6,n_th] <- apply(sens.m,2,mean)
    
    print("Mean TNR:")
    print(apply(spec.m,2,mean))
    table_1[7:9,n_th] <- apply(spec.m,2,mean)
  
  }
  )
  
  print(table_1)





### 2, Replicate the simulation study with size of real Marfan data####

  # load functions
  sapply(list.files(pattern="[.]R$", path="./functions/", full.names=TRUE), source)
  
  set.seed(4)
  
  N = 37    # number of individuals
  k = 4      # number of datasets
  m = 2      # number of latent variables (LV's) per dataset
  # number of irrelevant variables per dataset
  p0=c(36000,18000,47,10)
  # number of variables associated with the LV's
  p1=c(100,100,80,2)
  
  
  #repeat simulation 10 times, Change this to 100 if you'd like to get 
  #the exact results of the manuscript
  nr_of_simulations <- 10
  start_pos <- 1
  
  sens.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations-start_pos)+1,ncol = 3)
  sens2.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations-start_pos)+1,ncol = 3)
  spec.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations-start_pos)+1,ncol = 3)
  ridge_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations-start_pos)+1,ncol = 4)
  lasso_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations-start_pos)+1,ncol = 4)
  
  for (i in start_pos:nr_of_simulations){
    
    print("nr of simulation")
    print(i)    
    
    # for replication
    nr_seed = runif(1, 10, 10^8)
    set.seed(nr_seed)
    
    multiblockdata <- multiset_data_generator(N, k,m,
                                              p0,
                                              p1)
    X1 <- multiblockdata$X1
    X2 <- multiblockdata$X2
    X3 <- multiblockdata$X3
    
    Data <- cbind(X1,X2,X3)
    EXPL_X = c(0,0,0)
    RESP_Y = c(1,0,0)
    EXPL_Z = c(1,1,0)
    path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z)
    
    # blocks of outer model
    blocks = list(1:dim(X1)[2], 
                  (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
                  (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]))
    
    
    modes = c("predictive","predictive", "predicted")
    
    #if the analysis takes too long, reduce the nonzero grid
    time_data <- system.time(
      s_satpls <- multi_sRDA(Data, path_matrix, blocks, modes,
                             scaled=T, penalization = "ust", nonzero =  c(150,100,80,50), 
                             lambda = Inf, maxiter = 100, cross_validate = T)
    )
    
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
  
  }
  
  print("TAVI_X1, TAVI_X2:")
  print(apply(sens.m,2,mean)[1:2])
  print("TPR_X1, TPR_X2:")
  print(apply(sens2.m,2,mean)[1:2])
  print("TNR_X1, TNR_X2:")
  print(apply(spec.m,2,mean)[1:2])
  
  str(s_satpls)




### 3, Replicate Table 2 from manuscript ####
  
  #The Marfan data we used in the manuscript is confidental and unfortunately not publicly available.
  #We provide the final model that was obtained from the analysis, as described in Section 4 in the manuscript.

  #The structure of the data object obtained from the multi_sRDA function is identical to the one that we obtained
  #in the previous step where we replicated the size of real Marfan data for the second part of the simulation study.

  #load the data
  load("./data/RealDataResults.RData")
  
  #extract  alpha weights for gene expression
  GENEXP_alphas <- s_satpls2$outer_model[s_satpls2$outer_model[,2]=="EXPRES",3]
  res2.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]
  res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
  position_GENEXP_alphas <- which(abs(s_satpls2$outer_model[s_satpls2$outer_model[,2]=="EXPRES",3])>0)
  nonzero_GENEXP_alphas <- res2.outer[position_GENEXP_alphas,]
  plot_EXPRESSION_alphas <- res2.outer[position_GENEXP_alphas,][,c(1,3)]
  plot_EXPRESSION_alphas[,2] <- abs(plot_EXPRESSION_alphas[,2])
  
  print(plot_EXPRESSION_alphas)

  
### 4, Replicate Table 3 from manuscript ####
  
  #The Marfan data we used in the manuscript is confidental and unfortunately not publicly available.
  #We provide the final model that was obtained from the analysis, as described in Section 4 in the manuscript.
  
  #load the data
  load("./data/RealDataResults.RData")
  
  #extract alpha weights for methylation
  METHYL_alphas <- s_satpls2$outer_model[s_satpls2$outer_model[,2]=="METHYL",3]
  res1.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
  
  position_METHYL_alphas <- which(abs(s_satpls2$outer_model[s_satpls2$outer_model[,2]=="METHYL",3])>0)
  nonzero_METHYL_alphas <- res1.outer[position_METHYL_alphas,]
  plot_METHYL_alphas <- nonzero_METHYL_alphas[,c(1,3)]
  plot_METHYL_alphas[,2] <- abs(plot_METHYL_alphas[,2])
  
  print(plot_METHYL_alphas)
  
### 5, Replicate Table 4 from manuscript ####
  
  #The Marfan data we used in the manuscript is confidental and unfortunately not publicly available.
  #We provide the final model that was obtained from the analysis, as described in Section 4 in the manuscript.
  
  #load the data
  load("./data/RealDataResults.RData")
  
  #extract the beta expression weights that are higher than 0.66
  res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
  highest_methyl_betas <- (abs(res2.inner[,"METHYL"])>0.66)
  
  highest_EXPRESSION_beta <- res2.inner[highest_methyl_betas,]
  plot_EXPRESSION_betas <- highest_EXPRESSION_beta[,c(1,3)]
  plot_EXPRESSION_betas[,2] <- abs(plot_EXPRESSION_betas[,2])
  
  
  #extract the beta cytokine weights that are higher than 0.366
  res3.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="CYTO"),]
  highest_cyto_betas <- abs(res3.inner[,"CYTO"])>0.366
  plot_CYTO_betas <- res3.inner[highest_cyto_betas,][,c(1,5)]
  
  plot_CYTO_betas[,2] <- abs(plot_CYTO_betas[,2])
  
  print(plot_EXPRESSION_betas)
  print(plot_CYTO_betas)
  
  
  
### 6, Replicate Figure 2 from manuscript ####
  
  #The Marfan data we used in the manuscript is confidental and unfortunately not publicly available.
  #We provide the final model that was obtained from the analysis, as described in Section 4 in the manuscript.
  
  #load the sum of absolute correlations stored from the 10-fold cross validation 
  load("./data/sum_abs_crors.m.RData")

  #set graphical parameters for plot
  par(mar=c(4, 4, 1, 4) + 0.1)
  
  # get the maxmimum values from the sum of absolute correlations
  max_sumabsZ <- sum_abs_crors.m[sum_abs_crors.m[,2] == max(sum_abs_crors.m[,2]),]
  max_sumabsY <- sum_abs_crors.m[sum_abs_crors.m[,1] == max(sum_abs_crors.m[,1]),]
  
  
  plot(x = sum_abs_crors.m[order(sum_abs_crors.m[,3]),3], 
       y = sum_abs_crors.m[order(sum_abs_crors.m[,3]),1], 
       ylim = c(4500,7000),
       xlab = expression(paste(lambda[1], " penalty")),
       ylab="",
       axes = F,
       pch = 16, las = 1,
       cex.lab = 1)
  
  axis(4, ylim=c(0,10),col="black",las=1)  ## las=1 makes horizontal labels
  axis(1, ylim=c(0,10),col="black",las=1,line=0)
  
  mtext("Sum of absolute correlations",
        side=2,
        line=2.5)
  
  lines(smooth.spline(sum_abs_crors.m[order(sum_abs_crors.m[,3]),3], 
                      sum_abs_crors.m[order(sum_abs_crors.m[,3]),1], 
                      df = 50), lty = 1, col = "grey")
  
  abline(v = max_sumabsY[3], col = "grey")
  axis(1, at = max_sumabsY[3], labels = max_sumabsY[3], 
       col="grey",col.axis="black",las=1, line = 1, tck=0.06)
  axis(1, at = max_sumabsY[3], labels = "", 
       col="grey",col.axis="black",las=1, line = 1, tck=-0.06)
  
  box()
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot the second plot and put axis scale on right
  plot(sum_abs_crors.m[order(sum_abs_crors.m[,3]),3], 
       sum_abs_crors.m[order(sum_abs_crors.m[,3]),2],
       xlab="", ylab="", ylim = c(6,10),
       axes=FALSE, col = "red",
       pch = 17)
  
  lines(smooth.spline(sum_abs_crors.m[order(sum_abs_crors.m[,3]),3], 
                      sum_abs_crors.m[order(sum_abs_crors.m[,3]),2], 
                      df = 70), lty = 1, col = "red")
  
  axis(2, ylim=c(0,10), col="red",col.axis="red",las=1)
  
  abline(v = max_sumabsZ[3], col = "red")
  axis(1, at = max_sumabsZ[3], labels = max_sumabsZ[3], 
       col="red",col.axis="black",las=1, line = 1, tck=0.06)
  axis(1, at = max_sumabsZ[3], labels = "", 
       col="red",col.axis="black",las=1, line = 1, tck=-0.06)

  
  
### 7, Replicate Figure 3 from manuscript ####
  
  #The Marfan data we used in the manuscript is confidental and unfortunately not publicly available.
  #We provide the final model that was obtained from the analysis, as described in Section 4 in the manuscript.
  
  #load the data; the final model of the real data analysis is stored in RealDataResults.RData, 
  #aggregate.m.RData stores the sum of absolulte correlations obtained through the permutation study,
  #and Bootstrap_run_100.RData contains the sum of absolute correlations from the bootstrapping, see Section 4
  load("./data/RealDataResults.RData")
  load("./data/aggregate.m.RData")
  load("./data/Bootstrap_run_100.RData")
  
  aggregate.m <- aggregate.m[1:100,]
  to_plot <- aggregate.m[,1]
  
  #the sum of absolute correlation between the selected alpha weights for the methylation markers and the gene
  #expression dataset can be calculated as 
  #
  #sum(abs(cor(s_satpls2$scores[,"METHYL"],Y)))
  #
  #but since we cannot supply the original data sets we give the value here
  real_value <-  6842.32
  bootstrap_data <- sum_abs_crors.m[,1]
  
  hist(to_plot,main = "",
       xlab = expression(paste("Sum of absolute correlations")),
       xlim = c(0,9800),breaks = 100,
       col="grey")
  
  points(bootstrap_data, y = sample(seq(4,5,by = 0.1),100,replace = T),col="red",pch=20)
  abline(v = real_value, col = "red", lwd = 3)
  axis(1, at = real_value, labels = real_value, 
       col="red",col.axis="black",las=1, line = 1, tck=0.1, lwd= 3)
  abline(v = quantile(bootstrap_data,probs = c(0.025,0.975)), col = "red", lty = 2, lwd = 3)
  

  
### 8, Replicate Figure 4 from manuscript ####
  
  # load functions
  sapply(list.files(pattern="[.]R$", path="./functions/", full.names=TRUE), source)

  #load the data
  load("./data/RealDataResults.RData")
  
  #extract  alpha weights for gene expression
  GENEXP_alphas <- s_satpls2$outer_model[s_satpls2$outer_model[,2]=="EXPRES",3]
  res2.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "EXPRES"),]
  res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
  position_GENEXP_alphas <- which(abs(s_satpls2$outer_model[s_satpls2$outer_model[,2]=="EXPRES",3])>0)
  nonzero_GENEXP_alphas <- res2.outer[position_GENEXP_alphas,]
  plot_EXPRESSION_alphas <- res2.outer[position_GENEXP_alphas,][,c(1,3)]

  #extract alpha weights for methylation
  METHYL_alphas <- s_satpls2$outer_model[s_satpls2$outer_model[,2]=="METHYL",3]
  res1.outer <- s_satpls2$outer_model[which(s_satpls2$outer_model[,2] == "METHYL"),]
  position_METHYL_alphas <- which(abs(s_satpls2$outer_model[s_satpls2$outer_model[,2]=="METHYL",3])>0)
  nonzero_METHYL_alphas <- res1.outer[position_METHYL_alphas,]
  plot_METHYL_alphas <- nonzero_METHYL_alphas[,c(1,3)]
  
  #extract the beta expression weights that are higher than 0.66
  res2.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="EXPRES"),]
  highest_methyl_betas <- (abs(res2.inner[,"METHYL"])>0.66)
  highest_EXPRESSION_beta <- res2.inner[highest_methyl_betas,]
  plot_EXPRESSION_betas <- highest_EXPRESSION_beta[,c(1,3)]
  
  #extract the beta cytokine weights that are higher than 0.366
  res3.inner <-  s_satpls2$crossloadings[which(s_satpls2$crossloadings[,2]=="CYTO"),]
  highest_cyto_betas <- abs(res3.inner[,"CYTO"])>0.366
  plot_CYTO_betas <- res3.inner[highest_cyto_betas,][,c(1,5)]
  
  #chose the weights above 0.0082
  plot_METHYL_alphas <- plot_METHYL_alphas[(abs(plot_METHYL_alphas$weight)>0.0082),]
  
  #chose all the betas
  plot_EXPRESSION_betas$METHYL <- abs(plot_EXPRESSION_betas$METHYL)
  
  #save the canonical vectors
  can.vecs = list(Methylation= t(t(plot_METHYL_alphas$weight)), 
                  Expression= t(t(plot_EXPRESSION_betas$METHYL)), 
                  Cytokine= t(t(plot_CYTO_betas$CYTO)))
  
  #assign the new vectors with the old names
  rownames(can.vecs[[1]]) = plot_METHYL_alphas[,1]
  rownames(can.vecs[[2]]) = plot_EXPRESSION_betas[,1]
  rownames(can.vecs[[3]]) = plot_CYTO_betas[,1]
  
  #correlation between latent variables
  cor_12 <- abs(cor(s_satpls2$scores[,1],s_satpls2$scores[,2]))
  cor_23 <- abs(cor(s_satpls2$scores[,2],s_satpls2$scores[,3]))
  cor_13 <- abs(cor(s_satpls2$scores[,1],s_satpls2$scores[,3]))
  
  #store the correlations in a correlation matrix
  correlations = diag(3)
  correlations[correlations ==0] = 0.9
  correlations[1,2] = cor_12
  correlations[2,1] = cor_12
  correlations[1,3] = cor_13
  correlations[3,1] = cor_13
  correlations[2,3] = cor_23
  correlations[3,2] = cor_23
  
  rownames(correlations) = c("Methylation", "Expression", "Cytokine")
  colnames(correlations) = c("Methylation", "Expression", "Cytokine")
  
  #scale the colors according the weights
  col1 <- rgb(188,206,224, maxColorValue = 255)
  par(bg = col1)
  a = c()
  for(i in 1:length(can.vecs))
  {	
    tmp = can.vecs[[i]]
    tmp = tmp[abs(tmp) >0,]
    a = rbind(a, data.frame(names=names(tmp), dataset=names(can.vecs)[i], value=tmp))
  }
  
  #create the angles according to the number of variables to be plotted
  size.plot=14
  a$angle = seq(0,360-(360/nrow(a)), length=nrow(a))
  a$pi = seq(0, 2*pi-(2*pi/nrow(a)), length=nrow(a))
  angle.inner = tapply(a$pi, a$dataset, mean)
  innercircle = t(sapply(angle.inner, function(x) cbind(4*sin(x), 4*cos(x))))
  
  #color matrix
  a$color = round(abs(a$value)*1000)
  a$color  = a$color - min(a$color) + 1
  graph.col = rainbow(max(a$color), s=0.7, v=0.9, start=0, end=0.4)

  #Create the canvas
  par(mar=c(0,0,0,0), xpd=TRUE)
  plot(c(-size.plot,size.plot+4), c(-size.plot,size.plot), type="n",axes=FALSE, xlab="", ylab="" )
  
  #plot the lines with colors
  for(j in 1:nrow(a))
  {
    z=(10+abs(a$value[j])*3)
    
    x1 = sin(a$pi[j])*z
    y1 = cos(a$pi[j])*z
    canvec.location = innercircle[as.numeric(a$dataset[j]),]
    lines(c(x1, canvec.location[1]), c(y1, canvec.location[2]), col=graph.col[a$color[j]])
  }
  
  #plot the weights with variable names
  for(j in 1:nrow(a))
  {
    z=(10+abs(a$value[j])*3)
    
    x1 = sin(a$pi[j])*z
    y1 = cos(a$pi[j])*z
    
    x2 = sin(a$pi[j])*(z+2)
    y2 = cos(a$pi[j])*(z+2)
    
    punt = points(x1,y1, pch=16, col=as.numeric(a$dataset[j]))
    text(x2,y2, a$names[j], srt=90-a$angle[j], cex=0.7, col=1)
  }	
  
  #plot the correlation between the latent variables
  for(i in 1:(nrow(innercircle)-1))
  {
    for(j in (i+1):nrow(innercircle))
    {
      lines(c(innercircle[i,][1],innercircle[j,][1]),  c(innercircle[i,][2],innercircle[j,][2]))  			
      
      x = (innercircle[j,][1] + innercircle[i,][1])/2
      y = (innercircle[j,][2] + innercircle[i,][2])/2
      
      boxed.labels(x,y,round(correlations[i,j],3), cex=0.5,bg=col1,border=NA) 
    }
  }
  
  #create placeholder circles for latent variables
  for(i in 1:nrow(innercircle))
  {
    points(innercircle[i,][1], innercircle[i,][2], pch=16, col=col1, cex=4)
    points(innercircle[i,][1], innercircle[i,][2], pch=1, cex=4)		
  }
  
  #create legend to indicate the correlation between the color scale and 
  maximum =  max(a$value)
  legend = round(seq(-maximum,maximum , length=5),2)
  color.legend(size.plot+2   ,-size.plot,size.plot+3 , size.plot, legend = legend , c(rev(graph.col), graph.col) ,gradient="y",align="rb")
  
  text(innercircle[1,][1], innercircle[1,][2], expression(xi[1]), cex=1,col="black")
  text(innercircle[2,][1], innercircle[2,][2], expression(xi[2]), cex=1,col="black") 
  text(innercircle[3,][1], innercircle[3,][2], expression(eta), cex=1,col="black") 
  
  