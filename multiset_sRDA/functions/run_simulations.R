# ################################################################################
# #
# #                                                                                     
# #   Filename:	  run_simulations.R    												  
# #                                                                                     
# #   Project :   BiomJ article "Multiset sparse redundancy analysis for high 
# #               dimensional omics data"
# #   Date    :   24-01-2018
# #
# #
# ################################################################################
# 
# # fintetune parameters ####
# ## data generation
# N = 50    # number of individuals
# k = 4      # number of datasets
# m = 2      # number of latent variables (LV's) per dataset
# # number of irrelevant variables per dataset
# p0=c(1000,500,200,10)
# # number of variables associated with the LV's
# p1=c(10,8,8,2)
# 
# 
# #repeat simulation ####
# 
# nr_of_simulations <- 1
# 
# sens.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
# sens2.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
# spec.m <- matrix(c(0,0,0),nrow = max(nr_of_simulations),ncol = 3)
# ridge_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations),ncol = 4)
# lasso_param <- matrix(c(0,0,0,0),nrow = max(nr_of_simulations),ncol = 4)
# 
# for (i in 1:nr_of_simulations){
#   
#   #pdf(paste(i, "resultsFULL.pdf", sep=""))
#   print("nr of simulation")
#   print(i)    
#   
#   # for replication
#   nr_seed = runif(1, 10, 10^8)
#   set.seed(nr_seed)
#   
#   multiblockdata <- multiset_data_generator(N, k,m,
#                                             p0,
#                                             p1)
#   X1 <- multiblockdata$X1
#   X2 <- multiblockdata$X2
#   X3 <- multiblockdata$X3
#   X4 <- multiblockdata$X4
#   
#   Data <- cbind(X1,X2,X3,X4)
#   EXPL_X = c(0,0,0,0)
#   RESP_Y = c(1,0,0,0)
#   EXPL_Z = c(0,1,0,0)
#   RESP_V = c(0,0,1,0)
#   path_matrix = rbind(EXPL_X, RESP_Y,EXPL_Z,RESP_V)
#   
#   # blocks of outer model
#   blocks = list(1:dim(X1)[2], (dim(X1)[2]+1):(dim(X1)[2]+dim(X2)[2]),
#                 (dim(X1)[2]+dim(X2)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]),
#                 (dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+1):(dim(X1)[2]+dim(X2)[2]+dim(X3)[2]+dim(X4)[2]))
#   
#   
#   modes = c("B","B","B","A")
#   
#   time_data <- system.time(
#     s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
#                        scaled=T, penalization = "ust", nonzero = c(5,10), 
#                        lambda = 1, maxiter = 100, cross_validate = T)
#   )
#   
#   #s_satpls$outer_model
#   print("s_satpls$model$iter")
#   print(s_satpls$model$iter)
#   
#   print("s_satpls$nonzero")
#   print(s_satpls$nonzero)
#   print("s_satpls$lambda")
#   print(s_satpls$lambda)
#   
#   ridge_param[i] <- s_satpls$lambda
#   lasso_param[i] <- s_satpls$nonzero
#   
#   # checking non zeros
#   # nzero_X_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])>0)
#   # nzero_Z_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])>0)
#   # nzero_Y_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])>0)
#   
#   # mode A response will not have nonzeros
#   #nzero_V_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_V",3])>0)
#   
#   # nzero_X_positiong
#   # nzero_Y_positiong
#   # nzero_Z_positiong
#   # c(dim(X1),dim(X2),dim(X3),dim(X4))
#   # 
#   # print("nr of associated variables per dataset")
#   # p1
#   
#   # s_satpls$outer_model[which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])>0),]
#   # s_satpls$outer_model[510+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])>0),]
#   # s_satpls$outer_model[818+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])>0),]
#   # s_satpls$outer_model[1026+which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_V",3])>0),]
#   
#   #calculate sens spec
#   
#   nzero_X_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])>0)
#   nzero_Y_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])>0)
#   nzero_Z_positiong <- which(abs(s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])>0)
#   
#   zero_X_positiong <- which((s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",3])==0)
#   zero_Y_positiong <- which((s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Y",3])==0)
#   zero_Z_positiong <- which((s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Z",3])==0)
#   
#   #X sensitivity
#   sensX = sum(nzero_X_positiong %in% p0[1]:(p1+p0)[1])/min(p1[1],length(nzero_X_positiong))
#   sensY = sum(nzero_Y_positiong %in% p0[2]:(p1+p0)[2])/min(p1[2],length(nzero_Y_positiong))
#   sensZ = sum(nzero_Z_positiong %in% p0[3]:(p1+p0)[3])/min(p1[3],length(nzero_Z_positiong))
#   
#   sens2X = sum(nzero_X_positiong %in% p0[1]:(p1+p0)[1])/p1[1]
#   sens2Y = sum(nzero_Y_positiong %in% p0[2]:(p1+p0)[2])/p1[2]
#   sens2Z = sum(nzero_Z_positiong %in% p0[3]:(p1+p0)[3])/p1[3]
#   
#   print("print(c(sensX,sensY,sensZ))")
#   print(c(sensX,sensY,sensZ))
#   
#   print("print(c(sens2X,sens2Y,sens2Z))")
#   print(c(sens2X,sens2Y,sens2Z))
#   
#   #specificiy
#   specX = sum(zero_X_positiong %in% 1:p0[1])/p0[1]
#   specY = sum(zero_Y_positiong %in% 1:p0[2])/p0[2]
#   specZ = sum(zero_Z_positiong %in% 1:p0[3])/p0[3]
#   
#   print("print(c(specX,specY,specZ))")
#   print(c(specX,specY,specZ))
#   
#   #sens = sum(as.numeric(a!=0)[1:highly_correlated])/highly_correlated
#   #spec = sum(as.numeric(a[(highly_correlated+1):length(a)]==0)) / (length(a)-highly_correlated)
#   
#   sens.m[i,]  =c(sensX,sensY,sensZ)
#   sens2.m[i,] =c(sens2X,sens2Y,sens2Z)
#   spec.m[i,]  =c(specX,specY,specZ)
#   
#   
#   iter    = s_satpls$model$iter
#   nonzero = s_satpls$nonzero
#   lambda  = s_satpls$lambda
#   
#   #save objects in RData file
#   save(nr_seed,
#        nzero_X_positiong,
#        nzero_Y_positiong,
#        nzero_Z_positiong,
#        zero_X_positiong,
#        zero_Y_positiong,
#        zero_Z_positiong,
#        iter,
#        nonzero,
#        lambda,
#        N,k,m,
#        p0,p1,
#        time_data,
#        file = paste(i, "_nrResults.RData", sep=""))
#   
#   #close pdf
#   #dev.off()
#   
#   print("matrix Sens:")
#   print(sens.m)
#   print(apply(sens.m,2,mean))
#   print("matrix Sens2:")
#   print(sens2.m)
#   print(apply(sens2.m,2,mean))
#   print("matrix spec:")
#   print(spec.m)
#   print(apply(spec.m,2,mean))
#   
# }
# 
# print("Mean Sens:")
# print(apply(sens.m,2,mean))
# print("Mean Sens2:")
# print(apply(sens2.m,2,mean))
# print("Mean spec:")
# print(apply(spec.m,2,mean))
# 
# 
# save(sens.m,
#      sens2.m,
#      spec.m,
#      ridge_param,
#      lasso_param,
#      file = paste(i, "sens_spec.RData", sep=""))
# 
# print("ridge")
# ridge_param
# print("lasso")
# lasso_param
# 
# print("We're done have a good1")
