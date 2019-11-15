get_non_parallel_cv <- function(X, 
                                Y, 
                                lambdas,
                                non_zeros, 
                                label,
                                penalization,
                                max_iterations,
                                tolerance){
  #non-parallel x-fold Cross validation for each alpha, for example alpha = 0, 0.1, ... , 0.9, 1.0
  #
  #input:     X             n*p matrix          - independent variables 
  #           Y             n*q matrix          - dependent variables 
  #           lambda        integer             - factor for Ridge penalty 
  #           labels         vector with n_subset unique elements with length of n - for crossvalidation
  # 
  #output:    abs_cor       vector of integers  - returns the summed abs correlation
  #           stime         int/time            - return algorithms running time
  
  nr_subsets      <-    length(unique(label))
  abs_cors        <-    c()
  iterations_m    <-    c()
  
  length(lambdas)
  length(non_zeros)
  
  #measure time
  stime <- system.time({
    for (l in 1:length(lambdas)){
      
      sub_abs_cor       <- c()
      sub_results       <- c()
      sub_iterations_m  <- c()
      
      for (nz in 1:length(non_zeros)){
        
        for (i in 1:nr_subsets){
          
          X.train   <- X[label!=i,]
          X.test    <- X[label==i,]
          
          dim(X.train);dim(X.test)
          
          Y.train   <- Y[label!=i,]
          Y.test    <- Y[label==i,]
          
          dim(Y.train);dim(Y.test)
          

          data_sets <- cbind(X.train,Y.train)
          
          # Only 2 DATASETS####
          EXPL_X = c(0,0)
          RESP_Y = c(1,0)
          path_matrix = rbind(EXPL_X, RESP_Y)
          
          # blocks of outer model
          blocks = list(1:dim(X)[2], dim(X)[2]+1:dim(Y)[2])
          
          modes = c("predictive","predicted")


          
          sub_results[[i]] <- multi_sRDA(data_sets, 
                                     path_matrix,
                                     blocks, 
                                     modes,
                                     penalization = penalization, 
                                     nonzero = non_zeros[nz],
                                     lambda = lambdas[l], 
                                     maxiter = max_iterations,
                                     tol = tolerance,
                                     warning_non_convergence = FALSE)
          
          ALPHA <- t(X.train)%*% sub_results[[i]]$scores[,1]%*%
            solve(t( sub_results[[i]]$scores[,1])%*%
                    sub_results[[i]]$scores[,1])
          
          XI.test = scale(X.test) %*% ALPHA
          
          #devide with dim(Y.train)[2]
          sub_abs_cor[[i]] <- sum((abs(cor(XI.test,Y.test))))/dim(Y.train)[2]
          
          #sub_iterations_m[[i]] <- sub_results[[i]]$Nr_iterations
          
          
        }#end of subset for loop
        
        abs_cors        <- cbind(abs_cors, sub_abs_cor)
        #iterations_m    <- cbind(iterations_m, sub_iterations_m)
        
      }#end of non_zeros loop
      
    }#end of lambda for loop 
  })[3]#end of measure time
  
  #Figure out lambdas and non-zeros columns in results
  labels_non_zeros  <- rep(non_zeros, dim(abs_cors)[2]/length(non_zeros))
  labels_non_zeros
  
  labels_lambdas    <- rep(lambdas, each=dim(abs_cors)[2]/length(lambdas))
  labels_lambdas
  length(labels_lambdas)
  
  all_abs_cors  <- rbind(labels_lambdas, labels_non_zeros, abs_cors)
  
  
  mean_abs_cors <- c()
  
  for (i in 1:length(labels_lambdas)){
    
    sub_result    <-  (c(labels_lambdas[i], labels_non_zeros[i], mean(abs_cors[,i])))
    mean_abs_cors <- rbind(mean_abs_cors, sub_result)
    
  }
  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Ridge Penalty",
                                 "Number of Nonzeros", "mean_abs_cors")
  mean_abs_cors
  
  
  #plot(mean_abs_cors[,1],mean_abs_cors[,3], pch=19, col=mean_abs_cors[,2])
  #text(mean_abs_cors[,1],mean_abs_cors[,3], labels=mean_abs_cors[,2], cex= 1.1, pos=2, pch=19, col=mean_abs_cors[,2])
  
  
  # #*********************************
  # 
  # plot2 <-
  #   ggplot(data=data.frame(mean_abs_cors),
  #          aes(x = factor(Lambda), y = mean_abs_cors,         
  #              group = factor(Number.of.Nonzeros),
  #              shape = factor(Number.of.Nonzeros),
  #              color = factor(Number.of.Nonzeros)))+ 
  #   geom_line() + 
  #   geom_point() +
  #   scale_x_discrete("Lambda") +
  #   scale_y_continuous("Mean absolute correlation") + 
  #   facet_grid(.~Number.of.Nonzeros )
  # 
  # #*********************************
  
  
  #Return section**********************   
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime
                           # plot2 = plot2,
                           #iterations_m = iterations_m
                           
  )
  
  result
}