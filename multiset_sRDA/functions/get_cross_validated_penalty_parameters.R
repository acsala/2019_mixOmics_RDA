#prepare for cross validation

get_cross_validated_penalty_parameters <- function(data_sets,
                                                   path_matrix,
                                                   blocks,
                                                   specs,
                                                   penalization,
                                                   ridge_penalty,
                                                   nonzero,
                                                   tolerance,
                                                   max_iterations,
                                                   nr_subsets = nr_subsets
                                                   ){

  best_ridge_penalties  <- c()
  best_nr_of_nonzeros   <- c()
  cv_results_list       <- c()
  
  # TESTER VALUES
  # ridge_penalty <- c(0.1,1)
  # nonzero <- c(20,40)
  # tolerance <- 0.001
  # max_iterations <- 100
  # penalization <- "enet"
  # data_sets <- Data
  # nr_subsets <- 10
  # specs <- s_satpls$model$specs
  # specs$modes <- c("A","A","A","B")
  
  
  #crossvaldate datasets if they are connected
  for (nr_sets in 1:length(specs$modes)){
    
    #nr of connections current dataset has
    current_set_path <- path_matrix[,nr_sets]
    
    #check if dataset itself mode predictive
    if (specs$modes[nr_sets] == "predictive"){
      
      #if dataset is itself mode predictive and first dataset
      if (sum(current_set_path)>=1 || nr_sets == length(specs$modes)){
        
        for (set_connected_checker in 1:length(current_set_path)){
          
          if(current_set_path[set_connected_checker]==1 || 
             (nr_sets == length(specs$modes) && 
              set_connected_checker == length(current_set_path))){
            
            #save current dataset in varaiable
            matrix_X <- data_sets[,blocks[[nr_sets]]]
            
            dim(matrix_X)
            
            matrix_Y <- data_sets[,blocks[[set_connected_checker]]]
            
            dim(matrix_Y)
            
                    print("dataset itself mode predictive and connected (nr set / connected to)")
                    print(nr_sets)
                    print(set_connected_checker)

                    shuffled <-  get_split_sets(X = matrix_X, Y = matrix_Y,
                                                nr_subsets = nr_subsets)

                    X.sampled     <-   shuffled$X.sampled
                    Y.sampled     <-   shuffled$Y.sampled
                    label         <-   shuffled$labels


                    cv_results <- get_non_parallel_cv(X = X.sampled,
                                                    Y = Y.sampled,
                                                    lambdas = ridge_penalty,
                                                    non_zeros = nonzero,
                                                    label = label,
                                                    penalization = penalization,
                                                    max_iterations = max_iterations,
                                                    tolerance = tolerance)

                    cv_results$abs_cors
                    cv_results$mean_abs_cors

                    a = cv_results$mean_abs_cors[,3]

                    best_values     <- cv_results$mean_abs_cors[which.max(a),]

                    best_ridge   <- best_values[1]
                    best_nonzero   <- best_values[2]
          }
        }
      }
    }else{
      #if connection to B dataset
      best_ridge    <- 0
      best_nonzero  <- 0
      cv_results    <- c()
    } 
    
    best_ridge_penalties[[nr_sets]]   <- best_ridge
    best_nr_of_nonzeros[[nr_sets]]    <- best_nonzero
    cv_results_list[[nr_sets]]        <- cv_results
    
  }
      
  
  cv_results_list
  best_ridge_penalties
  best_nr_of_nonzeros
  
  #**********************
  result <- list(cv_results = cv_results_list,
                 best_ridge = best_ridge_penalties,
                 best_nonzero = best_nr_of_nonzeros
                 )
  
  
  result
  
}