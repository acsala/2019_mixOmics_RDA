
get_weights_p <- function(X, path_matrix, blocks, specs, penalization, 
                          nonzero, lambda, warning_non_convergence = TRUE)
{
  
  # X <- Data
  # specs <- s_satpls$model$specs 
  # specs$modes <- c("B","B")
  # penalization <- "enet"
  # nonzero <- 40
  # lambda <- 1
  # warning_non_convergence = TRUE

  lambda <- rep(lambda,length(specs$modes))
  nonzero <- rep(nonzero,length(specs$modes))
  
  lvs = nrow(path_matrix)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)) / nrow(X))   # std.dev factor correction
  blockinds = indexify(blocks)

  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  a_b_weigths = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w_old = rowSums(a_b_weigths)    
  iter = 1
  
  repeat 
  {
    # estimate the latent variables
    Z = X %*% a_b_weigths
    Z = scale(Z) * sdv
    # calculate the cor/reg weights between latent variables 
    Theta <- get_path_scheme(path_matrix, Z)
    # internal estimation of LVs 'Z'
    W = Z %*% Theta  
    # computing outer weights 'w'
    for (j in 1:lvs)
    {
      if (specs$modes[j] == "predicted")
        a_b_weigths[blockinds==j,j] <- (1/nrow(X)) * W[,j] %*% X[,blockinds==j]
        a_b_weigths[blockinds==j,j] <- scale(a_b_weigths[blockinds==j,j],
                                           center = F,scale = T)
      
      if (specs$modes[j] == "predictive")
        a_b_weigths[blockinds==j,j] <- 
          switch(penalization,
                 "none" = t(solve(t(X[,blockinds==j]) %*% X[,blockinds==j])%*% 
                              (t(X[,blockinds==j]) %*% W[,j])),
                 "ust" = get_ust(X[,blockinds==j],W[,j],nonzero[j]),
                 "enet" = calculateVectorEnet(X[,blockinds==j], W[,j], 
                                              lambda[j], nonzero[j])
                 )
        a_b_weigths[blockinds==j,j] <- scale(a_b_weigths[blockinds==j,j],
                                             center = F,scale = T)
    }
    
    w_new = rowSums(a_b_weigths)                
    w_dif = sum((abs(w_old) - abs(w_new))^2) 
    if (w_dif < specs$tol || iter == specs$maxiter) break
    w_old = w_new
    iter = iter + 1
  } # end repeat   
  
  # preparing results
  if (iter == specs$maxiter && warning_non_convergence) {
    print(paste("Iterative process did not converge with 'maxiter'=", 
                specs$maxiter, ", 'tol'=", specs$tol, ", 'nonzero'=",
                nonzero, ", and 'lambda'=", lambda,
                ". Results might be suboptimal.",
                sep=""))
  }
  
  # preparing results
  a_b_weigths = a_b_weigths %*% diag(1/(apply(X %*% a_b_weigths, 2, sd)*sdv), lvs, lvs)
  w_new = rowSums(a_b_weigths)
  names(w_new) = colnames(X)
  dimnames(a_b_weigths) = list(colnames(X), rownames(path_matrix))
  results = list(w = w_new, W = a_b_weigths, ODM = ODM, iter = iter)
  
  results
}
