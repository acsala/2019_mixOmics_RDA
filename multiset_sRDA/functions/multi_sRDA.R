################################################################################
#
#                                                                                     
#   Filename:	  multi_sRDA.R    												  
#                                                                                     
#   Project :   BiomJ article "Multiset sparse redundancy analysis for high 
#               dimensional omics data"
#   Date    :   24-01-2018
#
################################################################################

library(tester)
library(turner)
library(diagram)
library(shape)
library(amap)
library(plotrix)

multi_sRDA <- function(Data, 
                   path_matrix, 
                   blocks, 
                   modes = NULL, 
                   scaling = NULL, 
                   scaled = TRUE, 
                   tol = 0.000001, 
                   maxiter = 100,
                   plscomp = NULL, 
                   boot.val = FALSE, 
                   br = NULL,
                   dataset = TRUE,
                   penalization = "none",
                   nonzero = 0,
                   lambda = 0,
                   cross_validate = FALSE,
                   nr_subsets = 10,
                   warning_non_convergence = TRUE){
  # =======================================================================
  # checking arguments
  # =======================================================================
  valid = check_args(Data=Data, path_matrix=path_matrix, blocks=blocks,
                     scaling=scaling, modes=modes, scheme="path",
                     scaled=scaled, tol=tol, maxiter=maxiter,
                     plscomp=plscomp, boot.val=boot.val, br=br,
                     dataset=dataset)

  Data = valid$Data
  path_matrix = valid$path_matrix
  blocks = valid$blocks
  specs = valid$specs
  boot.val = valid$boot.val
  br = valid$br
  dataset = valid$dataset
  
  # =======================================================================
  # Preparing data and blocks indexification
  # =======================================================================
  # building data matrix 'MV'
  MV = get_manifests(Data, blocks)
  #check_MV = test_manifest_scaling(MV, specs$scaling)
  # generals about obs, mvs, lvs
  gens = get_generals(MV, path_matrix)
  # blocks indexing
  names(blocks) = gens$lvs_names
  block_sizes = lengths(blocks)
  blockinds = indexify(blocks)
  
  # transform to numeric if there are factors in MV
  if (test_factors(MV)) {
    numeric_levels = get_numerics(MV)
    MV = numeric_levels$MV
    categories = numeric_levels$categories
  }
  # apply corresponding treatment (centering, reducing, ranking)
  X = get_treated_data(MV, specs)
  
  # =======================================================================
  # Cross validation for best optimal penalty parameters
  # =======================================================================
  
  if (cross_validate){
    
    cv_results <- get_cross_validated_penalty_parameters(data_sets = Data,
                                                         path_matrix = path_matrix,
                                                         blocks = blocks,
                                                         specs = specs,
                                                         penalization = penalization,
                                                         ridge_penalty = lambda,
                                                         nonzero = nonzero,
                                                         tolerance = tol,
                                                         max_iterations = maxiter,
                                                         nr_subsets = nr_subsets)
    lambda = cv_results$best_ridge
    nonzero = cv_results$best_nonzero
    
  }
  
  # =======================================================================
  # Outer weights and LV scores
  # =======================================================================
  metric = get_metric(specs$scaling)
  if (metric) {
    # object 'weights' contains outer w's, W, ODM, iter
    weights = get_weights_p(X, path_matrix, blocks, specs, penalization,
                            nonzero, lambda, 
                            warning_non_convergence = warning_non_convergence)
    #ok_weights = test_null_weights(weights, specs)
    outer_weights = weights$w
    LV = get_scores(X, weights$W)
  } else {
    # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
    weights = get_weights_nonmetric(X, path_matrix, blocks, specs)
    ok_weights = test_null_weights(weights, specs)
    outer_weights = weights$w
    LV = weights$Y
    X = weights$QQ  # quantified MVs
    colnames(X) = gens$mvs_names
  }
  
  # =======================================================================
  # Path coefficients
  # =======================================================================
  inner_results = get_paths(path_matrix, LV)
  #inner_model = inner_results[[1]]
  Path = inner_results[[2]]
  #R2 = inner_results[[3]]
  
  # =======================================================================
  # Outer model: loadings, communalities, redundancy, crossloadings
  # =======================================================================
  xloads = cor(X, LV, use = 'pairwise.complete.obs')
  loadings = rowSums(xloads * weights$ODM)
  communality = loadings^2
  #R2_aux = rowSums(weights$ODM %*% diag(R2, gens$lvs, gens$lvs))
  #redundancy = communality * R2_aux
  crossloadings = data.frame(xloads, row.names=1:gens$mvs)
  crossloadings$name = factor(gens$mvs_names, levels = unique(gens$mvs_names))
  crossloadings$block = factor(rep(gens$lvs_names, block_sizes),
                               levels = gens$lvs_names)
  crossloadings = crossloadings[,c('name','block',colnames(xloads))]
  
  # outer model data frame
  outer_model = data.frame(
    name = factor(gens$mvs_names, levels = unique(gens$mvs_names)),
    block = factor(rep(gens$lvs_names, block_sizes),
                   levels = gens$lvs_names),
    weight = outer_weights,
    loading = loadings,
    communality = communality,
    #redundancy = redundancy,
    row.names = 1:gens$mvs)
  
  # Summary Inner model
  inner_summary = get_inner_summary(path_matrix, blocks, specs$modes,
                                    communality)
  
  # =======================================================================
  # Results
  # =======================================================================
  # deliver dataset?
  if (dataset) data = MV else data = NULL
  # deliver bootstrap validation results?
  bootstrap = FALSE
  if (boot.val)
  {
    if (nrow(X) <= 10) {
      warning("Bootstrapping stopped: very few cases.")
    } else {
      bootstrap = get_boots(MV, path_matrix, blocks, specs, br)
    }
  }
  
  # list with model specifications
  model = list(IDM=path_matrix, blocks=blocks, specs=specs,
               iter=weights$iter, boot.val=boot.val, br=br, gens=gens)
  
  ## output
  res = list(weights = outer_model,
             path_coefs = Path,
             scores = LV,
             crossloadings = crossloadings,
             data = data,
             manifests = X,
             model = model,
             lambda2 = lambda,
             nonzero = nonzero)
  class(res) = "multi_sRDA"
  return(res)
}


#'@S3method print multi_sRDA
print.multi_sRDA <- function(x, ...)
{
  cat("Multiset sparse redundancy analysis", "\n")
  cat(rep("-", 45), sep="")
  cat("\n   NAME            ", "DESCRIPTION")
  cat("\n1  $weights        ", "Alpha and Beta weights")
  cat("\n3  $path_coefs     ", "Path coefficients matrix")
  cat("\n4  $scores         ", "Latent variables")
  cat("\n5  $lambda2        ", "Ridge penalty")
  cat("\n6  $nonzero        ", "Number of nonzeros")
  cat("\n7  $crossloadings  ", "cross-loadings")
  cat("\n8  $data           ", "Predictive and predicted datasets")
  cat("\n")
  cat(rep("-", 45), sep="")
  cat("\nYou can also use the function 'summary'", "\n\n")
  invisible(x)
}
