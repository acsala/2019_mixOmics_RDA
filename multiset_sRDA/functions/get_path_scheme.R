
get_path_scheme <- function(path_matrix, LV)
{
  # matrix for inner weights
  E = path_matrix
  
  for (k in 1:ncol(path_matrix)) 
  {
    # followers
    follow <- path_matrix[k,] == 1
    if (sum(follow) > 0)
      E[follow,k] <- lm(LV[,k] ~ LV[,follow] - 1)$coef
    # predecesors
    predec <- path_matrix[,k] == 1
    if (sum(predec) > 0)
      E[predec,k] <- cor(LV[,k], LV[,predec])
  } 
  
  # output
  E
}
