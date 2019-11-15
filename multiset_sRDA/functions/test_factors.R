
test_factors <- function(DF)
{
  contains_factors = FALSE
  # to be used when DF is a data.frame containing factors
  if (is.data.frame(DF)) {
    mvs_class = sapply(DF, class)
    mvs_as_factors <- mvs_class == "factor"
    # tell me if there are MVs as factors
    if (sum(mvs_as_factors) > 0)
      contains_factors = TRUE
  }
  contains_factors
}



get_numerics <- function(MV)
{  
  mvs_class = sapply(MV, class)
  mvs_as_factors <- mvs_class == "factor"
  categorical = which(mvs_as_factors)
  categories = vector(mode="list", ncol(MV))
  
  # only keep levels of categorical mvs in 'factor' format
  for (j in seq_along(categorical)) {
    # keep original levels
    categories[[categorical[j]]] = levels(MV[,categorical[j]])
    # convert to numeric
    MV[,categorical[j]] = as.numeric(MV[,categorical[j]])
  }    
  
  # output
  list(MV=MV, categories=categories)
}
