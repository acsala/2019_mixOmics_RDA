

get_treated_data <- function(MV, specs)
{
  metric = get_metric(specs$scaling)
  
  if (metric) {
    # standardize if all numeric
    if (specs$scaled) {
      sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
      X = scale(MV, scale=sd.X)
    } else {
      # center if all raw
      X = scale(MV, scale=FALSE)
    }
  } else {
    # standardize
#    sd.X = sqrt((nrow(MV)-1)/nrow(MV)) * apply(MV, 2, sd)
#    X = scale(MV, scale=sd.X)
    X = scale(MV) / sqrt((nrow(MV)-1)/nrow(MV))
    # all variables quantified as "ord"/"nom" are pre-treated with get_rank
    scaling = unlist(specs$scaling)
    nominal_ordinal = which(scaling %in% c("ord", "nom"))

    for (j in nominal_ordinal) {
      X[,j] = get_rank(X[,j])
    }
  }
  
  # add names
  dimnames(X) = dimnames(MV)
  # output
  X
}


