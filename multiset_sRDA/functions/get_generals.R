
get_generals <- function(MV, path_matrix)
{
  list(obs = nrow(MV),
       obs_names = rownames(MV),
       mvs = ncol(MV),
       mvs_names = colnames(MV),
       lvs = nrow(path_matrix),
       lvs_names = rownames(path_matrix))
}
