
get_inner_summary <- 
function(path_matrix, blocks, modes, communality, redundancy=1, R2=1)
{
  blocklist = indexify(blocks)  
  exo_endo = rep("Exogenous", nrow(path_matrix))
  exo_endo[rowSums(path_matrix) != 0] = "Endogenous"
  avg_comu = rep(0, nrow(path_matrix))
  avg_redu = rep(0, nrow(path_matrix))
  AVE = rep(0, nrow(path_matrix))
  
  for (k in 1:nrow(path_matrix))
  {
    avg_comu[k] = mean(communality[blocklist == k])
    avg_redu[k] = mean(redundancy[blocklist == k])
    if (modes[k] == "A")
    {
      ave_num = sum(communality[blocklist == k])
      ave_denom = ave_num + sum(1 - (communality[blocklist == k]))
      AVE[k] = ave_num / ave_denom
    }
  }
  
  # output
  data.frame(Type = exo_endo,
             AVE = AVE,
             row.names = rownames(path_matrix))
}
