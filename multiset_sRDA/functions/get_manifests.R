
get_manifests <- function(Data, blocks)
{
  # building data matrix 'MV'
  MV = Data[,unlist(blocks)]
  
  # add row and column names
  mvs_names = colnames(Data)[unlist(blocks)]
  dimnames(MV) = list(rownames(Data), mvs_names)

  # output
  MV
}
