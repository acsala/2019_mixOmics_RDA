
get_scores <- function(X, W) 
{
  lvs = ncol(W)
  # correlations between MVs and LVs
  LV = X %*% W
  cor.XY = cor(X, LV)
  # sign ambiguity
  ODM = W
  ODM[W != 0] = 1
  w_sign = sign(colSums(sign((cor.XY * ODM))))
  if (any(w_sign <= 0)) {
    w_sign[w_sign == 0] = -1
    # scores
    LV = LV %*% diag(w_sign, lvs, lvs)    
  }
  colnames(LV) = colnames(W)
  LV
}
