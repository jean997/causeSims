# This function LD prunes using the correlation matrix rather than the two column table used by the cause function.
# It is useful for simulations

ld_prune_cormat <- function(R, snp_names, p_vals,  p_val_thresh, r2_thresh){
  stopifnot(nrow(R) == length(snp_names))
  stopifnot(length(snp_names) == length(p_vals))

  ix <- which(p_vals < p_val_thresh)

  if(length(ix) == 0) return(c())

  snp_names <- snp_names[ix]
  R <- R[ix, ix, drop=FALSE]
  p_vals <- p_vals[ix]
  o_c <- order(p_vals, decreasing=FALSE)
  keep <- c()
  while(length(o_c) > 0){
    keep <- c(keep, snp_names[o_c[1]])
    myld <- R[o_c, o_c, drop=FALSE]
    remove_ix <- which(myld[,1]^2 > r2_thresh)
    o_c <- o_c[-remove_ix]
  }
  return(keep)
}
