
#'@export
mrmix <- function(dat, n1, n2, p_val_thresh=5e-8, no_ld = FALSE){
  if(no_ld) dat <- process_dat_nold(dat)
  dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
  if(nrow(dat) < 3){
    R <- list("z"=NA, "p" = NA)
    return(R)
  }
  dat_std <- with(dat, standardize(beta_hat_1, beta_hat_2,
                         seb1, seb2, xtype = "continuous",
                         ytype = "continuous", n1, n2, MAF = NULL))
  res <- with(dat_std, MRMix(betahat_x = betahat_x_std, betahat_y = betahat_y_std,
                         sx=sx_std, sy = sy_std))
  R <- list(z = res$zstat_theta, p = res$pvalue_theta, est = res$theta)
  return(R)
}
