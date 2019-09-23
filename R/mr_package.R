
#'@export
ivw_re_nonome_MR <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
  if(no_ld) dat <- process_dat_nold(dat)
  dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
  if(nrow(dat) < 2){
    R <- list("z"=NA, "p" = NA)
    return(R)
  }
  dat <- with(dat, mr_input(bx = beta_hat_1, bxse = seb1,
                                    by = beta_hat_2, byse = seb2,
                                    snps = snp))
  f <- mr_ivw(dat, model="random", weights="delta")
  R <- list(z = f@Estimate/f@StdError, p = f@Pvalue)
  return(R)
}

#'@export
egger_re_MR <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
  if(no_ld) dat <- process_dat_nold(dat)
  dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
  if(nrow(dat) < 2){
    R <- list("z"=NA, "p" = NA)
    return(R)
  }
  dat <- with(dat, mr_input(bx = beta_hat_1, bxse = seb1,
                            by = beta_hat_2, byse = seb2,
                            snps = snp))
  f <- mr_egger(dat)
  R <- list(z = f@Estimate/f@StdError.Est, p = f@Pvalue.Est)
  return(R)
}

#'@export
wtd_median_MR <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
  if(no_ld) dat <- process_dat_nold(dat)
  dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
  if(nrow(dat) < 2){
    R <- list("z"=NA, "p" = NA)
    return(R)
  }
  dat <- with(dat, mr_input(bx = beta_hat_1, bxse = seb1,
                            by = beta_hat_2, byse = seb2,
                            snps = snp))
  f <- mr_median(dat)
  R <- list(z = f@Estimate/f@StdError, p = f@Pvalue)
  return(R)
}

#'@export
mbe_MR <- function(dat, p_val_thresh=5e-8, no_ld = FALSE,
                   weighting="weighted", no_NOME=TRUE, phi=1){
  if(no_ld) dat <- process_dat_nold(dat)
  dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
  if(nrow(dat) < 2){
    R <- list("z"=NA, "p" = NA)
    return(R)
  }
  dat <- with(dat, mr_input(bx = beta_hat_1, bxse = seb1,
                            by = beta_hat_2, byse = seb2,
                            snps = snp))
  if(no_NOME) stderror <- "delta"
    else stderror <- "simple"
  f <- mr_mbe(dat, weighting=weighting, stderror=stderror, phi=phi)
  R <- list(z = f@Estimate/f@StdError, p = f@Pvalue)
  return(R)
}


