
#'@title Run GSMR on simulated data
#'@param dat Data frame of simulated data
#'@param evd_list List of eigen value decompositions for LD regions
#'@param p_value_thresh p-value threshold
#'@param no_ld Use no LD data
#'@export
gsmr_sims <- function(dat, evd_list, p_val_thresh  = 5e-8, no_ld = FALSE){

    ldrho <- get_ld_gsmr(dat, evd_list, p_val_thresh, no_ld)

    if(no_ld) dat <- process_dat_nold(dat)

    ix <- with(dat, which(p_value < p_val_thresh & ld_prune == TRUE))
    if(length(ix) < 3) return(NULL)

    ldrho <- data.frame(ldrho)
    names(ldrho) <- seq_along(ix)
    rownames(ldrho) <- seq_along(ix)

    dat <- dat[ix,] %>%
           mutate(p_value2 = 2*pnorm(-abs(beta_hat_2/seb2)))

    res <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value,
                           beta_hat_2, seb2, bzy_pval = p_value2,
                           ldrho=ldrho, snpid=seq_along(ix),
                           n_ref = 1, nsnps_thresh=1, gwas_thresh=p_val_thresh)))
    if(class(res)=="try-error") return(NULL)
    return(res)
}


