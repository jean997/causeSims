
#'@title Estimate CAUSE parameters for simulated data
#'@param dat  A simulated data frame created with sum_stats
#'@param null_wt Null weight in dirichlet prior on mixing parameters
#'@param no_ld Run with the nold data (T/F)
#'@return
#'@export
cause_params_sims <- function(dat, null_wt = 10, no_ld=FALSE){

    if(no_ld){
      dat <- dat %>% select(-beta_hat_1, -beta_hat_2) %>%
            rename(beta_hat_1 = beta_hat_1_nold, beta_hat_2 = beta_hat_2_nold)
    }

    X <- dat  %>%
         select(snp, beta_hat_1, seb1, beta_hat_2, seb2) %>%
         new_cause_data(.)

    params <- est_cause_params(X, X$snp, null_wt = null_wt)
    return(params)
}

