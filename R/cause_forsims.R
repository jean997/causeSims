
cause_forsims <- function(dat, params, sigma_g, qalpha=1, qbeta=10,
                         ascertained_col="ld_prune", nold = FALSE,
                         thresh = 1e-3){

    stopifnot(ascertained_col %in% names(dat))

    if(no_ld){
      dat <- select(-beta_hat_1, -beta_hat_2, -p_value, -ld_prune) %>%
             rename(beta_hat_1 = beta_hat_1_nold,
                    beta_hat_2 = beta_hat_2_nold,
                    p_value = p_value_nold,
                    ld_prune = ld_prune_nold)

    }
    vars <- filter(dat, ld_prune == TRUE & p_value < thresh) %>% with(., snp)
    X <- new_cause_data(dat)

    #get sigma_g from data
    if(missing(sigma_g)) sigma_g <- eta_gamma_prior(X, vars)
    if(is.na(sigma_g)) sigma_g <- eta_gamma_prior(X, vars)

    res <- cause(X=X, variants = vars, param_ests = params, sigma_g = sigma_g, qalpha = qalpha, qbeta = qbeta, force=TRUE)
    return(res)
}


get_sigma <- function(effect, q){
    if(is.na(q) | q == -1) return(NA)
    fct <- function(s, q){
        abs(q - pnorm(effect, 0, s))}

    o <- optimize(fct, maximum = FALSE, q = q, lower = 0, upper = 40)
    return(o$minimum)
}
