
#'@title Run CAUSE on simulated data
#'@param dat A simulated data frame created with sum_stats
#'@param params A cause params object. You can make this with cause_params_sims
#'@param sigma_g Prior variance for gamma and eta. If left missing will be estimated.
#'@param qalpha,qbeta parameters defining prior distribution on q. q~Beta(qalpha,qbeta)
#'@param no_ld Run with the nold data (T/F)
#'@param thresh p-value threshold
#'@export
cause_sims <- function(dat, params, sigma_g, qalpha=1, qbeta=10,
                         no_ld = FALSE,
                         thresh = 1e-3){

  if(no_ld) dat <- process_dat_nold(dat)
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
