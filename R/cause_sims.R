
#'@title Run CAUSE on simulated data
#'@param dat A simulated data frame created with sum_stats
#'@param params A cause params object. You can make this with cause_params_sims
#'@param sigma_g Prior variance for gamma and eta. If left missing will be estimated.
#'@param qalpha,qbeta parameters defining prior distribution on q. q~Beta(qalpha,qbeta)
#'@param no_ld Run with the nold data (T/F)
#'@param thresh p-value threshold
#'@export
cause_sims <- function(dat, param_ests, sigma_g, qalpha=1, qbeta=10,
                         no_ld = FALSE,
                         thresh = 1e-3){

  if(no_ld) dat <- process_dat_nold(dat)

  vars <- filter(dat, ld_prune == TRUE & p_value < thresh) %>% with(., snp)
  X <- new_cause_data(dat)

  #get sigma_g from data
  if(missing(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)
  if(is.na(sigma_g)) sigma_g <- cause:::eta_gamma_prior(X, vars)

  res <- cause::cause(X=X, variants = vars, param_ests = param_ests, sigma_g = sigma_g,
               qalpha = qalpha, qbeta = qbeta, force=TRUE)
  return(res)
}

#this is a helper function for switching from with ld to no ld data
process_dat_nold <- function(dat){
  dat <- dat %>%  select(-beta_hat_1, -beta_hat_2, -p_value, -ld_prune) %>%
    rename(beta_hat_1 = beta_hat_1_nold,
           beta_hat_2 = beta_hat_2_nold,
           p_value = p_value_nold) %>%
    mutate(ld_prune = TRUE)
  return(dat)
}

28

#'@export
#'@title Helper function for simulations testing prior on gamma and eta
#'@param effect
#'@param q
#'@description Given an effect and a quantile, q, find the value of sigma_g
#'such that pnorm(effect, 0, sigma_g) = q
get_sigma <- function(effect, q){
  if(is.na(q) | q == -1) return(NA)
  fct <- function(s, q){
  abs(q - pnorm(effect, 0, s))}
  o <- optimize(fct, maximum = FALSE, q = q, lower = 0, upper = 40)
  return(o$minimum)
}
