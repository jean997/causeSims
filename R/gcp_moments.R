
#'@export
gcp_moments <- function(dat, h1, h2){
  M <- nrow(dat)
  alpha_1 <- with(dat, (b1*sqrt(2*AF*(1-AF))) * sqrt(M) / sqrt(h1)) # effect size normalization
  alpha_2 <- with(dat, (b2*sqrt(2*AF*(1-AF))) * sqrt(M) / sqrt(h2))

  rho <- cor(alpha_1, alpha_2) # genetic correlation
  kappa_1 <- mean(alpha_1^3 * alpha_2) - 3 * rho
  kappa_2 <- mean(alpha_1 * alpha_2^3) - 3 * rho
  x <- seq(0, 1, 0.02)  # search for gcp value using a grid
  diff <- (kappa_1 - kappa_2)/sqrt(kappa_1^2 + kappa_2^2)-(1-rho^(2*x))/ sqrt(1+rho^(4*x))  # objective function
  idx <- which.min(abs(diff))
  return(list(gcp = x[idx], obj = abs(diff[idx])))
}
