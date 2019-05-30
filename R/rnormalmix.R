#'@title Simulate from a normal mixture distribution
#'@param n Number of points to simulate
#'@param g Object of class normalmix (from ashr package). g has elements pi, mean and sd each of which has class K.
#'@param return.Z if TRUE, also return a vector of indicators indicating which of the K classes each sample belongs to
#'@return If return.Z=TRUE, returns a list with elements beta (samples) and Z (indicators). Otherwise returns a length n vector of samples.
#'@export
rnormalmix <- function(n, g, return.Z=FALSE){
  stopifnot(class(g)=="normalmix")
  K <- length(g$pi)
  Z <- sample(1:K, size=n, replace = TRUE, prob=g$pi)
  beta <- rep(NA, n)
  for(k in 1:K){
    nk <- sum(Z==k)
    if(nk==0) next
    beta[Z==k] <- rnorm(n=nk , mean=g$mean[k], sd = g$sd[k])
  }
  if(return.Z) return(list("beta"=beta, "Z"=Z))
  return(beta)
}
