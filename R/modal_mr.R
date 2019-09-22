
#'@title Modal MR
#'@description Implementation of the Modal MR method from Hartwig et al 2017.
#'This code is copied from the supplement of that paper.
#'@param BetaXG vector of instrument-exposure regression coefficients
#'@param BetaYG vector of instrument-outcome regression coefficients
#'@param seBetaXG vector of instrument-exposure standard errors (SEs)
#'@param seBetaYG vector of instrument-outcome SEs
#'@param phi tunning parameter (e.g., 1=default bandwidth; 0.5=half of the default bandwidth)
#'@param n_boot number of bootstrap iterations
#'@param alpha alpha level of the confidence intervals (e.g., alpha=0.05 corresponds to 1-0.05=95\% confidence intervals)
#'@export
MBE <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, phi=c(1,0.5,0.25), n_boot=1e4, alpha=0.05) {

  #--------------------------------------#
  #Function to compute the point estimate#
  #--------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  beta <- function(BetaIV.in, seBetaIV.in) {

    #Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)

    #Standardised weights
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)

    beta <- NULL

    for(cur_phi in phi) {
      #Define the actual bandwidth
      h <- s*cur_phi
      #Compute the smoothed empirical density function
      densityIV <- density(BetaIV.in, weights=weights, bw=h)
      #Extract the point with the highest density as the point estimate
      beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
    }
    return(beta)
  }

  #------------------------------------------#
  #Function to estimate SEs through bootstrap#
  #------------------------------------------#
  #BetaIV.in: ratio estimates
  #seBetaIV.in: standard errors of ratio estimates
  #beta_MBE.in: point causal effect estimates
  boot <- function(BetaIV.in, seBetaIV.in, beta_MBE.in) {

    #Set up a matrix to store the results from each bootstrap iteration
    beta.boot <- matrix(nrow=n_boot, ncol=length(beta_MBE.in))

    for(i in 1:n_boot) {
      #Re-sample each ratio estimate using SEs derived not assuming NOME
      BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
      #Re-sample each ratio estimate using SEs derived under NOME
      BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])

      #Simple MBE, not assuming NOME
      beta.boot[i,1:length(phi)]                     <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)))
      #Weighted MBE, not assuming NOME
      beta.boot[i,(length(phi)+1):(2*length(phi))]   <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1])
      #Simple MBE, assuming NOME
      beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)))
      #Weighted MBE, assuming NOME
      beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2])
    }
    return(beta.boot)
  }

  #Ratio estimates
  BetaIV   <- BetaYG/BetaXG
  #SEs of ratio estimates
  seBetaIV <- cbind(sqrt((seBetaYG^2)/(BetaXG^2) + ((BetaYG^2)*(seBetaXG^2))/(BetaXG^4)), #SEs NOT assuming NOME
                    seBetaYG/abs(BetaXG))                                                 #SEs ASSUMING NOME

  #Point causal effect estimate using the simple MBE
  beta_SimpleMBE        <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)))
  #Point causal effect estimate using the weighted MBE (not asusming NOME)
  beta_WeightedMBE      <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1])
  #Point causal effect estimate using the weighted MBE (asusming NOME)
  beta_WeightedMBE_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2])
  #Combine all point effect estimates in a single vector
  beta_MBE <- rep(c(beta_SimpleMBE, beta_WeightedMBE,
                    beta_SimpleMBE, beta_WeightedMBE_NOME))

  #Compute SEs, confidence intervals and P-value
  beta_MBE.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_MBE.in=beta_MBE)
  se_MBE <- apply(beta_MBE.boot, 2, mad)

  CIlow_MBE <- beta_MBE-qnorm(1-alpha/2)*se_MBE
  CIupp_MBE <- beta_MBE+qnorm(1-alpha/2)*se_MBE

  P_MBE <- pt(abs(beta_MBE/se_MBE), df=length(BetaXG)-1, lower.tail=F)*2

  #Vector to indicate the method referring to each row
  Method <- rep(c('Simple', 'Weighted', 'Simple (NOME)', 'Weighted (NOME)'), each=length(phi))

  #Return a data frame containing the results
  Results <- data.frame(Method, phi, beta_MBE, se_MBE, CIlow_MBE, CIupp_MBE, P_MBE)
  colnames(Results) <- c('Method', 'phi', 'Estimate', 'SE', 'CI_low', 'CI_upp', 'P')

  return(Results)
}
