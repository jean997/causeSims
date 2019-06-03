

#'@title Simulate summary statistics, nold
#'@param p Total number of snps
#'@param n1 Sample size in study 1
#'@param n2 sample size in study 2
#'@param h1 Expected heritability of trait 1
#'@param h2 Expected heritability of trait 2
#'@param neffect1 Average number of trait 1 effect snps
#'@param neffect2 Average number of trait 2 effect snps. Note that this is in addition to SNPs that effect Y through M.
#'@param gamma gamma
#'@param eta eta
#'@param q q
#'@param tau tau = h1*(gamma^2)/h2 Proportion of Y heritability explained by M through causal mechanism
#'@param omega omega = h1*(eta)^2/h2 omega*q^2 is proportion of Y heritability explained by M through confounder mechanism
#'@export
sum_stats_nold <- function(p, n1, n2, h1, h2,
                    neffect1, neffect2,
                    gamma, eta, q, tau, omega){


    if(!missing(gamma) & !missing(tau)) stop("Please provide only one of gamma or tau")
    if(missing(gamma) & missing(tau)) stop("Please provide only one of gamma or tau")


    if(missing(gamma)){
        gamma <- sqrt(tau*sum(h2)/sum(h1))
    }
    if(missing(eta)){
        eta <- sqrt(abs(omega)*sum(h2)/sum(h1))
        if(omega < 0) eta <- -1*eta
    }

    q = as.numeric(q)
    gamma = as.numeric(gamma)
    eta=as.numeric(eta)

    ## We will assume all SNPs have been normalized to have variance 1 and
    ## draw normalized effects (beta*sqrt(2*f*(1-f))) from a spike and slab
    sigma_1 <-sqrt( h1/neffect1)
    p1 <- neffect1/p
    g1 <- normalmix(pi=c(1-p1, p1),
                    mean=rep(0, 2),
                    sd=c(0, sigma_1))
    #generate trait 1 effects
    # sigma_1 = sqrt(h^2/(2*E[p*(1-p)]))
    #e_snp_var <- with(snps, mean(2*AF*(1-AF)))
    #sigma_1 <-sqrt( h1/(e_snp_var*neffect1))



    #Trait 2 effects
    #How much of trait 2 heritability is already explained by trait 1
    h2_from_h1 <- gamma^2*h1 + q*(eta^2) * h1
    stopifnot(sum(h2_from_h1) < h2)
    h2_remaining <- h2 - h2_from_h1
    p2 <- neffect2/p
    sigma_2 <-sqrt(h2_remaining/neffect2)
    #sigma_2 <-sqrt( h2_remaineder/(e_snp_var*neffect2))
    g2 <- normalmix(pi=c(1-p2, p2),
                    mean=rep(0, 2),
                    sd=c(0, sigma_2))


    b1 <- rnormalmix(p, g1)
    Z <- rbinom(n=p, size=1, prob=q)
    b2 <- rnormalmix(p, g2) + gamma*b1 + Z*eta*b1

    seb1 <- sqrt(1/n1)
    seb2 <- sqrt(1/n2)
    beta_hat_1 <- rnorm(n=p, mean = b1, sd = seb1)
    beta_hat_2 <- rnorm(n=p, mean = b2, sd = seb2)

    dat <- data.frame(beta_hat_1, seb1, beta_hat_2, seb2, snp = seq(p)) %>%
            mutate(p_value_M = 2*pnorm(-abs(beta_hat_1/seb1)),
                   p_value_Y = 2*pnorm(-abs(beta_hat_2/seb2)))

    return(dat)
}

