

#'@title Simulate summary statistics
#'@param n_copies How many times to replicate chromosome 19
#'@param n1 Sample size in study 1
#'@param n2 sample size in study 2
#'@param h1 Expected heritability of trait 1
#'@param h2 Expected heritability of trait 2
#'@param neffect1 Average number of trait 1 effect snps
#'@param neffect2 Average number of trait 2 effect snps
#'@param gamma gamma
#'@param eta eta
#'@param q q
#'@param tau tau = h1*(gamma^2)/h2 Proportion of Y heritability explained by M through causal mechanism
#'@param omega omega = h1*(eta)^2/h2 omega*q^2 is proportion of Y heritability explained by M through confounder mechanism
#'@export
ld_data <- function(snps, evd_list,
                    n_copies = 30,
                    n1, n2, h1, h2,
                    neffect1, neffect2,
                    gamma, eta, q, tau, omega,
                    cores, ld_prune_pval_thresh = 1e-3, r2_thresh = 0.1){

    plan(multiprocess, workers=cores)

    stopifnot(inherits(snps, "data.frame"))

    if(!all(c("AF", "SNP", "region_id") %in% names(snps))){
      stop("snps data frame should have at least collumns AF, SNP, and region_id.\n")
    }
    if(!all(snps$region_id %in% seq(length(evd_list)))){
      stop("region_id should be in 1:length(evd_list).\n")
    }
    #check dimensions
    region_ids <- unique(snps$region_id)
    for(r in region_ids){
      p <- sum(snps$region_id==r)
      stopifnot(p == nrow(evd_list[[r]]$vectors))
    }


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
    params <- c(q=q, gamma=gamma, eta=eta)
    n_snps <- nrow(snps)*n_copies

    ## We will assume all SNPs have been normalized to have variance 1 and
    ## draw normalized effects (beta*sqrt(2*f*(1-f))) from a spike and slab
    sigma_1 <-sqrt( h1/neffect1)
    p1 <- neffect1/n_snps
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
    p2 <- neffect2/n_snps
    sigma_2 <-sqrt(h2_remaining/neffect2)
    #sigma_2 <-sqrt( h2_remaineder/(e_snp_var*neffect2))
    g2 <- normalmix(pi=c(1-p2, p2),
                    mean=rep(0, 2),
                    sd=c(0, sigma_2))

    dat <- future_map_dfr(region_ids, function(r){
    #dflist <- lapply(region_ids, function(r){
    #dat <- future_map(region_ids, function(r){
                    #evd <- readRDS(paste0("data/EVD_hapmap_only/EVD_", r, ".RDS"))
                    R <- with(evd_list[[r]], crossprod(sqrt(values)*t(vectors)))

                    p <- sum(snps$region_id==r)
                    b1 <- replicate(n=n_copies, rnormalmix(p, g1))
                    Z <- replicate(n=n_copies, rbinom(n=p, size=1, prob=q))
                    b2 <- replicate(n=n_copies, rnormalmix(p, g2)) + gamma*b1 + Z*eta*b1

                    #LD transformed effects
                    ld_z1 <- apply(b1, 2, function(x){
                               sqrt(n1)* (R %*% x)})
                    ld_z2 <- apply(b2, 2, function(x){
                      sqrt(n2)* (R %*% x)})

                    z_hat_1 <- apply(ld_z1, 2, function(x){
                      mvrnorm_eig1(n=1, mu=x, eS=evd_list[[r]])
                    })

                    z_hat_2 <- apply(ld_z2, 2, function(x){
                      mvrnorm_eig1(n=1, mu=x, eS=evd_list[[r]])
                    })

                    ### Convert everything from normalized to non-normalized effects
                    f <- snps$AF[snps$region_id==r]
                    seb1 <- sqrt(1/(2*n1*f*(1-f)))
                    seb2 <- sqrt(1/(2*n2*f*(1-f)))
                    b1 <- b1*seb1*sqrt(n1)
                    b2 <- b2*seb2*sqrt(n2)
                    ld_b1 <- ld_z1*seb1
                    ld_b2 <- ld_z2*seb2
                    beta_hat_1 <- z_hat_1*seb1
                    beta_hat_2 <- z_hat_2*seb2


                    df <- data.frame(beta_hat_1 = as.vector(beta_hat_1),
                                     beta_hat_2 = as.vector(beta_hat_2),
                                     seb1 = rep(seb1, n_copies), seb2 = rep(seb2, n_copies),
                                     ld_b1 = as.vector(ld_b1), ld_b2 = as.vector(ld_b2),
                                     b1 = as.vector(b1), b2 = as.vector(b2))

                    ### Make also some no ld data with the same effects
                    df$beta_hat_1_nold <- with(df, rnorm(n=nrow(df), mean = b1, sd = seb1))
                    df$beta_hat_2_nold <- with(df, rnorm(n=nrow(df), mean = b1, sd = seb2))


                    # add snp data
                    snpdat <- filter(snps, region_id == r)
                    #df <- cbind(snpdat, df)
                    df$snp <- paste0(rep(snpdat$SNP, n_copies), ".", rep(seq(n_copies), each=p))
                    df$region_id <- rep(snpdat$region_id, n_copies)
                    df$AF <- rep(snpdat$AF, n_copies)
                    df$ldscore_hm3 <- rep(snpdat$ldscore_hm3, n_copies)
                    df$rep <- rep(seq(n_copies), each=p)

                    #LD pruning subset
                    df <- df %>% mutate(p_value = 2*pnorm(-abs(beta_hat_1/seb1)),
                                        p_value_nold = 2*pnorm(-abs(beta_hat_1_nold/seb1)))
                    keep <- sapply(seq(n_copies), function(i){
                      strt <- ((i-1)*p) + 1
                      stp <- i*p
                      ld_prune_cormat(R, df$SNP[strt:stp], df$p_value[strt:stp],  ld_prune_pval_thresh, r2_thresh)
                    }) %>% unlist()

                    keep_nold <- sapply(seq(n_copies), function(i){
                      strt <- ((i-1)*p) + 1
                      stp <- i*p
                      ld_prune_cormat(R, df$SNP[strt:stp], df$p_value_nold[strt:stp],  ld_prune_pval_thresh, r2_thresh)
                    }) %>% unlist()

                    df <- df %>% mutate(ld_prune = case_when(!SNP %in% keep ~ FALSE,
                                                                 TRUE ~ TRUE),
                                        ld_prune_nold = case_when(!SNP %in% keep_nold ~ FALSE,
                                                            TRUE ~ TRUE)
                                        )
                    #cat(apply(df, 2, class))
                    return(df)
                })

    return(dat)
    #return(dflist)
}

