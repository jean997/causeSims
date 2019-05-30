
process_dat_nold <- function(dat){
    dat <- dat %>%  select(-beta_hat_1, -beta_hat_2, -p_value, -ld_prune) %>%
      rename(beta_hat_1 = beta_hat_1_nold,
             beta_hat_2 = beta_hat_2_nold,
             p_value = p_value_nold) %>%
      mutate(ld_prune = TRUE)
    return(dat)
}


#'@export
ivw <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    if(nrow(dat) < 2){
       R <- list("z"=NA, "p" = NA)
       return(R)
    }
    f1 <- lm(beta_hat_2 ~ beta_hat_1-1, data=dat, weights=1/(seb2^2))
    z_mr <- as.numeric(f1$coefficients/sqrt(vcov(f1)))
    p_mr <- pt(abs(z_mr), df=nrow(dat)-1, lower.tail=F)*2
    R <- list(z = z_mr, p = p_mr)
    return(R)
}

#'@export
egger <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    if(nrow(dat) < 3){
       R <- list("z"=NA, "p" = NA)
       return(R)
    }
    f2 <- lm(beta_hat_2 ~ beta_hat_1, data=dat, weights=1/(seb2^2))
    z_mre <- as.numeric(f2$coefficients["beta_hat_1"]/sqrt(vcov(f2)[2,2]))
    p_mre <- pt(abs(z_mre), df=nrow(dat)-2, lower.tail=F)*2
    R <- list(z = z_mre, p = p_mre)
    return(R)
}


#'@export
mrpresso <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    x <- try(mr_presso(BetaOutcome = "beta_hat_2", BetaExposure = "beta_hat_1",
              SdOutcome = "seb2", SdExposure = "seb1", OUTLIERtest = TRUE,
              DISTORTIONtest = TRUE, data = dat,
              NbDistribution = 1000,  SignifThreshold = 0.05), silent=TRUE)
    if(class(x) == "try-error"){
        return(list(z = NA, p = NA))
    }
    t_mrp <- x["Main MR results"]$`Main MR results`[2,5]
    p_mrp <- x["Main MR results"]$`Main MR results`[2,6]
    if(is.na(t_mrp)){
        t_mrp <- x["Main MR results"]$`Main MR results`[1,5]
        p_mrp <- x["Main MR results"]$`Main MR results`[1,6]
    }
    R <- list(z = t_mrp, p = p_mrp)
    return(R)

}

#Gusev et al 2016
#'@export
twas <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    if(nrow(dat) < 2){
       R <- list("z"=NA, "p" = NA)
       return(R)
    }
    zg <- with(dat, beta_hat_2/seb2)
    numerator <- with(dat, sum(beta_hat_1*zg))
    denom <- with(dat, sqrt(sum(beta_hat_1^2)))
    z_twas <- numerator/denom
    p_twas <- 2*pnorm(abs(z_twas), lower.tail=FALSE)
    R <- list(z = z_twas, p = p_twas)
    return(R)
}


    #Code coppied from appendix 2 of Bowden et al Gen Epi
weighted.median <- function(betaIV.in, weights.in) {
    betaIV.order = betaIV.in[order(betaIV.in)]
    weights.order = weights.in[order(betaIV.in)]
    weights.sum = cumsum(weights.order)-0.5*weights.order
    weights.sum = weights.sum/sum(weights.order)
    below = max(which(weights.sum<0.5))
    weighted.est = betaIV.order[below] +
    (betaIV.order[below+1]-betaIV.order[below])* (0.5-weights.sum[below])/
                                                (weights.sum[below+1]-weights.sum[below])
    return(weighted.est)
}
weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
   med = NULL
   for(i in 1:1000){
        betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
        betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
        betaIV.boot = betaYG.boot/betaXG.boot
        med[i] = weighted.median(betaIV.boot, weights.in)
   }
   return(sd(med))
}

#'@export
weighted_median <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    if(nrow(dat)< 3){
        ret <- list(z = NA,  p = NA)
        return(ret)
    }
    betaYG <- dat$beta_hat_2
    betaXG <- dat$beta_hat_1
    sebetaYG <- dat$seb2
    sebetaXG <- dat$seb1
    betaIV= betaYG/betaXG # ratio estimates
    weights = (sebetaYG/betaXG)^-2 # inverse-variance weights
    betaIVW = sum(betaYG*betaXG*sebetaYG^-2)/sum(betaXG^2*sebetaYG^-2) # IVW estimate

    betaWM = weighted.median(betaIV, weights) # weighted median estimate
    sebetaWM  = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights) # standard error
    z_wm <- betaWM/sebetaWM
    ret <- list(z = z_wm,  p = 2*pnorm(abs(z_wm), lower.tail=FALSE))
    return(ret)
}

#'@export
pweighted_median <- function(dat, p_val_thresh=5e-8, no_ld = FALSE){
    if(no_ld) dat <- process_dat_nold(dat)
    dat <- dat %>% filter(p_value < p_val_thresh & ld_prune == TRUE)
    if(nrow(dat)< 3){
        ret <- list(z = NA,  p = NA)
        return(ret)
    }

    betaYG <- dat$beta_hat_2
    betaXG <- dat$beta_hat_1
    sebetaYG <- dat$seb2
    sebetaXG <- dat$seb1
    betaIV= betaYG/betaXG # ratio estimates
    weights = (sebetaYG/betaXG)^-2 # inverse-variance weights
    betaIVW = sum(betaYG*betaXG*sebetaYG^-2)/sum(betaXG^2*sebetaYG^-2) # IVW estimate

    penalty = pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
    pen.weights = weights*pmin(1, penalty*20) # penalized weights

    betaPWM= weighted.median(betaIV, pen.weights) # penalized weighted median estimate
    sebetaPWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, pen.weights) # standard error
    z_pwm <- betaPWM/sebetaPWM

    ret <- list(z=z_pwm,  p = 2*pnorm(abs(z_pwm), lower.tail=FALSE))
    return(ret)
}
