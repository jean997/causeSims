library(dplyr)
library(gsmr)

gsmr_res <- function(dat, pval_thresh=5e-8, p_val_column="p_value"){
    dat <- dat %>%
            filter(get(p_val_column) < pval_thresh)
    dat$p1 <- dat[[p_val_column]]
    if(nrow(dat) < 3) return(NULL)    
    ldrho <- data.frame(diag(rep(1, nrow(dat))))
    names(ldrho) <- seq(nrow(dat))
    rownames(ldrho) <- seq(nrow(dat))
    res <-  try(with(dat, gsmr(beta_hat_1, seb1, p1, 
                           beta_hat_2, seb2, ldrho=ldrho, snpid=seq(nrow(dat)), 
                           n_ref = 1, nsnps_thresh=1, gwas_thresh=pval_thresh)))
    if(class(res)=="try-error") return(NULL)
    return(res)
}


