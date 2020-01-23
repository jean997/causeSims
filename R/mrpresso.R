
#'@title Run MR-PRESSO on simulated data
#'@param dat A simulated data frame created with sum_stats
#'@param p_val_thresh p-value threshold
#'@param no_ld Run with the nold data (T/F)
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
    t_mrp <- x$`Main MR results`[2,5]
    p_mrp <- x$`Main MR results`[2,6]
    e_mrp <- x$`Main MR results`[2,3]
    if(is.na(t_mrp)){
        t_mrp <- x$`Main MR results`[1,5]
        p_mrp <- x$`Main MR results`[1,6]
        e_mrp <- x$`Main MR results`[1,3]
    }
    R <- list(z = t_mrp, p = p_mrp, est = e_mrp)
    return(R)

}
