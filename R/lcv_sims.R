
#'@title Run LCV with simulated data
#'@param dat Simulated data set
#'@param no_ld Use no LD data
#'@param intercpets use crosstrait.interecpt=1 and ldsc.intercept=1 in RunLCV
#'@param sig.thresh Signicance threshold for LCV to estimate crosstrait.intercept and ldsc.intercept
#'@export
lcv_sims <- function(dat, no_ld = FALSE, intercepts = !no_ld, sig.thresh = 30){

  if(no_ld){
    dat <- process_dat_nold(dat) %>%
      mutate(ldscore_hm3 = 1)
  }
  if(!intercepts){
    crosstrait.intercept=0
    ldsc.intercept=0
  }else{
    crosstrait.intercept=1
    ldsc.intercept=1
  }


  lcv_res <- with(dat,
                   RunLCV(ldscore_hm3,
                          beta_hat_1/seb1,beta_hat_2/seb2,
                          no.blocks=100,crosstrait.intercept,ldsc.intercept,
                          sig.threshold=sig.thresh,n.1=1,n.2=1,intercept.12=0)
                   )
  return(lcv_res)
}
