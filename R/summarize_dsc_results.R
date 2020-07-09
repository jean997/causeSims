
#'@export
summarize_dsc_results <- function(causedf, lcvdf, mrdf, direction = c("both", "pos", "neg"),
                              p_thresh = 0.05, gcp_mean_thresh = 0.6){

    direction <- match.arg(direction)
    if(direction=="both"){
      sgn <- c(-1, 1)
    }else if(direction=="pos"){
      sgn <- c(1)
    }else{
      sgn <- c(-1)
    }
    cause_summ <- causedf %>%
                    group_by(q, tau, omega,  n1, n2, cause.qbeta) %>%
                    summarize(n_sig = sum(cause.p < p_thresh & sign(cause.gamma_3_med) %in% sgn, na.rm=T),
                              missing = sum(is.na(cause.p)),
                              gamma = first(gamma),
                              eta = first(eta)) %>%
                    mutate(analysis = paste0("cause_", cause.qbeta)) %>%
                    select(-cause.qbeta) %>% ungroup()

    lcv_summ <- lcvdf %>%
                    group_by(q, tau, omega,  n1, n2) %>%
                    summarize(n_sig = sum(LCV.p < p_thresh & LCV.gcp_mean > 0, na.rm=T),
                              missing = sum(is.na(LCV.p)),
                              gamma = first(gamma),
                              eta = first(eta)) %>%
                    mutate(analysis = "lcv_p")   %>% ungroup()

    lcv_summ2 <- lcvdf %>%
                    group_by(q, tau, omega,  n1, n2) %>%
                    summarize(n_sig = sum(LCV.gcp_mean > gcp_mean_thresh, na.rm=T),
                              missing = sum(is.na(LCV.p)),
                              gamma = first(gamma),
                              eta = first(eta)) %>%
                    mutate(analysis = "lcv_mean") %>% ungroup()

    mr_summ <- mrdf %>%
                    group_by(q, tau, omega,  n1, n2, mr) %>%
                    summarize(n_sig = sum(mr.p < p_thresh & sign(mr.z) %in% sgn, na.rm=T),
                              missing = sum(is.na(mr.p)),
                              gamma = first(gamma),
                              eta = first(eta)) %>%
                    rename(analysis = "mr") %>% ungroup()
    res <- bind_rows(cause_summ, lcv_summ, lcv_summ2, mr_summ)


   # res <- res %>%
   #      mutate(n1 = factor(paste0("N[M]==", n1)),
   #             n2 = factor(paste0("N[Y]==", n2)))



    return(res)
}

#'@export
roc_data <- function(causedf, lcvdf, mrdf){
  ###### ROC Curves
  x1 <- causedf %>%
        mutate(analysis = paste0("cause_", qbeta), stat = -log10(cause.p)) %>%
        select(tag, q, omega, tau, n1, n2, analysis, stat)

  x2 <- lcvdf %>%
        mutate(lcv_p = pmax(0, sign(LCV.gcp_mean)*(-log10(LCV.p))),
               lcv_mean = LCV.gcp_mean) %>%
        select(tag, q, omega, tau, n1, n2, lcv_p, lcv_mean) %>%
        gather("analysis", "stat", -tag, -q, -omega, -tau, -n1, -n2) %>%
        select(tag, q, omega, tau, n1, n2, analysis, stat)

  x3 <- mrdf %>%
        mutate(stat = -log10(mr.p)) %>%
        rename(analysis = mr) %>%
        select(tag, q, omega, tau, n1, n2, analysis, stat)

  x <- bind_rows(x1, x2, x3)
  return(x)
}
