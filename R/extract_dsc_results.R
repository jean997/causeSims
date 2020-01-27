
#'@export
extract_dsc_results <- function(dir, extract_cause=TRUE, extract_mr = TRUE,
                                extract_lcv = TRUE, sigma_g = FALSE, version1 = FALSE){

  if(version1){
    paramdf <- dscquery(dsc.outdir=dir,
                        targets=c("simulate",
                                  "gw_sig.tau","gw_sig.gamma", "gw_sig.q",
                                  "gw_sig.omega", "gw_sig.eta", "gw_sig.h1", "gw_sig.h2",
                                  "gw_sig.n1", "gw_sig.n2", "gw_sig.neffect1",
                                  "gw_sig.neffect2", "gw_sig.lcv_q1", "gw_sig.lcv_q2",
                                  "gw_sig.lcv_gcp",  "gw_sig.m_sig") ,
                        module.output.files = c("simulate"),
                        return.type="data.frame")
    paramdf <- paramdf %>% rename(q = gw_sig.q, tau = gw_sig.tau,
                                  eta = gw_sig.eta, gamma = gw_sig.gamma,
                                  omega = gw_sig.omega, h1 = gw_sig.h1,
                                  h2 = gw_sig.h2, n1 = gw_sig.n1, n2 = gw_sig.n2,
                                  neffect1 = gw_sig.neffect1, neffect2 = gw_sig.neffect2,
                                  m_sig = gw_sig.m_sig,
                                  lcv_q1 = gw_sig.lcv_q1, lcv_q2 = gw_sig.lcv_q2,
                                  gcp = gw_sig.lcv_gcp) %>%
      mutate(params = paste0("(", tau, ",", omega, ",", q, ")"))
    res <- list(params = paramdf)
  }else{
    paramdf <- dscquery(dsc.outdir=dir,
                     targets=c("simulate",
                               "data_summ.tau","data_summ.gamma", "data_summ.q",
                               "data_summ.omega", "data_summ.eta", "data_summ.h1", "data_summ.h2",
                               "data_summ.n1", "data_summ.n2", "data_summ.neffect1",
                               "data_summ.neffect2", "data_summ.lcv_q1", "data_summ.lcv_q2",
                               "data_summ.lcv_gcp",  "data_summ.m_sig") ,
                     module.output.files = c("simulate"),
                     return.type="data.frame")
    paramdf <- paramdf %>% rename(q = data_summ.q, tau = data_summ.tau,
                                    eta = data_summ.eta, gamma = data_summ.gamma,
                                    omega = data_summ.omega, h1 = data_summ.h1,
                                h2 = data_summ.h2, n1 = data_summ.n1, n2 = data_summ.n2,
                                neffect1 = data_summ.neffect1, neffect2 = data_summ.neffect2,
                                m_sig = data_summ.m_sig,
                                lcv_q1 = data_summ.lcv_q1, lcv_q2 = data_summ.lcv_q2,
                                gcp = data_summ.lcv_gcp) %>%
                  mutate(params = paste0("(", tau, ",", omega, ",", q, ")"))
    res <- list(params = paramdf)
  }
  #saveRDS(paramdf, file=paste0(lab, "_paramdf.RDS"))



  #cause
  if(extract_cause & ! sigma_g){
    causedf <- dscquery(dsc.outdir=dir,
                   targets=c("simulate", "cause", "cause.qbeta", "cause.sigma_g",
                             "cause.p", "cause.eta_2_med", "cause.q_2_med",
                             "cause.eta_2_lower", "cause.q_2_lower",
                             "cause.eta_2_upper", "cause.q_2_upper",
                             "cause.eta_3_med", "cause.q_3_med", "cause.gamma_3_med",
                             "cause.eta_3_lower", "cause.q_3_lower", "cause.gamma_3_lower",
                             "cause.eta_3_upper", "cause.q_3_upper", "cause.gamma_3_upper"),
                   module.output.files = c("simulate", "cause"),
                   return.type="data.frame")

      causedf <- causedf %>%
                 full_join(paramdf, ., by=c("DSC", "simulate.output.file")) %>%
                 mutate(tag = paste0(lab, "_", simulate.output.file))
      res[["cause"]] <- causedf
  }else if(extract_case & sigma_g){
    causedf <- dscquery(dsc.outdir=dir,
                        targets=c("simulate", "cause", "cause.qbeta", "cause.quant", "cause.sigma_g",
                                  "cause.p", "cause.eta_2_med", "cause.q_2_med",
                                  "cause.eta_2_lower", "cause.q_2_lower",
                                  "cause.eta_2_upper", "cause.q_2_upper",
                                  "cause.eta_3_med", "cause.q_3_med", "cause.gamma_3_med",
                                  "cause.eta_3_lower", "cause.q_3_lower", "cause.gamma_3_lower",
                                  "cause.eta_3_upper", "cause.q_3_upper", "cause.gamma_3_upper"),
                        module.output.files = c("simulate", "cause"),
                        return.type="data.frame")

     causedf <- causedf %>%
                full_join(paramdf, ., by=c("DSC", "simulate.output.file")) %>%
                mutate(tag = paste0(lab, "_", simulate.output.file))
     res[["cause"]] <- causedf
  }



# mr
  if(extract_mr){
    mrdf <- dscquery(dsc.outdir=dir,
                   targets=c("simulate", "mr",
                             "mr.phi", "mr.p",
                             "mr.z"),
                   module.output.files = c("simulate", "mr"),
                   return.type="data.frame")

    mrdf <- mrdf %>%
        full_join(paramdf, ., by=c("DSC", "simulate.output.file")) %>%
        mutate(mr = case_when(!is.na(mr.phi) ~ paste0(mr, "_", mr.phi),
                                   TRUE ~ mr)) %>%
             select(-mr.phi)
    res[["mr"]] <- mrdf
  }


  if(extrac_lcv){
    lcvdf <- dscquery(dsc.outdir=dir,
         targets=c("simulate",
                   "LCV.p", "LCV.gcp_mean", "LCV.gcp_pse"),
                   module.output.files = c("simulate"),
                   return.type="data.frame")

    lcvdf <- lcvdf %>%
             full_join(paramdf, ., by=c("DSC", "simulate.output.file")) %>%
             mutate(tag = paste0(lab, "_", simulate.output.file))
    res[["lcv"]] <- lcvdf
  }
  return(res)

}
