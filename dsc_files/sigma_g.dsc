#These simulations are for supplementary note section 5.2

DSC:
    run: simulate*gw_sig,
         simulate*cause_params*cause_sigma_g
    replicate: 100
    output: sigma_g 

simulate: R(library(causeSims);
            snps <- readRDS("data/chr19_snpdata_hm3only.RDS");
            evd_list <- readRDS("data/evd_list_chr19_hm3.RDS");
            dat <- sum_stats(snps, evd_list,
                  n_copies = 30,
                  n1 = n1, n2=n2,
                  h1=h1, h2=h1,
                  neffect1 = neffect1, neffect2 = neffect2,
                  tau = qot["tau"], omega = qot["omega"], q = qot["q"],
                  cores = 4, ld_prune_pval_thresh = 0.01,
                  r2_thresh = 0.1))

     qot: c(q =0, omega = 0, tau = 0),
          c(q = 0, omega = 0, tau = 0.05),
          c(q = 0.3, omega = 0.05, tau = 0)
     n1: 40000
     n2: 40000
     neffect1: 1000
     neffect2: 1000
     h1: 0.25
     h2: 0.25
     $sim_params: c(qot,  h1 = h1, h2 = h2, n1 = n1, n2 = n2, neffect1 = neffect1, neffect2 =neffect2)
     $dat: dat


# Count how many variants are genome-wide significant for M. Also collect parameters,
# it will be faster to get them from this module than from simulate because the objects are smaller
gw_sig: R(library(causeSims);
           m_sig_nold <- with($(dat), sum(p_value_nold < thresh ));
           y_sig_nold <- with($(dat), sum(2*pnorm(-abs(beta_hat_2_nold/seb2)) < thresh ));
           m_sig <- with($(dat), sum(p_value < thresh & ld_prune==TRUE));

          params <- $(sim_params);
          tau <- as.numeric(params["tau"]);
          omega <- as.numeric(params["omega"]);
          q <- as.numeric(params["q"]);
          h1 <- as.numeric(params["h1"]);
          h2 <- as.numeric(params["h2"]);
          neffect1 <- as.numeric(params["neffect1"]);
          neffect2 <- as.numeric(params["neffect2"]);
          gamma <- sqrt(tau*sum(h2)/sum(h1))*sign(tau);
          eta <- sqrt(abs(omega)*sum(h2)/sum(h1))*sign(omega);

          #LCV parameters
          if(q == 0 & tau == 0){
            q.1 <- q.2 <- 0;
            gcp <- 0;
          }else if(q == 0){
            q.1 <- 1;
            q.2 <- 0;
            gcp <- 1;
          }else{
            q.1 <- sqrt(q);
            q.2 <- eta * sqrt(q) * sqrt(h1) / sqrt(h2);
            gcp = (log(q.2^2) - log(q.1^2)) / (log(q.2^2) + log(q.1^2));
          };
          gcp_mom <- gcp_moments($(dat), h1, h2)
          )
    thresh: 5e-8
    $m_sig: m_sig
    $m_sig_nold: m_sig_nold
    $y_sig_nold: y_sig_nold
    $q: q
    $omega: omega
    $tau: tau
    $gamma: gamma
    $eta: eta
    $h1: h1
    $h2: h2
    $n1: as.numeric(params["n1"])
    $n2: as.numeric(params["n2"])
    $neffect1: as.numeric(params["neffect1"])
    $neffect2: as.numeric(params["neffect2"])
    $lcv_q1: q.1
    $lcv_q2: q.2
    $lcv_gcp: gcp
    $lcv_gcp_mom: gcp_mom$gcp

cause_params: R(library(causeSims);
                p1 <- cause_params_sims($(dat), null_wt = null_wt, no_ld=FALSE);
                )
    null_wt: 10
    $cause_params_ld: p1

cause_sigma_g: R(library(causeSims);
                 library(cause);
                 effect <- sqrt(0.05);
                 sigma_g <- get_sigma(effect, quant);
                 cause_res <- cause_sims($(dat), $(cause_params_ld), sigma_g = sigma_g,
                                          no_ld = FALSE);
                 z <- -1*summary(cause_res)$z;
                 p <- pnorm(-z);
                 sigma_g <- cause_res$sigma_g;
                 quants <- summary(cause_res)$quants;
                 eta_med_2 <-  quants[[1]][1,2];
                 q_med_2 <- quants[[1]][1,3];
                 eta_med_3 <- quants[[2]][1,2];
                 gamma_med_3 <- quants[[2]][1,1];
                 q_med_3 <- quants[[2]][1,3];
          )
    quant: 0.51, 0.65, 0.8
    $cause_res: cause_res
    $sigma_g: sigma_g
    $eta_med_2: eta_med_2
    $q_med_2: q_med_2
    $eta_med_3: eta_med_3
    $gamma_med_3: gamma_med_3
    $q_med_3: q_med_3
    $z: z


