#nohup dsc --replicate 100 --host config.yml -c 4 sigma_g.dsc > sg.out &
#
#These simulations are for supplementary note section SN2

DSC:
    run: simulate*data_summ,
         simulate*cause_params*cause
    replicate: 100
    output: sigma_g 

# Simulate data
# We consider three scnearios
# a) gamma = tau = q = 0
# b) gamma = sqrt(0.05), q = 0, eta = 0
# c) gamma = 0, q = 0.3, eta = sqrt(0.05)
# We are interested in knowing how the prior on gamma and eta (we use the same prior for both)
# affects estimation of gamma and eta in each scenario
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
     $sim_params: c(q = qot["q"], omega = qot["omega"], tau = qot["tau"],  h1 = h1, h2 = h2, n1 = n1, n2 = n2, neffect1 = neffect1, neffect2 =neffect2)
     $dat: dat


# This module computes some summaries about the data
data_summ: R(library(causeSims);
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

### CAUSE
cause_params: R(library(causeSims);
                p1 <- cause_params_sims($(dat), null_wt = null_wt, no_ld=FALSE);
                )
    null_wt: 10
    $cause_params_ld: p


# In this module "quant" determines the value of sigma_{gamma eta} (variable name sigma_g) 
# which defines the prior on gamma and eta. 
# The prior distributions on gamma and eta are N(0, sigma_g).
# sigma_g is chosen so that p_prior(gamma < sqrt(0.05)) = quant
cause: R(library(causeSims);
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
    $z: z
    $p: p
    $eta_2_upper: qs[[1]][3,2]
    $q_2_upper: qs[[1]][3,3]
    $eta_2_lower: qs[[1]][2,2]
    $q_2_lower: qs[[1]][2,3]
    $eta_2_med: qs[[1]][1,2]
    $q_2_med: qs[[1]][1,3]
    $eta_3_upper: qs[[2]][3,2]
    $q_3_upper: qs[[2]][3,3]
    $gamma_3_upper: qs[[2]][3,1]
    $eta_3_lower: qs[[2]][2,2]
    $q_3_lower: qs[[2]][2,3]
    $gamma_3_lower: qs[[2]][2,1]
    $eta_3_med: qs[[2]][1,2]
    $q_3_med: qs[[2]][1,3]
    $gamma_3_med: qs[[2]][1,1]


