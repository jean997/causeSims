#nohup dsc --replicate 100 --host config.yml -c 4 mixture2.dsc > mix2.out &
#
# These simulations correspond to Supplementary Note Section SN6.3
DSC:
    define:
      mr: ivw_MR, egger_MR, wm_MR, mbe_MR, mrp, gsmr
    run: simulate*data_summ,
         simulate*LCV,
         simulate*mr,
         simulate*cause_params*cause
    replicate: 100
    output: mix2 


# Simulate data
# Each simulated data frame includes data both with and without LD
# In this DSC file we only use the "with LD" data
# The "no_ld" flag seen in lots of the functions controls which data we use
# So here we always have no_ld=FALSE
simulate: R(library(causeSims);
            snps <- readRDS("data/chr19_snpdata_hm3only.RDS");
            evd_list <- readRDS("data/evd_list_chr19_hm3.RDS");
            dat <- sum_stats(snps, evd_list,
                  n_copies = 30,
                  n1 = n1, n2=n2,
                  h1=h1, h2=h1,
                  neffect1 = neffect1, neffect2 = neffect2,
                  tau = tau, omega = omega, q = q,
                  cores = 4, ld_prune_pval_thresh = 0.01,
                  r2_thresh = 0.1))
     q: 0.1, 0.3
     omega: -0.02
     tau: 0.02, 0.08
     n1: 12000, 40000
     n2: 12000, 40000
     neffect1: 1000
     neffect2: 1000
     h1: 0.25
     h2: 0.25
     $sim_params: c(q = q, omega = omega, tau = tau,  h1 = h1, h2 = h2, n1 = n1, n2 = n2, neffect1 = neffect1, neffect2 =neffect2)
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


######### CAUSE #############
# Parameter estimation is separate from analysis for flexibility
# null_wt determines the dirichlet prior on the mixing parameters
cause_params: R(library(causeSims);
                p <- cause_params_sims($(dat), null_wt = null_wt, no_ld=FALSE);
                )
    null_wt: 10
    $cause_params_ld: p

# This mdule runs CAUSE with three different priors on q 
# The module returns the z-score and p-value for CAUSE
# as well as posterior medians and 95% credible intervals for 
# q, eta, and gamma in both models (gamma is only estimated in the causal model)
# In the variable names, the model is specified by a number
# 2 = sharing model (2 free parameters), 3 = causal model (three free parameters)
cause: R(library(causeSims);
         library(cause);
         cause_res <- cause_sims($(dat), $(cause_params_ld), no_ld = FALSE, 
                                qalpha = 1, qbeta = qbeta);
         z <- -1*summary(cause_res)$z;
         p <- pnorm(-z);
         qs <- summary(cause_res)$quants)
    qbeta: 10
    $cause_res: cause_res
    $sigma_g: cause_res$sigma_g
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


## Other MR methods

# MR-PRESSO
mrp: R(library(causeSims);
       res <- mrpresso($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $est: res$est
    $z: res$z
    $p: res$p


# GSMR
# GSMR uses LD information
gsmr: R(library(causeSims);
        evd_list <- readRDS("data/evd_list_chr19_hm3.RDS");
        res <- gsmr_sims($(dat), evd_list, p_val_thresh  = 5e-8, no_ld = FALSE);
        if(!is.null(res)){
           z <- res$bxy/res$bxy_se;
           est <- res$bxy;
           p <- res$bxy_pval;
        }else{
           z <- est <- p <-  NA;
        }
      )
    thresh: 5e-8
    $z: z
    $p: p
    $est: est
    $gsmr: res


# LCV
LCV: R(library(causeSims);
       res <- try(lcv_sims($(dat),no_ld = no_ld, sig.thresh = thresh));
       if(class(res) == "try-error"){
            res <- list(pval.gcpzero.2tailed=NA, gcp.pm = NA, gcp.pse = NA);
        };
       )
    thresh: 30
    no_ld: FALSE
    $p: res$pval.gcpzero.2tailed
    $gcp_mean: res$gcp.pm
    $gcp_pse: res$gcp.pse
    $gcp_obj: res


# These modules all use the MendelianRandomization package
# (inside wrappers from the causeSims package that make it easy to run with our data
ivw_MR: R(library(causeSims);
       res <- ivw_re_nonome_MR($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $est: res$est
    $z: res$z
    $p: res$p

egger_MR: R(library(causeSims);
       res <- egger_re_MR($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $est: res$est
    $z: res$z
    $p: res$p

wm_MR: R(library(causeSims);
       res <- wtd_median_MR($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $est: res$est
    $z: res$z
    $p: res$p

# This one is the weighted mode (Hartwig et al. 2017)
# We run with three different values of phi
#
mbe_MR: R(library(causeSims);
       res <-mbe_MR($(dat), p_val_thresh=thresh, no_ld = no_ld, 
                weighting=weighting, no_NOME=no_NOME, phi=phi);
       )
    thresh: 5e-8
    no_ld: FALSE
    weighting: "weighted"
    no_NOME: TRUE
    phi: 1, 0.5, 0.25
    $est: res$est
    $z: res$z
    $p: res$p

