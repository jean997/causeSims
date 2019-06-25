#nohup dsc --replicate 100 --host config-broadwl.yml -c 4 false_positives.dsc > fp_dsc.out &

DSC:
    define:
        mr: ivw, egger, mrp, wm, twas, wm, gsmr
    run: simulate*gw_sig,
         simulate*LCV,
         simulate*mr,
         simulate*cause_params*cause
    replicate: 20
    output: fp
    exec_path: R


######## Simulate
simulate: R(library(causeSims);
            snps <- readRDS("data/chr19_snpdata_hm3only.RDS");
            evd_list <- readRDS("data/evd_list_chr19_hm3.RDS");
            neffect1 <- effects[1];
            neffect2 <- effects[2];
            dat <- sum_stats(snps, evd_list,
                  n_copies = 30,
                  n1 = n1, n2=n2,
                  h1=h1, h2=h1,
                  neffect1 = neffect1, neffect2 = neffect2,
                  tau = tau, omega = omega, q = q,
                  cores = 4, ld_prune_pval_thresh = 0.01,
                  r2_thresh = 0.1))

     q: 0.1, 0.2, 0.3, 0.4, 0.5
     omega: 0.02, 0.05
     tau: 0
     n1: 12000, 40000
     n2: 12000, 40000
     effects: c(1000, 1000)
     h1: 0.25
     h2: 0.25
     $sim_params: c(q = q, omega = omega, tau = tau,  h1 = h1, h2 = h2, n1 = n1, n2 = n2, neffect1 = neffect1, neffect2 =neffect2)
     $dat: dat


# Count how many variants are genome-wide significant for M. Also collect parameters,
# it will be faster to get them from this module than from simulate because the objects are smaller
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
          gamma <- sqrt(tau*sum(h2)/sum(h1));
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


######### CAUSE #############
# Parameter estimation is separate from analysis for flexibility

cause_params: R(library(causeSims);
                p1 <- cause_params_sims($(dat), null_wt = null_wt, no_ld=FALSE);
                )
                #p2 <- cause_params_sims($(dat), null_wt = null_wt, no_ld=TRUE))
    null_wt: 10
    $cause_params_ld: p1

cause: R(library(causeSims);
         library(cause);
         cause_res <- cause_sims($(dat), $(cause_params_ld), no_ld = FALSE);
         z <- -1*summary(cause_res)$z;
         p <- pnorm(-z);
         quants <- summary(cause_res)$quants)
    $cause_res: cause_res
    $sigma_g: cause_res$sigma_g
    $eta_med_2: quants[[1]][1,2]
    $q_med_2: quants[[1]][1,3]
    $eta_med_3: quants[[2]][1,2]
    $gamma_med_3: quants[[2]][1,1]
    $q_med_3: quants[[2]][1,3]
    $z: z
    $p: p

## Other MR methods
ivw: R(library(causeSims);
       res <- ivw($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $z: res$z
    $p: res$p

egger: R(library(causeSims);
       res <- egger($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $z: res$z
    $p: res$p

mrp: R(library(causeSims);
       res <- mrpresso($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $z: res$z
    $p: res$p

twas: R(library(causeSims);
       res <- twas($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $z: res$z
    $p: res$p

wm: R(library(causeSims);
       res <- weighted_median($(dat), p_val_thresh=thresh, no_ld = no_ld);
       )
    thresh: 5e-8
    no_ld: FALSE
    $z: res$z
    $p: res$p

gsmr: R(library(causeSims);
        evd_list <- readRDS("data/evd_list_chr19_hm3.RDS");
        n <- with($(dat), sum(p_value < 5e-8));
        if(n < 3){
            res <- NULL;
        }else{
            res <- gsmr_sims($(dat), evd_list, p_val_thresh  = 5e-8, no_ld = FALSE);
        };
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
    $est_gsmr: est
    $gsmr: res


##LCV
LCV: R(library(causeSims);
       res <- try(lcv_sims($(dat),no_ld = no_ld, sig.thresh = thresh));
       if(class(res) == "try-error"){
            res <- list(pval.gcpzero.2tailed=NA, gcp.pm = NA, gcp.pse = NA);
        };
       )
    thresh: 30
    no_ld: FALSE
    $p: res$pval.gcpzero.2tailed
    $gcp_med: res$gcp.pm
    $gcp_pse: res$gcp.pse
    $gcp_obj: res


