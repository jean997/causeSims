#'@title Reverse a simulated data set
#'@param dat Data set simulated with sum_stats
#'@param evd_list List of eigen value decompositions for LD regions
#'@param ld_prune_pval_thresh Minimum p-value to be included in LD pruning
#'@param r2_thres r^2 threshold for LD pruning
#'@export
reverse_data <- function(dat, evd_list, ld_prune_pval_thresh = 1e-3, r2_thresh = 0.1){
  dat <- dat %>%
         rename(bh1 = beta_hat_2, beta_hat_2 = beta_hat_1,
                s1 = seb2, seb2 = seb1,
                l1 = ld_b2, ld_b2 = ld_b1,
                bb1 = b2, b2 = b1,
                bh1n = beta_hat_2_nold, beta_hat_2_nold = beta_hat_1_nold) %>%
          rename(beta_hat_1 = bh1, seb1 = s1, ld_b1 = l1, b1 = bb1,
                 beta_hat_1_nold = bh1n) %>%
          mutate(p_value = 2*pnorm(-abs(beta_hat_1/seb1)),
                 p_value_nold = 2*pnorm(-abs(beta_hat_1_nold/seb1)))
  # New LD pruning
  # Relying on the fact that dat was made by sum_stats so the variants are in order of
  # replicates of block 1, replicates of block 2, etc
  n_snp <- with(dat, sum(rep==1))
  n_copies <- nrow(dat)/n_snp
  region_ids <- unique(dat$region_id)
  for(r in region_ids){
    R <- with(evd_list[[r]], crossprod(sqrt(values)*t(vectors)))
    keep <- sapply(seq(n_copies), function(i){
      strt <- with(dat, min(which(region_id==r & rep==i)))
      stp <- with(dat, max(which(region_id==r & rep==i)))
      causeSims:::ld_prune_cormat(R, dat$snp[strt:stp], dat$p_value[strt:stp],  ld_prune_pval_thresh, r2_thresh)
    }) %>% unlist()
    strt <- min(which(dat$region_id==r ))
    stp <- max(which(dat$region_id==r ))
    dat$ld_prune[strt:stp] <- dat$snp[strt:stp] %in% keep
  }
  return(dat)
}
