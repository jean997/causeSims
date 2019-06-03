
#'@export
update_ld_prune <- function(dat, evd_list, n_copies, ld_prune_pval_thresh = 1e-3, r2_thresh = 0.1){

  keep <-sapply(seq(length(evd_list)), function(r){
                R <- with(evd_list[[r]], crossprod(sqrt(values)*t(vectors)))
                p <- nrow(R)
                df <- filter(dat, region_id == r)
                k <- sapply(seq(n_copies), function(i){
                  strt <- ((i-1)*p) + 1
                  stp <- i*p
                  causeSims:::ld_prune_cormat(R, df$snp[strt:stp], df$p_value[strt:stp],  ld_prune_pval_thresh, r2_thresh)
                }) %>% unlist()
                return(k)
                }) %>% unlist()

  return(keep)

}
