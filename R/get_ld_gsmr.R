
#'@title Get LD matrix for running GSMR
#'@export
get_ld_gsmr <- function(dat, evd_list, p_val_thresh  = 5e-8, no_ld = FALSE){

    if(no_ld){
      n <- with(dat, sum(p_value_nold < p_val_thresh))
      return(diag(rep(1, n)))
    }

    ix <- with(dat, which(p_value < p_val_thresh & ld_prune == TRUE))
    if(length(ix) < 1) return(NULL)
    #This is a data frame of all the region/replicate combos with a significant snp
    regions <- filter(dat, p_value < p_val_thresh & ld_prune ==TRUE) %>%
               select(region_id, rep) %>%
               unique()
    regions$n <- apply(regions, 1, function(x){
      with(dat, sum(p_value < p_val_thresh & ld_prune==T & region_id == x[1] & rep == x[2]))
    })
    ldmats <- list()
    for(i in seq(nrow(regions))){
        if(regions$n[i] == 1){
            ldmats[[i]] <- matrix(c(1), nrow=1)
        }else{
            R <- with(evd_list[[regions$region_id[i]]],  crossprod(sqrt(values)*t(vectors)))
            ix <- filter(dat, region_id == regions$region_id[i] & rep == regions$rep[i]) %>%
                 with(., which(p_value < p_val_thresh & ld_prune==TRUE))
            ldmats[[i]] <- R[ix, ix]
        }
    }
    ldrho <- as.matrix(bdiag(ldmats))
    return(ldrho)
}
