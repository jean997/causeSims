library(dplyr)
library(gsmr)
library(Matrix)
library(intervals)

get_ld_gsmr <- function(dat, ix){
    stopifnot(all(diff(ix) > 0))
    x <- rle(dat$region_id)
    n <- length(x$lengths)
    starts <- cumsum(c(1, x$lengths))[seq(n)]
    stops <- starts + x$lengths -1
    labels <- dat$region_id[starts]

    I <- Intervals(cbind(starts, stops))
    blocks <- unlist(interval_overlap(ix, I))
    block_rle <- rle(blocks)

    ldmats <- list()
    for(i in seq_along(block_rle$lengths)){
        if(block_rle$lengths[i] == 1){
            ldmats[[i]] <- matrix(c(1), nrow=1)
        }else{
            sr <- sum(c(1, block_rle$lengths)[seq(i)])
            sp <- sr + block_rle$lengths[i] -1
            bl <- labels[block_rle$values[i]]
            vars <- dat$region_ix[ix[sr:sp]]
            v <- readRDS(paste0("/project2/mstephens/sherlock2/chr19_EVD/th5/EVD/Q_", bl, ".RDS"))
            d <- readRDS(paste0("/project2/mstephens/sherlock2/chr19_EVD/th5/EVD/D_", bl, ".RDS"))
            R <- crossprod(sqrt(d)*t(v))
            ldmats[[i]] <- R[vars, vars]
        }
    }
    ldrho <- as.matrix(bdiag(ldmats))
    return(ldrho)
}

gsmr_res <- function(dat, ix, pval_thresh){

    if(length(ix) < 3) return(NULL)

    ldrho <- get_ld_gsmr(dat, ix)

    ldrho <- data.frame(ldrho)
    names(ldrho) <- seq_along(ix)
    rownames(ldrho) <- seq_along(ix)

    dat <- dat[ix,]
    dat$p_value <- with(dat, 2*pnorm(-abs(beta_hat_1/seb1)))
    res <-  try(with(dat, gsmr(beta_hat_1, seb1, p_value, 
                           beta_hat_2, seb2, ldrho=ldrho, snpid=seq_along(ix), 
                           n_ref = 1, nsnps_thresh=1, gwas_thresh=pval_thresh)))
    if(class(res)=="try-error") return(NULL)
    return(res)
}


