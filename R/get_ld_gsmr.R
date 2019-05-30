library(intervals)
library(Matrix)

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
