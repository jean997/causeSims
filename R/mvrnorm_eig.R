
mvrnorm_eig1 <- function (n = 1, mu, eS, tol = 1e-06){
    p <- length(mu)
    if (length(eS$values) != p) stop("incompatible arguments")
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L])))
          stop("'Sigma' is not positive definite")
    X <- matrix(zrnorm(p * n), n)
    X <- drop(mu) + eS$vectors %*% (diag(sqrt(pmax(ev, 0)), p) %*% t(X))
    if (n == 1) drop(X)
       else t(X)
}

mvrnorm_eig <- function(n=1, mu, eS, tol=1e-6){
    p <- length(mu)
    if(class(eS)=="list"){
        nchunk <- length(eS)
        chunksize <- sapply(eS, function(x){length(x$values)})
        stopifnot(sum(chunksize)==p)
        strt <- cumsum(c(1, chunksize[-nchunk]))
        stp <- cumsum(chunksize)
        X <- sapply(1:nchunk, function(i){
                    #for(i in 1:nchunk){
                    #    cat(i, " ")
                    mvrnorm_eig1(n, mu[strt[i]:stp[i]], eS[[i]], tol)
                   })
        if(n==1) X <- unlist(X)
            else X <- do.call(cbind, X)
        return(X)
     }else{
        return(mvrnorm_eig1(n, mu, eS, tol))
     }
}

