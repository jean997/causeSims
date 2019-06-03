

## Simulate summary statistics from LCV model, with no LD. See ExampleScript.R for usage.
#'@export
SimulateLCV <- function(M,N.1,N.2,h2.1,h2.2,q.1,q.2,p.pi,p.g1,p.g2){

    pi<-sample(c( rnorm(M*p.pi), rep(0,M*(1-p.pi)) ))
    if(p.pi > 0) pi<-pi / sqrt(var(pi))

    gamma.1<-sample(c( rnorm(M*p.g1), rep(0,M*(1-p.g1)) ))
    if(p.g1 > 0) gamma.1 <- gamma.1 / sqrt(var(gamma.1)) * sqrt(1-q.1^2)

    gamma.2<-sample(c( rnorm(M*p.g2), rep(0,M*(1-p.g2)) ))
    if(p.g2 > 0) gamma.2 <- gamma.2 / sqrt(var(gamma.2))* sqrt(1-q.2^2)

    b1<-q.1*pi + gamma.1
    b1<-b1 / sqrt(t(b1)%*%b1) * sqrt(h2.1)
    b2<-q.2*pi + gamma.2
    b2<-b2 / sqrt(t(b2)%*%b2) * sqrt(h2.2)

    a1<-b1+rnorm(M,0,sqrt(1/N.1))
    a2<-b2+rnorm(M,0,sqrt(1/N.2))

    alpha<-data.frame(a1,b1,a2,b2)


    dat <- alpha %>% rename(beta_hat_1 = a1,  beta_hat_2 = a2) %>%
      mutate(seb1 = 1/sqrt(N.1), seb2 = 1/sqrt(N.2),
             pval1 = 2*pnorm(-abs(beta_hat_1/seb1)),
             snp = row_number())

    return(dat)
}
