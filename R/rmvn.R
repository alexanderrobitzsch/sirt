## File Name: rmvn.R
## File Version: 0.04
## File Last Change: 2022-04-17

rmvn <- function(N, mu, Sigma, exact=TRUE)
{
    Sigma_svd <- svd(Sigma)
    D <- ncol(Sigma)
    dat0 <- matrix( stats::rnorm(N*D, mean=0, sd=1), ncol=D)
    dat00 <- dat0

    #-- compute data with exact zero means and identity covariance matrix
    if (exact){
        c1 <- stats::cov.wt(dat0, method="ML")
        dat0 <- dat0 - matrix( c1$center, nrow=N, ncol=D, byrow=TRUE)
        COV0 <- c1$cov
        c0 <- svd(COV0)
        c00 <- t(c0$u) %*% diag(1/sqrt(c0$d))
        c00 <- diag(1/sqrt(c0$d)) %*% t(c0$u)
        dat00 <- dat0 %*% t(c00)
    }

    #-- compute data with prescribed distribution
    c11 <- Sigma_svd$u %*% diag( sqrt(Sigma_svd$d) )
    dat1 <- dat00 %*% t(c11)
    dat1 <- dat1 + matrix( mu, nrow=N, ncol=D, byrow=TRUE)

    #--- output
    return(dat1)
}
