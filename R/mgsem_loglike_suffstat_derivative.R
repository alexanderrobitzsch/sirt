## File Name: mgsem_loglike_suffstat_derivative.R
## File Version: 0.172


mgsem_loglike_suffstat_derivative <- function(suffstat, Mu, Sigma )
{
    N <- suffstat$N
    M <- suffstat$M
    S <- suffstat$S
    if (missing(Sigma)){
        res <- Mu
        Mu <- res$Mu
        Sigma <- res$Sigma
    }
    S1 <- mgsem_ginv(X=Sigma)
    p <- length(Mu)
    m1 <- M-Mu

    #*** mean
    dermean <- as.vector( N*( crossprod(m1, S1 )))

    #*** covariance matrix
    y <- S1 %*% m1
    S2 <- S %*% S1
    S3 <- S1 %*% S2

    #-- output
    res <- list(dermean=dermean, y=y, S1=S1, S2=S2, S3=S3)
    return(res)
}
