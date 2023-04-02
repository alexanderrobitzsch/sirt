## File Name: mgsem_loglike_suffstat.R
## File Version: 0.136

mgsem_loglike_suffstat <- function(suffstat, Mu, Sigma, output_all=FALSE )
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

    #*** mean structure
    # t1 <- ( t(M-Mu) %*% S1 %*% (M-Mu) )[1,1]
    # t1 <- ( crossprod(m1, S1) %*% m1 )[1,1]
    t1 <- sirt_rcpp_mgsem_quadform(y=M-Mu, A=S1)

    #*** covariance structure
    # t2b <- sum( diag(S %*% S1 ) )
    t2b <- sirt_rcpp_mgsem_trace_product(A=S, B=S1)
    # t2 <- log( det(Sigma) ) + t2b
    t2 <- sirt_logdet(x=Sigma) + t2b

    # normalization
    t3 <- p*log(2*pi)
    # total log likelihood
    res <- -N/2*(t1+t2+t3)

    #-- output_all
    if (output_all){
        res <- list(loglike=res, Sigma=Sigma, S1=S1, mean_residual=m1)
    }
    return(res)
}
