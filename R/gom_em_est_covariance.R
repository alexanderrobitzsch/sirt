## File Name: gom_em_est_covariance.R
## File Version: 0.05


#-- GOM: estimation of covariance
gom_em_est_covariance <- function( f.qk.yi, Sigma, theta.kM, N  )
{
    Sigma.cov <- Sigma
    delta.theta <- 1
    hwt <- f.qk.yi
    theta.k <- theta.kM
    thetabar <- hwt %*% theta.k
    # calculation of mu
    mu <- colSums( thetabar ) / N
    mu[1] <- 0
    # calculation of the covariance matrix
    theta.k.adj <- theta.k - matrix( mu, nrow=nrow(theta.k),
                                            ncol=ncol(theta.k), byrow=TRUE)
    D <- 2
    for (dd1 in 1L:D){
        for (dd2 in dd1:D){
            tk <- theta.k.adj[,dd1]*theta.k.adj[,dd2]
            h1 <- ( hwt %*% tk ) * delta.theta
            Sigma.cov[dd1,dd2] <- sum( h1 ) / N
            if (dd1 < dd2 ){ Sigma.cov[dd2,dd1] <- Sigma.cov[dd1,dd2] }
        }
    }
    diag(Sigma.cov) <- diag(Sigma.cov) + 10^(-10)
    pi.k <- sirt_dmvnorm_discrete( theta.k, mean=mu, sigma=Sigma.cov, as_matrix=TRUE )
    res <- list(mu=mu, Sigma=Sigma.cov, pi.k=pi.k )
    return(res)
}

