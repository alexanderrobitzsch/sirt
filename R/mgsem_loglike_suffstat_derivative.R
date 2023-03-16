## File Name: mgsem_loglike_suffstat_derivative.R
## File Version: 0.16


mgsem_loglike_suffstat_derivative <- function(suffstat, Mu, Sigma )
{
    requireNamespace("MASS")
    N <- suffstat$N
    M <- suffstat$M
    S <- suffstat$S
    if (missing(Sigma)){
        res <- Mu
        Mu <- res$Mu
        Sigma <- res$Sigma
    }
    S1 <- MASS::ginv(Sigma)
    p <- length(Mu)
    m1 <- M-Mu

    #*** mean
    dermean <- as.vector( N*( crossprod(m1, S1 )))

    #*** covariance matrix
    # dercov <- dercov0 <- matrix(0, nrow=p, ncol=p)
    y <- S1 %*% m1
    S2 <- S %*% S1
    S3 <- S1 %*% S2
        #    for (ii in 1:p){
        #        for (jj in ii:p){
        #            # Sigma_der <- add_increment(dercov0, 1, ii,jj, symm=TRUE)
        #            # t2 <- -( t(y) %*% Sigma_der %*% y )[1,1]
        #            t2 <- -ifelse(ii==jj, y[ii]^2, y[ii]*y[jj]+y[jj]*y[ii])
        #            # t3 <- - sum( diag( S3 %*% Sigma_der ) )
        #            t3 <- -ifelse(ii==jj, S3[ii,jj], 2*S3[ii,jj] )
        #            # t4 <- sum(diag(S1 %*% Sigma_der))
        #            t4 <- ifelse(ii==jj, S1[ii,jj], 2*S1[ii,jj] )
        #            dercov[ii,jj] <- dercov[jj,ii] <- -N/2*(t2+t3+t4)
        #        }
        #    }
    #-- output
    res <- list(dermean=dermean, y=y, S1=S1, S2=S2, S3=S3)
    return(res)
}
