## File Name: brm.sim.R
## File Version: 0.124


#-- brm.sim
brm.sim <- function( theta, delta, tau, K=NULL)
{
    I <- length(delta)
    N <- length(theta)
    dat <- matrix( 0, nrow=N, ncol=I )
    colnames(dat) <- paste0( 'I', 1L:9 )
    if ( ! is.null(K) ){
        br <- seq( 0, 1, len=K+1 )
    }
    for (ii in 1L:I){
        m1 <- exp( ( theta - delta[ii] + tau[ii] ) / 2 )
        n1 <- exp( ( - theta + delta[ii] + tau[ii] ) / 2 )
        dat[,ii] <- stats::rbeta( N, shape1=m1, shape2=n1 )
        if ( ! is.null(K) ){
            dat[,ii] <- as.numeric(cut( dat[,ii], breaks=br )) - 1
        }
    }
    return(dat)
}
