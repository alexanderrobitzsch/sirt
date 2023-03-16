## File Name: sim.rasch.dep.R
## File Version: 0.22
## File Last Change: 2018-12-30


#**** simulation of Rasch copula model
sim.rasch.dep <- function( theta, b, itemcluster, rho )
{
    probmat <- stats::plogis( outer( theta, b, "-" ) )
    I <- length(b)
    n <- length(theta)
    # covariance matrix of dependencies
    cov.dep <- diag(1, I)
    clusters <- unique(itemcluster[ itemcluster > 0 ] )
    CC <- length(clusters)
    for (cc in 1:CC){
        v1 <- which( itemcluster==cc )
        for (ii in v1){
            for (jj in v1) {
                if ( ii!=jj ){
                    cov.dep[ii,jj] <- rho[cc]
                }
            }
        }
    }
    random.gen <- stats::pnorm( sirt_rmvnorm( n, mean=rep(0,I), sigma=cov.dep ) )
    dat <- 1 * ( probmat > random.gen )
    colnames(dat) <- paste0( "I", substring(100+1:I,2))
    return(dat)
}

