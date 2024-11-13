## File Name: truescore.irt.R
## File Version: 0.252


#--- true score item response function
truescore.irt <- function( A, B, c=NULL, d=NULL, theta=seq(-3,3,len=21),
        error=NULL, pid=NULL, h=.001 )
{
    if ( is.vector(B) ){ B <- matrix( B, ncol=1 ) }
    if ( is.vector(A) ){ A <- matrix( A, ncol=1 ) }
    nB <- ncol(B)
    maxK <- nB - rowSums( is.na( B ))
    I <- nrow(B)
    if ( is.null(c) ){ c <- rep(0,I ) }
    if ( is.null(d) ){ d <- rep(1,I ) }
    if ( is.null(pid) ){ pid <- 1L:(length(theta)) }
    # true score
    truescore <- truescore_irt_irf( A=A, B=B, c=c, d=d, theta=theta )
    # calculate standard error of true score
    if ( is.null( error ) ){
        truescore.error <- NULL } else {
        truescore1 <- truescore_irt_irf( A=A, B=B, c=c, d=d, theta=theta + h )
        truescore.error <- sqrt( ( ( truescore1 - truescore ) /  h )^2 ) * error
    }
    percscore <- truescore / sum( maxK )
    percscore.error <- NULL
    if ( ! is.null( error ) ){
        percscore.error <- truescore.error / sum( maxK )
    }
    # optimization function values
    ind <- intersect( which( !is.na( theta ) ), which( !is.na(percscore) ) )
    x0 <- theta[ind]
    y0 <- percscore[ind]
    minf <- sum( ifelse( maxK==1, c, 0 ) ) / sum(maxK)
    maxf <- sum( ifelse( maxK==1, d, maxK ) ) / sum(maxK)
    critfct <- function(x) {
        a <- x[1]
        b <- x[2]
        # define fit function
        sum( ( y0 - minf - (maxf - minf) * stats::plogis( a*x0 + b ) )^2 )
    }
    h1 <- stats::optim( c( .5, 0 ), critfct )

    #--- OUTPUT
    res <- data.frame( pid=pid, truescore=truescore )
    res$truescore.error <- truescore.error
    res2 <- data.frame( percscore=percscore )
    res2$percscore.error <- percscore.error
    res3 <- data.frame( lower=minf, upper=maxf, a=h1$par[1], b=h1$par[2] )
    res <- cbind( res, res2, res3 )
    return(res)
}

