## File Name: dirichlet.mle.R
## File Version: 0.13



#--- Maximum likelihood estimation of distribution parameters
dirichlet.mle <- function( x, weights=NULL, eps=10^(-5), convcrit=.00001,
        maxit=1000, oldfac=.3, progress=FALSE)
{
    N <- nrow(x)
    K <- ncol(x)
    # compute log pbar
    x <- ( x+eps ) / ( 1 + 2*eps )
    x <- x / rowSums(x)
    N <- nrow(x)
    if ( is.null(weights) ){
        weights <- rep(1,N)
    }
    weights <- N * weights / sum( weights )
    log.pbar <- colMeans( weights * log( x ) )
    # compute inits
    alphaprob <- colMeans( x * weights )
    p2 <- mean( x[,1]^2 * weights )
    xsi <- ( alphaprob[1] - p2 ) / ( p2 - ( alphaprob[1] )^2 )
    alpha <- xsi * alphaprob
    K1 <- matrix(1,K,K)
    conv <- 1
    iter <- 1

    #--- BEGIN iterations
    while( ( conv > convcrit ) & (iter < maxit) ){
        alpha0 <- alpha
        g <- N * digamma( sum(alpha ) ) - N * digamma(alpha) + N * log.pbar
        z <- N * sirt_digamma1( sum(alpha ))
        H <- diag( -N * sirt_digamma1( alpha ) ) + z
        alpha <- alpha0 - solve(H, g )
        alpha[ alpha < 0 ] <- 1e-10
        alpha <- alpha0 + oldfac*( alpha - alpha0 )
        conv <- max( abs( alpha0 - alpha ) )
        if (progress){
            print( paste( iter, sum(alpha), conv) )
            utils::flush.console()
        }
        iter <- iter+1
    }
    alpha0 <- sum(alpha)
    xsi <- alpha / alpha0
    res <- list( alpha=alpha, alpha0=alpha0, xsi=xsi )
    return(res)
}

