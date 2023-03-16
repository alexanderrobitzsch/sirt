## File Name: marginal.truescore.reliability.R
## File Version: 0.171
## File Last Change: 2023-03-15


marginal.truescore.reliability <- function(  b, a=1+0*b, c=0*b, d=1+0*b,
        mean.trait=0, sd.trait=1, theta.k=seq( -6, 6, len=200) )
{
    TT <- length(theta.k)
    I <- length(b)
    phi.k <- sirt_dnorm_discrete( theta.k, mean=mean.trait, sd=sd.trait )
    phi.kM <- matrix( phi.k, nrow=TT, ncol=I )

    aM <- sirt_matrix2( a, nrow=TT)
    bM <- sirt_matrix2( b, nrow=TT)
    cM <- sirt_matrix2( c, nrow=TT)
    dM <- sirt_matrix2( d, nrow=TT)

    # item response functions
    theta.kM <- matrix( theta.k, nrow=TT, ncol=I )
    icc.theta <- cM + (dM-cM)* stats::plogis( aM*(theta.kM - bM ) )

    # calculate pi_i (predicted probabilities)
    pi.i <- colSums( icc.theta * phi.kM )

    #  expected number correct
    mu <- sum(pi.i)

    # compute error variance
    sig2.error <- colSums( icc.theta * ( 1 - icc.theta ) * phi.k )

    # compute total error variance
    sig2.total.error <- sum( sig2.error )

    # compute total true score variance
    iccmean <- rowMeans( icc.theta )
    sig2.total.tau <- sum( ( I*iccmean )^2 * phi.k ) - ( sum( I*iccmean*phi.k ) )^2

    # item level true score variance
    sig2.tau <- pi.i * ( 1 - pi.i ) - sig2.error

    # collect all formulas for item
    item <- data.frame( item=1:I, pi=pi.i, sig2.tau=sig2.tau, sig2.error=sig2.error )
    item$rel.item <- item$sig2.tau /  ( item$sig2.tau + item$sig2.error )

    rel.test <- sig2.total.tau / ( sig2.total.tau + sig2.total.error )
    cat('Reliability', '=', round( rel.test,3 ),'\n')
    # Formula (15)
    # sum( sqrt( outer( sig2.tau, sig2.tau ) ) )
    res <- list(rel.test=rel.test, item=item, pi=mu/I, sig2.tau=sig2.total.tau,
                    sig2.error=sig2.total.error )
    invisible(res)
}
