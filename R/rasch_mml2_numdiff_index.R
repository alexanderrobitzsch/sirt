## File Name: rasch_mml2_numdiff_index.R
## File Version: 1.127



#** general function for numerical differentiation
#** diffindex aggregates across super items
rasch_mml2_numdiff_index <- function( pjk, pjk1, pjk2, n.ik, diffindex,
        max.increment, numdiff.parm, eps=1e-16, shortcut=TRUE )
{
    eps2 <- 1e-10
    h <- numdiff.parm
    an.ik <- n.ik

    ll0 <- rowSums( an.ik * log(pjk+eps) )
    ll0 <- rowsum(ll0, diffindex)[,1]

    if (! shortcut){

        ll1 <- rowSums( an.ik * log(pjk1+eps) )
        ll2 <- rowSums( an.ik * log(pjk2+eps) )
        # ll0 <- stats::aggregate( ll0, list(diffindex), sum )[,2]
        # ll1 <- stats::aggregate( ll1, list(diffindex), sum )[,2]
        # ll2 <- stats::aggregate( ll2, list(diffindex), sum )[,2]
        ll1 <- rowsum(ll1, diffindex)[,1]
        ll2 <- rowsum(ll2, diffindex)[,1]

        d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?

        # second order derivative
        # f(x+h)+f(x-h)=2*f(x) + f''(x)*h^2
        d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2

    }

    if (shortcut){

        p1 <- ( pjk1 - pjk2 ) / (2*h)
        ll1a <- rowSums( an.ik * p1 / pjk )
        ll1a <- rowsum(ll1a, diffindex)[,1]
        d1 <- ll1a

        p2 <- ( pjk1 + pjk2 - 2*pjk) / h^2
        ll1a <- rowSums( an.ik * (p2*pjk-p1^2)/pjk^2 )
        ll1a <- rowsum(ll1a, diffindex)[,1]
        d2 <- ll1a

    }

    # change in item difficulty
    d2[ abs(d2) < eps2 ] <- eps2
    increment <- - d1 / d2
    ci <- ceiling( abs(increment) / ( abs( max.increment) + eps2 ) )
    increment <- ifelse( abs( increment) > abs(max.increment),
                                increment/(2*ci), increment )

    # handle issues
    ind <- which( ! is.finite(increment) )
    if (length(ind)>0){
        increment[ind] <- 0
        d2[ind] <- 1e10
    }

    max_increment <- max( abs( increment ) )
    se <- sqrt( abs( 1 / d2 ) )
    increment <- increment[ match(diffindex, sort(unique(diffindex))) ]
    increment[ is.na(increment) ] <- 0
    #-- output
    res <- list(increment=increment, d2=d2, ll0=ll0, max_increment=max_increment,
                    se=se)
    return(res)
}


.mml2.numdiff.index <- rasch_mml2_numdiff_index
