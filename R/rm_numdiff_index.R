## File Name: rm_numdiff_index.R
## File Version: 0.24
## File Last Change: 2019-05-03


####################################################################
# general function for numerical differentiation
# diffindex aggregates across super items
rm_numdiff_index <- function( pjk, pjk1, pjk2, n.ik, diffindex,
        max.increment, numdiff.parm, eps=1E-80, eps2=1E-10, prior=NULL, value=NULL )
{
    h <- numdiff.parm
    an.ik <- aperm( n.ik, c(2,3,1) )
    # str("pjk")
    # [items, categories, nodes]

    #--- evaluate expected likelihood
    ll0 <- rm_grouped_expected_likelihood(pjk=pjk, n.ik=an.ik, diffindex=diffindex, eps=eps)
    ll1 <- rm_grouped_expected_likelihood(pjk=pjk1, n.ik=an.ik, diffindex=diffindex, eps=eps)
    ll2 <- rm_grouped_expected_likelihood(pjk=pjk2, n.ik=an.ik, diffindex=diffindex, eps=eps)

    if (! is.null(prior)){
        if (length(value) > length(ll0)){
            value <- stats::aggregate( matrix(value,ncol=1), list(diffindex), mean )[,2]
        }
        M <- prior[1]
        SD <- sqrt(prior[2])
        ll0 <- ll0 + sirt_dnorm( value, mean=M, sd=SD, log=TRUE)
        ll1 <- ll1 + sirt_dnorm( value+h, mean=M, sd=SD, log=TRUE)
        ll2 <- ll2 + sirt_dnorm( value-h, mean=M, sd=SD, log=TRUE)
    }

    #-- discrete differences
    res <- rm_numdiff_discrete_differences(ll0=ll0, ll1=ll1, ll2=ll2, h=h)
    d1 <- res$d1
    d2 <- res$d2

    #-- compute incremental change in item parameters
    d2[ abs(d2) < eps2 ] <- eps2
    increment <- - d1 / d2


    #-- trim increment
    increment <- rm_numdiff_trim_increment( increment=increment, max.increment=max.increment, eps2=eps2 )
    #--- output
    res <- list(increment=increment, d1=d1, d2=d2, ll0=ll0, h=h)
    return(res)
}


.rm.numdiff.index <- rm_numdiff_index
