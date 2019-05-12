## File Name: gom_em_numdiff_index.R
## File Version: 0.13


#--- general function for numerical differentiation
#--- diffindex aggregates across super items
gom_em_numdiff_index <- function( pjk, pjk1, pjk2, an.ik, diffindex,
        max.increment, numdiff.parm, eps=1e-20 )
{
    h <- numdiff.parm
    ll0 <- rowSums( an.ik * log(pjk+eps) )
    ll1 <- rowSums( an.ik * log(pjk1+eps) )
    ll2 <- rowSums( an.ik * log(pjk2+eps) )

    #- derivatives
    res <- rasch_mml2_difference_quotient(ll0=ll0, ll1=ll1, ll2=ll2, h=numdiff.parm)
    d1 <- res$d1
    d2 <- res$d2
    increment <- - d1 / d2
    increment <- sirt_trim_increment(increment=increment, max_increment=max.increment)
    #-- output
    res <- list(increment=increment, d2=d2, ll0=ll0)
    return(res)
}


.gom.numdiff.index <- gom_em_numdiff_index
