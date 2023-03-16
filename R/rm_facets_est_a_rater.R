## File Name: rm_facets_est_a_rater.R
## File Version: 0.14
## File Last Change: 2018-12-30


#####################################################
# estimation of slope parameter for rater
rm_facets_est_a_rater <- function( b.item, b.rater, Qmatrix, tau.item,
        VV, K, I, TP, a.item, a.rater, item.index, rater.index,
        n.ik, numdiff.parm=.001, max.b.increment=1,theta.k, msteps,
        mstepconv, a.rater.center, a.rater.fixed, a_lower, a_upper )
{
    h <- numdiff.parm
    diffindex <- rater.index
    RR <- length(b.rater)
    cat("  M steps a.rater parameter  |")
    it <- 0
    conv1 <- 1000
    #--- args clacprobs
    args <- list( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP,
                    a.item=a.item, a.rater=a.rater, item.index=item.index, rater.index=rater.index,
                    theta.k=theta.k, RR=RR )
    #--- begin M-steps
    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        a.rater0 <- a.rater
        args$a.rater <- a.rater
        pjk <- do.call( what=rm_facets_calcprobs, args=args)
        args$a.rater <- a.rater + h
        pjk1 <- do.call( what=rm_facets_calcprobs, args=args)
        args$a.rater <- a.rater - h
        pjk2 <- do.call( what=rm_facets_calcprobs, args=args)
        #-- increments
        res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex,
                    max.increment=max.b.increment, numdiff.parm=numdiff.parm )
        a.rater <- a.rater + res$increment
        #-- bound parameter estimates
        a.rater <- rm_squeeze(x=a.rater, lower=a_lower, upper=a_upper )
        if ( ! is.null( a.rater.fixed) ){
            ind <- which( ! is.na( a.rater.fixed  ) )
            a.rater[ind] <- a.rater.fixed[ind]
        res$d2[ ind ] <- -1E10
        }
        #-- center rater discriminations
        a.rater <- rm_center_vector( vec=a.rater, center_type=a.rater.center, do_log=TRUE)
        conv1 <- max( abs( a.rater - a.rater0 ) )
        it <- it+1
        cat("-")
    }
    cat(" ", it, "Step(s) \n")
    res <- list(a.rater=a.rater, se.a.rater=sqrt( abs(-1/res$d2 )), ll=sum(res$ll0) )
    return(res)
}

.rm.facets.est.a.rater <- rm_facets_est_a_rater
