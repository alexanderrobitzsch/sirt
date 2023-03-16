## File Name: rm_facets_est_tau_item.R
## File Version: 0.20
## File Last Change: 2018-12-30


#####################################################
# estimation of tau.item parameters
rm_facets_est_tau_item <- function( b.item, b.rater, Qmatrix, tau.item,
        VV, K, I, TP, a.item, a.rater, item.index, rater.index,
        n.ik, numdiff.parm=.001, max.b.increment=1,theta.k, msteps,
        mstepconv, tau.item.fixed, tau.item.fixed_val )
{

    h <- numdiff.parm
    diffindex <- item.index
    RR <- length(b.rater)
    Q0 <- matrix(0,nrow=VV, ncol=K)
    se.tau.item <- Q0
    cat("  M steps tau.item parameter |")
    it <- 0
    conv1 <- 1000
    #--- input rm_facets_calcprobs
    args <- list( b.item=b.item, b.rater=b.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K,
                    I=I, TP=TP, a.item=a.item, a.rater=a.rater, item.index=item.index,
                    rater.index=rater.index, theta.k=theta.k, RR=RR )
    #--- begin M-steps
    while( ( it < msteps ) & ( conv1 > mstepconv ) ){
        tau.item11 <- tau.item0 <- tau.item
        for (kk in 1:K){
            Q1 <- Q0
            Q1[,kk] <- 1
            #-- compute expected likelihood
            args$tau.item <- tau.item11
            pjk <- do.call( what=rm_facets_calcprobs, args=args)
            args$tau.item <- tau.item11 + h*Q1
            pjk1 <- do.call( what=rm_facets_calcprobs, args=args)
            args$tau.item <- tau.item11 - h*Q1
            pjk2 <- do.call( what=rm_facets_calcprobs, args=args)
            #-- compute increments
            res <- rm_numdiff_index( pjk=pjk, pjk1=pjk1, pjk2=pjk2, n.ik=n.ik, diffindex=diffindex,
                        max.increment=max.b.increment, numdiff.parm=numdiff.parm )
            increment <- Q1*matrix( res$increment, nrow=VV, ncol=K)
            tau.item <- tau.item + increment
            se.tau.item[,kk] <- sqrt(abs(-1/res$d2)    )
        }

        if ( ! is.null( tau.item.fixed_val ) ){
            MK <- ncol( tau.item.fixed_val )
            for ( kk in 1:MK){
                ind <- which( ! is.na( tau.item.fixed_val[,kk]) )
                if ( length(ind) > 0 ){
                    tau.item[ ind, kk] <- tau.item.fixed_val[ ind, kk]
                }
            }
        }
        conv1 <- max( abs( tau.item - tau.item0 ) )
        it <- it+1
        cat("-")
        if ( ! is.null(tau.item.fixed) ){
            tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
        }
    }
    cat(" ", it, "Step(s) \n")
    #-- output
    res <- list(tau.item=tau.item, se.tau.item=se.tau.item, ll=sum(res$ll0) )
    return(res)
}

.rm.facets.est.tau.item <- rm_facets_est_tau_item
