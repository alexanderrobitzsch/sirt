## File Name: rasch_mml2_mstep_raschtype.R
## File Version: 0.154
## File Last Change: 2019-10-27


#--- M-Step Rasch Type Model
rasch_mml2_mstep_raschtype <- function( theta.k, b, n.k, n, n.jk, r.jk, pi.k, I,
                    conv1, constraints, mitermax, pure.rasch, trait.weights, fixed.a,
                    fixed.c, fixed.d, alpha1, alpha2, h=.0025, designmatrix, group,
                    numdiff.parm, Qmatrix=NULL, old_increment, est.b, center.b,
                    min.b, max.b, prior.b=NULL )
{
    abs.change <- 1
    miter <- 0
    # group estimation
    G <- ncol(n.k)
    b00 <- b
    # number of subjects within groups
    n <- colSums(n.k)
    if ( is.null(Qmatrix) ){
            NT <- length(theta.k)
    } else {
        if ( is.matrix(theta.k) ){
            NT <- nrow(theta.k)
        } else {
            NT <- length(theta.k)
        }
    }

    #*** begin loop
    eps <- numdiff.parm
    h <- h1 <- eps
    h2 <- 1 + 2*eps

    #-- start iterations M-step
    while( abs.change > conv1 & miter < mitermax ){

        b0 <- b

        probs_args <- list( theta.k=theta.k, b=b,
                    fixed.a=fixed.a, fixed.c=fixed.c, fixed.d=fixed.d, alpha1=alpha1,
                    alpha2=alpha2, Qmatrix=Qmatrix, h=0, h1=h1, h2=h2, incr="b" )
        res <- do.call(what=rasch_mml2_prob_genlogis_4pl_evaluate, args=probs_args )
        pjk.M <- res$pjk.M
        qjk.M <- res$qjk.M

        probs_args$h <- h
        res <- do.call(what=rasch_mml2_prob_genlogis_4pl_evaluate, args=probs_args )
        pjk1.M <- res$pjk.M
        qjk1.M <- res$qjk.M

        probs_args$h <- -h
        res <- do.call(what=rasch_mml2_prob_genlogis_4pl_evaluate, args=probs_args )
        pjk2.M <- res$pjk.M
        qjk2.M <- res$qjk.M

        # update item parameters
        ll0 <- ll1 <- ll2 <- matrix( 0, nrow=nrow(n.jk), ncol=G)
        for ( gg in 1:G){
            exp_r <- r.jk[,,gg]
            exp_n <- n.jk[,,gg]
            ll0[,gg] <- rasch_mml2_mstep_calc_loglike( exp_r=exp_r, prob1=pjk.M,
                                exp_n=exp_n, prob0=qjk.M)
            ll1[,gg] <- rasch_mml2_mstep_calc_loglike( exp_r=exp_r, prob1=pjk1.M,
                                exp_n=exp_n, prob0=qjk1.M)
            ll2[,gg] <- rasch_mml2_mstep_calc_loglike( exp_r=exp_r, prob1=pjk2.M,
                                exp_n=exp_n, prob0=qjk2.M)
        }
        # sum across all groups
        ll0 <- rowSums(ll0)
        ll1 <- rowSums(ll1)
        ll2 <- rowSums(ll2)

        # evaluate prior distribution
        parm <- b
        prior_fct <- stats::dnorm
        prior <- prior.b
        res <- rasch_mml2_raschtype_mstep_parameter_group_evaluate_prior(parm=parm,
                h=h, prior=prior, ll0=ll0, ll1=ll1, ll2=ll2, prior_fct=prior_fct )
        ll0 <- res$ll0
        ll1 <- res$ll1
        ll2 <- res$ll2

        #- derivatives
        res <- rasch_mml2_difference_quotient(ll0=ll0, ll1=ll1, ll2=ll2, h=h)
        d1 <- res$d1
        d2 <- res$d2

        #-- sum over contributions
        if ( ! is.null(est.b) ){
            d1 <- stats::aggregate( d1, list(est.b), sum )
            i1 <- d1[,1]
            d1 <- d1[,2]
            d2 <- stats::aggregate( d2, list(est.b), sum )[,2]
        }
        increment <- - d1 / d2
        increment <- sirt_trim_increment(increment=increment,
                            max_increment=max(abs(old_increment)))
        # define old_increment here
        old_increment <- increment
        if ( ! is.null(est.b) ){
            increment <- increment[ match( est.b, i1 ) ]
        }
        b <- b + increment
        # linear parameter constraints
        if ( ! is.null( designmatrix ) & is.null(est.b) ){
            mod <- stats::lm( b ~ 0 + designmatrix  )
            b <- stats::fitted(mod)
        }
        # last item is the sum of all other item difficulties
        center <- is.null(constraints)
        if ( !is.null( constraints) ){
            b[ constraints[,1] ] <- constraints[,2]
        }
        abs.change <- max( abs( b0 - b ) )
        miter <- miter+1
    }
    #*** end iteration loop

    #*** squeeze parameter estimates
    b <- squeeze.mml2( b, c( min.b, max.b ) )

    #-- center b
    if ( center.b ){
        D <- ncol(Qmatrix)
        for (dd in 1:D){
            ind.dd <- which( Qmatrix[,dd] > 0 )
            b[ind.dd] <- b[ind.dd] - sum( Qmatrix[ind.dd,dd] * b[ ind.dd] ) /
                    sum( Qmatrix[ind.dd,dd]  )
        }
    }
    #-- recompute probabilities
    pjk <- prob_genlogis_4pl(theta=theta.k, b=b, a=fixed.a, c=fixed.c, d=fixed.d,
                alpha1=alpha1, alpha2=alpha2, Qmatrix=Qmatrix)

    # calculate log likelihood
    ll <- sapply( 1:G, FUN=function(gg){
            sum( rasch_mml2_mstep_calc_loglike( exp_r=r.jk[,,gg], prob1=t(pjk), exp_n=n.jk[,,gg]) )
            } )
    ll <- sum(ll)
    #-- output
    res <- list( b=b, pi.k=pi.k, ll=ll, miter=miter, center=center, G=G,
                    old_increment=increment, se.b=sqrt( 1 /abs(d2) ) )
    return(res)
}


.m.step.raschtype <- rasch_mml2_mstep_raschtype
