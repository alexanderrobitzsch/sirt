## File Name: rasch_mml2_mstep_missing1.R
## File Version: 0.159
## File Last Change: 2022-10-18


#*** M-step for missing data model
rasch_mml2_mstep_missing1 <- function( theta.k, n.ik, mitermax, conv1,
        b, beta, delta.miss, pjk, numdiff.parm,
        constraints, est.delta, min.beta, est_delta, fixed.a=NULL, est.a=NULL,
        min.delta=-200, max.delta=200, irtmodel="missing1")
{
    I <- length(b)
    if (is.null(fixed.a)){
        fixed.a <- rep(1,fixed.a)
    }
    if (is.null(est.a)){
        est.a <- rep(0,I)
    }
    est_a <- sum(est.a) > 0
    se.a <- rep(0,I)

    h <- numdiff.parm
    se.delta <- 0
    miter <- 0
    diffindex <- est.b <- 1:I
    diffindex1 <- est.delta

    prob_fun <- rasch_mml2_calcprob_missing1

    max_incr_b <- 1
    max_incr_beta <- 1
    max_incr_delta <- 1
    max_incr_a <- 1
    miterchange <- 1000
    a0 <- fixed.a
    #---
    while( ( miterchange > conv1 ) & ( miter < mitermax ) ){

        #--- update b
        b0 <- b

        max.increment <- max_incr_b
        entry <- "b"
        diffindex <- est.b

        args0 <- list(theta.k=theta.k, b=b, beta=beta, delta.miss=delta.miss,
                    pjk=pjk, fixed.a=fixed.a, irtmodel=irtmodel)
        res <- rasch_mml2_mstep_one_step(args0=args0, prob_fun=prob_fun, entry=entry,
                        n.ik=n.ik, diffindex=diffindex, max.increment=max.increment,
                        numdiff.parm=numdiff.parm)
        args0 <- res$args0
        b <- args0[[entry]]
        max_incr_b <- res$max_increment
        se.b <- res$se
        if ( ! is.null(constraints) ){
            b[ constraints[,1] ] <- constraints[,2]
        }

        #--- update beta
        beta0 <- beta

        max.increment <- max_incr_beta
        entry <- "beta"
        diffindex <- est.b
        res <- rasch_mml2_mstep_one_step(args0=args0, prob_fun=prob_fun, entry=entry,
                        n.ik=n.ik, diffindex=diffindex, max.increment=max.increment,
                        numdiff.parm=numdiff.parm)
        args0 <- res$args0
        beta <- sirt_squeeze(args0[[entry]], lower=min.beta)
        max_incr_beta <- res$max_increment
        se.beta <- res$se

        #--- update delta
        if (est_delta){
            delta0 <- delta.miss
            max.increment <- max_incr_delta
            entry <- "delta.miss"
            diffindex <- est.delta
            res <- rasch_mml2_mstep_one_step(args0=args0, prob_fun=prob_fun, entry=entry,
                            n.ik=n.ik, diffindex=diffindex, max.increment=max.increment,
                            numdiff.parm=numdiff.parm)
            args0 <- res$args0
            delta.miss <- sirt_squeeze(args0[[entry]], lower=min.delta, upper=max.delta)
            args0[[entry]] <- delta.miss
            max_incr_delta <- res$max_increment
            se.delta <- res$se
        }
        #--- update a
        if (est_a){
            a0 <- fixed.a
            max.increment <- max_incr_a
            entry <- "fixed.a"
            diffindex <- est.a
            res <- rasch_mml2_mstep_one_step(args0=args0, prob_fun=prob_fun, entry=entry,
                            n.ik=n.ik, diffindex=diffindex, max.increment=max.increment,
                            numdiff.parm=numdiff.parm)
            args0 <- res$args0
            fixed.a <- args0[[entry]]
            max_incr_a <- res$max_increment
            se.a <- res$se

        }

        miter <- miter + 1
    }
    a_change <- max(abs(fixed.a-a0))

    #-- output
    res <- list(b=b, se.b=b, beta=beta, se.beta=se.beta,
            delta.miss=delta.miss, se.delta=se.delta, fixed.a=fixed.a,
            max_incr_a=max_incr_a, se.a=se.a, a_change=a_change)
    return(res)
}


.mstep.mml.missing1 <- rasch_mml2_mstep_missing1
