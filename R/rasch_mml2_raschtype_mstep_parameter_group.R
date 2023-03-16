## File Name: rasch_mml2_raschtype_mstep_parameter_group.R
## File Version: 0.160
## File Last Change: 2019-10-27


#--- estimation of c parameter
rasch_mml2_raschtype_mstep_parameter_group <- function( theta.k, b, fixed.a,
            fixed.c, fixed.d, pjk, alpha1, alpha2, h, G, I, r.jk, n.jk, est_val,
            min_val, max_val, iter, old_increment, Qmatrix, parameter, prior=NULL)
{
    numdiff.parm <- h
    old_increment.c <- old_increment
    min.c <- min_val
    max.c <- max_val
    est.c <- est_val

    h1 <- h / 2
    h2 <- 1 + 2*h
    probs_args <- list( theta.k=theta.k, b=b, fixed.a=fixed.a, fixed.c=fixed.c,
                    fixed.d=fixed.d, alpha1=alpha1, alpha2=alpha2, Qmatrix=Qmatrix,
                    h=0, h1=h1, h2=h2, incr=parameter )
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

    # first order derivative
    # f(x+h) - f(x-h)=2* f'(x) * h
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
    ll0 <- rowSums(ll0)
    ll1 <- rowSums(ll1)
    ll2 <- rowSums(ll2)
    # aggregate with respect to estimation of a
    a1 <- stats::aggregate( cbind( ll0, ll1, ll2 ), list(est.c), sum, na.rm=TRUE)
    a1 <- a1[ a1[,1] > 0,]
    ll0 <- a1[,2]
    ll1 <- a1[,3]
    ll2 <- a1[,4]

    #** evaluate prior distribution
    prior_fct <- NULL
    parm <- NULL
    if (parameter=="a"){
        parm <- fixed.a
        prior_fct <- stats::dnorm
    }
    if (parameter=="c"){
        parm <- fixed.c
        parm[ parm < h] <- 2*h
        prior_fct <- stats::dbeta
    }
    if (parameter=="d"){
        parm <- fixed.d
        parm[ parm > 1-h] <- 1-2*h
        prior_fct <- stats::dbeta
    }
    res <- rasch_mml2_raschtype_mstep_parameter_group_evaluate_prior(parm=parm,
                h=h, prior=prior, ll0=ll0, ll1=ll1, ll2=ll2, prior_fct=prior_fct )
    ll0 <- res$ll0
    ll1 <- res$ll1
    ll2 <- res$ll2

    #- derivatives
    res <- rasch_mml2_difference_quotient(ll0=ll0, ll1=ll1, ll2=ll2, h=h)
    d1 <- res$d1
    d2 <- res$d2
    # change in item difficulty
    parm_change <- - d1 / d2
    parm_change <- sirt_trim_increment(increment=parm_change,
                            max_increment=max(abs(old_increment.c)))
    parm_change <- parm_change[ match( est.c, a1[,1] ) ]

    if ( any(est.c==0) ){
        parm_change[est.c==0] <- 0
    }
    if (parameter=="a"){ parm <- fixed.a }
    if (parameter=="b"){ parm <- b }
    if (parameter=="c"){ parm <- fixed.c }
    if (parameter=="d"){ parm <- fixed.d }

    parm <- parm + parm_change
    if (parameter=="c"){
        parm[ parm < 2*numdiff.parm ] <- 2*numdiff.parm
    }
    if (parameter=="d"){
        parm[ parm > 1 - 2 * numdiff.parm ] <- 1 - 2*numdiff.parm
    }
    parm[ parm > max.c ] <- max.c
    parm[ parm < min.c ] <- min.c
    #--- output
    res <- list("parm"=parm, "se"=sqrt( 1 /abs(d2) ) )
    return(res)
}
