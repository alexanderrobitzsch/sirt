## File Name: rm_sdt_pem_acceleration.R
## File Version: 0.12

rm_sdt_pem_acceleration <- function( iter, pem_parameter_index, pem_parameter_sequence,
        c.rater, Qmatrix, tau.item, VV, K, I, TP, a.item, d.rater, item.index, rater.index, theta.k, RR,
        dat2, dat2.resp, pi.k, dat2.ind.resp, center_log_a, ll, mu, sigma, pem_pars, a_center_type,
        PEM_itermax )
{
    PEM <- TRUE
    #-- collect all parameters in a list
    parmlist <- list( mu=mu, tau.item=tau.item, a.item=a.item, c.rater=c.rater, d.rater=d.rater, sigma=sigma)
    #-- transform into a vector
    pem_parm <- sirt_pem_collect_parameters( parmlist=parmlist, pem_parameter_index=pem_parameter_index )
    #-- collect parameters in initial iterations
    pem_parameter_sequence <- sirt_pem_parameter_sequence_initial_iterations( pem_parm=pem_parm,
                                            pem_parameter_sequence=pem_parameter_sequence, iter=iter )

    if ( ( iter %% 2==0 ) & ( iter > 0 ) & ( iter < PEM_itermax ) ){
        pem_parameter_sequence$P2 <- pem_parm

        #-- extract parameters
        ll_args <- list( c.rater=c.rater, Qmatrix=Qmatrix, tau.item=tau.item, VV=VV, K=K, I=I, TP=TP, a.item=a.item,
                            d.rater=d.rater, item.index=item.index, rater.index=rater.index, theta.k=theta.k, RR=RR,
                            dat2=dat2, dat2.resp=dat2.resp, pi.k=pi.k, dat2.ind.resp=dat2.ind.resp, center_log_a=center_log_a,
                            a_center_type=a_center_type)

        #** baseline likelihood
        ll_args <- sirt_pem_include_ll_args( ll_args=ll_args, pem_parm=pem_parm, pem_pars=pem_pars,
                            pem_parameter_index=pem_parameter_index )
        res <- do.call( what=rm_sdt_calc_loglikelihood, args=ll_args )
        ll0 <- ll <- res$ll

        P0 <- pem_parameter_sequence$P0
        P1 <- pem_parameter_sequence$P1
        P2 <- pem_parameter_sequence$P2
        iterate <- TRUE
        ii <- 0
        #--- begin PEM iterations
        while (iterate){
            ll_args0 <- ll_args
            res0 <- res
            ll0 <- ll
            tt <- sirt_pem_algorithm_compute_t( i=ii )
            Pnew <- sirt_pem_algorithm_compute_Pnew( tt=tt, P0=P0, P1=P1, P2=P2 )
            ll_args <- sirt_pem_include_ll_args( ll_args=ll_args, pem_parm=Pnew, pem_pars=pem_pars,
                            pem_parameter_index=pem_parameter_index )
            res <- do.call( what=rm_sdt_calc_loglikelihood, args=ll_args )
            ll <- res$ll
            if (is.na(ll)){
                ll <- -Inf
                iterate <- FALSE
            }
            if ( ll < ll0 ){
                iterate <- FALSE
            }
            ii <- ii + 1
        }
        c.rater <- ll_args0$c.rater
        d.rater <- ll_args0$d.rater
        tau.item <- ll_args0$tau.item
        mu <- ll_args0$mu
        sigma <- ll_args0$sigma
        pi.k <- res0$pi.k
        a.item <- res0$a.item
        ll <- res0$ll
        pem_parameter_sequence$P0 <- P1
        pem_parameter_sequence$P1 <- P2
    }
    if (iter > PEM_itermax){
        PEM <- FALSE
    }
    #--- output
    res <- list(ll=ll, pem_parameter_sequence=pem_parameter_sequence, c.rater=c.rater, d.rater=d.rater, a.item=a.item,
                    tau.item=tau.item, pi.k=pi.k, mu=mu, sigma=sigma, PEM=PEM)
    return(res)
}

