## File Name: xxirt_em_algorithm.R
## File Version: 0.070

xxirt_em_algorithm <- function(maxit, verbose1, verbose2, verbose3, disp, item_list,
            items, Theta, ncat, partable, partable_index, dat, resp_index,
            dat_resp, dat_resp_bool, dat1, dat1_resp, group, customTheta, G, par0, maxK,
            group_index, weights, mstep_iter, eps, mstep_reltol, mstep_method,
            item_index, h, use_grad, penalty_fun_item, par1, globconv, conv,
            verbose_index)
{
    iter <- 1
    dev <- 1E100
    converged <- FALSE
    while (
            ( iter < ( maxit + 1 ) ) & ( ! converged )
                ){

        if (verbose1){
            cat(disp)
            cat('Iteration', iter, '   ', paste( Sys.time() ), '\n' )
        }
        dev0 <- dev

        #*** item probabilities
        probs_items <- xxirt_compute_itemprobs( item_list=item_list,
                            items=items, Theta=Theta, ncat=ncat,
                            partable=partable, partable_index=partable_index )

        #*** compute individual likelihood
        p.xi.aj <- xxirt_compute_likelihood( probs_items=probs_items, dat=dat,
                             resp_index=resp_index, dat_resp_bool=dat_resp_bool )

        #*** compute prior distribution
        prior_Theta <- xxirt_compute_priorDistribution( Theta=Theta,
                              customTheta=customTheta, G=G )

        #*** compute posterior distribution and expected counts
        res <- xxirt_compute_posterior( prior_Theta=prior_Theta, p.xi.aj=p.xi.aj,
                        group=group, G=G, weights=weights, dat1=dat1,
                        dat_resp=dat_resp, maxK=maxK, group_index=group_index,
                        dat1_resp=dat1_resp )
        n.ik <- res$n.ik
        p.aj.xi <- res$p.aj.xi
        N.ik <- res$N.ik
        N.k <- res$N.k
        post_unnorm <- res$post_unnorm

        #*** M-step item parameters
        par00 <- par0
        res <- xxirt_mstep_itemParameters( partable=partable, item_list=item_list,
                        items=items, Theta=Theta, ncat=ncat,
                        partable_index=partable_index, N.ik=N.ik,
                        mstep_iter=mstep_iter, par0=par0, eps=eps,
                        mstep_reltol=mstep_reltol, mstep_method=mstep_method,
                        item_index=item_index, h=h, use_grad=use_grad,
                        penalty_fun_item=penalty_fun_item)
        ll1 <- res$ll1
        partable <- res$partable
        par0 <- res$par0
        pen_val <- res$pen_val
        prior_par <- res$prior_par

        #*** M-step theta distribution
        par10 <- par1
        res <- xxirt_mstep_ThetaParameters( customTheta=customTheta, G=G, eps=eps,
                        mstep_iter=mstep_iter, N.k=N.k, par1=par1,
                        mstep_reltol=mstep_reltol, Theta=Theta )
        ll2 <- res$ll2
        customTheta <- res$customTheta
        par1 <- res$par1

        #*** compute deviance
        ll_case <- log( rowSums( post_unnorm ) )
        dev <- - 2 * sum( weights * ll_case )
        dev00 <- dev
        dev <- dev + pen_val
        opt_fun_value <- dev00/2 + pen_val + prior_par
        globconv_temp <- abs( ( - dev + dev0 ) / dev0 )

        conv0 <- 0
        if ( length(par0) > 0){
            conv0 <- max( abs(par0-par00))
        }
        conv1 <- 0
        if ( !is.null(par1) ){
            conv1 <- max( abs(par10-par1))
        }
        conv_temp <- max( conv0, conv1)
        converged <- ( globconv_temp < globconv ) & ( conv_temp < conv )

        #-- print progress
        res <- xxirt_print_progress( opt_fun_value=opt_fun_value, dev=dev, dev0=dev0,
                    dev00=dev00, pen_val=pen_val, conv0=conv0, conv1=conv1, iter=iter,
                    verbose1=verbose1, verbose2=verbose2, verbose_index=verbose_index,
                    verbose3=verbose3, prior_par=prior_par)
        iter <- iter + 1

    }
    #---- output
    res <- list(iter=iter, converged=converged, probs_items=probs_items,
                    p.xi.aj=p.xi.aj, prior_Theta=prior_Theta, n.ik=n.ik,
                    p.aj.xi=p.aj.xi, N.ik=N.ik, N.k=N.k,
                    post_unnorm=post_unnorm, ll1=ll1,
                    partable=partable, par0=par0, pen_val=pen_val,
                    prior_par=prior_par, ll2=ll2, customTheta=customTheta, par1=par1,
                    ll_case=ll_case, dev=dev, dev00=dev00)
    return(res)
}
