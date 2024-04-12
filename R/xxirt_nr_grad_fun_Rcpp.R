## File Name: xxirt_nr_grad_fun_Rcpp.R
## File Version: 0.175

xxirt_nr_grad_fun_Rcpp <- function(x, em_args, eps=1e-100)
{
    NP <- em_args$NP
    NPI <- em_args$NPI
    NPT <- em_args$NPT
    free_pars_design <- em_args$free_pars_design
    h <- em_args$h
    eps2 <- 1e-300

    grad <- 0*x

    #*** compute prior distribution
    prior_Theta0 <- xxirt_compute_prior_Theta_from_x(x=x, em_args=em_args)

    #* compute item response probabilities
    probs_items0 <- xxirt_compute_prob_item_from_x(x=x, em_args=em_args)
    p.xi.aj0 <- xxirt_compute_likelihood( probs_items=probs_items0,
                            dat=em_args$dat, dat_resp_bool=em_args$dat_resp_bool )
    ll_case0 <- xxirt_compute_casewise_likelihood(prior_Theta=prior_Theta0,
                            group=em_args$group, p.xi.aj=p.xi.aj0)
    ll0 <- sum( em_args$weights*log(ll_case0+eps2) )

    partable <- xxirt_partable_include_freeParameters( partable=em_args$partable, x=x )

    #---
    # derivatives of probabilities
    MIGC <- em_args$MIGC
    probs_items_der <- list()
    ratio_list <- list()

    for (mm in seq_len(MIGC) ){
        free_pars_design_mm <- free_pars_design[ free_pars_design$item_group_comp==mm, ]
        x1 <- x[ em_args$parindex_items ]
        x1 <- sirt_add_increment(x=x1, pos=free_pars_design_mm$pid, value=h)
        partable1 <- xxirt_partable_include_freeParameters( partable=em_args$partable,
                                x=x1 )
        probs_temp <- xxirt_compute_itemprobs( item_list=em_args$item_list,
                                items=em_args$items, Theta=em_args$Theta,
                                ncat=em_args$ncat, partable=partable1,
                                partable_index=em_args$partable_index )
        probs_items_der[[mm]] <- (probs_temp-probs_items0)/h
        ratio_list[[mm]] <- probs_items_der[[mm]] / ( probs_items0 + eps )
    }

    #--- loop over item parameters
    for (pp in seq_len(NPI) ){

        free_pars_design_pp <- free_pars_design[pp,]
        partable_pp <- partable[ partable$parlabel==free_pars_design_pp$parlabel, ]
        item_pp <- free_pars_design_pp$itemnr
        item_group_comp_pp <- free_pars_design_pp$item_group_comp

        #- shortcut for item parameter per one item
        if (free_pars_design_pp$one_item){
            ratio <- (ratio_list[[item_group_comp_pp]])[item_pp,,]
            grad[pp] <- sirt_rcpp_xxirt_newton_raphson_derivative_par(
                                dat=em_args$dat, dat_resp_bool=em_args$dat_resp_bool,
                                ratio=ratio, p_xi_aj=p.xi.aj0, item=item_pp,
                                prior_Theta=prior_Theta0, group0=em_args$group0,
                                weights=em_args$weights, ll_case0=ll_case0, eps=eps)
        } else {
            x1 <- x[ em_args$parindex_items ]
            x1 <- sirt_add_increment(x=x1, pos=pp, value=h)
            probs_items_temp <- xxirt_compute_prob_item_from_x(x=x1, em_args=em_args)
            p.xi.aj.temp <- xxirt_compute_likelihood( probs_items=probs_items_temp,
                                dat=em_args$dat, dat_resp_bool=em_args$dat_resp_bool)
            ll_case_temp <- xxirt_compute_casewise_likelihood(prior_Theta=prior_Theta0,
                                group=em_args$group, p.xi.aj=p.xi.aj.temp)
            ll_temp <- sum( em_args$weights*log(ll_case_temp+eps2) )
            grad[pp] <- -( ll_temp - ll0 ) / h

        }

    }        # end pp

    #*** prior distributions
    if (NPI > 0){
        index_items <- 1L:NPI
        x1 <- x[ em_args$parindex_items ]
        partable1 <- xxirt_partable_include_freeParameters( partable=em_args$partable,
                            x=x1)
        par_prior <- xxirt_mstep_itemParameters_evalPrior(partable1, h=0)
        par_prior1 <- xxirt_mstep_itemParameters_evalPrior(partable1, h=h)
        par_prior_der <- ( par_prior1 - par_prior) / h
        grad[index_items] <- grad[index_items] + par_prior_der
    }

    #*** penalty function
    if (NPI>0){
        if (!is.null(em_args$penalty_fun_item)){
            pen0 <- em_args$penalty_fun_item(x=x1)
            pen1 <- rep(0,NPI)
            for (pp in 1L:NPI){
                xh <- sirt_add_increment(x=x1, pos=pp, value=h)
                pen1[pp] <- em_args$penalty_fun_item(x=xh)
            }
            grad[index_items] <- grad[index_items] + (pen1-pen0)/h
        }
    }

    #--- loop over distribution parameters
    if (NPT>0){
        for (pp in em_args$parindex_Theta){
            x1 <- sirt_add_increment(x=x, pos=pp, value=h)
            prior_Theta <- xxirt_compute_prior_Theta_from_x(x=x1, em_args=em_args)
            prior_Theta_der <- ( prior_Theta - prior_Theta0 ) / h
            ll_case <- xxirt_compute_casewise_likelihood(prior_Theta=prior_Theta_der,
                                group=em_args$group, p.xi.aj=p.xi.aj0)
            grad[pp] <- -sum( em_args$weights*ll_case/(ll_case0 + eps) )
        }
    }

    return(grad)
}
