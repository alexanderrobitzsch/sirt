## File Name: xxirt_nr_optim_fun.R
## File Version: 0.145

xxirt_nr_optim_fun <- function(x, em_args, output_all=FALSE)
{
    pen_val <- 0

    customTheta <- em_args$customTheta

    #*** compute prior distribution
    prior_Theta <- xxirt_compute_prior_Theta_from_x(x=x, em_args=em_args)

    #* compute item response probabilities
    probs_items <- xxirt_compute_prob_item_from_x(x=x, em_args=em_args)

    #*** individual likelihood
    p.xi.aj <- xxirt_compute_likelihood( probs_items=probs_items,
                                dat=em_args$dat, dat_resp_bool=em_args$dat_resp_bool )

    #*** compute likelihood function
    ll_case <- xxirt_compute_casewise_likelihood(prior_Theta=prior_Theta,
                                group=em_args$group, p.xi.aj=p.xi.aj,
                                customTheta=customTheta)
    ll0 <- ll <- -sum( em_args$weights*log(ll_case) )

    #*** add prior distributions
    x1 <- x[ em_args$parindex_items ]
    partable1 <- xxirt_partable_include_freeParameters( partable=em_args$partable, x=x1)
    par_prior <- sum( xxirt_mstep_itemParameters_evalPrior(partable1, h=0)    )
    ll <- ll+par_prior

    #*** add penalty function
    if (!is.null(em_args$penalty_fun_item)){
        pen_val <- em_args$penalty_fun_item(x=x1)
        ll <- ll + pen_val
    }

    x2 <- x[ em_args$parindex_Theta ]
    if (!is.null(em_args$penalty_fun_theta)){
        pen2 <- em_args$penalty_fun_theta(x=x2)
        ll <- ll + pen2
        pen_val <- pen_val + pen2
    }

    res <- ll
    if (output_all){
        res <- list(opt_fun=ll, pen_val=pen_val, par_prior=par_prior, ll=ll0)
    }
    return(res)
}

