## File Name: rm_sdt_mstep_item_function_value.R
## File Version: 0.06


rm_sdt_mstep_item_function_value <- function(x, par_index, partable_item, nik.item,
    Qmatrix, theta.k, VV, K, TP, eps=1E-10)
{
    res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable_item,
                            parm0=x, type="item" )
    partable1 <- res$partable
    tau.item <- res$parm_list$tau.item
    a.item <- res$parm_list$a.item
    #-- compute log-likelihood
    ll <- rm_sdt_calc_expected_likelihood_item( nik.item=nik.item,
                    a.item=a.item, tau.item=tau.item, Qmatrix=Qmatrix, theta.k=theta.k,
                    VV=VV, K=K, TP=TP, eps=eps )
    ll <- -sum(ll)
    # add prior
    prior <- rm_sdt_evaluate_prior(partable=partable1)
    post <- ll + sum(prior)
    return(post)
}
