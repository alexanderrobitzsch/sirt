## File Name: rm_sdt_mstep_type_function_value.R
## File Version: 0.16


rm_sdt_mstep_type_function_value <- function(x, par_index, partable, type,
    probs_args, probs_fun, probs_dim, update_probs_args, nik)
{
    res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable,
                                parm0=x, type=type )
    partable1 <- res$partable
    parm_list <- res$parm_list
    probs_args <- rm_sdt_mstep_include_probs_args(probs_args=probs_args,
                    parm_list=parm_list, update_probs_args=update_probs_args )
    probs <- do.call( what=probs_fun, args=probs_args )
    ll <- sirt_rcpp_rm_sdt_calc_gradient_likelihood_item_llgrad( logprob_D1=probs,
                        prob_D1_dim=probs_dim, nik_item=nik )
    ll <- -sum(ll)
    prior <- rm_sdt_evaluate_prior(partable=partable1)
    post <- ll + sum(prior)
    return(post)
}
