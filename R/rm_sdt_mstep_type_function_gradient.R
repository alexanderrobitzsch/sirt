## File Name: rm_sdt_mstep_type_function_gradient.R
## File Version: 0.11



rm_sdt_mstep_type_function_gradient <- function(x, par_index, partable, type,
    pargroup_type, probs_args, probs_fun, probs_dim, update_probs_args, nik,
    numdiff.parm)
{
    res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable,
                                parm0=x, type=type )
    partable1 <- res$partable
    probs_nr <- as.list(1:2)
    grad_out <- rep(NA, pargroup_type$np)
    numdiff_fac <- c(1,-1)
    K <- probs_dim[2] - 1
    for (pp in seq_len(pargroup_type$max_pargroup) ){
        type_pp <- paste( pargroup_type$pargroup_type[[pp]] )
        #- x0 + h and x0 - h
        for ( vv in 1:2){
            probs_args1 <- probs_args
            parm1 <- rm_sdt_extract_par_from_partable_add_increment(partable=partable1,
                            pargroup=pp, increment=numdiff.parm*numdiff_fac[vv] )
            res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable,
                                        parm0=parm1, type=type )
            probs_args1 <- rm_sdt_mstep_include_probs_args(probs_args=probs_args1,
                    parm_list=res$parm_list, update_probs_args=update_probs_args )
            probs_nr[[vv]] <- do.call( what=probs_fun, args=probs_args1)
        }
        logprobs <- ( probs_nr[[1]] - probs_nr[[2]] ) / ( 2 * numdiff.parm)
        ll_grad <- rm_sdt_calc_gradient_likelihood_item_llgrad2( logprob_D1=logprobs,
                            nik.item=nik, diffindex=pargroup_type$pargroup_diffindex[[pp]], K=K,
                            prob_D1_dim=probs_dim)
        grad_out[ pargroup_type$pargroup_index[[pp]] ] <- ll_grad
    }
    grad_ll <- - grad_out
    grad_prior <- rm_sdt_evaluate_prior_derivative(partable=partable1, h=numdiff.parm)
    grad_post <- grad_ll + grad_prior
    grad_post <- ifelse(is.na(grad_post), 0, grad_post )
    return(grad_post)
}
