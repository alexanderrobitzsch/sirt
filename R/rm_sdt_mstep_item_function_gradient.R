## File Name: rm_sdt_mstep_item_function_gradient.R
## File Version: 0.16
## File Last Change: 2019-07-21


rm_sdt_mstep_item_function_gradient <- function(x, par_index, partable_item, Qmatrix, theta.k,
    VV, K, TP, pargroup_item, nik.item, numdiff.parm, eps=1E-10, conv_array=TRUE )
{
    res <- rm_sdt_fill_par_to_partable( par_index=par_index, partable=partable_item,
                parm0=x, type="item" )
    partable1 <- res$partable
    tau.item <- res$parm_list$tau.item
    a.item <- res$parm_list$a.item
    args <- list(a.item=a.item, tau.item=tau.item, Qmatrix=Qmatrix, theta.k=theta.k, VV=VV, K=K,
                TP=TP, eps=eps, use_log=FALSE )
    prob0 <- do.call( what=rm_sdt_calc_probs_gpcm_rcpp, args=args )
    args$use_log <- TRUE
    prob_D1_dim <- c(VV, K+1, TP)
    prob_D1 <- array(0, dim=prob_D1_dim )

    #*** loop over pargroups
    grad_out <- rep(NA, pargroup_item$np)
    for (pp in seq_len(pargroup_item$max_pargroup) ){
        pt_pp <- partable1[ partable1$pargroup==pp, ][1,]
        type_pp <- pt_pp$type
        #-- tau parameter
        if ( type_pp=="tau"){
            cat_pp <- pt_pp$col
            prob_D1 <- sirt_rcpp_rm_sdt_calc_gradient_item_deriv_tau( prob0=as.vector(prob0),
                                    prob_D1_dim=prob_D1_dim, cat_pp=cat_pp)
            prob_D1 <- array(prob_D1, dim=prob_D1_dim)
        }
        #-- a parameter
        if (type_pp=="a" ){
            prob_D1 <- sirt_rcpp_rm_sdt_calc_gradient_item_deriv_a( prob0=as.vector(prob0),
                            prob_D1_dim=prob_D1_dim, theta_k=theta.k )
            prob_D1 <- array(prob_D1, dim=prob_D1_dim)
        }
        #*** compute gradient
        ll_grad <- rm_sdt_calc_gradient_likelihood_item_llgrad2( logprob_D1=prob_D1,
                        nik.item=nik.item, diffindex=pargroup_item$pargroup_diffindex[[pp]], K=K,
                        prob_D1_dim=prob_D1_dim)
        grad_out[ pargroup_item$pargroup_index[[pp]] ] <- ll_grad
    }
    grad_ll <- - grad_out
    grad_prior <- rm_sdt_evaluate_prior_derivative(partable=partable1, h=numdiff.parm)
    grad_post <- grad_ll + grad_prior
    return(grad_post)
}
