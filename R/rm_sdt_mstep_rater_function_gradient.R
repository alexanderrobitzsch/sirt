## File Name: rm_sdt_mstep_rater_function_gradient.R
## File Version: 0.09


rm_sdt_mstep_rater_function_gradient <- function(x, par_index, partable_rater,
        pargroup_rater, I, K, numdiff.parm, nik.rater, eps=1E-10)
{
    probs_fun <- rm_sdt_calc_probs_grm_rcpp
    K1 <- K+1
    probs_dim <- c(I, K1, K1)
    probs_args <- list( I=I, K=K, eps=eps, use_log=TRUE )
    update_probs_args <- c('c.rater', 'd.rater')
    grad_post <- rm_sdt_mstep_type_function_gradient( x=x,
            par_index=par_index, partable=partable_rater, type='rater',
            pargroup_type=pargroup_rater, probs_args=probs_args, probs_fun=probs_fun,
            probs_dim=probs_dim, update_probs_args=update_probs_args, nik=nik.rater,
            numdiff.parm=numdiff.parm )
    return(grad_post)
}
