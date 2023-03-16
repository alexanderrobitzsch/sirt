## File Name: rm_sdt_mstep_rater_function_value.R
## File Version: 0.09


rm_sdt_mstep_rater_function_value <- function(x, par_index, partable_rater, I, K,
        nik_rater, eps=1E-10)
{
    K1 <- K + 1
    probs_dim <- c(I, K1, K1)
    probs_fun <- rm_sdt_calc_probs_grm_rcpp
    update_probs_args <- c('c.rater', 'd.rater')
    probs_args <- list( I=I, K=K, eps=eps, use_log=TRUE )
    post <- rm_sdt_mstep_type_function_value( x=x, par_index=par_index,
                partable=partable_rater, type='rater', probs_args=probs_args,
                probs_fun=probs_fun, probs_dim=probs_dim,
                update_probs_args=update_probs_args, nik=nik_rater )
    return(post)
}
