## File Name: rasch_mml2_raschtype_mstep_parameter_group_evaluate_prior.R
## File Version: 0.06
## File Last Change: 2019-10-27


rasch_mml2_raschtype_mstep_parameter_group_evaluate_prior <- function(parm,
    h, prior, ll0, ll1, ll2, prior_fct )
{
    if (!is.null(prior)){
        prior_dens_args <- list(parm, prior[1], prior[2], log=TRUE)
        m0 <- do.call(prior_fct, args=prior_dens_args)
        prior_dens_args[[1]] <- parm + h
        m1 <- do.call(prior_fct, args=prior_dens_args)
        prior_dens_args[[1]] <- parm - h
        m2 <- do.call(prior_fct, args=prior_dens_args)
        ll0 <- ll0 + m0
        ll1 <- ll1 + m1
        ll2 <- ll2 + m2
    }
    res <- list(ll0=ll0, ll1=ll1, ll2=ll2)
    return(res)
}
