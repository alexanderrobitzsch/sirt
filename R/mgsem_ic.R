## File Name: mgsem_ic.R
## File Version: 0.12

mgsem_ic <- function(opt_fun_output, opt_fun_args, partable, technical)
{

    ic <- list(estimator=opt_fun_args$estimator)
    ic$deviance <- -2*opt_fun_output$loglike
    ic$loglike <- opt_fun_output$loglike
    ic$pen <- -opt_fun_output$pen_all
    ic$n <- opt_fun_args$N

    ic$np_all <- max(partable$index)

    #-- extract free and penalized parameters
    partable <- partable[ partable$unique==1, ]
    ic$np_penal <- 0
    eps <- technical$eps_count_penal
    ic$np_penal <- sum( ( partable$pen_lp > 0 ) * ( abs(partable$est) < eps ) )
    ic$np <- ic$np_all-ic$np_penal

    # compute information criteria
    ic <- xxirt_ic_compute_criteria(ic=ic, compute_np=FALSE)

    #-- output
    return(ic)
}
