## File Name: linking_haberman_itempars_convert.R
## File Version: 0.092


linking_haberman_itempars_convert <- function(itempars=NULL, lambda=NULL,
    nu=NULL, a=NULL, b=NULL)
{
    if (!is.null(nu)){
        b <- -nu/lambda
        a <- lambda
    }
    if (!is.null(b)){
        if (is.null(a)){
            a <- 1+0*b
        }
        itempars <- linking_haberman_itempars_prepare(b=t(b), a=t(a))
    }

    if (! is.null(itempars)){
        itempars_list <- linking_proc_itempars(itempars=itempars)
    }

    #- reconvert into a and b
    a <- linking_haberman_itempars_convert_process_matrices(
                mat=itempars_list$aM, est_pars=itempars_list$est_pars)
    b <- linking_haberman_itempars_convert_process_matrices(
                mat=itempars_list$bM, est_pars=itempars_list$est_pars)

    #- reconvert into lambda and mu
    lambda <- t(a)
    nu <- -t(b/a)

    #-- output
    res <- list(itempars=itempars, lambda=lambda, nu=nu, a=a, b=b,
                    itempars_list=itempars_list)
    return(res)
}
