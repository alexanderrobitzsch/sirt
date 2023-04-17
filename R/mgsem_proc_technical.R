## File Name: mgsem_proc_technical.R
## File Version: 0.180

mgsem_proc_technical <- function(technical, control, p_me, p_pen, eps_approx, suffstat,
            estimator, diffpar_pen=NULL, cd_control=NULL)
{
    #*** defaults technical; update list
    technical0 <- list(maxiter=1000, h=1e-5, eps_zero=1e-12, optimizer='nlminb',
                        use_deriv=TRUE,    num_approx=FALSE,
                        eps_count_penal=1e-2,
                        approx_method='lp', use_rcpp_penalty=TRUE)
    use_approx_method <- ! is.null(technical$approx_method)
    technical <- mgsem_update_list_entries(add_list=technical, output_list=technical0)

    #*** update list cd_control
    cd_control0 <- list(maxiter=20, tol=5*1e-4, interval_length=0.05,
                            method="exact")
    cd_control <- mgsem_update_list_entries(add_list=cd_control, output_list=cd_control0)

    if ((p_pen==2) & (p_me==2) & ( ! use_approx_method) ){
        technical$approx_method    <- 'l2'
    }

    if (technical$approx_method=='l2' & ( (p_pen!=2) | ( p_me!=2) ) ){
        technical$approx_method    <- 'lp'
    }
    # if (technical$approx_method=='l2'){
    #     eps_approx <- 0
    # }

    #- optimizer
    if (technical$optimizer=='nlminb'){
        control$iter.max <- technical$maxiter
    }
    if (technical$optimizer=='optim'){
        control$maxit <- technical$maxiter
    }
    if (technical$optimizer=='bobyqa'){
        fac <- max(15, 3*nrow(suffstat[[1]][[2]]))
        technical$maxiter <- technical$maxiter*fac
        control$maxfun <- technical$maxiter
    }
    if (technical$optimizer=='nloptr'){
        control$maxeval <- technical$maxiter
    }

    p <- p_pen
    if (estimator=='ME'){
        p <- p_me
    }

    #--- output
    res <- list(control=control, technical=technical, eps_approx=eps_approx, p=p,
                    cd_control=cd_control)
    return(res)
}
