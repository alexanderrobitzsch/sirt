## File Name: mgsem_proc_technical.R
## File Version: 0.166

mgsem_proc_technical <- function(technical, control, p_me, p_pen, eps_approx, suffstat)
{
    #*** defaults technical; update list
    technical0 <- list(maxiter=1000, h=1e-5, eps_zero=1e-12, optimizer="nlminb",
                        use_deriv=TRUE,    num_approx=FALSE,
                        eps_count_penal=1e-2,
                        approx_method="lp", use_rcpp_penalty=TRUE)
    use_approx_method <- ! is.null(technical$approx_method)
    technical <- mgsem_update_list_entries(add_list=technical, output_list=technical0)

    if ((p_pen==2) & (p_me==2) & ( ! use_approx_method) ){
        technical$approx_method    <- "l2"
    }

    if (technical$approx_method=="l2" & ( (p_pen!=2) | ( p_me!=2) ) ){
        technical$approx_method    <- "lp"
    }
    # if (technical$approx_method=="l2"){
    #     eps_approx <- 0
    # }

    #- optimizer
    if (technical$optimizer=="nlminb"){
        control$iter.max <- technical$maxiter
    }
    if (technical$optimizer=="optim"){
        control$maxit <- technical$maxiter
    }
    if (technical$optimizer=="bobyqa"){
        fac <- max(15, 3*nrow(suffstat[[1]][[2]]))
        technical$maxiter <- technical$maxiter*fac
        control$maxfun <- technical$maxiter
    }
    if (technical$optimizer=="nloptr"){
        control$maxeval <- technical$maxiter
    }

    #--- output
    res <- list(control=control, technical=technical, eps_approx=eps_approx)
    return(res)
}
