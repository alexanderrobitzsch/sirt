## File Name: mgsem_evaluate_penalties.R
## File Version: 0.332


mgsem_evaluate_penalties <- function(x, partable, prior_list, technical,
        h, p, eps_approx, deriv=FALSE, difflp_info=NULL, loop_parms,
        pen_type="lasso", a_scad=3.7)
{
    ND <- nrow(partable)
    if (is.null(loop_parms)){
        loop_parms <- (1:ND)[ partable$unique==1]
    }
    NP <- max(partable$index)
    if (!deriv){
        NP <- 1
    }
    pen_prior <- pen_l2 <- pen_lp <- pen_difflp <- rep(0,NP)

    res <- list( x=x, pen_prior=pen_prior, pen_l2=pen_l2,
                pen_lp=pen_lp, pen_difflp=pen_difflp, loop_parms=loop_parms,
                h=h, p=p, eps_approx=eps_approx)

    #- loop over parameters
    for (dd in loop_parms){
        index <- partable[dd,'index']
        res <- mgsem_evaluate_penalties_evaluate_entry(
                    x=x, res=res, dd=dd, index=index,
                    partable=partable, technical=technical,
                    h=h, p=p, eps_approx=eps_approx,
                    prior_list=prior_list, deriv=deriv,
                    difflp_info=difflp_info)
    } # end dd

    partable2 <- partable[ loop_parms, ]
    n <- partable2$N_group

    #*** L2 penalty
    if (technical$is_pen_l2){
        lambda <- partable2$pen_l2
        if (deriv){
            pen_l2 <- 2*lambda*n*x
        } else {
            pen_l2 <- sum(lambda*n*x^2)
        }
        res$pen_l2 <- pen_l2
    }

    use_rcpp_penalty <- technical$use_rcpp_penalty

    #*** Lp penalty
    if (technical$is_pen_lp){
        fac <- partable2$pen_lp
        args_pen <- list(x=x, p=p, n=n, fac=fac, eps=eps_approx, deriv=deriv,
                        pen_type=pen_type, a=a_scad, h=h)

        if (!use_rcpp_penalty){
            fun_pen <- 'mgsem_eval_lp_penalty_vector'
        } else {
            fun_pen <- 'sirt_rcpp_mgsem_eval_lp_penalty'
        }
        val <- do.call(what=fun_pen, args=args_pen)
        if (!deriv){
            val <- sum(val)
        }
        res$pen_lp <- val
    }

    #*** diffLp penalty
    if (technical$is_pen_difflp){
        indices <- difflp_info$difflp_indices
        x1 <- x[indices]
        fac <- difflp_info$lpdiff_facmat
        fac_logical <- difflp_info$lpdiff_facmat_logical
        n <- difflp_info$lpdiff_n

        if (!deriv){
            #--- no derivative
            args_pen <- list(x=x1, fac=fac, p=p, n=n, h=h, eps_approx=eps_approx,
                            a_scad=a_scad, pen_type=pen_type)
            if (! use_rcpp_penalty){
                fun_pen <- 'mgsem_eval_lp_penalty_matrix'
            } else {
                args_pen$h <- NULL
                args_pen$fac_logical <- fac_logical
                fun_pen <- 'sirt_rcpp_mgsem_eval_lpdiff_penalty'
            }
            z <- do.call( what=fun_pen, args=args_pen)

        } else {
            #--- derivative
            z <- rep(0, length(x))
            args_pen <- list(fac=fac, p=p, eps_approx=eps_approx, h=h,
                            a_scad=a_scad, pen_type=pen_type, n=n)
            if (! use_rcpp_penalty){
                args_pen$par <- x1
                args_pen$FUN <- mgsem_eval_lp_penalty_matrix
                args_pen$gradient <- TRUE
                args_pen$hessian <- FALSE
                fun_pen <- CDM::numerical_Hessian
            } else {
                args_pen$x <- x1
                args_pen$fac_logical <- fac_logical
                fun_pen <- 'sirt_rcpp_mgsem_eval_lpdiff_penalty_deriv'
            }
            val <- do.call(what=fun_pen, args=args_pen)
            z[indices] <- val
        }
        res$pen_difflp <- z
    }

    #*** penalty for diffpar
    res$pen_diffpar_lp <- 0
    if (technical$is_diffpar_pen){
        diffpar_pen_list_entries <- technical$diffpar_pen$diffpar_pen_list_entries
        NDP <- nrow(diffpar_pen_list_entries)
        p <- technical$diffpar_pen$p
        n <- partable2$N_group
        # vector of differences of parameter
        z <- x[ diffpar_pen_list_entries$index1 ] - x[ diffpar_pen_list_entries$index2 ]
        n2 <- sqrt(n[ diffpar_pen_list_entries$index1 ] *
                            n[ diffpar_pen_list_entries$index2 ] )
        args_pen <- list(x=z, fac=diffpar_pen_list_entries$W, n=n2, p=p,
                        eps_approx=eps_approx, pen_type=pen_type, h=h, deriv=deriv)
        fun_pen <- 'mgsem_eval_lp_penalty_vector'
        val <- do.call(what=fun_pen, args=args_pen)

        #* no derivative
        if (!deriv){
            val <- sum(val)
        }

        #* derivative
        if (deriv){
            der_z <- val
            NP <- length(x)
            val <- rep(0,NP)
            for (hh in 1:NDP){
                i1 <- diffpar_pen_list_entries$index1[hh]
                val[i1] <- val[i1] + der_z[hh]
                i2 <- diffpar_pen_list_entries$index2[hh]
                val[i2] <- val[i2] - der_z[hh]
            }
        }
        res$pen_diffpar_lp <- val
    }

    #*** sum all penalties
    res$pen_all <- res$pen_prior - res$pen_l2 - res$pen_lp - res$pen_difflp -
                            res$pen_diffpar_lp

    #--- output
    return(res)
}
