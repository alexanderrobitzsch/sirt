## File Name: mgsem_grad_fun.R
## File Version: 0.169


mgsem_grad_fun <- function(x, opt_fun_args, output_all=FALSE)
{

    estimator <- opt_fun_args$estimator
    NP <- opt_fun_args$NP
    ND <- opt_fun_args$ND
    G <- opt_fun_args$G
    grad_contrib <- matrix(0, nrow=ND, ncol=G)
    grad <- rep(0,NP)

    partable <- mgsem_coef2partable(coef=x, partable=opt_fun_args$partable)
    model <- mgsem_partable2model(partable=partable,
                            model=opt_fun_args$model, index=FALSE)

    #- group-wise calculations
    res <- mgsem_list_elements_est_total_implied(model=model, is_B=opt_fun_args$is_B)
    implied0 <- res$implied0
    est_total0 <- res$est_total0

    dermoments0 <- list()
    for (gg in 1:G){
        grad_suffstat_fun_args <- list(suffstat=opt_fun_args$suffstat[[gg]],
                                        Mu=implied0[[gg]])
        if (estimator=="ME"){
            grad_suffstat_fun_args$p <- opt_fun_args$p
            grad_suffstat_fun_args$eps <- opt_fun_args$eps_approx
            grad_suffstat_fun_args$deriv <- TRUE
            grad_suffstat_fun_args$only_deriv <- TRUE
            grad_suffstat_fun_args$approx_method <- opt_fun_args$technical$approx_method
        }
        dermoments0[[gg]] <- do.call(what=opt_fun_args$grad_suffstat_fun,
                                        args=grad_suffstat_fun_args)
    }

    for (dd in 1:ND){

        group_dd <- partable$group[dd]
        i1_dd <- partable$i1[dd]
        i2_dd <- partable$i2[dd]
        type_dd <- paste(partable$type[dd])
        index_dd <- partable$index[dd]
        recycle_dd <- partable$recycle[dd]
        recycling <- recycle_dd > 0

        gg <- 1
        sel_groups <- seq_len(G)
        if (group_dd>0){
            sel_groups <- group_dd
        }

        for (gg in sel_groups){
            if (! recycling){
                implied <- implied0[[gg]]
                dermoments <- dermoments0[[gg]]
                grad_param_fun_args <- list(est=est_total0[[gg]],
                            dermoments=dermoments, suffstat=opt_fun_args$suffstat[[gg]],
                            type=type_dd, i1=i1_dd, i2=i2_dd, h=opt_fun_args$technical$h,
                            is_B=opt_fun_args$is_B,    eps=opt_fun_args$technical$eps_zero,
                            num_approx=opt_fun_args$technical$num_approx )
                gr1 <- do.call(what=opt_fun_args$grad_param_fun,
                                args=grad_param_fun_args)
            } else {
                gr1 <- grad_contrib[ recycle_dd, gg]
            }
            grad_contrib[dd,gg] <- gr1
            grad[index_dd] <- grad[index_dd] + gr1
        }
    }
    grad0 <- grad

    #-- penalty function
    if (opt_fun_args$use_penalty){
        res <- mgsem_evaluate_penalties(x=x, partable=partable,
                            prior_list=opt_fun_args$prior_list,
                            technical=opt_fun_args$technical,
                            h=opt_fun_args$technical$h, p=opt_fun_args$p_pen,
                            eps_approx=opt_fun_args$eps_approx,
                            deriv=TRUE, difflp_info=opt_fun_args$difflp_info,
                            loop_parms=opt_fun_args$loop_parms,
                            pen_type=opt_fun_args$pen_type, a_scad=opt_fun_args$a_scad)
        pen_grad_output <- res
        grad <- grad + res$pen_all
    }


    #-- negative function (minimization problem)
    grad <- -grad

    #-- output_list
    res <- grad
    if (output_all){
        res <- list(grad_fun=grad, grad=grad0, pen=pen_grad_output$pen_all,
                    dermoments=dermoments0, implied=implied0,
                    pen_grad_output=pen_grad_output, estimator=estimator,
                    N=opt_fun_args$N, N_group=opt_fun_args$N_group)
    }

    #-- output
    return(res)

}
