## File Name: mgsem_opt_fun.R
## File Version: 0.216


mgsem_opt_fun <- function(x, opt_fun_args, output_all=FALSE)
{

    estimator <- opt_fun_args$estimator
    partable <- mgsem_coef2partable(coef=x, partable=opt_fun_args$partable)
    model <- mgsem_partable2model(partable=partable,
                        model=opt_fun_args$model, index=FALSE)
    G <- opt_fun_args$G
    eval_pen_res <- NULL
    if (output_all){
        implied_list <- list()
        est_tot_list <- list()
        S1_list <- list()
        mean_residual_list <- list()
    }

    ll <- 0
    for (gg in seq_len(G) ){
        est0 <- model[[1]]$est
        est_gg <- model[[gg+1]]$est
        est_tot_gg <- mgsem_add_list_entries(list1=est0, add_list=est_gg, output_list=est0)
        implied <- mgsem_compute_model_implied_moments(est=est_tot_gg,
                            is_B=opt_fun_args$is_B, calc_Sigma=TRUE, calc_Mu=TRUE)
        if (output_all){
            implied_list[[gg]] <- implied
            est_tot_list[[gg]] <- est_tot_gg
        }

        #- function evaluation
        eval_args <- list(suffstat=opt_fun_args$suffstat[[gg]],
                            Mu=implied$Mu, Sigma=implied$Sigma,
                            output_all=output_all)
        if (estimator=="ME"){
            eval_args$p <- opt_fun_args$p_me
            eval_args$eps <- opt_fun_args$eps_approx
            eval_args$deriv <- FALSE
            eval_args$approx_method <- opt_fun_args$technical$approx_method
        }
        ll_gg <- do.call(what=opt_fun_args$eval_fun, args=eval_args)
        if (output_all & (estimator=="ML") ){
            S1_list[[gg]] <- ll_gg$S1
            mean_residual_list[[gg]] <- ll_gg$mean_residual
            ll_gg <- ll_gg$loglike
        }
        ll <- ll+ll_gg
    }
    ll0 <- ll

    #-- penalty function
    if (opt_fun_args$use_penalty){
        res <- mgsem_evaluate_penalties(x=x, partable=partable,
                        prior_list=opt_fun_args$prior_list,
                        technical=opt_fun_args$technical,
                        h=opt_fun_args$technical$h,
                        p=opt_fun_args$p_pen,
                        eps_approx=opt_fun_args$eps_approx,
                        deriv=FALSE, difflp_info=opt_fun_args$difflp_info,
                        loop_parms=opt_fun_args$loop_parms,
                        pen_type=opt_fun_args$pen_type, a_scad=opt_fun_args$a_scad)
        eval_pen_res <- res
        ll <- ll + res$pen_all
    }


    #-- negative function (minimization problem)
    ll <- - ll

    #- whole output
    res <- ll
    if (output_all){
        res <- list(loglike=ll0, eval_pen_res=eval_pen_res, opt_val=ll,
                pen_all=eval_pen_res$pen_all,
                implied=implied_list, est_tot=est_tot_list,
                S1=S1_list,
                suffstat=opt_fun_args$suffstat, mean_residual=mean_residual_list,
                G=G, estimator=estimator )
    }

    #-- output
    return(res)
}
