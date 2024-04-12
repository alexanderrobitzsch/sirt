## File Name: mgsem_opt_fun.R
## File Version: 0.275


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
    suffstat <- opt_fun_args$suffstat

    ll <- 0
    chisq <- 0
    N <- 0
    for (gg in seq_len(G) ){
        est0 <- model[[1]]$est
        est_gg <- model[[gg+1]]$est
        est_tot_gg <- mgsem_add_list_entries(list1=est0, add_list=est_gg,
                            output_list=est0)
        implied <- mgsem_compute_model_implied_moments(est=est_tot_gg,
                            is_B=opt_fun_args$is_B, calc_Sigma=TRUE, calc_Mu=TRUE)
        if (output_all){
            implied_list[[gg]] <- implied
            est_tot_list[[gg]] <- est_tot_gg
        }
        S_gg <- suffstat[[gg]]$S
        p <- nrow(S_gg)

        #- function evaluation
        eval_args <- list(suffstat=opt_fun_args$suffstat[[gg]],
                            Mu=implied$Mu, Sigma=implied$Sigma,
                            output_all=output_all)
        if (estimator=='ME'){
            eval_args$p <- opt_fun_args$p_me
            eval_args$eps <- opt_fun_args$eps_approx
            eval_args$deriv <- FALSE
            eval_args$approx_method <- opt_fun_args$technical$approx_method
        }
        ll_gg <- do.call(what=opt_fun_args$eval_fun, args=eval_args)
        if (output_all & (estimator=='ML') ){
            ll0_gg <- ll_gg
            S1_list[[gg]] <- ll_gg$S1
            mean_residual_list[[gg]] <- ll_gg$mean_residual
            ll_gg <- ll_gg$loglike
            N_gg <- suffstat[[gg]]$N
            S_gg <- suffstat[[gg]]$S
            p <- nrow(S_gg)
            ll_gg_adj <- ll_gg + N_gg/2*p*log(2*pi)
            chi_gg <- -2*( (N_gg/2)*( sirt_logdet(x=S_gg) + p ) + ll_gg_adj )
            chisq <- chisq + chi_gg
            N <- N + N_gg
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
        # chi square statistic and RMSEA
        p_mu <- 0
        for (gg in 1L:G){
            mu1 <- suffstat[[gg]]$M
            p_mu <- p_mu + sum(abs(mu1)>1e-14)
        }
        chisq_df <- p_mu + G*p*(p+1)/2 - max(partable$index)
        chisq_p <- 1-stats::pchisq(chisq, df=chisq_df)
        lambda <- max( chisq - chisq_df, 0 )
        fac <- chisq_df * N
        RMSEA <- sqrt( lambda / fac )
        # fit proportion
        pen_all <- abs(eval_pen_res$pen_all)
        pen_prop <- pen_all / ( chisq/2 + pen_all )

        #-- output
        res <- list(loglike=ll0, eval_pen_res=eval_pen_res, opt_val=ll,
                pen_all=pen_all, implied=implied_list, est_tot=est_tot_list,
                S1=S1_list, suffstat=opt_fun_args$suffstat,
                mean_residual=mean_residual_list, G=G, estimator=estimator,
                chisq=chisq, chisq_df=chisq_df, chisq_p=chisq_p,
                RMSEA=RMSEA, pen_prop=pen_prop)
    }

    #-- output
    return(res)
}
