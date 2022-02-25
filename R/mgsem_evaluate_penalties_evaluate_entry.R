## File Name: mgsem_evaluate_penalties_evaluate_entry.R
## File Version: 0.186


mgsem_evaluate_penalties_evaluate_entry <- function(x, res, dd, index, partable, technical,
    h, p, eps_approx, prior_list, deriv=FALSE, difflp_info=NULL)
{
    x0 <- x
    x <- x0[index]

    if (deriv){
        ind_add <- index
        null_vec <- rep(0, length(x))
    } else {
        ind_add <- 1
        null_vec <- 0
    }

    # res[["pen_l2"]] <- null_vec

    # fun_eval="none" produces a value of zero

    #*** priors
    if (technical$is_prior){
        pen_entry <- "pen_prior"
        prior <- partable[dd,"prior"]
        if (prior!="none"){
            fun_eval <- prior_list[[prior]]
            args_eval <- list()
            val <- mgsem_evaluate_penalties_evaluate_entry_fun_eval(x=x,
                        args_eval=args_eval, fun_eval=fun_eval, h=h, deriv=deriv)
            res[[pen_entry]][ind_add] <- res[[pen_entry]][ind_add] + val
        }
    }  # end is_prior

        ##
        ##      #*** l2 and lp
        ##      # if (technical$is_pen_l2){
        ##      if (FALSE){
        ##          pen_entry <- "pen_l2"
        ##          fun_eval <- "mgsem_power_fun_differentiable_approx"
        ##          approx_method <- "l2"
        ##          # args_eval <- list(approx_method=approx_method, deriv=FALSE,
        ##          #                   p=2, eps=0)
        ##          # val <- mgsem_evaluate_penalties_evaluate_entry_fun_eval(x=x,
        ##          #           args_eval=args_eval, fun_eval=fun_eval, h=h, deriv=deriv)
        ##          val <- mgsem_power_fun_differentiable_approx(x=x, p=2,
        ##                              eps=0, deriv=deriv, approx_method="l2")
        ##          val <- partable[dd,pen_entry]*val
        ##          res[[pen_entry]][ind_add] <- res[[pen_entry]][ind_add] + val
        ##      }  # end is_l2
        ##
        ##      # if (technical$is_pen_lp){
        ##      if (FALSE){
        ##          pen_entry <- "pen_lp"
        ##          fun_eval <- "mgsem_power_fun_differentiable_approx"
        ##          approx_method <- "lp"
        ##          if (partable[dd,pen_entry]>0){
        ##              # args_eval <- list(approx_method=approx_method, deriv=FALSE,
        ##              #                   p=p, eps=eps_approx)
        ##              # val <- mgsem_evaluate_penalties_evaluate_entry_fun_eval(x=x,
        ##              #           args_eval=args_eval, fun_eval=fun_eval, h=h, deriv=deriv)
        ##              val <- mgsem_power_fun_differentiable_approx(x=x, p=p,
        ##                              eps=eps_approx, deriv=deriv, approx_method="lp")
        ##              val <- partable[dd,pen_entry]*val
        ##              res[[pen_entry]][ind_add] <- res[[pen_entry]][ind_add] + val
        ##          }
        ##      }  # end is_lp
        ##
        ##      #if (technical$is_pen_difflp){
        ##      if (FALSE){
        ##          pen_entry <- "pen_difflp"
        ##          fun_eval <- "mgsem_power_fun_differentiable_approx"
        ##          approx_method <- "lp"
        ##          if (partable[dd,pen_entry]>0){
        ##              # ind <- which( partable$type==partable$type[dd] & partable$i1==partable$i1[dd] &
        ##              #   partable$i2==partable$i2[dd] & (partable$group > 0) & (partable$rid!=dd) )
        ##              # partable2 <- partable[ ind,, drop=FALSE ]
        ##              # y <- x0[ partable2$index ]
        ##              y <- x0[ difflp_info$lpdiff_diff_indices[[dd]] ]
        ##              x2 <- x-y
        ##              # args_eval <- list(approx_method=approx_method, deriv=FALSE,
        ##              #                   p=p, eps=eps_approx)
        ##              # val <- mgsem_evaluate_penalties_evaluate_entry_fun_eval(x=x2,
        ##              #           args_eval=args_eval, fun_eval=fun_eval, h=h, deriv=deriv)
        ##              val1 <- mgsem_power_fun_differentiable_approx(x=x2, p=p,
        ##                              eps=eps_approx, deriv=deriv, approx_method="lp")
        ##              # val <- val1*sqrt(partable[dd,pen_entry]*partable2[,pen_entry])
        ##              val <- val1*difflp_info$pen_difflp_factor[[dd]]
        ##              res[[pen_entry]][ind_add] <- res[[pen_entry]][ind_add] + sum(val)
        ##          }
        ##      }  # end is_lp

    # --- output
    return(res)
}
