## File Name: mgsem_evaluate_penalties_evaluate_entry.R
## File Version: 0.188
## File Last Change: 2022-05-16


mgsem_evaluate_penalties_evaluate_entry <- function(x, res, dd, index, partable,
            technical, h, p, eps_approx, prior_list, deriv=FALSE, difflp_info=NULL)
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

    # --- output
    return(res)
}
