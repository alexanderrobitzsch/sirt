## File Name: xxirt_newton_raphson.R
## File Version: 0.194


xxirt_newton_raphson <- function(em_out, em_args, maxit_nr, optimizer_nr,
            control_nr, verbose=TRUE)
{
    partable <- em_out$partable
    customTheta <- em_out$customTheta

    G <- em_args$G
    Theta <- em_args$Theta
    items <- em_args$items
    item_list <- em_args$item_list
    ncat <- em_args$ncat
    partable_index <- em_args$partable_index
    dat <- em_args$dat
    dat_resp_bool <- em_args$dat_resp_bool
    resp_index <- em_args$resp_index
    group <- em_args$group
    weights <- em_args$weights

    par1 <- par_items <- xxirt_partable_extract_freeParameters( partable=partable )
    par2 <- par_Theta <- xxirt_parTheta_extract_freeParameters( customTheta=customTheta )
    par0 <- par <- c(par1, par2)

    NPI <- length(par_items)
    NPT <- length(par_Theta)
    NP <- NPI+NPT

    em_args$parindex_items <- 1:NPI
    em_args$parindex_Theta <- (NPI+1):(NPI+NPT)
    x <- par0

    # define parameter tables for estimation
    I <- ncol(dat)
    partable$item_group_comp <- 0
    for (ii in 1:I){
        ind_ii <- which(partable$itemnr==ii)
        partable$item_group_comp[ ind_ii ] <- 1:length(ind_ii)
    }
    partable$item_group_comp[ partable$parfree==0 ] <- 0
    t2 <- table(partable$parindex)
    i3 <- as.numeric( names(t2)[ t2 > 1] )
    partable$item_group_comp[ partable$parindex %in% i3 ] <- 0

    MIGC <- max(partable$item_group_comp)

    em_args$partable <- partable
    em_args$MIGC <- MIGC
    em_args$customTheta <- customTheta
    em_args$NPI <- NPI
    em_args$NPT <- NPT
    em_args$NP <- NP

    i1 <- stats::aggregate(partable$itemnr, list(partable$parindex), min )
    i2 <- stats::aggregate(partable$itemnr, list(partable$parindex), max )

    free_pars_design <- data.frame( pid=1:NPI, type='item', parlabel=names(par1) )
    free_pars_design$one_item <- i1[,2]==i2[,2]
    free_pars_design$itemnr <- ifelse(free_pars_design$one_item, i1[,2], -9 )
    partable_free <- partable[ partable$parfree==1 & partable$est, ]
    free_pars_design$item_group_comp <- partable_free$item_group_comp

    em_args$free_pars_design <- free_pars_design
    em_args$group0 <- em_args$group - 1

    test <- TRUE
    test <- FALSE
    if (test){
        ll1 <- xxirt_nr_optim_fun(x=x, em_args=em_args)
        grad1 <- xxirt_nr_grad_fun_numapprox(x=x, em_args=em_args)
        grad2 <- xxirt_nr_grad_fun_Rcpp(x=x, em_args=em_args)
        # Revalpr_round('grad1-grad2',5)
    }

    #-- optimize
    control_nr$maxit <- maxit_nr
    partable_free <- partable[ partable$parfree==1, ]
    lower <- c( partable_free$lower, customTheta$lower )
    upper <- c( partable_free$upper, customTheta$upper )
    names(upper) <- names(lower) <- names(x)

    if (verbose){
        cat( paste0('****** Newton-Raphson Optimization ********\n'))
        utils::flush.console()
    }
    res_opt_nr <- sirt_optimizer(optimizer=optimizer_nr,
                            par=x, fn=xxirt_nr_optim_fun,
                            grad=xxirt_nr_grad_fun_Rcpp, hessian=FALSE,
                            control=control_nr, em_args=em_args,
                            lower=lower, upper=upper )

    #-- collect output
    x <- res_opt_nr$par
    partable <- xxirt_partable_include_freeParameters( partable=em_args$partable,
                                        x=x[ em_args$parindex_items ] )
    customTheta <- xxirt_parTheta_include_freeParameters(
                                customTheta=em_args$customTheta,
                                x=x[ em_args$parindex_Theta ])

    #-- optimization function values
    opt_values <- xxirt_nr_optim_fun(x=x, em_args=em_args, output_all=TRUE)

    #--- output
    res <- list(res_opt_nr=res_opt_nr, partable=partable, customTheta=customTheta,
                    opt_values=opt_values)
    return(res)
}
