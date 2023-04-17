## File Name: xxirt_mstep_itemParameters.R
## File Version: 0.385


#--- M-step item parameters
xxirt_mstep_itemParameters <- function( partable, item_list, items, Theta,
            ncat, partable_index, N.ik, mstep_iter, par0, eps,
            mstep_reltol, mstep_method, item_index, h, use_grad,
            penalty_fun_item=NULL)
{
    #-------------------------------------------------
    #**** define likelihood function
    like_items1 <- function( x, ... ){
        partable <- xxirt_partable_include_freeParameters( partable=partable, x=x )
        probs1 <- xxirt_compute_itemprobs( item_list=item_list, items=items,
                            Theta=Theta, ncat=ncat, partable=partable,
                            partable_index=partable_index )
        log_probs1 <- log( probs1 + eps )
        # log likelihood value
        ll <- - sum( N.ik * log_probs1 )
        # add prior distributions
        pen <- xxirt_mstep_itemParameters_evalPrior(partable=partable, h=0)
        # add penalty
        ll <- ll + sum(pen)
        # add penalty function
        if (!is.null(penalty_fun_item)){
            ll <- ll + penalty_fun_item(x=x)
        }
        return(2*ll)
    }
    #**** end definition likelihood
    #-------------------------------------------------

    #-------------------------------------------------
    #*** define gradient here
    #### x <- partable[ partable$parfree==1, 'value']
    grad_items1 <- function(x, ...){
        partable0 <- xxirt_partable_include_freeParameters( partable=partable, x=x )
        probs1 <- xxirt_compute_itemprobs( item_list=item_list, items=items,
                            Theta=Theta, ncat=ncat, partable=partable0,
                            partable_index=partable_index )
        NP <- sum( partable$parfree==1)
        pen1 <- grad1 <- rep(0,NP)
        pen0 <- 0
        for (pp in 1:NP){
            item_index_pp <- item_index[[pp]]
            xh <- x
            xh[pp] <- xh[pp] + h
            N.ik_pp <- N.ik[item_index_pp,,, drop=FALSE ]
            partable_h <- xxirt_partable_include_freeParameters( partable=partable0,
                                    x=xh )
            probs1h <- xxirt_compute_itemprobs( item_list=item_list, items=items,
                                    Theta=Theta, ncat=ncat, partable=partable_h,
                                    partable_index=partable_index,
                                    item_index=item_index_pp)
            ll0 <- - sum( N.ik_pp * log( probs1[ item_index_pp,,,drop=FALSE] + eps ) )
            ll1 <- - sum( N.ik_pp * log( probs1h + eps ) )
            grad1[pp] <- ( ll1 - ll0 )
        }
        #--- prior distributions
        pen <- xxirt_mstep_itemParameters_evalPrior(partable=partable0, h=0)
        pen_h <- xxirt_mstep_itemParameters_evalPrior(partable=partable0, h=h)
        #-- penalty function
        if (!is.null(penalty_fun_item)){
            pen0 <- penalty_fun_item(x=x)
            for (pp in 1:NP){
                xh <- x
                xh[pp] <- xh[pp] + h
                pen1[pp] <- penalty_fun_item(x=xh)
            }
        }
        grad1 <- ( grad1 +  pen_h - pen + pen1 - pen0 )/h
        return(grad1)
    }
    #**** end definition gradient
    #-------------------------------------------------

    opt_fun <- stats::optim

    like_items <- function(x)
    {
        like_items1(x=x)
    }
    grad_items <- function(x)
    {
        grad_items1(x=x)
    }

    # mstep_method <- 'nloptr'
    use_nloptr <- mstep_method=='nloptr'
    if (use_nloptr){
        requireNamespace('nloptr')
        opt_fun <- nloptr::nloptr
    }

    if (use_nloptr){
        arg_control <- list(maxeval=mstep_iter)
    } else {
        arg_control <- list(maxit=mstep_iter)
    }
    if (! use_nloptr){
        if ( mstep_method=='BFGS'){
            arg_control$reltol <- mstep_reltol
        }
        if ( mstep_method=='L-BFGS-B'){
            arg_control$factr <- mstep_reltol
        }
    } else {
        arg_control <- list(ftol_rel=mstep_reltol, algorithm='NLOPT_LD_LBFGS')
    }
    # method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
    if (use_nloptr){
        arg_list <- list(x0=par0, eval_f=like_items, opts=arg_control)
    } else {
        arg_list <- list( par=par0, fn=like_items, method=mstep_method,
                        control=arg_control)
    }
    if (use_grad){
        if (use_nloptr){
            arg_list$eval_grad_f <- grad_items
        } else {
            arg_list$gr <- grad_items
        }
    }
    lb <- partable[ partable$parfree==1, 'lower']
    ub <- partable[ partable$parfree==1, 'upper']
    if ( mstep_method=='L-BFGS-B'){
        arg_list$lower <- lb
        arg_list$upper <- ub
    }
    if (use_nloptr){
        arg_list$lb <- lb
        arg_list$ub <- ub
    }


    mod <- do.call( what=opt_fun, args=arg_list )

    partable <- xxirt_partable_include_freeParameters( partable=partable, x=mod$par )
    pen_val <- 0
    if (!is.null(penalty_fun_item)){
        pen_val <- penalty_fun_item(x=mod$par)
    }

    #-- output
    res <- list( ll1=mod$value, partable=partable, par0=mod$par, pen_val=pen_val )
    return(res)
}
