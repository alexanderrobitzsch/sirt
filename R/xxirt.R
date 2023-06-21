## File Name: xxirt.R
## File Version: 1.024


#--- user specified item response model
xxirt <- function( dat, Theta=NULL, itemtype=NULL, customItems=NULL,
                partable=NULL, customTheta=NULL, group=NULL, weights=NULL,
                globconv=1E-6, conv=1E-4, maxit=1000, mstep_iter=4,
                mstep_reltol=1E-6, h=1E-4, use_grad=TRUE,
                verbose=TRUE, penalty_fun_item=NULL,
                np_fun_item=NULL, verbose_index=NULL,
                cv_kfold=0, cv_maxit=10)
{
    #*** preliminaries
    CALL <- match.call()
    s1 <- Sys.time()

    #*** some data processing of dat
    res <- xxirt_data_proc(dat=dat, group=group, weights=weights )
    N <- res$N
    W <- res$W
    G <- res$G
    group <- res$group
    items <- res$items
    group0 <- res$group0
    groups_unique <- res$groups_unique
    I <- res$I
    maxK <- res$maxK
    ncat <- res$ncat
    weights <- res$weights
    group_index <- res$group_index
    dat_resp <- res$dat_resp
    resp_index <- res$resp_index
    dat1 <- res$dat1
    dat_resp_bool <- res$dat_resp_bool

    #*** default Theta
    if ( is.null(Theta) ){
        Theta <- matrix( seq(-6,6,length=21), ncol=1 )
    }
    Theta <- as.matrix(Theta)
    TP <- nrow(Theta)

    #*** eps - handle numerical instabilities
    eps <- 1E-8

    # create partable if not provided
    if ( is.null(partable) ){
        partable <- xxirt_createParTable( dat=dat, itemtype=itemtype,
                            customItems=customItems )
    }

    # process partable and itemtype
    res <- xxirt_proc_ParTable( itemtype=itemtype, partable=partable, items=items )
    itemtype <- res$itemtype
    partable <- res$partable
    partable_index <- res$partable_index
    ncat <- res$ncat
    maxK <- res$maxK
    mstep_method <- res$mstep_method
    item_index <- res$item_index
    dat <- as.matrix(dat)

    # create item list
    item_list <- xxirt_createItemList( customItems=customItems, itemtype=itemtype,
                        items=items, partable=partable )

    # shortcut for calculating expected counts
    dat1_resp <- xxirt_prepare_response_data(G=G, group_index=group_index,
                        weights=weights, dat1=dat1, dat_resp=dat_resp, maxK=maxK )

    #*** starting values item parameters
    par0 <- xxirt_partable_extract_freeParameters( partable=partable )
    par1 <- xxirt_ThetaDistribution_extract_freeParameters( customTheta=customTheta )

    #*** verbose
    verbose1 <- verbose==1
    verbose2 <- verbose==2

    disp <- '...........................................................\n'
    iter <- 1
    dev <- 1E100
    converged <- FALSE

    #--- create list with arguments for EM algorithm
    em_args <- list( maxit=maxit, verbose1=verbose1, disp=disp,
                item_list=item_list, items=items, Theta=Theta, ncat=ncat,
                partable=partable, partable_index=partable_index, dat=dat,
                resp_index=resp_index, dat_resp=dat_resp, dat_resp_bool=dat_resp_bool,
                dat1=dat1, dat1_resp=dat1_resp, customTheta=customTheta, G=G, par0=par0,
                maxK=maxK, group_index=group_index, weights=weights,
                mstep_iter=mstep_iter, eps=eps, mstep_reltol=mstep_reltol,
                mstep_method=mstep_method, item_index=item_index, h=h, use_grad=use_grad,
                penalty_fun_item=penalty_fun_item, group=group, par1=par1,
                globconv=globconv, conv=conv, verbose2=verbose2, verbose3=FALSE,
                verbose_index=verbose_index)
                
    #--- run EM algorithm
    res <- do.call(what=xxirt_em_algorithm, args=em_args)
    #--- collect EM output
    iter <- res$iter
    converged <- res$converged
    probs_items <- res$probs_items
    p.xi.aj <- res$p.xi.aj
    prior_Theta <- res$prior_Theta
    n.ik <- res$n.ik
    p.aj.xi <- res$p.aj.xi
    N.ik <- res$N.ik
    N.k <- res$N.k
    post_unnorm <- res$post_unnorm
    ll1 <- res$ll1
    partable <- res$partable
    par0 <- res$par0
    pen_val <- res$pen_val
    ll2 <- res$ll2
    customTheta <- res$customTheta
    par1 <- res$par1
    ll_case <- res$ll_case
    dev <- res$dev
    dev00 <- res$dev00

    #*** optional cross-validation
    cv_loglike <- NA
    if (cv_kfold>0){
        N <- nrow(dat)
        cv_zone <- ( ( 0:(N-1) ) %% cv_kfold ) + 1
        cv_loglike <- 0
        em_args$partable <- partable
        em_args$customTheta <- customTheta
        for (jj in 1:cv_kfold){
            weights1 <- weights*(cv_zone!=jj )
            weights2 <- weights*(cv_zone==jj )
            em_args$weights <- weights1
            em_args$maxit <- cv_maxit
            em_args$verbose1 <- FALSE
            em_args$verbose2 <- FALSE
            em_args$verbose3 <- TRUE
            cat(paste0('\nCross-validation sample ', jj ), ' |' )
            utils::flush.console()
            res <- do.call(what=xxirt_em_algorithm, args=em_args)
            # compute cv log likelihood
            ll_jj <- sum( weights2 * res$ll_case )
            cv_loglike <- cv_loglike + ll_jj
        }
        cat('\n')
    }

    #**** post processing
    opt_val <- dev
    dev <- dev00

    #-- collect parameters
    res <- xxirt_postproc_parameters( partable=partable, customTheta=customTheta,
                    items=items, probs_items=probs_items, np_fun_item=np_fun_item )
    par_items <- res$par_items
    par_Theta <- res$par_Theta
    probs_items <- res$probs_items
    par_items_summary <- res$par_items_summary
    par_items_bounds <- res$par_items_bounds
    np_item <- res$np_item

    #-- information criteria
    ic <- xxirt_ic( dev=dev, N=W, par_items=par_items,
                par_Theta=par_Theta, I=I, par_items_bounds=par_items_bounds,
                np_item=np_item)

    #-- compute EAP
    EAP <- xxirt_EAP(p.aj.xi=p.aj.xi, Theta=Theta )

    #--- output
    s2 <- Sys.time()
    res <- list( partable=partable, par_items=par_items,
                par_items_summary=par_items_summary, par_items_bounds=par_items_bounds,
                par_Theta=par_Theta, Theta=Theta, probs_items=probs_items,
                probs_Theta=prior_Theta, deviance=dev, loglike=-dev/2,
                opt_val=opt_val, pen_val=pen_val, ic=ic,
                item_list=item_list, customItems=customItems,
                customTheta=customTheta, cv_loglike=cv_loglike,
                p.xi.aj=p.xi.aj, p.aj.xi=p.aj.xi, ll_case=ll_case,
                n.ik=n.ik, EAP=EAP, dat=dat, dat_resp=dat_resp,
                weights=weights, item_index=item_index, G=G, group=group,
                group_orig=group0, ncat=ncat, mstepItem_method=mstep_method,
                partable_index=partable_index, items=items, dat1=dat1,
                dat1_resp=dat1_resp, dat_resp_bool=dat_resp_bool,
                group_index=group_index, maxK=maxK,
                resp_index=resp_index, converged=converged, iter=iter-1,
                CALL=CALL, s1=s1, s2=s2 )
    class(res) <- 'xxirt'
    return(res)
}
