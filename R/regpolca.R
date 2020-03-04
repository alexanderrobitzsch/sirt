## File Name: regpolca.R
## File Version: 0.116


#- Regularized polytomous latent class analysis
regpolca <- function(dat, nclasses, weights=NULL, group=NULL,
    regular_type="scad", regular_lam=0, sd_noise_init=1,
    par_item_init=NULL,    par_classprobs_init=NULL, random_starts=1, random_iter=20,
    conv=1e-5, h=1e-4, mstep_iter=10, maxit=1000, verbose=TRUE, par_item_max=10,
    set_equal=.01)
{
    #*** preliminaries
    CALL <- match.call()
    s1 <- Sys.time()

    #** analyze response patterns
    res <- regpolca_proc_data(dat=dat, group=group)
    ncats <- res$ncats
    lca_dich <- res$lca_dich
    I <- res$I
    N <- res$N
    group <- res$group
    groups <- res$groups
    G <- res$G
    Ni <- res$Ni

    #- define theta class distribution
    K <- nclasses
    par_Theta <- xxirt_classprobs_lca_init_par(K=K, G=G, random_sd=0)
    customTheta  <- xxirt_createThetaDistribution( par=par_Theta,
                            est=rep(TRUE,G*(K-1)), P=xxirt_classprobs_lca)
    Theta <- diag(K)

    #- define item response functions
    res <- regpolca_define_customItems( ncats=ncats, K=K, dat=dat,
                    par_item_max=par_item_max )
    customItems <- res$customItems
    partable <- res$partable
    itemtype <- res$itemtype

    #-- include penalty function
    penalty_fun_item <- NULL
    if (regular_lam>0){
        combis <- t( utils::combn(K,2) )
        if (regular_type=="scad"){ penalty_used <- penalty_D1_scad }
        if (regular_type=="mcp"){ penalty_used <- penalty_D1_mcp }
        if (regular_type=="lasso"){ penalty_used <- penalty_D1_lasso }
        penalty_fun_item <- function(x, ...){
            pen <- 0
            #* fused probabilities among classes
            for (ii in 1:I){
                x_ii <- x[ partable$itemnr==ii ]
                diff_ii <- x_ii[ combis[,1] ] - x_ii[ combis[,2] ]
                a1 <- penalty_used( x=diff_ii, lambda=regular_lam, eps=eps )
                pen <- pen + Ni[ii]*sum(a1)
            }
            return(pen)
        }
    }

    #-- create argument list for xxirt
    args <- list( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta=customTheta, maxit=random_iter, mstep_iter=mstep_iter,
                penalty_fun_item=penalty_fun_item, h=h, use_grad=TRUE, verbose=2 )

    #-- random starts if required
    args <- regpolca_run_xxirt_random_starts( args=args, random_starts=random_starts,
                    sd_noise_init=sd_noise_init )

    #-- arguments for final xxirt model
    args$verbose <- TRUE
    args$maxit <- maxit

    #-- run xxirt in a final model
    res <- do.call(what=xxirt, args=args)
    res$iter <- res$iter + random_iter*(random_starts>0)

    #- process output
    res$probs_Theta <- regpolca_postproc_prob_Theta(probs_Theta=res$probs_Theta)
    item <- regpolca_postproc_irf(probs_items=res$probs_items, dat=dat,
                    lca_dich=lca_dich)
    res0 <- regpolca_postproc_count_regularized_parameters(item=item,
                    set_equal=set_equal, lca_dich=lca_dich, probs_items=res$probs_items)
    item1_index <- res0$item1_index
    n_reg <- res0$n_reg
    item <- res0$item
    res$probs_items <- res0$probs_items

    #- adapt information criteria
    res$ic <- regpolca_postproc_ic(ic=res$ic, n_reg=n_reg)

    #-- arrange output
    res$CALL <- CALL
    res2 <- list(s1=s1, s2=Sys.time(), lca_dich=lca_dich, nclasses=nclasses,
                    item=item, regular_lam=regular_lam, regular_type=regular_type,
                    item1_index=item1_index, n_reg=n_reg)
    res <- sirt_add_list_elements(res=res, res2=res2)

    class(res) <- "regpolca"
    return(res)
}
