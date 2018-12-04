## File Name: invariance_alignment_find_parameter_constraints.R
## File Version: 0.24

invariance_alignment_find_parameter_constraints <- function(parm,
    parm_tol, wgt=NULL, maxiter=10, conv=1E-4)
{
    parm_est <- parm
    G <- nrow(parm)
    I <- ncol(parm)
    parm_names <- colnames(parm)
    parm_joint <- sirt_colMedian(x=parm)
    if (is.null(wgt)){
        wgt <- 1+0*parm
    }
    iter <- 1
    iterate <- TRUE

    while(iterate){
        parm_est_old <- parm_est
        parm_est <- TAM::tam_matrix2( parm_joint, nrow=G)
        ind <- abs(parm-parm_est) > parm_tol
        ind0 <- ! ind
        parm_est[ind] <- parm[ind]
        parm_joint <- colSums(parm_est*ind0*wgt) / colSums(ind0*wgt)
        change <- max(abs(parm_est-parm_est_old))
        if ((change < conv) | (iter >=maxiter)){
            iterate <- FALSE
        }
        iter <- iter + 1
    }
    parm_est <- sirt_matrix_names(x=parm_est, extract_names=parm)

    #--- count number of estimated parameters
    pj <- TAM::tam_matrix2(x=parm_joint, nrow=G)
    N_unique_parm_items <- colSums( parm_est !=pj )
    names(N_unique_parm_items) <- parm_names
    N_parm_items <- ( colSums( parm_est==pj ) > 0 ) + N_unique_parm_items
    names(N_parm_items) <- parm_names
    N_parm_all <- sum(N_parm_items)
    N_unique_parm_groups <- rowSums( parm_est !=pj )
    names(N_unique_parm_groups) <- rownames(parm)
    N_total <- G*I
    prop_noninvariance <- 100*sum(N_unique_parm_groups)/(G*I)

    #--- output
    res <- list( parm_est=parm_est, parm_joint=parm_joint, parm=parm, parm_tol=parm_tol,
                    N_parm_items=N_parm_items, N_parm_all=N_parm_all,
                    N_unique_parm_items=N_unique_parm_items, N_unique_parm_groups=N_unique_parm_groups,
                    prop_noninvariance=prop_noninvariance, iter=iter, maxiter=maxiter,
                    G=G, I=I, wgt=wgt, N_total=N_total)
    return(res)
}
