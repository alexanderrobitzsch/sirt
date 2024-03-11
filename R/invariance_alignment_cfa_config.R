## File Name: invariance_alignment_cfa_config.R
## File Version: 0.261


invariance_alignment_cfa_config <- function(dat, group, weights=NULL,
    model="2PM", verbose=FALSE, ...)
{
    CALL <- match.call()
    is_data <- sirt_is_data(dat=dat)
    if (! is_data){
        mu_list <- dat[[1]]
        Sigma_list <- dat[[2]]
        N_list <- dat[[3]]
        I <- length(mu_list[[1]])
        group <- seq(1,length(mu_list))
        is_data <- FALSE
    }
    if (is_data){
        ind <- order(group)
        dat <- dat[ ind, ]
        group <- group[ind]
        I <- ncol(dat)
        if (!is.null(weights)){
            weights <- weights[ind]
        }
    }
    groups <- unique(group)
    G <- length(groups)
    items <- colnames(dat)
    N <- rep(NA, G)
    names(N) <- groups
    nu <- matrix(NA, nrow=G, ncol=I)
    rownames(nu) <- groups
    colnames(nu) <- items
    lambda <- nu
    err_var <- nu
    weights_gg <- NULL
    vcov <- matrix(0, nrow=2*I*G, ncol=2*I*G)
    names1 <- c()
    for (gg in 1L:G){
        names2 <- c( paste0(items, '_lam_Gr',gg), paste0(items, '_nu_Gr',gg) )
        names1 <- c( names1, names2)
    }
    colnames(vcov) <- rownames(vcov) <- names1

    vcov_ind0 <- 0
    for (gg in 1L:G){
        if (is_data){
            dat_gg <- dat[ group==groups[gg], ]
            dat_gg <- dat_gg[, colMeans(is.na(dat_gg)) < 1 ]
            items_gg <- colnames(dat_gg)
            ind_gg <- match(items_gg, items)
            if (!is.null(weights)){
                weights_gg <- weights[group==groups[gg]]
            }
            args <- list(dat_gg=dat_gg, weights_gg=weights_gg, model=model,
                            N=nrow(dat_gg),...)
        }
        if (!is_data){
            dat_gg <- list(mu=mu_list[[gg]], Sigma=Sigma_list[[gg]], N=N_list[[gg]])
            args <- list(dat_gg=dat_gg, weights_gg=NULL, model=model)
            args$N <- dat_gg$N
            ind_gg <- 1L:I
        }
        cat( paste0('Compute CFA for group ', gg, ' | model ', model, '\n') )
        utils::flush.console()
        res <- do.call(what='invariance_alignment_cfa_config_estimate', args=args)
        nu[gg, ind_gg] <- res$nu
        lambda[gg, ind_gg] <- res$lambda
        err_var[gg, ind_gg] <- res$err_var
        N[gg] <- res$nobs
        ind_gg1 <- vcov_ind0 + 1L:(2*I)
        vcov[ind_gg1, ind_gg1] <- res$vcov
        vcov_ind0 <- vcov_ind0+2*I
    }
    #-- output
    res <- list(nu=nu, lambda=lambda, err_var=err_var, vcov=vcov, N=N, G=G, I=I,
                    items=items, groups=groups, CALL=CALL)
    return(res)
}
