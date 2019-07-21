## File Name: invariance_alignment_cfa_config.R
## File Version: 0.224


invariance_alignment_cfa_config <- function(dat, group, weights=NULL,
    verbose=FALSE, ...)
{
    CALL <- match.call()
    groups <- unique(group)
    G <- length(groups)
    I <- ncol(dat)
    items <- colnames(dat)
    N <- rep(NA, G)
    names(N) <- groups
    nu <- matrix(NA, nrow=G, ncol=I)
    rownames(nu) <- groups
    colnames(nu) <- items
    lambda <- nu
    err_var <- nu
    weights_gg <- NULL
    for (gg in 1:G){
        dat_gg <- dat[ group==groups[gg], ]
        dat_gg <- dat_gg[, colMeans(is.na(dat_gg)) < 1 ]
        items_gg <- colnames(dat_gg)
        ind_gg <- match(items_gg, items)
        if (!is.null(weights)){
            weights_gg <- weights[group==groups[gg]]
        }
        cat( paste0("Compute CFA for group ", gg, "\n") )
        utils::flush.console()
        res <- invariance_alignment_cfa_config_estimate(dat_gg=dat_gg,
                        weights_gg=weights_gg, ...)
        nu[gg, ind_gg] <- res$nu
        lambda[gg, ind_gg] <- res$lambda
        err_var[gg, ind_gg] <- res$err_var
        N[gg] <- res$nobs
    }
    #-- output
    res <- list(nu=nu, lambda=lambda, err_var=err_var, N=N, G=G, I=I,
                    items=items, groups=groups, CALL=CALL)
    return(res)
}
