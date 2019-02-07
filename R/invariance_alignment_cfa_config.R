## File Name: invariance_alignment_cfa_config.R
## File Version: 0.12


invariance_alignment_cfa_config <- function(dat, group, ...)
{
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
    for (gg in 1:G){
        dat_gg <- dat[ group==groups[gg], ]
        dat_gg <- dat_gg[, colMeans(is.na(dat_gg)) < 1 ]
        items_gg <- colnames(dat_gg)
        ind_gg <- match(items_gg, items)
        res <- invariance_alignment_cfa_config_estimate(dat_gg=dat_gg, ...)
        nu[gg, ind_gg] <- res$nu
        lambda[gg, ind_gg] <- res$lambda
        err_var[gg, ind_gg] <- res$err_var
        N[gg] <- res$nobs
    }
    #-- output
    res <- list(nu=nu, lambda=lambda, err_var=err_var, N=N, G=G, I=I,
                    items=items, groups=groups)
    return(res)
}
