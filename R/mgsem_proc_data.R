## File Name: mgsem_proc_data.R
## File Version: 0.04

mgsem_proc_data <- function(data, group, weights)
{
    groups <- unique(group)
    G <- length(groups)
    N <- nrow(data)
    if (is.null(weights)){
        weights <- rep(1,N)
    }
    idgroup <- match(group, groups)
    suffstat <- list()
    for (gg in 1:G){
        ind_gg <- which(group==groups[gg])
        dat_gg <- data[ ind_gg, ]
        w <- weights[ind_gg]
        res <- stats::cov.wt( x=dat_gg, w=w, method="ML")
        s_gg <- list(M=res$center, S=res$cov, N=sum(w) )
        suffstat[[gg]] <- s_gg
    }
    #-- output
    res <- list(suffstat=suffstat, weights=weights, groups=groups,
                    G=G, N=N, idgroup=idgroup, data=data)
    return(res)
}
