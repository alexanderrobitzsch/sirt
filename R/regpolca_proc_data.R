## File Name: regpolca_proc_data.R
## File Version: 0.03
## File Last Change: 2020-02-24

regpolca_proc_data <- function(dat, group)
{
    ncats <- apply(dat, 2, max, na.rm=TRUE)+1
    lca_dich <- max(ncats)==2
    I <- ncol(dat)
    N <- nrow(dat)
    if (is.null(group)){
        group <- rep(1, N)
    }
    groups <- unique(sort(group))
    group <- match(group, groups)
    G <- length(groups)
    Ni <- colSums(1-is.na(dat))
    #- output
    res <- list(ncats=ncats, lca_dich=lca_dich, I=I, N=N, group=group,
                groups=groups, G=G, Ni=Ni)
    return(res)
}
