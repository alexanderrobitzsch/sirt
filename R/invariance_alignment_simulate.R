## File Name: invariance_alignment_simulate.R
## File Version: 0.05

invariance_alignment_simulate <- function(nu, lambda, err_var, mu, sigma, N)
{
    N_tot <- sum(N)
    G <- nrow(nu)
    I <- ncol(nu)
    items <- colnames(nu)
    if (is.null(items)){
        items <- paste0("I",1:I)
    }
    n_end <- cumsum(N)
    n_start <- c(1,n_end+1)[-c(G+1)]
    #* simulate data
    group <- rep(1:G, N)
    dat <- matrix(NA, nrow=N_tot, ncol=I+1)
    colnames(dat) <- c("group",items)
    dat <- as.data.frame(dat)
    dat$group <- group
    for (gg in 1:G){
        N_gg <- N[gg]
        ind_gg <- seq(n_start[gg], n_end[gg])
        fac <- stats::rnorm(N_gg, mean=mu[gg], sd=sigma[gg])
        for (ii in 1:I){
            err_ii <- stats::rnorm(N_gg, mean=0, sd=sqrt(err_var[gg,ii]) )
            dat[ind_gg, ii+1] <- nu[gg,ii] + lambda[gg,ii]*fac + err_ii
        }
    }
    #--- output
    return(dat)
}
