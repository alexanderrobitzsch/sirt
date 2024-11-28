## File Name: mcmc_Rhat.R
## File Version: 0.121

mcmc_Rhat <- function( mcmc_object, n_splits=3 )
{
    n_samples <- nrow(mcmc_object)
    n_pars <- ncol(mcmc_object)
    n_within <- floor( n_samples / n_splits )
    rhat_vec <- rep(NA, n_pars)
    names(rhat_vec) <- colnames(mcmc_object)
    for (pp in 1L:n_pars){
        matr <- matrix( NA, nrow=n_within, ncol=n_splits)
        for (ss in 1L:n_splits){
            matr[,ss] <- mcmc_object[ (ss-1)* n_within + 1L:n_within, pp ]
        }
        rhat_vec[pp] <- Rhat1(matr)
    }
    return(rhat_vec)
}
