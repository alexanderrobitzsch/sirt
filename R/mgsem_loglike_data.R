## File Name: mgsem_loglike_data.R
## File Version: 0.02
## File Last Change: 2022-02-25


mgsem_loglike_data <- function(dat, Mu, Sigma)
{
    requireNamespace("mvtnorm")
    ll <- sum(mvtnorm::dmvnorm(dat, mean=Mu, sigma=Sigma, log=TRUE))
    return(ll)
}
