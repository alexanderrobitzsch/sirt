## File Name: mgsem_loglike_data.R
## File Version: 0.03


mgsem_loglike_data <- function(dat, Mu, Sigma)
{
    requireNamespace('mvtnorm')
    ll <- sum(mvtnorm::dmvnorm(dat, mean=Mu, sigma=Sigma, log=TRUE))
    return(ll)
}
