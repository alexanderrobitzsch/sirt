## File Name: xxirt_hessian_compute_loglike.R
## File Version: 0.03
## File Last Change: 2022-10-23


xxirt_hessian_compute_loglike <- function(p.xi.aj, prior1, weights)
{
    post_unnorm <- prior1 * p.xi.aj
    # compute log-likelihood
    dev <- sum( weights * log( rowSums( post_unnorm ) ) )
    return(dev)
}
