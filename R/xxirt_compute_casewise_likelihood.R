## File Name: xxirt_compute_casewise_likelihood.R
## File Version: 0.01

xxirt_compute_casewise_likelihood <- function(prior_Theta, group, p.xi.aj)
{
    prior1 <- t( prior_Theta[, group ] )
    p.aj.xi <- prior1 * p.xi.aj
    ll_case <- rowSums(p.aj.xi)
    return(ll_case)
}
