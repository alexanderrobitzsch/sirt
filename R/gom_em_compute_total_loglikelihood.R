## File Name: gom_em_compute_total_loglikelihood.R
## File Version: 0.02

gom_em_compute_total_loglikelihood <- function(f.yi.qk, pi.k, weights)
{
    N <- nrow(f.yi.qk)
    ll <- sum( weights*log( rowSums( f.yi.qk * sirt_matrix2( x=pi.k, nrow=N ) ) ) )
    return(ll)
}
