## File Name: rm_facets_pp_mle_calc_ll_theta.R
## File Version: 0.11
## File Last Change: 2019-01-02



rm_facets_pp_mle_calc_ll_theta <- function( data, a, b, theta )
{
    N <- length(theta)
    I <- ncol(data)
    ll0 <- rep(0,N)
    for (ii in 1L:I){
        probs.ii <- rm_facets_pp_mle_calc_pcm( theta=theta, a=a, b=b, ii=ii )
        res <- rm_facets_pp_mle_calc_ll( probs=probs.ii, data=data, ii=ii )
        ll0 <- ll0 + res
    }
    return(ll0)
}


