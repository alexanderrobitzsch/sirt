## File Name: rm_facets_pp_mle_calc_pcm.R
## File Version: 0.13


#--- calculate item response probabilities
rm_facets_pp_mle_calc_pcm <- function( theta, a, b, ii )
{
    K <- ncol(b)
    N <- length(theta)
    matrK <- sirt_matrix2( x=0:K, nrow=N)
    eta <- a[ii] * theta * matrK - sirt_matrix2( x=c(0,b[ii,]), nrow=N)
    eta <- exp(eta)
    probs <- eta / rowSums(eta, na.rm=TRUE)
    return(probs)
}
