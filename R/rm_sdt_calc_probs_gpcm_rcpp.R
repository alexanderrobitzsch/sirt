## File Name: rm_sdt_calc_probs_gpcm_rcpp.R
## File Version: 0.12
## File Last Change: 2018-12-30

rm_sdt_calc_probs_gpcm_rcpp <- function(a.item, tau.item, Qmatrix, theta.k, VV, K, TP,
    eps=0, use_log=FALSE, as_vector=FALSE)
{
    K1 <- K+1
    prob_dim <- c(VV, K1, TP)
    res <- sirt_rcpp_rm_sdt_calc_probs_gpcm( a=a.item, tau=tau.item,
                theta_k=theta.k, VV=VV, K1=K1, TP=TP, eps=eps, use_log=use_log )
    if ( ! as_vector ){
        res <- array(res, dim=prob_dim)
    }
    return(res)
}
