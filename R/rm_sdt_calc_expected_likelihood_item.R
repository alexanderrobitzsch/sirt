## File Name: rm_sdt_calc_expected_likelihood_item.R
## File Version: 0.16
## File Last Change: 2018-12-30


rm_sdt_calc_expected_likelihood_item <- function( nik.item, a.item, tau.item, Qmatrix,
        theta.k, VV, K, TP, eps=1E-10)
{
    prob_dim <- c(VV, K+1, TP)
    nik_item <- as.vector(nik.item)
    probs <- rm_sdt_calc_probs_gpcm_rcpp( a.item=a.item, tau.item=tau.item,
                Qmatrix=Qmatrix, theta.k=theta.k, VV=VV, K=K, TP=TP, eps=eps, use_log=TRUE,
                as_vector=TRUE)
    ll <- sirt_rcpp_rm_sdt_calc_gradient_likelihood_item_llgrad( logprob_D1=probs,
                    prob_D1_dim=prob_dim, nik_item=nik_item )
    return(ll)
}
