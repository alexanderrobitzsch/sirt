## File Name: rm_sdt_calc_gradient_likelihood_item_llgrad2.R
## File Version: 0.16

rm_sdt_calc_gradient_likelihood_item_llgrad2 <- function(logprob_D1, nik.item, diffindex, K,
    prob_D1_dim )
{
    nik_item <- as.vector(nik.item)
    logprob_D1 <- as.vector(logprob_D1)
    ll_grad <- sirt_rcpp_rm_sdt_calc_gradient_likelihood_item_llgrad( logprob_D1=logprob_D1,
                    prob_D1_dim=prob_D1_dim, nik_item=nik_item )
    ll_grad <- rowsum(ll_grad, diffindex )
    ll_grad <- ll_grad[ rownames(ll_grad) > 0, 1 ]
    return(ll_grad)
}
