## File Name: rm_sdt_calc_gradient_likelihood_item_llgrad.R
## File Version: 0.06
## File Last Change: 2018-12-30

rm_sdt_calc_gradient_likelihood_item_llgrad <- function(logprob_D1, nik.item, diffindex, K)
{
    ll <- rowSums( logprob_D1[,1,] * nik.item[,1,] )
    for (kk in 2:(K+1) ){
        ll <- ll + rowSums( logprob_D1[,kk,] * nik.item[,kk,] )
    }
    ll_grad <- ll
    ll_grad <- rowsum(ll_grad, diffindex )
    ll_grad <- ll_grad[ rownames(ll_grad) > 0, 1 ]
    return(ll_grad)
}
