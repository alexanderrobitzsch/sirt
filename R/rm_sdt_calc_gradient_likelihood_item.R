## File Name: rm_sdt_calc_gradient_likelihood_item.R
## File Version: 0.06

rm_sdt_calc_gradient_likelihood_item <- function(logprob1, logprob2,
    numdiff.parm, nik.item, diffindex)
{
    K <- dim(nik.item)[2] - 1
    logprob_D1 <- ( logprob1 - logprob2 ) / (2*numdiff.parm)
    ll <- rowSums( logprob_D1[,1,] * nik.item[,1,] )
    for (kk in 2:(K+1) ){
        ll <- ll + rowSums( logprob_D1[,kk,] * nik.item[,kk,] )
    }
    ll_grad <- ll
    ll_grad <- rowsum(ll_grad, diffindex )
    ll_grad <- ll_grad[ rownames(ll_grad) > 0, 1 ]
    return(ll_grad)
}
