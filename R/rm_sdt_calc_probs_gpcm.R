## File Name: rm_sdt_calc_probs_gpcm.R
## File Version: 0.07

rm_sdt_calc_probs_gpcm <- function(a.item, tau.item, Qmatrix, theta.k, VV, K, TP,
    eps=0, use_log=FALSE)
{
    a <- a.item
    b <- tau.item
    res <- .rm.pcm.calcprobs( a=a, b=b, Qmatrix=Qmatrix, theta.k=theta.k, I=VV, K=K, TP=TP )
    if (eps > 0){
        res[ res < eps ] <- eps
    }
    if (use_log){
        res <- log(res)
    }
    return(res)
}
