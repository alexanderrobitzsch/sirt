## File Name: rm_sdt_mstep_numdiff_diffindex.R
## File Version: 0.02


rm_sdt_mstep_numdiff_diffindex <- function(ll1, ll2, numdiff.parm,
        diffindex )
{
    ll_grad <- ( ll1 - ll2 ) / (2*numdiff.parm)
    ll_grad <- rowsum(ll_grad, diffindex )
    ll_grad <- ll_grad[ rownames(ll_grad) > 0, 1 ]
    return(ll_grad)
}
