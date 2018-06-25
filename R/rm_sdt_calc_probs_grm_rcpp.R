## File Name: rm_sdt_calc_probs_grm_rcpp.R
## File Version: 0.03


rm_sdt_calc_probs_grm_rcpp <- function(c.rater, d.rater, I, K, eps=0, use_log=FALSE)
{
    res <- sirt_rcpp_rm_sdt_calc_probs_grm_rater( c_rater=c.rater, d_rater=d.rater, I=I, K=K,
                    eps=eps, use_log=use_log)
    K1 <- K + 1
    res <- array(res, dim=c(I,K1,K1))
    return(res)
}
