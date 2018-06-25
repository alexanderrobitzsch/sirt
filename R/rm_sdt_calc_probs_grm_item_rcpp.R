## File Name: rm_sdt_calc_probs_grm_item_rcpp.R
## File Version: 0.03


rm_sdt_calc_probs_grm_item_rcpp <- function( tau.item, a.item, theta.k, VV,
    K, TP, eps=0, use_log=FALSE)
{
    prob_item <- sirt_rcpp_rm_sdt_calc_probs_grm_item( tau_item=tau.item, a_item=a.item,
                    theta_k=theta.k, VV=VV, K=K, TP=TP, eps=eps, use_log=use_log)
    prob.item <- array(prob_item, dim=c(VV,K+1,TP) )
    return(prob.item)
}
