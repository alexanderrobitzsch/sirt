## File Name: rm_hrm_est_c_rater_order_parameters.R
## File Version: 0.01

rm_hrm_est_c_rater_order_parameters <- function(c.rater, K)
{
    for (kk in 1:K){
        if (kk>1){
            ind <- which( c.rater[,kk] < c.rater[,kk-1] )
            if (length(ind)>0){
                l1 <- c.rater[ind,kk-1]
                c.rater[ind, kk-1] <- c.rater[ind,kk]
                c.rater[ ind , kk ] <- l1
            }
        }
    }
    return(c.rater)
}
    
