## File Name: linking_haberman_bisquare_weight.R
## File Version: 0.02

linking_haberman_bisquare_weight <- function(x, cutoff)
{
    wgt_adj <- ( 1 - ( x / cutoff )^2 )^2
    wgt_adj <- (abs(x)<=cutoff)*wgt_adj
    return(wgt_adj)
}
