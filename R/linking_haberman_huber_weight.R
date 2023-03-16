## File Name: linking_haberman_huber_weight.R
## File Version: 0.04
## File Last Change: 2019-01-14

linking_haberman_huber_weight <- function(x, cutoff)
{
    eps <- 1e-10
    wgt_adj <- (abs(x) >=cutoff)*cutoff / abs(x) + (abs(x) < cutoff)*1
    return(wgt_adj)
}
