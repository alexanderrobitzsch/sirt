## File Name: rasch_mml2_mstep_calc_loglike.R
## File Version: 0.02
## File Last Change: 2019-05-10


rasch_mml2_mstep_calc_loglike <- function( exp_r, prob1, exp_n, prob0=NULL)
{
    if (is.null(prob0)){
        prob0 <- 1 - prob1
    }
    y <- rowSums( exp_r * log(prob1) + (exp_n - exp_r) * log(prob0) )
    return(y)
}
