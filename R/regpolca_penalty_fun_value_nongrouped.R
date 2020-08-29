## File Name: regpolca_penalty_fun_value_nongrouped.R
## File Version: 0.02

regpolca_penalty_fun_value_nongrouped <- function(x_ii, combis_ii,
    regular_lam, eps, penalty_used)
{
    diff_ii <- x_ii[ combis_ii[,1] ] - x_ii[ combis_ii[,2] ]
    a1 <- penalty_used( x=diff_ii, lambda=regular_lam, eps=eps )
    return(a1)
}
