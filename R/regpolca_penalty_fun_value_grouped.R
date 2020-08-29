## File Name: regpolca_penalty_fun_value_grouped.R
## File Version: 0.01

regpolca_penalty_fun_value_grouped <- function(x_ii, combis_ii,
    regular_lam, eps, penalty_used)
{
    diff_ii <- x_ii[ combis_ii[,1] ] - x_ii[ combis_ii[,2] ]
    norm_ii <- regpolca_grouped_norm(x=diff_ii)
    a1 <- penalty_used( x=norm_ii, lambda=regular_lam, eps=eps )
    return(a1)
}
