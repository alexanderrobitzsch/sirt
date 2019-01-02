## File Name: invariance_alignment_define_parameters.R
## File Version: 0.03

invariance_alignment_define_parameters <- function(x, ind_alpha, ind_psi)
{
    alpha0 <- c(0, x[ind_alpha])
    psi0 <- c(1, x[ind_psi])
    #--- output
    res <- list(alpha0=alpha0, psi0=psi0)
    return(res)
}
