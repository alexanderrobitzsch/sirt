## File Name: noharm_sirt_implied_cov.R
## File Version: 0.01

noharm_sirt_implied_cov <- function(Fmat, Pmat, Psimat)
{
    gamma_val <- (Fmat %*% Pmat) %*% t(Fmat) + Psimat
    return(gamma_val)
}
