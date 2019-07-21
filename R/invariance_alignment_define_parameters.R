## File Name: invariance_alignment_define_parameters.R
## File Version: 0.09

invariance_alignment_define_parameters <- function(x, ind_alpha, ind_psi,
        fix_first_psi=TRUE, reparam=FALSE)
{
    alpha0 <- c(0, x[ind_alpha])
    if (fix_first_psi){
        psi0 <- c(1, x[ind_psi])
    } else {
        psi0 <- c(x[ind_psi])
    }
    if (reparam){
        prod_psi <- prod(psi0)
        NX <- length(psi0)
        psi0 <- psi0 / ( prod_psi^(1/NX) )
    }
    #--- output
    res <- list(alpha0=alpha0, psi0=psi0)
    return(res)
}
