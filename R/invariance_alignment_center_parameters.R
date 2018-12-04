## File Name: invariance_alignment_center_parameters.R
## File Version: 0.01

invariance_alignment_center_parameters <- function(alpha0, psi0, center)
{
    if (center){
        alpha0 <- alpha0 - mean(alpha0)
        log_psi <- log(psi0)
        psi0 <- exp(log_psi - mean(log_psi))
    }
    #--- output
    res <- list(alpha0=alpha0, psi0=psi0)
    return(res)
}
