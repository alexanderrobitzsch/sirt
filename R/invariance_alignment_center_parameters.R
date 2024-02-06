## File Name: invariance_alignment_center_parameters.R
## File Version: 0.130

invariance_alignment_center_parameters <- function(alpha0, psi0, center,
    reparam=FALSE, convert=FALSE, meth=1)
{
    if (reparam & convert){
        alpha0 <- alpha0 * psi0
    }
    if (center){
        alpha0 <- alpha0 - mean(alpha0)
        log_psi <- log(psi0)
        psi0 <- exp(log_psi - mean(log_psi))
    }
    if (meth %in% c(0,0.5) ){
        fac <- 1/psi0[1]
        psi0 <- fac*psi0
        alpha0 <- fac*alpha0
        if (meth==0){
            alpha0 <- alpha0 - alpha0[1]
        }
    }
    #--- output
    res <- list(alpha0=alpha0, psi0=psi0)
    return(res)
}
