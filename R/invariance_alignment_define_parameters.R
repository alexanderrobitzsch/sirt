## File Name: invariance_alignment_define_parameters.R
## File Version: 0.215

invariance_alignment_define_parameters <- function(x, ind_alpha, ind_psi,
        fix_first_psi=TRUE, reparam=FALSE, constraint="prod", overparam=FALSE,
        meth=1)
{
    # meth=0=> free optimization with ridge regularization
    alpha0 <- c(0, x[ind_alpha])
    if (meth==0){
        alpha0 <- alpha0[-c(1)]
    }
    if (fix_first_psi){
        psi0 <- c(1, x[ind_psi])
    } else {
        psi0 <- c(x[ind_psi])
    }
    if (( reparam & ( constraint=='prod') ) | (meth %in% c(0,0.5) ) ){
        prod_psi <- prod(psi0)
        NX <- length(psi0)
        log_psi0 <- log(psi0)
        psi0 <- exp( log_psi0 - mean(log_psi0) )
        # psi0 <- psi0 / ( prod_psi^(1/NX) )
    }
    if (reparam & ( constraint=='sum') ){
        psi0 <- psi0 - mean(psi0) + 1
    }
    if (meth==0){
        # alpha0[1] <- -sum(alpha0[-c(1)])
        # alpha0 <- alpha0 - mean(alpha0)
        # Implementing the constraints leads to bias
    }
    #--- output
    res <- list(alpha0=alpha0, psi0=psi0)
    return(res)
}
