## File Name: linking_haberman_lq_pw_le_hess_delta.R
## File Version: 0.06


linking_haberman_lq_pw_le_hess_delta <- function(par_delta, par_gamma, des, h=1e-4)
{
    I <- des$I
    G <- des$G
    NPD <- length(par_delta)
    hess <- matrix(0, nrow=NPD, ncol=NPD)
    rownames(hess) <- colnames(hess) <- names(par_delta)
    grad_fun <- linking_haberman_lq_pw_le_grad
    for (pp in 1L:NPD){
        par_delta2 <- mgsem_add_increment(x=par_delta, h=h, i1=pp)
        args <- list(par_delta=par_delta2, par_gamma=par_gamma, des=des,h=h)
        opt1 <- do.call(what=grad_fun, args=args)
        args$par_delta <- mgsem_add_increment(x=par_delta, h=-h, i1=pp)
        opt2 <- do.call(what=grad_fun, args=args)
        hess[,pp] <- (colSums(opt1)-colSums(opt2))/(2*h)
    }
    return(hess)
}
