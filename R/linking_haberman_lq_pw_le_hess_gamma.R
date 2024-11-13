## File Name: linking_haberman_lq_pw_le_hess_gamma.R
## File Version: 0.04


linking_haberman_lq_pw_le_hess_gamma <- function(par_delta, par_gamma, des, h=1e-4)
{
    I <- des$I
    G <- des$G
    NPD <- length(par_delta)
    NPG <- length(par_gamma)
    hess <- matrix(0, nrow=NPD, ncol=NPG)
    rownames(hess) <- names(par_delta)
    colnames(hess) <- names(par_gamma)
    grad_fun <- linking_haberman_lq_pw_le_grad
    pgi <- is.na(par_gamma)
    if (sum(pgi)>0){
        par_gamma[ pgi ] <- 0
    }

    for (gg in 1L:G){
        for (ss in 1L:2){
            ind_gg <- 2*I*(gg-1)+seq(ss,2*I,by=2)
            par_gamma2 <- mgsem_add_increment(x=par_gamma, h=h, i1=ind_gg)
            args <- list(par_delta=par_delta, par_gamma=par_gamma2, des=des,h=h)
            opt1 <- do.call(what=grad_fun, args=args)
            args$par_gamma <- mgsem_add_increment(x=par_gamma, h=-h, i1=ind_gg)
            opt2 <- do.call(what=grad_fun, args=args)
            hess1 <- (opt1-opt2)/(2*h)
            hess[,ind_gg] <- t(hess1)
        }
    }

    return(hess)
}
