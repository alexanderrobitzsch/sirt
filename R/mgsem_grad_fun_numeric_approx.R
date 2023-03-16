## File Name: mgsem_grad_fun_numeric_approx.R
## File Version: 0.06
## File Last Change: 2022-02-25


mgsem_grad_fun_numeric_approx <- function(x, opt_fun_args)
{
    h <- opt_fun_args$technical$h
    NP <- opt_fun_args$NP
    grad <- rep(0,NP)
    for (pp in seq_len(NP) ){
        coef1 <- mgsem_add_increment(x=x,h=h, i1=pp)
        coef2 <- mgsem_add_increment(x=x,h=-h, i1=pp)
        ll1 <- mgsem_opt_fun(x=coef1, opt_fun_args=opt_fun_args)
        ll2 <- mgsem_opt_fun(x=coef2, opt_fun_args=opt_fun_args)
        D1 <- (ll1-ll2)/(2*h)
        grad[pp] <- grad[pp] + D1
    }
    return(grad)
}
