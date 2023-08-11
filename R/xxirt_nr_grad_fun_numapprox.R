## File Name: xxirt_nr_grad_fun_numapprox.R
## File Version: 0.04

xxirt_nr_grad_fun_numapprox <- function(x, em_args)
{
    grad <- 0*x
    ll1 <- xxirt_nr_optim_fun(x=x, em_args=em_args)
    h <- em_args$h
    for (pp in 1:em_args$NP){
        x1 <- sirt_add_increment(x=x, pos=pp, value=h)
        ll2 <- xxirt_nr_optim_fun(x=x1, em_args=em_args)
        grad[pp] <- (ll2-ll1) / h
    }
    return(grad)
}
