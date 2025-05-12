## File Name: xxirt_nr_grad_fun_numapprox.R
## File Version: 0.092

xxirt_nr_grad_fun_numapprox <- function(x, em_args, opt_fun)
{
    grad <- 0*x
    ll1 <- opt_fun(x=x, em_args=em_args)
    h <- em_args$h
    for (pp in 1L:em_args$NP){
        x1 <- sirt_add_increment(x=x, pos=pp, value=h)
        ll2 <- opt_fun(x=x1, em_args=em_args)
        x1 <- sirt_add_increment(x=x, pos=pp, value=-h)
        ll2b <- opt_fun(x=x1, em_args=em_args)
        grad[pp] <- (ll2-ll2b) / (2*h)
    }
    return(grad)
}
