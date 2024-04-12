## File Name: xxirt_nr_grad_fun_pml_casewise.R
## File Version: 0.05


#-- case-wise optimization function for PML
xxirt_nr_grad_fun_pml_casewise <- function(x, args, opt_fun, h=1e-4)
{
    NP <- length(x)
    f0 <- do.call(what=opt_fun, args=args)
    N <- length(f0)
    grad <- matrix(0, nrow=N, ncol=NP)
    colnames(grad) <- names(x)
    for (pp in 1L:NP){
        args1 <- args
        args1$x <- sirt_add_increment(x=x, pos=pp, value=h)
        f1 <- do.call(what=opt_fun, args=args1)
        grad[,pp] <- (f1-f0)/h
    }
    #- output
    return(grad)
}
