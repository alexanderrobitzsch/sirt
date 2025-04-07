## File Name: linking_2groups_numerical_gradient.R
## File Version: 0.03


linking_2groups_numerical_gradient <- function(fun, args, h=1e-4)
{
    NP <- length(args$x)
    grad <- rep(0, NP)
    args2 <- args
    for (pp in 1L:NP){
        args2$x <- mgsem_add_increment(x=args$x, h=h, i1=pp)
        f1 <- do.call( what=fun, args=args2)
        args2$x <- mgsem_add_increment(x=args$x, h=-h, i1=pp)
        f2 <- do.call( what=fun, args=args2)
        grad[pp] <- (f1-f2)/(2*h)
    }
    return(grad)
}
