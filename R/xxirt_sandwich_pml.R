## File Name: xxirt_sandwich_pml.R
## File Version: 0.118

xxirt_sandwich_pml <- function(object, h=1e-4)
{
    requireNamespace('MASS')
    weights <- object$weights
    partable <- object$partable
    customTheta <- object$customTheta
    par1 <- xxirt_partable_extract_freeParameters( partable=partable )
    par2 <- xxirt_parTheta_extract_freeParameters( customTheta=customTheta )
    par0 <- par <- x <- c(par1, par2)
    NP <- length(x)
    parnames <- names(x)

    # compute gradient function
    opt_fun <- xxirt_nr_opt_fun_pml_casewise
    args <- list(x=x, object=object, eps=1e-14)
    grad <- xxirt_nr_grad_fun_pml_casewise(x=x, args=args, opt_fun=opt_fun, h=h)

    #- compute variance matrix
    G1 <- grad*sqrt(weights)
    A <- crossprod(G1)
    rownames(A) <- colnames(A) <- parnames

    #- compute hessian
    hess <- matrix(NA, nrow=NP, ncol=NP, dimnames=list(parnames,parnames))

    pp <- 1
    for (pp in 1L:NP){
        x1 <- sirt_add_increment(x=x, pos=pp, value=h)
        grad1 <- xxirt_nr_grad_fun_pml_casewise(x=x1, args=args, opt_fun=opt_fun, h=h)
        hess_pp <- (grad1-grad)/h
        hess[,pp] <- colSums(hess_pp*weights)
    }
    hess <- ( hess + t(hess) ) / 2
    B1 <- MASS::ginv(hess)

    V <- B1 %*% A %*% t(B1)
    se <- sqrt(diag(V))
    names(se) <- parnames

    #-- output
    res <- list(V=V, A=A, hess=hess, se=se)
    return(res)
}
