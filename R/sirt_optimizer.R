## File Name: sirt_optimizer.R
## File Version: 0.365

sirt_optimizer <- function(optimizer, par, fn, grad=NULL, method="L-BFGS-B",
        hessian=TRUE, control=list(), ...)
{
    s1 <- Sys.time()
    h <- 1e-5

    #-- optim
    if (optimizer=="optim"){
        res <- stats::optim(par=par, fn=fn, gr=grad, method=method,
                        hessian=hessian, control=control, ...)
        res$iter <- res$counts['function']
    }

    #-- Rvmmin
    if (optimizer=="Rvmmin"){
        requireNamespace("optimx")
        args <- list(...)
        NP <- length(par)
        if ( ! ( "lower" %in% names(args) ) ){
            args$lower <- rep(-Inf,NP)
        }
        if ( ! ( "upper" %in% names(args) ) ){
            args$upper <- rep(Inf,NP)
        }
        # args$bdmsk <- 1*( ( args$lower > - Inf ) + ( args$upper < Inf )==0 )
        args$par <- par
        args$fn <- fn
        args$gr <- grad
        args$method <- "Rvmmin"
        args$control <- control
        optimx_fun <- optimx::optimx
        res0 <- do.call(what=optimx_fun, args=args)
        # res <- optimx::Rvmmin(par=par, fn=fn, gr=grad, control=control, ...)
        res <- list()
        m1 <- as.vector(as.numeric(res0[1,1:NP]))
        res$par <- m1
        res$convergence <- res0$convcode[1]
        res$iter <- res0$feval[1]
        res$value <- res0$value[1]
    }

    #-- nlminb
    if (optimizer=="nlminb"){
        res <- stats::nlminb(start=par, objective=fn, gradient=grad,
                        control=control, ...)
        res$value <- res$objective
        res$iter <- res$iterations
    }

    #-- bobyqa
    if (optimizer=="bobyqa"){
        requireNamespace("minqa")
        # control=list(iprint=2, maxfun=10000, ...)
        res <- minqa::bobyqa(par=par, fn=fn, control=control, ...)
        res$value <- res$fval
        res$iter <- res$feval
        res$convergence <- res$ierr
    }

    #-- nloptr
    if (optimizer=="nloptr"){
        requireNamespace("nloptr")
        # print_level, maxevel
        # control$algorithm <- "NLOPT_LD_LBFGS"
        control$algorithm <- "NLOPT_LD_MMA"
        args <- list(...)
        args <- sirt_rename_list_names(x=args, old="lower", new="lb")
        args <- sirt_rename_list_names(x=args, old="upper", new="ub")
        args$x0 <- par
        args$eval_f <- fn
        args$eval_grad_f <- grad
        args$opts <- control
        res <- do.call( what=nloptr::nloptr, args=args )
        res$par <- res$solution
        res$value <- res$objective
        res$iter <- res$iterations
        res$convergence <- 1 - ( res$status %in% c(4) )
    }

    #*****
    #-- compute Hessian
    comp_hess_optimizers <- c("nlminb", "bobyqa")
    if (hessian & ( optimizer %in% comp_hess_optimizers ) ){
        res <- sirt_optimizer_hessian(res=res, fn=fn, grad=grad, h=h, ...)
    }

    #-- names for Hessian
    if (!is.null(res$hessian) & (!is.null(names(par))) ){
        rownames(res$hessian) <- colnames(res$hessian) <- names(par)
    }

    s2 <- Sys.time()
    res$time_diff <- s2 - s1
    res$optimizer <- optimizer
    res$converged <- res$convergence==0
    #--- output
    return(res)
}
