## File Name: locpolycor_est_thresh_item.R
## File Version: 0.189


locpolycor_est_thresh_item <- function(y, data.mod, x0, w, model, par_init=NULL,
        eps=1e-10, optimizer="optim", use_deriv=TRUE)
{
    K <- max(y)
    N <- length(y)
    x1 <- data.mod
    # starting values
    par0 <- stats::qnorm( seq(1/(K+1), K/(K+1), length=K) )
    if (model=='lin'){
        par1 <- 0*par0
        names(par1) <- paste0('lin', 1L:K)
        par0 <- c(par0, par1)
    }
    if (! is.null(par_init) ){
        par0 <- par_init
    }
    #- optimization
    if (use_deriv){
        grad <- locpolycor_est_thresh_grad_fun
    } else {
        grad <- NULL
    }

    res_optim <- sirt_optimizer(optimizer=optimizer, par=par0,
                        fn=locpolycor_est_thresh_opt_fun,
                        grad=locpolycor_est_thresh_grad_fun,
                        y=y, x1=x1, x0=x0, w=w, model=model, method='L-BFGS-B',
                        K=K, eps=eps, hessian=FALSE )

    #- local thresholds across all persons at x0
    thresh <- res_optim$par
    if (model=='lin'){
        thresh <- thresh[1L:K]
    }
    #- individual threshold predictions
    thresh_ind <- matrix(NA, nrow=N, ncol=K)
    colnames(thresh_ind) <- paste0('Cat',1L:K)
    for (kk in 1L:K){
        if (model=='const'){
            thresh_ind[,kk] <- thresh[kk]
        }
        if (model=='lin'){
            thresh_ind[,kk] <- thresh[kk] + res_optim$par[K+kk] * (x1-x0)
        }
    }

    #*** output
    res <- list(thresh=thresh, thresh_ind=thresh_ind, x0=x0, res_optim=res_optim, N=N,
                        W=sum(w), K=K, model=model)
    return(res)
}
