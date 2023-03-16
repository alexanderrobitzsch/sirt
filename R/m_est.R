## File Name: m_est.R
## File Version: 0.150

m_est <- function(data, par, optfun_case=NULL, gradfun_case=NULL,
                bread=TRUE, optimizer="optim",
                method="L-BFGS-B", control=list(), ... )
{
    requireNamespace('MASS')

    is_optfun_case <- ! is.null(optfun_case)
    is_gradfun_case <- ! is.null(gradfun_case)
    optfun_case0 <- optfun_case
    gradfun_case0 <- gradfun_case
    N <- nrow(data)

    #-- define case-wise gradient
    if (! is_gradfun_case ){
        gradfun_case <- function(par, data, h=1e-4)
        {
            what_optfun <- 'optfun_case'
            args <- list(par=par, data=data)
            val1 <- do.call(what=what_optfun, args=args)
            N <- length(val1)
            NP <- length(par)
            grad_case <- matrix(NA, nrow=N, ncol=NP)
            colnames(grad_case) <- names(par)
            args1 <- args
            for (pp in 1:NP){
                par1 <- m_est_add_increment(x=par, pos=pp, h=h)
                args1$par <- par1
                val2 <- do.call(what=what_optfun, args=args1)
                grad_case[,pp] <- (val2-val1)/h
            }
            return(grad_case)
        }
    }

    #-- define gradient
    gradfun <- function(x, data)
    {
        args <- list(par=x, data=data)
        res <- do.call(what='gradfun_case', args=args)
        res <- colSums(res)
        return(res)
    }

    #-- define optimization function
    if (is_optfun_case){
        optfun <- function(x, data)
        {
            args <- list(par=x, data=data)
            val <- do.call(what='optfun_case', args=args)
            val <- sum(val)
            return(val)
        }
    } else {
        optfun <- function(x, data)
        {
            args <- list(par=x, data=data)
            val <- do.call(what='gradfun_case', args=args)
            val <- sum( colSums(val)^2 )/2
            return(val)
        }
        gradfun <- NULL
    }

    #***** optimization
    mod <- sirt_optimizer(optimizer=optimizer, par=par, fn=optfun, data=data,
                    grad=gradfun, method=method, hessian=TRUE, control=control, ... )
    coef <- mod$par
    parnames <- names(par)
    NP <- length(parnames)

    #-- vcov based on observed information
    V_obs <- MASS::ginv(X=mod$hessian)
    V_obs <- rowcolnames(x=V_obs, names=parnames)
    B <- mod$hessian / N

    #-- compute scores and sandwich vcov
    scores <- -gradfun_case(data=data, par=coef)
    A <- crossprod(scores)
    if (bread){
        NP <- ncol(scores)
        hess <- array(0, dim=c(N,NP,NP))
        h <- 1e-4
        pp <- 1
        for (pp in 1:NP){
            par_pp <- m_est_add_increment(x=coef, pos=pp, h=h)
            scores2 <- -gradfun_case(data=data, par=par_pp)
            hess[,pp,] <- (scores2-scores)/h
        }
        B <- matrix(0, nrow=NP, ncol=NP)
        for (pp in 1:NP){
            for (hh in 1:NP){
                B[pp,hh] <- mean( hess[,pp,hh] )
            }
        }
        B <- -B
        # B <- ( B + t(B) ) / 2
        A1 <- A / N^2
        res <- dmlavaan_sandwich_formula(A=A1, B=B, parnames=parnames)
        V_sw <- res$V
        V_obs <- NULL
    } else {
        V_sw <- V_obs %*% A %*% V_obs
    }

    #-- create parameter table
    partable <- data.frame(id=1:NP, parname=parnames, est=coef)
    if (!bread){
        partable$se_obs <- sqrt_diag(x=V_obs)
    }
    partable$se_sw <- sqrt_diag(x=V_sw)
    rownames(partable) <- NULL

    #--- output
    res <- list(coef=coef, partable=partable, vcov=V_sw, parnames=parnames, mod=mod,
                    scores=scores, V_obs=V_obs, V_sw=V_sw, A=A, B=B, NP=NP,
                    is_optfun_case=is_optfun_case, is_gradfun_case=is_gradfun_case,
                    optfun_case=optfun_case0, gradfun_case=gradfun_case0,
                    optimizer=optimizer, method=method,
                    control=control, data=data, N=N)
    class(res) <- 'm_est'
    return(res)
}
