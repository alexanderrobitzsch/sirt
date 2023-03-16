## File Name: lq_fit.R
## File Version: 0.152
## File Last Change: 2021-03-07

lq_fit <- function(y, X, w=NULL, pow=2, eps=1e-3, beta_init=NULL,
        est_pow=FALSE, optimizer="optim", eps_vec=10^seq(0,-10, by=-.5),
        conv=1e-4, miter=20, lower_pow=.1, upper_pow=5)
{
    N <- length(y)
    if (is.null(w)){
        w <- rep(1,N)
    }
    ind <- ! is.na(y)
    y <- y[ind]
    X <- X[ind,,drop=FALSE]
    w <- w[ind]
    if (is.null(beta_init)){
        mod <- stats::lm.wfit(y=y, x=X, w=w)
        beta_init <- mod$coefficients
    }

    #-- define optimization functions
    Xs <- sirt_rcpp_lq_fit_analyze_matrix(X=X)

    fct_optim <- function(x, y, X, pow, eps, w, Xs=NULL)
    {
        beta <- x
        # e <- y - X %*% beta
        # e <- sirt_rcpp_lq_fit_matrix_mult( Z=Xs, y=y, beta=beta)
        # val <- sum( w*(e^2 + eps )^(pow/2) )
        val <- sirt_rcpp_lq_fit_fct_optim(Z=Xs, y=y, beta=beta, pow=pow, w=w,
                    eps=eps)
        return(val)
    }

    grad_optim <- function(x, y, X, Xs, pow, eps, w)
    {
        beta <- x
        # e <- ( y - X %*% beta )[,1]
        e <- sirt_rcpp_lq_fit_matrix_mult( Z=Xs, y=y, beta=beta)
        pow2 <- pow/2-1
        h1 <- pow*exp(pow2 * log( e^2 + eps ))*e*w
        # der <- - colSums(X*h1)
        px <- ncol(X)
        der <- sirt_rcpp_lq_fit_sparse_matrix_derivative(Z=Xs, h1=h1, px=px)
        return(der)
    }

    #- define epsilon sequence
    eps_vec <- sirt_define_eps_sequence(eps=eps, eps_vec=eps_vec)
    pow00 <- pow

    #-- optimize
    for (eps in eps_vec){
        iterate_powers <- TRUE
        iter <- 0
        pow <- pow00
        while (iterate_powers){
            pow0 <- pow
            res_optim <- sirt_optimizer(optimizer=optimizer, par=beta_init, fn=fct_optim,
                            grad=grad_optim, y=y, X=X, Xs=Xs, pow=pow, eps=eps,w=w)
            beta <- res_optim$par
            converged <- res_optim$converged
            e <- ( y - X %*% beta )[,1]
            beta_init <- beta
            if (est_pow){
                pow <- lq_fit_estimate_power(e=e, pow_init=pow, lower_pow=lower_pow,
                            upper_pow=upper_pow)
            }
            pow_change <- abs(pow-pow0)
            if (( pow_change < conv ) | ( iter > miter ) ){
                iterate_powers <- FALSE
            }
            iter <- iter + 1
            # cat( paste0("eps=",eps, " | pow=",pow,"\n"))
        }
    }

    #-- output
    res <- list(coefficients=beta, residuals=e, pow=pow, eps=eps,
                converged=converged, res_optim=res_optim, eps_vec=eps_vec)
    return(res)
}
