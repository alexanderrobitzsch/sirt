## File Name: lq_fit.R
## File Version: 0.094

lq_fit <- function(y, X, w=NULL, pow=2, eps=1e-3, beta_init=NULL,
        optimizer="optim", eps_vec=10^seq(0,-10, by=-.5) )
{
    N <- length(y)
    if (is.null(w)){
        w <- rep(1,N)
    }
    ind <- ! is.na(y)
    y <- y[ind]
    X <- X[ind,]
    w <- w[ind]
    if (is.null(beta_init)){
        mod <- stats::lm.wfit(y=y, x=X, w=w)
        beta_init <- mod$coefficients
    }

    #-- define optimization functions
    fct_optim <- function(x, y, X, pow, eps,w)
    {
        beta <- x
        e <- y - X %*% beta
        val <- sum( w*(e^2 + eps )^(pow/2) )
        return(val)
    }
    grad_optim <- function(x, y, X, pow, eps,w)
    {
        beta <- x
        e <- ( y - X %*% beta )[,1]
        h1 <- pow*(e^2 + eps)^(pow/2-1)*e*w
        der <- - colSums(X*h1)
        return(der)
    }

    #- define epsilon sequence
    eps_vec <- sirt_define_eps_sequence(eps=eps, eps_vec=eps_vec)

    #-- optimize
    for (eps in eps_vec){
        res_optim <- sirt_optimizer(optimizer=optimizer, par=beta_init, fn=fct_optim,
                        gr=grad_optim, y=y, X=X, pow=pow, eps=eps,w=w)
        beta <- res_optim$par
        converged <- res_optim$converged
        e <- ( y - X %*% beta )[,1]
        beta_init <- beta
    }

    #-- output
    res <- list(coefficients=beta, residuals=e, pow=pow, eps=eps,
                converged=converged, res_optim=res_optim, eps_vec=eps_vec)
    return(res)
}
