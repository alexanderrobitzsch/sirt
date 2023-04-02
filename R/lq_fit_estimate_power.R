## File Name: lq_fit_estimate_power.R
## File Version: 0.15


lq_fit_estimate_power <- function(e, pow_init=2, lower_pow=.1, upper_pow=10)
{
    pow0 <- pow_init
    me <- mean(e)
    n <- length(e)
    est_pow <- TRUE
    if (stats::var(e)<1e-5){
        est_pow <- FALSE
        p <- pow0
    }
    if (est_pow){
        I <- sum( abs(e-me) / n ) / sqrt( sum( (e-me)^2 ) / (n-1) )
        vi <- 1/I
        # copied part from normalp::estimatep function
        fvi <- function(p){
            (vi - sqrt(gamma(1/p) * gamma(3/p))/gamma(2/p))^2
        }
        res <- stats::optim(par=pow0, fn=fvi, lower=lower_pow, upper=upper_pow,
                            method='L-BFGS-B')
        p <- res$par
    }
    return(p)
}
