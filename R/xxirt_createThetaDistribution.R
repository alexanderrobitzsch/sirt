## File Name: xxirt_createThetaDistribution.R
## File Version: 0.169

xxirt_createThetaDistribution <- function( par, est, P, prior=NULL,
        prior_par1=NULL, prior_par2=NULL, lower=NULL, upper=NULL )
{
    res <- list()
    res$par <- par
    res$est <- est
    res$P <- P
    NP <- length(par)
    res$prior <- prior
    res$prior_par1 <- prior_par1
    res$prior_par2 <- prior_par2

    NPT <- sum(est)
    np1 <- which(est)
    if (is.null(lower)){
        lower <- rep(-Inf,NPT)
    } else {
        lower <- lower[np1]
    }
    if (is.null(upper)){
        upper <- rep(Inf,NPT)
    } else {
        upper <- upper[np1]
    }
    names(lower) <- names(upper) <- np1
    res$lower <- lower
    res$upper <- upper
    res$some_bound <- ( sum(res$lower>-Inf)+sum(res$upper<Inf) )>0
    class(res) <- 'ThetaDistribution'
    return(res)
}
