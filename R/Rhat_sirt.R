## File Name: Rhat_sirt.R
## File Version: 1.12


####################################
# auxiliary functions for Rhat statistic
    ############################################################################
    # Code from rube package
    # Source: http://www.stat.cmu.edu/~hseltman/rube/rube0.2-16/R/Rhat.R
    # Inference from Iterative Simulation Using Multiple Sequences
    # Author(s): Andrew Gelman and Donald B. Rubin
    # Source: Statistical Science, Vol. 7, No. 4 (Nov., 1992), pp. 457-472
    # Stable URL: http://www.jstor.org/stable/2246093
    ## Matches gelman.diag() from package "coda", but not WinBUGS() "summary" component.
    ## Better than gelman.diag() because multivariate stat is not
    ## bothered to be calculated
Rhat1 <- function(mat)
{
    m <- ncol(mat)
    n <- nrow(mat)
    b <- apply(mat,2,mean)
    B <- sum((b-mean(mat))^2)*n/(m-1)
    w <- apply(mat,2,stats::var)
    W <- mean(w)
    s2hat <- (n-1)/n*W + B/n
    Vhat <- s2hat + B/m/n
    covWB <- n /m * (stats::cov(w,b^2)-2*mean(b)*stats::cov(w,b))
    varV <- (n-1)^2 / n^2 * stats::var(w)/m +
                (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
                2 * (m-1)*(n-1)/m/n^2 * covWB
    df <- 2 * Vhat^2 / varV
    R <- sqrt((df+3) * Vhat / (df+1) / W)
    return(R)
}


Rhat <- function(arr)
{
    dm <- dim(arr)
    if (length(dm)==2) return(Rhat1(arr))
    if (dm[2]==1) return(NULL)
    if (dm[3]==1) return(Rhat1(arr[,,1]))
    return(apply(arr,3,Rhat1))
}
