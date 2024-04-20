## File Name: linking_haberman_compute_lts_mean.R
## File Version: 0.091

linking_haberman_compute_lts_mean <- function(x, w, lts_prop, maxiter=10)
{
    res <- linking_haberman_remove_missings_vector(x=x,w=w)
    x <- res$x
    w <- res$w
    n <- length(x)
    k <- max(2, ceiling(n*lts_prop))
    if (n < k){
        k <- n
    }
    index <- 1L:k
    m1 <- stats::weighted.mean(x=x, w=w)
    iter <- 0
    iterate <- TRUE
    dfr <- data.frame(index=1L:n, x=x, e=x-m1, w=w)
    while(iterate){
        m0 <- m1
        dfr <- dfr[ order(abs(dfr$e)), ]
        m1 <- stats::weighted.mean(x=dfr$x[1L:k], w=dfr$w[1L:k])
        dfr$e <- dfr$x - m1
        iter <- iter + 1
        if (iter>maxiter){ iterate <- FALSE }
        if (abs(m0-m1) < 1e-4 ){ iterate <- FALSE }
    }
    res <- m1
    return(res)
}
