## File Name: rasch_pairwise_iterations.R
## File Version: 0.10


rasch_pairwise_iterations <- function(eps, y.ij, delta.ij, conv, maxiter,
    zerosum, progress, fix_first=FALSE)
{
    change <- 1
    iter <- 0

    #* start estimation algorithm
    while( change > conv & iter < maxiter ){
        eps0 <- eps
        b0 <- b
        eps <- sqrt( rowSums( y.ij * eps * delta.ij ) / colSums( y.ij / eps ) )
        if (zerosum){
            eps <- rasch_pairwise_zerosum(eps=eps)
        }
        if (fix_first){
            eps[1] <- 1
        }
        b <- -log(eps)
        change <- max( abs( eps0 - eps ) )
        change_b <- max(abs( -b0 + b))
        iter <- iter + 1
        if ( progress ){
            cat( "PL Iter.", iter, ": max. parm. change=", round( change_b, 6 ), "\n")
            utils::flush.console()
        }
    } #* end estimation algorithm
    if (zerosum){
        eps <- rasch_pairwise_zerosum(eps=eps)
        b <- -log(eps)
    }
    #--- output
    res <- list(eps=eps, b=b, iter=iter)
    return(res)
}
