## File Name: lsdm_est_logist_rasch.R
## File Version: 0.283


lsdm_est_logist_rasch <- function( y, theta, wgt_theta, method="L-BFGS-B" )
{
    y <- as.numeric(y)

    #-- define optimization function
    dist_irf_val <- function(x){
        irf_est <- stats::plogis(theta-x)
        diff_irf <- y - irf_est
        val <- sum( diff_irf^2 * wgt_theta )
        return(val)
    }
    dist_irf_grad <- function(x){
        irf_est <- stats::plogis(theta-x)
        diff_irf <- y - irf_est
        grad1 <- 2*diff_irf*irf_est*(1-irf_est)
        grad <- sum(grad1*wgt_theta)
        return(grad)
    }

    #* optimization
    par0 <- theta[ which.min( abs(y-0.5) ) ]
    res <- stats::optim(par=par0, fn=dist_irf_val, gr=dist_irf_grad, method=method)

    #--- output
    res <- c( res$par, sqrt(res$value) )
    return(res)
}
