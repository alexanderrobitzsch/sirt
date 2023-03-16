## File Name: lsdm_est_logist_2pl.R
## File Version: 0.225


lsdm_est_logist_2pl <- function( y, theta, wgt_theta, method="L-BFGS-B",
        b0=NULL, a0=NULL)
{
    y <- as.numeric(y)

    #-- define optimization function
    dist_irf_val <- function(x){
        a <- x[2]
        b <- x[1]
        irf_est <- stats::plogis(a*(theta-b))
        diff_irf <- y - irf_est
        val <- sum( diff_irf^2 * wgt_theta )
        return(val)
    }
    dist_irf_grad <- function(x){
        grad <- rep(0,2)
        a <- x[2]
        b <- x[1]
        irf_est <- stats::plogis(a*(theta-b))
        diff_irf <- y - irf_est
        der1 <- 2*diff_irf*irf_est*(1-irf_est)*wgt_theta
        grad[1] <- sum( der1*a )
        grad[2] <- -sum( der1*(theta-b) )
        return(grad)
    }

    #* perform optimization
    if (is.null(b0)){
        b0 <- theta[ which.min( abs(y-0.5) ) ]
        par0 <- c(b0, 1)
    } else {
        par0 <- c( b0, a0)
    }
    res <- stats::optim( par=par0, fn=dist_irf_val, gr=dist_irf_grad, method=method)

    #--- output
    res <- c( res$par, sqrt(res$value) )
    return(res)
}
