## File Name: rm_facets_pp_mle.R
## File Version: 0.183


#*** person parameter estimation in partial credit model
rm_facets_pp_mle <- function( data, a, b, theta, WLE=FALSE,
    maxiter=30, maxincr=3, h=.001, convP=.0001, maxval=9.99,
    progress=TRUE )
{

    theta0 <- theta
    N <- length(theta)
    I <- ncol(data)
    iter <- 1
    conv <- 1E5
    
    args <- list(data=data, a=a, b=b, theta=theta)    
    args_change <- function(args, theta){
        args$theta <- theta
        return(args)
    }
    
    #-------- begin algorithm
    while ( ( conv > convP ) & (  iter <=maxiter ) ){

        theta0 <- theta        
        ll0 <- do.call(rm_facets_pp_mle_calc_ll_theta, args_change(args, theta))
        llP1 <- do.call(rm_facets_pp_mle_calc_ll_theta, args_change(args, theta+h))
        llM1 <- do.call(rm_facets_pp_mle_calc_ll_theta, args_change(args, theta-h))
        if (WLE){
            llP2 <- do.call(rm_facets_pp_mle_calc_ll_theta, args_change(args, theta+2*h))
            llM2 <- do.call(rm_facets_pp_mle_calc_ll_theta, args_change(args, theta-2*h))
        }

        ll1 <- (llP1 - llM1) / (2*h)
        ll2 <- (llP1 - 2*ll0 + llM1) / h^2
        incr <- - ll1 / ll2
        if (WLE){
            ll3 <- ( llP2 - 2*llP1 + 2*llM1 - llM2 ) / (2*h^3 )
            incr <- - ll1 / ll2 - ll3 / ( 2*ll2^2 )
        }

        maxincr <- maxincr / 1.05
        incr <- rm_numdiff_trim_increment( increment=incr, max.increment=maxincr, eps2=0 )        
        theta <- theta + incr
        theta <- ifelse( abs(theta) > maxval, sign(theta)*maxval, theta )
        conv <- max( abs( theta - theta0) )
        if (progress){
            cat("* Iteration", iter, ":", "maximum parameter change", round( conv, 5), "\n")
            utils::flush.console()
        }
        iter <- iter + 1
    }

    #-- output
    se <- sqrt( abs( - 1/ll2 ) )
    se <- ifelse( abs(theta)==maxval, NA, se )
    theta <- ifelse( theta==maxval, Inf, theta )
    theta <- ifelse( theta==- maxval, -Inf, theta )
    res <- data.frame( est=theta, se=se )
    return(res)
}


