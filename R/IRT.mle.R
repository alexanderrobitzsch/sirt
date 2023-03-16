## File Name: IRT.mle.R
## File Version: 0.34
## File Last Change: 2023-03-08


#--- IRT.mle
IRT.mle <- function(data, irffct, arg.list, theta=rep(0,nrow(data)), type="MLE",
        mu=0, sigma=1, maxiter=20, maxincr=3, h=0.001, convP=1e-04, maxval=9,
        progress=TRUE)
{
    N <- length(theta)
    I <- ncol(data)
    iter <- 1
    MAP <- WLE <- FALSE
    conv <- 1E3
    eps <- 1E-6

    if ( type=="WLE" ){ WLE <- TRUE }
    if ( type=="MAP" ){ MAP <- TRUE }

    #-------- begin algorithm
    while ( ( conv > convP ) & (  iter <=maxiter  ) ){
        theta0 <- theta
        ll0 <- ll1 <- ll2 <- rep(0,N)
        llP2 <- llP1 <- llP0 <- llM1 <- llM2 <- rep(0,N)
        for (ii in 1:I){
            # ii <- 2
            arg.list$ii <- ii
            # theta
            arg.list$theta <- theta
            llP0 <- llP0 + calc.ll( do.call( irffct, arg.list ), data, ii )
            # theta + h
            arg.list$theta <- theta+h
            llP1 <- llP1 + calc.ll( do.call( irffct, arg.list ), data, ii )
            # theta - h
            arg.list$theta <- theta-h
            llM1 <- llM1 + calc.ll( do.call( irffct, arg.list ), data, ii )
            if (WLE){
                # theta + 2*h
                arg.list$theta <- theta+2*h
                llP2 <- llP2 + calc.ll( do.call( irffct, arg.list ), data, ii )
                # theta - 2*h
                arg.list$theta <- theta-2*h
                llM2 <- llM2 + calc.ll( do.call( irffct, arg.list ), data, ii )
            }
        }# end item ii
        #-----
        # prior distribution
        if (MAP){
           llP0 <- llP0 + sirt_dnorm( theta, mean=mu,sd=sigma, log=TRUE)
           llP1 <- llP1 + sirt_dnorm( theta+h, mean=mu,sd=sigma, log=TRUE)
           llP2 <- llP2 + sirt_dnorm( theta+2*h, mean=mu,sd=sigma, log=TRUE)
           llM1 <- llM1 + sirt_dnorm( theta-h, mean=mu,sd=sigma, log=TRUE)
           llM2 <- llM2 + sirt_dnorm( theta-2*h, mean=mu,sd=sigma, log=TRUE)
        }
        #-----
        #*** likelihood
        ll0 <- llP0
        ll1 <- ( llP1 - llP0 ) / h
        ll2 <- ( llP1 - 2*llP0 + llM1 ) / (h^2)
        # ability increments
        M1 <- ll1
        ll2 <- ll2 + eps
        incr <- - ll1 / ll2
        # WLE bias correction term
        if (WLE){
            ll3 <- ( llP2 - 2*llP1 + 2*llM1 - llM2 ) / (2*h^3 )
            incr <- - ll1 / ll2 - ll3 / ( 2*ll2^2 )
        }

        maxincr <- maxincr / 1.05
        incr <- ifelse( abs(incr) > maxincr, sign(incr)*maxincr, incr )
        theta <- theta + incr
        theta <- ifelse( abs(theta) > maxval, sign(theta)*maxval, theta )
        conv <- max( abs( theta - theta0) )
        if (progress){
            cat("* Iteration", iter, ":", "maximum parameter change",
                round( conv, 5 ), "\n")
            utils::flush.console()
        }
        iter <- iter + 1
    }
    #------------- end algorithm

    #--- output
    se <- sqrt( abs( - 1 / ll2 ) )
    se <-  ifelse( abs(theta)==maxval, NA, se )
    theta <- ifelse( theta==maxval, Inf, theta )
    theta <- ifelse( theta==- maxval, -Inf, theta )
    res <- data.frame( "est"=theta, "se"=se )

    #--- calculate reliability
    attr(res, "type") <- type
    if ( type %in% c("MLE","WLE") ){
        attr(res, "reliability") <- mle.reliability( res$est, res$se )
    }
    if ( type %in% c("MAP") ){
        attr(res, "reliability") <- eap.reliability( res$est, res$se )
    }
    #--- output
    return(res)
}


#--- calculate individual likelihood for item ii
calc.ll <- function( probs, data, ii, eps=1e-50 )
{
    N <- nrow(data)
    probs <- log(probs + eps)
    m1 <- matrix(1:N, nrow=N, ncol=2)
    m1[,2] <- data[,ii] + 1
    h1 <- probs[ m1 ]
    h1[ is.na(h1) ] <- 0
    return(h1)
}

