## File Name: linking.robust.R
## File Version: 1.263


#*** Robust linking
linking.robust <- function(  itempars )
{
    itempars0 <- itempars
    itempars <- na.omit(itempars)
    pars <- itempars[,2] - itempars[,3]
    items <- paste(itempars[,1])
    names(pars) <- items
    I <- length(pars)
    x <- sort(pars)
    kvec <- seq(1, floor( (I-1)/2  ) )
    KK <- length(kvec)
    se <- meanpars <- rep(NA, KK )
    # define trimming factor
    for (kk in 1L:KK){
        # arrange calculations
        N <- length(x)
        k <- kk
        indkk <- seq( k+1,  N-k,1)
        x0 <- x[ indkk ]
        # compute winsorized mean
        trim.mean <- mean( x0 )
        swk2 <- k * ( x[ k] - trim.mean )^2 + sum( ( x0 - trim.mean )^2 )
                                    + k * ( x[ N - k + 1] - trim.mean )^2
        # standard error
        se.trimmean <- sqrt( swk2 ) / sqrt( (N-2*k) * ( N - 2*k - 1 ) )
        # output
        meanpars[kk] <- trim.mean
        se[kk] <- se.trimmean
    }

    v1 <- paste0('k', 0:KK)
    meanpars <- c( mean(x), meanpars )
    se <- c( stats::sd(x) / sqrt(I), se )
    names(meanpars) <- v1
    names(se) <- v1

    # arrange output
    res1 <- list()
    res1$ind.kopt <- ind.kopt <- which.min( se )
    res1$kopt <- kvec[ ind.kopt ] - 1
    res1$meanpars.kopt <- meanpars[ ind.kopt ]
    res1$se.kopt <- se[ ind.kopt ]
    res1$meanpars <- meanpars
    res1$se <- se
    res1$sd <- stats::sd(x)
    res1$mad <- stats::mad(x)
    res1$k.robust <- c(0,kvec)
    res1$I <- I
    res1$itempars <- itempars0
    class(res1) <- 'linking.robust'
    return(res1)
}

