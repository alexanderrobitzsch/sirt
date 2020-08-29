## File Name: detect.index.R
## File Version: 0.34



detect.index <- function( ccovtable, itemcluster )
{
    ccovtable <- ccovtable$ccov.table
    ccovtable$delta <- ifelse( itemcluster[ ccovtable$item1ID ]==
                            itemcluster[ ccovtable$item2ID ], 1, -1 )

    #--- calculate unweighted and weighted indizes
    ccov <- ccovtable$ccov            # conditional covariance
    delta <- ccovtable$delta        # indicator for partition
    N <- ccovtable$N
    sqrt_N <- sqrt(N)
    sign_ccov <- sign(ccov)
    abs_ccov <- abs(ccov)
    # number of parameters
    np <- 5
    parnames <- weighted.indizes <- indizes <- rep(NA,np)

    #--- DETECT
    ii <- 1
    indizes[ii] <- 100*mean(ccov*delta)
    weighted.indizes[ii] <- 100*stats::weighted.mean( ccov * delta, sqrt_N )
    parnames[ii] <- "DETECT"
    #--- ASSI
    ii <- 2
    indizes[ii] <- mean( sign_ccov * delta )
    weighted.indizes[ii] <- stats::weighted.mean( sign_ccov * delta, sqrt_N )
    parnames[ii] <- "ASSI"
    #--- RATIO
    ii <- 3
    indizes[ii] <- sum( ccov * delta ) / sum( abs_ccov )
    weighted.indizes[ii] <- sum( ccov * delta * sqrt_N ) / sum( abs_ccov * sqrt_N )
    parnames[ii] <- "RATIO"
    #--- MADCOV
    ii <- 4
    indizes[ii] <- 100 * mean( abs_ccov )
    weighted.indizes[ii] <- 100* stats::weighted.mean( abs_ccov, sqrt_N )
    parnames[ii] <- "MADCOV100"
    #--- MCOV
    ii <- 5
    indizes[ii] <- 100 * mean( ccov )
    weighted.indizes[ii] <- 100* stats::weighted.mean( ccov, sqrt_N )
    parnames[ii] <- "MCOV100"

    #--- output
    res <- data.frame( "unweighted"=indizes, "weighted"=weighted.indizes )
    rownames(res) <- parnames
    return(res)
}


