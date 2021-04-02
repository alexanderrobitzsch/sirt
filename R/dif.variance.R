## File Name: dif.variance.R
## File Version: 0.193



#*** Routine for calculating DIF variance   (Camilli & Penfield, 1997)      #
dif.variance <- function( dif, se.dif, items=paste("item",1:length(dif),sep="") )
{
    # calculating the weights
    wi <-  1/se.dif^2
    # mean DIF
    md <- mean(dif)
    # weighted tau2 estimator
    # formula (5) of PENFIELD & ALGINA (2006)
    weighted.tauq <- sum( wi^2 * ( dif - md )^2 - wi  ) / sum( wi^2 )
    # unweighted tau2 estimator
    # formula (5) of PENFIELD & ALGINA (2006)
    unweighted.tauq <- sum( ( dif - md )^2 - se.dif^2 ) / ( length(items) - 1 )
    # calculation of variances v_i
    vi <- weighted.tauq + se.dif^2
    weighted.tauq[ weighted.tauq < 0 ] <- 0
    unweighted.tauq[ unweighted.tauq < 0 ] <- 0
    # Empirical Bayes DIF estimate
    lambda.i <- weighted.tauq / vi
    eb.dif <- lambda.i * ( dif - md ) + ( 1 - lambda.i) *md

    #* output
    res <- list( weighted.DIFSD=sqrt(weighted.tauq),
                unweighted.DIFSD=sqrt(unweighted.tauq),
                mean.se.dif=sqrt( mean( se.dif^2 ) ), eb.dif=eb.dif  )
    return(res)
}

