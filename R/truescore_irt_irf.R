## File Name: truescore_irt_irf.R
## File Version: 0.03


#---- true score IRT - item response function
truescore_irt_irf <- function( A, B, c, d, theta )
{
    TP <- length(theta)
    truescore <- rep(0, length(theta))
    # extract maximum item categories from B parameters
    nB <- ncol(B)
    maxK <- nB - rowSums( is.na( B ))
    I <- nrow(B)
    scoreM <- matrix( 1:max(maxK), nrow=TP, ncol=max(maxK), byrow=TRUE )
    for (ii in 1:I){
        prob.ii <- matrix( NA, TP, maxK[ii] )
        for (kk in seq(1, maxK[ii] ) ){
            prob.ii[,kk] <- exp( theta * A[ii,kk] + B[ii,kk] )
        }
        prob.ii <- prob.ii / ( rowSums( prob.ii ) + 1 )
        if ( ( maxK[ii]==1 ) & ( abs(c[ii])+abs(1-d[ii]) > 0 ) ){
            prob.ii <- c[ii] + ( d[ii] - c[ii] ) * prob.ii
        }
        truescore <- truescore + rowSums( prob.ii * scoreM[, seq( 1, maxK[ii] ),
                            drop=FALSE] )
    }
    return(truescore)
}

