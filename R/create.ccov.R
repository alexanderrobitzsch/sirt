## File Name: create.ccov.R
## File Version: 1.04


#*********************************************************
# auxiliary function for creating a covariance matrix
.create.ccov <- function( cc, data ){
    ccc <- cc$ccov.table
    I <- max( ccc$item1ID, ccc$item2ID )
    ccov.matrix <- matrix( 0, I, I)
    rownames(ccov.matrix) <- colnames(ccov.matrix) <- colnames(data)
    LL <- nrow(ccc)
    for (ll in 1:LL){
        ccov.matrix[ ccc$item1ID[ll], ccc$item2ID[ll] ] <- ccc$ccov[ll]
        ccov.matrix[ ccc$item2ID[ll], ccc$item1ID[ll] ] <-
                        ccov.matrix[ ccc$item1ID[ll], ccc$item2ID[ll] ]
    }
    return( ccov.matrix)
}
#*********************************************************
