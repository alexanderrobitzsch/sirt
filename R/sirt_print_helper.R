## File Name: sirt_print_helper.R
## File Version: 2.145
## File Last Change: 2019-01-04


# This is an auxiliary function which helps for printing some progress
sirt_print_helper <- function( object, digits)
{
    ow <- options()$warn
    if ( length(dim(object))==2 ){
        options( warn=-1 )
        if ( nrow(object) >=1 ){
            g1a <- apply( object, 2, as.numeric )
        } else {
            g1a <- object
        }
        g1a <- matrix(g1a, nrow=nrow(object), ncol=ncol(object))
        colnames(g1a) <- colnames(object)
        g1 <- colMeans( g1a )
        g1 <- which( ! is.na( g1 ) )
        options( warn=ow )
        object1 <- object
        object1[, g1 ] <- round( object1[, g1 ], digits )
        print( object1 )
    } else {
        print( round( object, digits ) )
    }
}


.pr <- sirt_print_helper
