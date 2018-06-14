## File Name: summary.linking.haberman.R
## File Version: 0.16


summary.linking.haberman <- function( object , digits = 3 , file=NULL , ... ){

    # open sink
    sirt_osink( file = file)

    cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
    cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )

    cat( object$description , "\n\n")

    cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")

    cat( "Date of Analysis:"  , paste( object$time ) , "\n" )
    cat("\n")

    cat("-----------------------------------------------------------------\n")
    cat("Transformation parameters (Haberman linking)\n")
    obji <- object$transf.pars
    obji[,-1] <- round( obji[,-1] , digits )
    print( obji )
    cat("\nLinear transformation for item parameters a and b\n")
    obji <- object$transf.itempars
    obji[,-1] <- round( obji[,-1] , digits )
    print( obji )
    cat("\nLinear transformation for person parameters theta\n")
    obji <- object$transf.personpars
    obji[,-1] <- round( obji[,-1] , digits )
    print( obji )

    cat("\n-----------------------------------------------------------------\n")
    cat("R-Squared Measures of Invariance (all items and unweighted)\n")
    obji <- object$es.invariance
    obji <- round( obji , digits )
    print( obji )
    cat("\nR-Squared Measures of Invariance (weighted)\n")
    obji <- object$es.robust
    obji <- round( obji , digits )
    print( obji )

    cat("\n-----------------------------------------------------------------\n")
    if ( ! object$linking_slopes ){
        cat("Estimated DIF effects in logarithms of item slopes \n")
        obji <- object$a.resid[ object$selitems , ]
        obji <- round( obji , digits )
        print( obji )
        cat("\n")
    }
    cat("Estimated DIF effects of item intercepts \n")
    obji <- object$b.resid[ object$selitems , ]
    obji <- round( obji , digits )
    print( obji )

    cat("\n-----------------------------------------------------------------\n")
    digits2 <- 1
    if ( ! object$linking_slopes ){
        cat("Used items in linking item slopes \n")
        obji <- object$a.item_stat
        obji[,-1] <- round( obji[,-1] , digits2 )
        print( obji )
        cat("\n")
    }
    cat("Used items in linking item intercepts \n")
    obji <- object$b.item_stat
    obji[,-1] <- round( obji[,-1] , digits2 )
    print( obji )

    # close sink
    sirt_csink( file = file )

}
