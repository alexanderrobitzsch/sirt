## File Name: print.xxirt.R
## File Version: 0.03


################################################################################
print.xxirt <- function(x, ... ){
    object <- x
    d1 <- utils::packageDescription("sirt")
    cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )
#    cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
#    cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")

    cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")

#    modeltype <- object$irtmodel
        cat( object$ic$n , "Cases, " , object$ic$I , "Items, " ,
                object$G , "Group(s)", # "," ,
                "\n")


    #*** parameters
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )
    cat( "  Number of estimated item parameters = " , object$ic$np.item ,
                "\n" )
    cat( "  Number of estimated distribution parameters = " , object$ic$np.Theta ,
                "\n\n" )

    #*** likelihood
    cat( paste0( "Log-Likelihood = " , round( x$loglike ,2 ) , "\n") )
    #*** information criteria
    cat( paste0( "AIC = " , round( x$ic$AIC ,0 ) , "\n") )
    cat( paste0( "BIC = " , round( x$ic$BIC ,0 ) , "\n") )
    invisible(x)
}
