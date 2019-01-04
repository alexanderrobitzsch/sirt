## File Name: summary.rasch.copula2.R
## File Version: 0.204


#** Summary for rasch.copula object
summary.rasch.copula2 <- function( object, file=NULL, digits=3, ... )
{
    is_copula2 <- class(object)=="rasch.copula2"
    is_copula3 <- class(object)=="rasch.copula3"

    # open sink for a file
    sirt_osink(file=file)

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat("Marginal Maximum Likelihood Estimation \n")
    cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
    cat( paste0("Function '", class(object), "'\n"))
    cat("alpha1","=", round(object$alpha1, digits)," alpha2", "=",
                        round(object$alpha2, digits), "\n")

    cat(display_string)
    cat( "Deviance", "=", round( object$deviance, 2 ), "\n" )
    cat( "Number of persons", "=", object$ic$n, " (", nrow(object$pattern), " Response Patterns)\n" )
    cat( "Number of estimated parameters", "=", object$ic$np, "\n" )
    cat( "Number of iterations", "=", object$iter, "\n" )
    cat( "AIC", "=", round( object$ic$AIC, 2 ), " | penalty", "=",
                    round( object$ic$AIC - object$ic$deviance,2 ), "\n" )
    cat( "AICc", "=", round( object$ic$AICc, 2 ), " | penalty", "=", round( object$ic$AICc - object$ic$deviance,2 ) )
    cat(" (bias corrected AIC)\n" )
    cat( "BIC", "=", round( object$ic$BIC, 2 ), " | penalty", "=", round( object$ic$BIC - object$ic$deviance,2 ), "\n" )
    cat( "CAIC", "=", round( object$ic$CAIC, 2 )," | penalty", "=", round( object$ic$CAIC - object$ic$deviance,2 ) )
    cat( " (consistent AIC) \n\n" )

    #--- trait distribution
    if (is_copula2){
        cat( "Trait Distribution (", length(object$theta.k), " Knots )\n",
                "Mean", "=", 0, " SD", "=", 1, "\n")
        cat(paste("\nEAP Reliability:", round( object$EAP.Rel, digits)),"\n\n")
    }
    if (is_copula3){
        cat( "Trait Distribution (", length(object$theta.k), " Knots )\n" )
        cat("\nMean Vector\n")
        print( round( object$mu, digits))
        cat("\nCorrelation Matrix\n")
        print( round( object$sigma, digits))
    }

    # item paramters
    cat(display_string)
    cat("Item Parameter \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=TRUE)

    # dependency parameters
    cat("\nDependency parameters\n")
    obji <- object$summary.delta
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=TRUE)

    sirt_csink(file=file)
}
