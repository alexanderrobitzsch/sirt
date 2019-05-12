## File Name: summary.gom.em.R
## File Version: 0.174


#--- summary for gom object
summary.gom <- function( object, file=NULL, ...)
{

    # open sink
    sirt_osink( file=file )

    cat("-----------------------------------------------------------------\n")
    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat("  Function 'gom.em' \n")

    if (object$model=="GOM"){
        cat("   Discrete Grade of Membership Model\n\n")
    }
    if (object$model=="GOMnormal"){
        cat("   Grade of Membership Model Based on Multivariate Normal Distribution\n\n")
    }
    if (object$model=="GOMRasch"){
        cat("   Rasch Grade of Membership Model\n\n")
    }
    modeltype <- object$irtmodel
    cat( "   ", object$ic$n, "Cases, ", ncol(object$dat2), "Items, ",
                object$K, "Classes", ",", object$TP, "Discrete Integration Points\n")

    cat("-----------------------------------------------------------------\n")
    cat( "Number of iterations", "=", object$iter, "\n" )
    cat( "Deviance", "=", round( object$deviance, 2 ), " | " )
    cat( "Log Likelihood", "=", round( -object$deviance/2, 2 ), "\n" )
    cat( "Number of persons", "=", object$ic$n, "\n" )

    cat( "Number of estimated parameters", "=", object$ic$np, "\n" )
    cat( "  Number of estimated item parameters", "=", object$ic$np.item, "\n" )
    cat( "  Number of estimated distribution parameters", "=", object$ic$np.trait, "\n\n" )

    #--- information criteria
    rm_summary_information_criteria(object=object)

    #--- trait distribution
    if (object$model %in% c("GOMRasch","GOMnormal") ){
        cat("-----------------------------------------------------------------\n")
        if (object$model=="GOMRasch"){
            cat("Trait Distribution (Location, Variability)\n")
        }
        if (object$model=="GOMnormal"){
            cat("Underlying multivariate normal distribution\n")
        }
        cat( " Means: ", round( object$mu, 3 ), "\n")
        cat( " Standard deviations: ", round( sqrt(diag(object$Sigma)), 3 ), "\n")
        if (object$model=="GOMRasch"){
            c1 <- stats::cov2cor(object$Sigma)
            cat( " Correlation ", round( c1[lower.tri(c1)], 3 ), "\n")
            cat("EAP Reliability: ", round(object$EAP.rel,3), "\n")
        }
        if (object$model=="GOMnormal"){
            c1 <- object$Sigma
            cat("Covariance matrix:\n")
            print(round(c1,3))
        }
    }

    #--- membership descriptives
    if ( ! ( object$plmat) ){
        cat("-----------------------------------------------------------------\n")
        cat("Membership Function Descriptives \n")
        obji <- object$classdesc
        sirt_summary_print_objects(obji=obji, digits=3, rownames_null=FALSE)
    }

    #-- item parameters
    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=3, from=2)

    # close sink
    sirt_csink( file=file )
}

