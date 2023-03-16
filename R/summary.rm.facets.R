## File Name: summary.rm.facets.R
## File Version: 0.26


# Summary for rm.facets object
summary.rm.facets <- function( object, file=NULL, ... ){

    # open sink for a file
    sirt_osink( file=file  )

    cat("-----------------------------------------------------------------\n")
    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    #-- number of parameters
    cat("Rater Facet Model with Item/Rater Intercepts and Slopes \n")
    cat("-----------------------------------------------------------------\n")
    cat( "Number of iterations","=", object$iter, "\n" )
    cat( "Deviance","=", round( object$deviance, 2 ), " | " )
    cat( "Log Likelihood","=", round( -object$deviance/2, 2 ), "\n" )
    cat( "Number of persons","=", object$ic$n, "\n" )
    cat( "Number of items","=", object$ic$VV, "\n" )
    cat( "Number of raters","=", object$ic$RR, "\n" )
    cat( "Number of estimated parameters","=", object$ic$np, "\n" )
    cat( "   Number of item parameters","=", object$ic$np.item, "\n" )
    cat( "   Number of rater parameters","=", object$ic$np.rater, "\n" )
    cat( "   Number of distribution parameters","=", object$ic$np.trait, "\n\n" )

    #--- information criteria
    rm_summary_information_criteria(object=object)

    #--- trait distribution
    rm_summary_trait_distribution(object=object)

    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=3, from=2)

    cat("-----------------------------------------------------------------\n")
    cat("Rater Parameters \n")
    obji <- object$rater
    sirt_summary_print_objects(obji=obji, digits=3, from=2)

    sirt_csink(file=file)
}
