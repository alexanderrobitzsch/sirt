## File Name: summary.rm.sdt.R
## File Version: 1.15


#*******************************************************
# Summary for rm.facets object                         *
summary.rm.sdt <- function( object , file=NULL, ...){

    # open sink for a file
    sirt_osink( file=file  )

    cat("-----------------------------------------------------------------\n")
    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)


    cat("Hierarchical Rater Model: Signal Detection Model \n")
    cat("-----------------------------------------------------------------\n")
    cat( "Iteration with minimal deviance =" , object$iter_opt , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )
    cat( "Number of persons = " , object$ic$n , "\n" )
    cat( "Number of items   = " , object$ic$VV , "\n" )
    cat( "Number of raters  = " , object$ic$RR , "\n" )
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )
    cat( "  Distribution parameters  = " , object$ic$np.skill , "\n" )
    cat( "  Item parameters  = " , object$ic$np.item , "\n" )
    cat( "  Rater parameters = " , object$ic$np.rater , "\n" )
    cat("\n")

    #--- information criteria
    rm_summary_information_criteria(object=object)

    cat("-----------------------------------------------------------------\n")
    cat( "Trait Distribution\n" ,
              "Mean=" , round(object$mu,3)  , " SD=" , round( object$sigma , 3) )

    if ( object$skillspace == "discrete" ){
        cat("\n\nDiscrete Skill Distribution\n")
        obji <- object$skill.distribution
        sirt_summary_print_objects(obji=obji, digits=4, from=1)
    }


    cat( "\n\nEAP Reliability = ")
    cat(round( object$EAP.rel,3 ) )
    cat( "\n")

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
#*******************************************************



