## File Name: summary.lsem.R
## File Version: 0.33

#############################################
# summary lsem
summary.lsem <- function( object, file=NULL, digits=3, ... )
{
    # open sink for a file
    sirt_osink( file=file  )

    cat("-----------------------------------------------------------------\n")
    cat("Local Structural Equation Model \n\n")

    #-- print packages
    packages <- c("sirt", "lavaan", "lavaan.survey")
    sirt_summary_print_packages(packages=packages)
    #-- print R session
    cat("\n")
    sirt_summary_print_rsession()

    cat(paste0("Function 'sirt::lsem.estimate', type='", object$type,"'"), "\n\n")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat( paste0( "Number of observations=", round(object$N,digits) ), "\n")
    if ( object$type=="LSEM"){
        cat( paste0( "Bandwidth factor=", round(object$h,digits) ), "\n")
        cat( paste0( "Bandwidth=", round(object$bw,digits) ), "\n")
        cat( paste0( "Number of focal points for moderator=",
                            length(object$moderator.grid ) ), "\n")
    }

    if ( object$type=="MGM"){
        cat( paste0( "Number of groups for moderator=",
                            length(object$moderator.grid ) ), "\n")
    }

    cat("\nlavaan Model\n")
    cat(object$lavmodel)

    cat("\n\n")
    cat("Parameter Estimate Summary\n\n")
    obji <- object$parameters_summary
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)

    cat("\n")
    cat("Distribution of Moderator: Density and Effective Sample Size\n\n")
    obji <- object$moderator.density
    sirt_summary_print_objects(obji=obji, digits=digits, from=1)

    cat("\n")
    obji <- object$moderator.stat
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)

    # close file
    sirt_csink(file)

}
#############################################
