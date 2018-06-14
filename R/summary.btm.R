## File Name: summary.btm.R
## File Version: 0.17



#*******************************************************
# summary.btm
summary.btm <- function( object , file=NULL , digits=4,... )
{
    # open sink
    sirt_osink( file = file )

    res <- object

    cat("------------------------------------------------------------\n")
    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)


    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    cat("Bradley-Terry Model with Ties and Home Advantage Parameters\n")

    cat("------------------------------------------------------------\n")
    cat( "Number of iterations =" , object$iter , "\n" )
    #    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    #    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )
    cat( "Number of individuals = " , object$ic$n , "\n" )
    cat( "Number of pairwise comparisons = " , object$ic$D , "\n" )
    #    cat( "Number of estimated parameters = " , object$ic$np , "\n" )

    cat("------------------------------------------------------------\n")
    cat("Ties and Home advantage parameters\n")
    obji <- res$pars
    sirt_summary_print_objects(obji=obji, digits=digits, from=3)

    cat("------------------------------------------------------------\n")
    cat("Summary of individual effects parameters\n")
    obji <- res$summary.effects
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=TRUE)

    cat("------------------------------------------------------------\n")
    cat("MLE reliability (separation reliability)\n")
    cat(paste0("MLE Rel = " , round( res$mle.rel , digits )  , "\n") )
    cat(paste0("Separation index = " , round( res$sepG , digits )  , "\n") )

    cat("------------------------------------------------------------\n")
    cat("Individual effects parameters\n")
    obji <- res$effects
    sirt_summary_print_objects(obji=obji, digits=digits, from=3)

    # close sink
    sirt_csink( file = file )
}
#*******************************************************



