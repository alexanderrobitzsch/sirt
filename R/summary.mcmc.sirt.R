## File Name: summary.mcmc.sirt.R
## File Version: 1.33
## File Last Change: 2018-12-30


# summary for MCMC item analysis in sirt
summary.mcmc.sirt <- function( object, digits=3, file=NULL, ... ){

    # open sink
    sirt_osink(file=file)

    cat("-----------------------------------------------------------------\n")

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- computation time
    sirt_summary_print_computation_time(object=object)

    #- information criteria
    mcmc_summary_print_information_criteria(object=object)

    #- sample characteristics
    cat( "Number of persons=", nrow(object$dat), "\n" )
    if ( object$model=="2pno.ml"){
        cat( "Number of groups=", object$ic$G, "\n")
        cat( "  Group sizes: M=", round(object$ic$M.n,3),
            " | SD=", round(object$ic$SD.n,3), "\n")
    }
    cat( "Number of items=", ncol(object$dat), "\n" )

    cat( "\nEAP Reliability=")
    cat(round( object$EAP.rel,3 ) )
    cat( "\n")

    cat("-----------------------------------------------------------------\n")
    cat("Item Parameters \n")
    obji <- object$summary.mcmcobj
    obji <- obji[ obji$parameter !="deviance", ]
    vars <- c("parameter", "Mean", "SD", "MAP", "Rhat", "effSize", "Q5", "Q95" )
    obji <- obji[, vars ]
    digits_vec <- sirt_vector_with_names(value=digits, names=vars)
    digits_vec["Rhat"] <- 2
    digits_vec["effSize"] <- 1
    sirt_summary_print_objects(obji=obji, digits=digits_vec, rownames_null=TRUE)

    #- close sink
    sirt_csink(file=file)
}
