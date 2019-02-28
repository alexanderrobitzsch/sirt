## File Name: summary.lsem.permutationTest.R
## File Version: 0.292
######################################################
summary.lsem.permutationTest <- function( object, file=NULL, digits=3, ... )
{

    # open sink for a file
    sirt_osink( file=file )

    cat("-----------------------------------------------------------------\n")
    cat("Permutation Test for Local Structural Equation Model \n\n")

    packages <- c("sirt", "lavaan")
    if (object$use_lavaan_survey){
        packages <- c(packages, "lavaan.survey")
    }
    sirt_summary_print_packages(packages=packages)

    cat("\nFunction 'sirt::lsem.permutationTest' \n\n")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    cat( "Date of Analysis:", paste( object$s2 ), "\n" )
    cat("Computation Time:", print(object$s2 - object$s1), "\n\n")

    cat( "Number of permutations", "=", object$B, "\n")
    cat( "Percentage of non-converged datasets",  "=", round(object$nonconverged_rate, digits=digits), "\n")
    cat( paste0( "Number of observations",  "=", round(object$N,digits) ), "\n")
    cat( paste0( "Bandwidth factor",  "=", round(object$h,digits) ), "\n")
    cat( paste0( "Bandwidth",  "=", round(object$bw,digits) ), "\n")
    cat( paste0( "Number of focal points for moderator",  "=",
                        length(object$moderator.grid ) ), "\n")

    cat("\nlavaan Model\n")
    cat(object$lavmodel)

    cat("\n\n")

    cat("Global Test Statistics\n\n")
    obji <- object$teststat
    VV <- ncol(obji)
    for (vv in 2:VV){
        obji[,vv] <- round( obji[,vv], digits )
    }
    print(obji)

    cat("\n")
    cat("Pointwise Test Statistics\n\n")
    obji <- object$parameters_pointwise_test
    vars <- c("est", "p")
    for (vv in vars){
        obji[,vv] <- round( obji[,vv], digits )
    }
    obji <- obji[, c("par", "parindex", "moderator", "est", "p") ]
    rownames(obji) <- NULL
    print(obji)

    # close file
    sirt_csink(file=file)

}
######################################################
