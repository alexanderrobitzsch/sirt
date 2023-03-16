## File Name: summary.rasch.jml.R
## File Version: 1.12
## File Last Change: 2018-12-30


summary.rasch.jml <- function( object, digits=3, ... )
{
    cat("-----------------------------------------------------------------\n")
    cat("Joint Maximum Likelihood Estimation \n")
    cat("Rasch Model \n\n")

    #-- print packages
    packages <- c("sirt")
    sirt_summary_print_package(pack="sirt")
    #-- print R session
    # sirt_summary_print_rsession()
    cat( TAM::tam_rsessinfo() )
    #- print call
    sirt_summary_print_call(CALL=object$CALL)
    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    #-- some informations about model specification
    sirt_summary_cat_label_equal_value("Number of persons", object$N, "\n")
    sirt_summary_cat_label_equal_value("Number of items", object$I, "\n\n")

    sirt_summary_cat_label_equal_value("method", object$method, "\n")
    sirt_summary_cat_label_equal_value("center", object$center, "\n")

    cat("\n-----------------------------------------------------------------\n")
    sirt_summary_cat_label_equal_value("Deviance", object$deviance, "\n", digits=2)
    sirt_summary_cat_label_equal_value("Number of JML iterations", object$iter, "\n", digits=2)
    cat("\n")

    #** person parameters
    cat("Person Parameters \n\n")
    obji <- object$theta_summary
    sirt_summary_print_objects(obji=obji, digits=digits, from=1)

    #** item parameters
    cat("\nItem Parameters \n\n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=1)
}
