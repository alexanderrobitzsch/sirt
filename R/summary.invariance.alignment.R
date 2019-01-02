## File Name: summary.invariance.alignment.R
## File Version: 0.282



summary.invariance.alignment <- function( object, digits=3, file=NULL, ...)
{
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

    #-- informations about optimization
    cat(display_string)
    invariance_alignment_summary_optimization(object=object, digits=digits)

    cat(display_string)
    cat("Number of items per group\n")
    print(object$numb_items)

    cat(display_string)
    cat("Effect Sizes of Approximate Invariance \n")
    obji <- object$es.invariance
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Group Means and Standard Deviations ")
    cat(paste0("(center=", object$center,")\n"))
    obji <- object$pars
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Summary Aligned Item Parameters Lambda \n")
    obji <- object$itempars.aligned
    sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                    rownames_null=FALSE, grep_string="lambda")

    cat(display_string)
    cat("Summary Aligned Item Parameters Nu \n")
    obji <- object$itempars.aligned
    sirt_summary_print_objects(obji=obji, digits=digits, from=1,
                    rownames_null=FALSE, grep_string="nu")

    cat(display_string)
    cat("Aligned Lambda Parameters \n")
    obji <- object$lambda.aligned
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Aligned Nu Parameters \n")
    obji <- object$nu.aligned
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Summary Absolute Residuals Loadings Lambda \n")
    obji <- as.vector(abs(object$lambda.resid))
    sirt_summary_print_vector_summary(obji=obji, digits=digits)

    cat(display_string)
    cat("Summary Absolute Residuals Intercepts Nu \n")
    obji <- as.vector(abs(object$nu.resid))
    sirt_summary_print_vector_summary(obji=obji, digits=digits)

    sirt_csink(file=file)
}
