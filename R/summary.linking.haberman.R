## File Name: summary.linking.haberman.R
## File Version: 0.284
## File Last Change: 2019-01-09


summary.linking.haberman <- function( object, digits=3, file=NULL, ... )
{
    # open sink
    sirt_osink(file=file)

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    cat("\n")
    cat(object$description, "")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object$time)

    cat(display_string)
    cat("Estimation information item slopes\n")
    linking_haberman_summary_estimation_information(res_opt=object$res_opt_slopes)

    cat(display_string)
    cat("Estimation information item intercepts\n")
    linking_haberman_summary_estimation_information(res_opt=object$res_opt_intercepts)

    cat(display_string)
    cat("Transformation parameters (Haberman linking)\n")
    obji <- object$transf.pars
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    cat("\nLinear transformation for item parameters a and b\n")
    obji <- object$transf.itempars
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    cat("\nLinear transformation for person parameters theta\n")
    obji <- object$transf.personpars
    sirt_summary_print_objects(obji=obji, digits=digits, from=2, rownames_null=FALSE)

    cat("\n")
    cat(display_string)
    cat("R-Squared Measures of Invariance (all items and unweighted)\n")
    obji <- object$es.invariance
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat("\nR-Squared Measures of Invariance (weighted)\n")
    obji <- object$es.robust
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat("\n-----------------------------------------------------------------\n")
    if ( ! object$linking_slopes ){
        cat("Estimated DIF effects in logarithms of item slopes \n")
        obji <- object$a.resid[ object$selitems, ]
        sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)
        cat("\n")
    }
    cat("Estimated DIF effects of item intercepts \n")
    obji <- object$b.resid[ object$selitems, ]
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat("\n-----------------------------------------------------------------\n")
    digits2 <- 1
    if ( ! object$linking_slopes ){
        cat("Used items in linking item slopes \n")
        obji <- object$a.item_stat
        sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)
        cat("\n")
    }
    cat("Used items in linking item intercepts \n")
    obji <- object$b.item_stat
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    # close sink
    sirt_csink( file=file )
}
