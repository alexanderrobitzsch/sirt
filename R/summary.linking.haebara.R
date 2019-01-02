## File Name: summary.linking.haebara.R
## File Version: 0.11



summary.linking.haebara <- function( object, digits=3, file=NULL, ...)
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
    linking_haebara_summary_optimization(object=object, digits=digits)

    cat(display_string)
    cat("Number of items per study\n")
    print(object$numb_items)

    cat(display_string)
    cat("Group Means and Standard Deviations ")
    cat(paste0("(center=", object$center,")\n"))
    obji <- object$pars
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Joint item parameters \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=2)

    cat(display_string)
    cat("Estimated DIF effects of item slopes \n")
    obji <- object$a.resid
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat(display_string)
    cat("Estimated DIF effects of item intercepts \n")
    obji <- object$b.resid
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    sirt_csink(file=file)
}
