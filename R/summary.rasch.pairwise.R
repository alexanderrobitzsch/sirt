## File Name: summary.rasch.pairwise.R
## File Version: 0.09
## File Last Change: 2018-12-30



# Summary for rasch.pairwise objects
summary.rasch.pairwise <- function( object, digits=3, file=NULL, ...)
{
    # open sink for a file
    sirt_osink(file=file)

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    #- package and R session
    sirt_summary_print_package_rsession(pack="sirt")

    #- print call
    sirt_summary_print_call(CALL=object$CALL)

    cat(paste0("Function '", object$fct, "'"), "\n\n")

    #-- print computation time
    sirt_summary_print_computation_time_s1(object=object)

    cat(display_string)
    cat("Pairwise likelihood estimation \n")
    cat("Rasch Model \n")
    cat(display_string)
    cat("Item Parameters \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    sirt_csink(file=file)
}

