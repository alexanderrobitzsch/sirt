## File Name: summary.invariance_alignment_constraints.R
## File Version: 0.08



summary.invariance_alignment_constraints <- function( object, digits=3, file=NULL, ...)
{
    # open sink for a file
    sirt_osink(file=file)

    display_string <- sirt_summary_print_display(symbol="-", len=65)
    cat(display_string)

    cat("Invariance Alignment with Post-Hoc Item Parameter Constraints\n\n")
    invariance_alignment_summary_optimization(object=object$model, digits=digits)

    cat(display_string)
    cat("Alignment Results Lambda Parameters \n\n")
    invariance_alignment_constraints_summary_print_item_summary(
            x=object$lambda_list, digits=digits)

    cat(display_string)
    cat("Alignment Results Nu Parameters \n\n")
    invariance_alignment_constraints_summary_print_item_summary(
            x=object$nu_list, digits=digits)

    sirt_csink(file=file)
}
