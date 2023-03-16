## File Name: mirt_summary.R
## File Version: 0.151


mirt_summary <- function(object, digits=4, file=NULL, ...)
{
    TAM::require_namespace_msg("mirt")

    sirt_osink(file=file)

    is_mg <- inherits(object,"MultipleGroupClass")

    #- print mirt object
    print(object)
    cat("\n\n")

    # print estimated coefficients
    res <- mirt.wrapper.coef(mirt.obj=object)
    cat("*** Item Parameters ***\n\n")
    sirt_summary_print_objects(obji=res$coef, digits=digits, from=2+is_mg)

    cat("\n*** Trait Distribution ***\n\n")
    G <- res$G
    if (G==1){
        sirt_summary_print_objects(obji=res$GroupPars, digits=digits)
    } else {
        for (gg in 1:G){
            groups <- res$groups
            cat(paste0("Group ", groups[gg],":\n"))
            sirt_summary_print_objects(obji=res$GroupPars[[gg]], digits=digits)
            cat("\n")
        }
    }

    sirt_csink(file=file)
}
