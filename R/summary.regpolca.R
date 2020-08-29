## File Name: summary.regpolca.R
## File Version: 0.09


#--- summary for regpolca object
summary.regpolca <- function( object, digits=3, file=NULL, ...)
{
    # open sink
    sirt_osink( file=file )
    len_disp <- 66

    # print summary based on xxirt object
    res <- xxirt_summary_parts(object=object, digits=digits, len_disp=len_disp)

    #-------- output regpolca specific -------
    sirt_display_function(length=len_disp)
    cat("Regularized polytomous latent class analysis (regpolca)\n\n")
    cat( "Regularization type","=", object$regular_type, "\n" )
    cat( "Group regularization","=", object$regular_grouped, "\n" )
    cat( "Regularization parameter(s)","=", object$regular_lam, "\n" )
    cat( "Number of regularized parameters","=", object$n_reg, "\n\n" )

    #* class probabilities
    sirt_display_function(length=len_disp)
    cat("Class Probabilities\n")
    obji <- object$probs_Theta
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    #* item response probabilities
    sirt_display_function(length=len_disp)
    cat("Item Response Functions \n")
    obji <- object$item
    sirt_summary_print_objects(obji=obji, digits=digits, from=3-object$lca_dich,
                rownames_null=TRUE)

    # close sink
    sirt_csink( file=file )
}
