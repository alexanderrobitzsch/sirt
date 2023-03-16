## File Name: invariance_alignment_constraints_summary_print_item_summary.R
## File Version: 0.08
## File Last Change: 2019-03-06

invariance_alignment_constraints_summary_print_item_summary <- function(x, digits)
{
    cat("Parameter tolerance value", "=", x$parm_tol, "\n")
    cat("Total number of items", "=", x$N_total, "\n")
    cat("Number of unique item parameters", "=", x$N_parm_all, "\n")
    cat("Percentage of noninvariance item parameters", "=",
                round(x$prop_noninvariance,1),"%\n")

    cat("\nUnique item parameters per item\n")
    print(x$N_unique_parm_items)

    cat("\nUnique item parameters per group\n")
    print(x$N_unique_parm_groups)

    cat("\nJoint item parameters\n")
    obji <- x$parm_joint
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat("\nEstimated item parameters\n")
    obji <- x$parm_est
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)

    cat("\nEstimated DIF effects\n")
    obji <- x$parm_dif
    sirt_summary_print_objects(obji=obji, digits=digits, from=1, rownames_null=FALSE)
}
