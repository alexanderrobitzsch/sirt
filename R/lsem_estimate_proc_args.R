## File Name: lsem_estimate_proc_args.R
## File Version: 0.14

lsem_estimate_proc_args <- function(lavaan.args, sufficient_statistics,
    pseudo_weights, lavmodel, data, use_lavaan_survey, est_joint=FALSE,
    par_invariant=NULL, par_linear=NULL, par_quadratic=NULL,
    partable_joint=NULL)
{
    use_pseudo_weights <- pseudo_weights > 0
    if ( sufficient_statistics | use_pseudo_weights ){
        use_lavaan_survey <- FALSE
    }
    if (use_pseudo_weights){
        use_lavaan_survey <- FALSE
    }

    #- variables in model
    partable <- sirt_import_lavaan_lavaanify(model=lavmodel)
    variables_model <- intersect( union( partable$lhs, partable$rhs ), colnames(data) )

    #- joint estimation
    par_vec <- union(union(par_invariant, par_linear), par_quadratic)
    if (length(par_vec)>0){
        cat("Use joint estimation\n")
        est_joint <- TRUE
        sufficient_statistics <- TRUE
    }
    if ( !is.null(partable_joint) ){
        est_joint <- TRUE
    }
    if (est_joint){
        sufficient_statistics <- TRUE
    }


    #- ordered variables
    variables_ordered <- NULL
    if ("ordered" %in% names(lavaan.args)){
        sufficient_statistics <- FALSE
        variables_ordered <- lavaan.args$ordered
    }

    data1 <- data[,variables_model,drop=FALSE]
    data_ordered <- rep(FALSE, ncol(data1))
    names(data_ordered) <- colnames(data1)
    NV <- ncol(data1)
    for (vv in 1:NV){
        data_ordered[vv] <- is.factor(data1[,vv])
    }
    if (any(data_ordered)){
        sufficient_statistics <- FALSE
        variables_ordered <- union(variables_ordered, names(data_ordered)[data_ordered] )
    }
    if (est_joint & ( ! sufficient_statistics)){
        stop("Cannot perform joint estimation for non-continuous data.")
    }

    #-- output
    res <- list(sufficient_statistics=sufficient_statistics, use_lavaan_survey=use_lavaan_survey,
                variables_model=variables_model, use_pseudo_weights=use_pseudo_weights,
                variables_ordered=variables_ordered, est_joint=est_joint)
    return(res)
}
