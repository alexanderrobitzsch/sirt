## File Name: lsem_estimate_proc_args.R
## File Version: 0.402

lsem_estimate_proc_args <- function(lavaan.args, sufficient_statistics,
    pseudo_weights, lavmodel, data, use_lavaan_survey, est_joint=FALSE,
    par_invariant=NULL, par_linear=NULL, par_quadratic=NULL,
    partable_joint=NULL, se=NULL, G=NULL, moderator.grid=NULL, verbose=TRUE)
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
    has_meanstructure <- "~" %in% partable$op

    #- joint estimation
    par_vec <- union(union(par_invariant, par_linear), par_quadratic)
    if (length(par_vec)>0){
        if (verbose){
            cat("Use joint estimation\n")
        }
        est_joint <- TRUE
        # sufficient_statistics <- TRUE
    }
    if ( ! is.null(partable_joint) ){
        est_joint <- TRUE
    }
    if (est_joint){
        sufficient_statistics <- TRUE
    }

    lavaan_args_names <- names(lavaan.args)
    if ( "missing" %in% lavaan_args_names){
        if ( lavaan.args[["missing"]]=="fiml" ){
            sufficient_statistics <- FALSE
        }
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
    some_ordinal <- FALSE
    if (any(data_ordered)){
        sufficient_statistics <- FALSE
        variables_ordered <- union(variables_ordered, names(data_ordered)[data_ordered] )
        some_ordinal <- TRUE
    }
    if (est_joint & ( ! sufficient_statistics)){
        if (pseudo_weights==0){
            G <- length(moderator.grid)
            pseudo_weights <- 10*G*nrow(data)
            use_pseudo_weights <- TRUE
        }
        if (some_ordinal){
            stop("Cannot perform joint estimation for non-continuous data.")
        }
    }
    if (is.null(se)){
        if (est_joint){
            se <- "none"
        } else {
            se <- "standard"

        }
    }
    compute_se <- se !="none"

    #-- output
    res <- list(sufficient_statistics=sufficient_statistics,
                use_lavaan_survey=use_lavaan_survey,
                variables_model=variables_model, use_pseudo_weights=use_pseudo_weights,
                variables_ordered=variables_ordered, est_joint=est_joint,
                partable=partable, has_meanstructure=has_meanstructure, se=se,
                compute_se=compute_se, pseudo_weights=pseudo_weights,
                some_ordinal=some_ordinal)
    return(res)
}
