## File Name: invariance_alignment_cfa_config_estimate.R
## File Version: 0.191

invariance_alignment_cfa_config_estimate <- function(dat_gg, N, weights_gg=NULL,
    model="2PM", ...)
{
    is_data <- sirt_is_data(dat=dat_gg)
    #-- create lavaan model
    if (is_data){
        I_gg <- ncol(dat_gg)
        items_gg <- colnames(dat_gg)
    } else {
        mu <- dat_gg[[1]]
        Sigma <- dat_gg[[2]]
        I_gg <- length(mu)
        items_gg <- names(mu)
        if (is.null(items_gg)){
            items_gg <- paste0("I",1:I)
        }
        names(mu) <- items_gg
        rownames(Sigma) <- items_gg
        colnames(Sigma) <- items_gg
    }
    lavmodel <- invariance_alignment_cfa_config_estimate_define_lavaan_model(
                    items_gg=items_gg, label_F="F", model=model)
    if (is.null(weights_gg)){
        weights_name <- NULL
    } else {
        dat_gg$weights <- weights_gg
        weights_name <- "weights"
    }

    #-- estimate lavaan model
    args <- list( model=lavmodel, std.lv=TRUE, meanstructure=TRUE, ...)
    if (is_data){
        args$data <- dat_gg
        args$sampling.weights <- weights_name
    } else {
        args$sample.cov <- Sigma
        args$sample.mean <- mu
        args$sample.nobs <- min(1e20,N)
    }
    mod <- do.call(what="sirt_import_lavaan_cfa", args=args)
    partable <- sirt_import_lavaan_parameterTable(object=mod)
    lambda <- partable[ partable$op=="=~", "est"]
    nu <- partable[ partable$op=="~1", "est"][1:I_gg]
    err_var <- partable[ partable$op=="~~", "est"][1:I_gg]
    nobs <- mod@Data@nobs[[1]]
    #--- output
    res <- list(lambda=lambda, nu=nu, err_var=err_var, nobs=nobs)
    return(res)
}
