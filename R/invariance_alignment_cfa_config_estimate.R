## File Name: invariance_alignment_cfa_config_estimate.R
## File Version: 0.08

invariance_alignment_cfa_config_estimate <- function(dat_gg, weights_gg=NULL, ...)
{
    I_gg <- ncol(dat_gg)
    items_gg <- colnames(dat_gg)
    lavmodel <- paste0("F=~", paste0(items_gg, collapse="+") )
    if (is.null(weights_gg)){
        weights_name <- NULL
    } else {
        dat_gg$weights <- weights_gg
        weights_name <- "weights"
    }
    mod <- lavaan::cfa(data=dat_gg, model=lavmodel, std.lv=TRUE,
                meanstructure=TRUE, sampling.weights=weights_name, ...)
    partable <- lavaan::parameterTable(object=mod)
    lambda <- partable[ partable$op=="=~", "est"]
    nu <- partable[ partable$op=="~1", "est"][1:I_gg]
    err_var <- partable[ partable$op=="~~", "est"][1:I_gg]
    nobs <- mod@Data@nobs[[1]]
    #--- output
    res <- list(lambda=lambda, nu=nu, err_var=err_var, nobs=nobs)
    return(res)
}
