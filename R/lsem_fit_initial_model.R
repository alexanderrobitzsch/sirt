## File Name: lsem_fit_initial_model.R
## File Version: 0.212

lsem_fit_initial_model <- function(lavmodel__, lavaan_est_fun, dat, variables_model,
    sampling_weights, has_meanstructure, sufficient_statistics, est_joint=FALSE,
    se="standard", use_lavaan_survey=FALSE, is_imputed=FALSE, Nimp=0, ...)
{
    if (est_joint){
        has_meanstructure <- TRUE
    }
    if (sufficient_statistics) {
        #- compute sufficient statistics
        res <- lsem_fit_initial_model_sufficient_statistics(dat=dat,
                    variables_model=variables_model, sampling_weights=sampling_weights,
                    has_meanstructure=has_meanstructure,
                    is_imputed=is_imputed, Nimp=Nimp)
        wmean <- res$wmean
        wcov <- res$wcov
        Nobs <- res$Nobs
        lavfit <- lavaan_est_fun(model=lavmodel__, sample.cov=wcov,
                    sample.nobs=Nobs, sample.mean=wmean, se=se, ...)
    } else {
        if (! use_lavaan_survey){
            lavfit <- lavaan_est_fun(model=lavmodel__, data=dat, se=se, ...)
        } else {
            lavfit <- lavaan_est_fun(model=lavmodel__, data=dat, ...)
        }
    }
    return(lavfit)
}
