## File Name: lsem_fit_initial_model.R
## File Version: 0.18

lsem_fit_initial_model <- function(lavmodel__, lavaan_est_fun, dat, variables_model,
    sampling_weights, has_meanstructure, sufficient_statistics, est_joint=FALSE,
    se="standard", use_lavaan_survey=FALSE, ...)
{
    if (sufficient_statistics) {
        if (est_joint){
            has_meanstructure <- TRUE
        }
        #- compute sufficient statistics
        res <- lsem_fit_initial_model_sufficient_statistics(dat=dat,
                    variables_model=variables_model, sampling_weights=sampling_weights,
                    has_meanstructure=has_meanstructure)
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
