## File Name: lsem_fitsem_sufficient_statistics_lavaan.R
## File Version: 0.04

lsem_fitsem_sufficient_statistics_lavaan <- function(gg, lavmodel, lavaan_est_fun,
    survey.fit, sample_stats, is_meanstructure, ...)
{
    wmean <- sample_stats$wmean
    wcov <- sample_stats$wcov
    Nobs <- sample_stats$Nobs
    if (gg==1){
        input_model <- lavmodel
    } else {
        input_model <- sirt_import_lavaan_parameterTable(survey.fit)
    }
    if (is_meanstructure){
        wmean_gg <- wmean[[gg]]
    } else {
        wmean_gg <- NULL
    }
    survey.fit <- lavaan_est_fun(input_model, sample.cov=wcov[[gg]],
                                sample.mean=wmean_gg, sample.nobs=Nobs[[gg]], ... )
    return(survey.fit)
}
