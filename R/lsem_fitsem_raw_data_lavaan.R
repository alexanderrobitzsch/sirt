## File Name: lsem_fitsem_raw_data_lavaan.R
## File Version: 0.136

lsem_fitsem_raw_data_lavaan <- function(dat, pseudo_weights, survey.fit,
        lavaan_est_fun, se, ...)
{
    #- define pseudo weights if requested
    res <- lsem_fitsem_raw_data_define_pseudo_weights(dat=dat,
                        pseudo_weights=pseudo_weights)
    dat1 <- res$dat
    sampling_weights <- res$sampling_weights  #=NULL for pseudo weights
    nobs_pseudo <- res$nobs_pseudo
    sum_weight <- res$sum_weight
    # use starting values
    partable <- sirt_import_lavaan_parameterTable(object=survey.fit)
    partable$start <- partable$est
    #- fit model with sampling weights
    survey.fit <- lavaan_est_fun(model=partable, data=dat1,
                        sampling.weights=sampling_weights, se=se, ... )

    #- adjust sample size
    survey.fit <- lavaan_object_adjust_sample_size(object=survey.fit, n_used=sum_weight)

    #-- output
    return(survey.fit)
}
