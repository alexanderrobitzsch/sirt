## File Name: lsem_fitsem_raw_data_lavaan_survey.R
## File Version: 0.095

lsem_fitsem_raw_data_lavaan_survey <- function(dat, lavmodel, lavfit)
{
    TAM::require_namespace_msg('lavaan.survey')
    TAM::require_namespace_msg('survey')
    datsvy <- survey::svydesign(id=~index, weights=~weight, data=dat)
    assign_args <- list( x='lavmodel__', value=lavmodel, pos=1)
    res0 <- do.call( what='assign', args=assign_args)
    what <- paste0('lavaan.survey','::', 'lavaan.survey')
    args <- list(lavaan.fit=lavfit,    survey.design=datsvy )
    survey.fit <- do.call(what=what, args=args)
    # survey.fit <- lavaan.survey::lavaan.survey(lavaan.fit=lavfit,
    #                            survey.design=datsvy )
    return(survey.fit)
}
