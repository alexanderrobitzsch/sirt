## File Name: sirt_import_lavaan_parameterEstimates.R
## File Version: 0.07

sirt_import_lavaan_parameterEstimates <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::parameterEstimates(...)
    return(res)
}
