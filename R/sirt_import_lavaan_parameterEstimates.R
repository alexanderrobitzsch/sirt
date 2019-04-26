## File Name: sirt_import_lavaan_parameterEstimates.R
## File Version: 0.05

sirt_import_lavaan_parameterEstimates <- function(...)
{
    res <- sirt_import_function_value(fun="parameterEstimates", pkg="lavaan", ...)
    return(res)
}
