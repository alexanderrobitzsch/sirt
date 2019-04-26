## File Name: sirt_import_lavaan_parameterTable.R
## File Version: 0.03

sirt_import_lavaan_parameterTable <- function(...)
{
    res <- sirt_import_function_value(fun="parameterTable", pkg="lavaan", ...)
    return(res)
}
