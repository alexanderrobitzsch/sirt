## File Name: sirt_import_lavaan_parameterTable.R
## File Version: 0.05

sirt_import_lavaan_parameterTable <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::parameterTable(...)
    return(res)
}
