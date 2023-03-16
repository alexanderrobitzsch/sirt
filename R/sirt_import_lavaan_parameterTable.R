## File Name: sirt_import_lavaan_parameterTable.R
## File Version: 0.05
## File Last Change: 2023-03-08

sirt_import_lavaan_parameterTable <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::parameterTable(...)
    return(res)
}
