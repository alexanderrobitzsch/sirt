## File Name: sirt_import_lavaan_standardizedSolution.R
## File Version: 0.05

sirt_import_lavaan_standardizedSolution <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::standardizedSolution(...)
    return(res)
}
