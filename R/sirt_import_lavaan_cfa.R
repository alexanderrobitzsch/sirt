## File Name: sirt_import_lavaan_cfa.R
## File Version: 0.05

sirt_import_lavaan_cfa <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::cfa(...)
    return(res)
}
