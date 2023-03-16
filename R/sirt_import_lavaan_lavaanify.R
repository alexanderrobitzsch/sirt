## File Name: sirt_import_lavaan_lavaanify.R
## File Version: 0.06
## File Last Change: 2023-03-08

sirt_import_lavaan_lavaanify <- function(...)
{
    TAM::require_namespace_msg('lavaan')
    res <- lavaan::lavaanify(...)
    return(res)
}
