## File Name: sirt_import_lavaan_lavaanify.R
## File Version: 0.05

sirt_import_lavaan_lavaanify <- function(...)
{
    TAM::require_namespace_msg("lavaan")
    res <- lavaan::lavaanify(...)
    return(res)
}
