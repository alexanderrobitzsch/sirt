## File Name: sirt_import_lavaan_fitMeasures.R
## File Version: 0.05

sirt_import_lavaan_fitMeasures <- function(...)
{
    TAM::require_namespace_msg("lavaan")
    res <- lavaan::fitMeasures(...)
    return(res)
}
