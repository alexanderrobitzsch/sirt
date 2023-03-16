## File Name: sirt_import_coda_effectiveSize.R
## File Version: 0.02
## File Last Change: 2023-03-08


sirt_import_coda_effectiveSize <- function(...)
{
    TAM::require_namespace_msg('coda')
    res <- coda::effectiveSize(...)
    return(res)
}
