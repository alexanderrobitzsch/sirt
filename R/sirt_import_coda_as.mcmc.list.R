## File Name: sirt_import_coda_as.mcmc.list.R
## File Version: 0.01

sirt_import_coda_as.mcmc.list <- function(...)
{
    TAM::require_namespace_msg("coda")
    res <- coda::as.mcmc.list(...)
    return(res)
}
