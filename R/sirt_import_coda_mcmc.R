## File Name: sirt_import_coda_mcmc.R
## File Version: 0.02

sirt_import_coda_mcmc <- function(...)
{
    TAM::require_namespace_msg('coda')
    res <- coda::mcmc(...)
    return(res)
}
