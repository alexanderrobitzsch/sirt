## File Name: sirt_import_coda_mcmc.R
## File Version: 0.02
## File Last Change: 2023-03-08

sirt_import_coda_mcmc <- function(...)
{
    TAM::require_namespace_msg('coda')
    res <- coda::mcmc(...)
    return(res)
}
