## File Name: sirt_import_psych_omega.R
## File Version: 0.02
## File Last Change: 2023-03-08

sirt_import_psych_omega <- function(m, nfactors, ...)
{
    TAM::require_namespace_msg('psych')
    y <- psych::omega( m=m, nfactors=nfactors, ...)
    return(y)
}
