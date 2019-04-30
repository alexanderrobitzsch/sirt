## File Name: sirt_import_psych_omega.R
## File Version: 0.01

sirt_import_psych_omega <- function(m, nfactors, ...)
{
    TAM::require_namespace_msg("psych")
    y <- psych::omega( m=m, nfactors=nfactors, ...)
    return(y)
}
