## File Name: sirt_import_psych_fa.R
## File Version: 0.03

sirt_import_psych_fa <- function(r, nfactors, ...)
{
    TAM::require_namespace_msg('psych')
    y <- psych::fa( r=r, nfactors=nfactors, ...)
    return(y)
}
