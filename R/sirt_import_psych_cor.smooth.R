## File Name: sirt_import_psych_cor.smooth.R
## File Version: 0.03

sirt_import_psych_cor.smooth <- function(x, ...)
{
    TAM::require_namespace_msg('psych')
    y <- psych::cor.smooth(x=x, ...)
    return(y)
}
