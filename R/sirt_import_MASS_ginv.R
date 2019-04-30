## File Name: sirt_import_MASS_ginv.R
## File Version: 0.12

sirt_import_MASS_ginv <- function(X,...)
{
    TAM::require_namespace_msg("MASS")
    y <- MASS::ginv(X=X, ...)
    return(y)
}
