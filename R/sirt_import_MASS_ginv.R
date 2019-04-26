## File Name: sirt_import_MASS_ginv.R
## File Version: 0.08

sirt_import_MASS_ginv <- function(...)
{
    res <- sirt_import_function_value(fun="ginv", pkg="MASS", ...)
    return(res)
}
