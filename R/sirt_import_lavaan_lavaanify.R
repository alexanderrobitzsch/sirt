## File Name: sirt_import_lavaan_lavaanify.R
## File Version: 0.04

sirt_import_lavaan_lavaanify <- function(...)
{
    res <- sirt_import_function_value(fun="lavaanify", pkg="lavaan", ...)
    return(res)
}
