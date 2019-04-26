## File Name: sirt_import_lavaan_cfa.R
## File Version: 0.03

sirt_import_lavaan_cfa <- function(...)
{
    res <- sirt_import_function_value(fun="cfa", pkg="lavaan", ...)
    return(res)
}
