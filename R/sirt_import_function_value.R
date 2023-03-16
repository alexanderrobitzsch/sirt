## File Name: sirt_import_function_value.R
## File Version: 0.091

sirt_import_function_value <- function(fun, pkg, ...)
{
    TAM::require_namespace_msg(pkg)
    fun1 <- NULL
    eval(parse(text=paste0('fun1 <- ', pkg,'::', fun)))
    res <- fun1(...)
    return(res)
}
