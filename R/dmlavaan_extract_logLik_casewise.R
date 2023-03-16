## File Name: dmlavaan_extract_logLik_casewise.R
## File Version: 0.03

dmlavaan_extract_logLik_casewise <- function(mod)
{
    requireNamespace('lavaan')
    res <- unlist( lavaan::lavInspect(mod, what='loglik.casewise', list.by.group=FALSE) )
    return(res)
}
