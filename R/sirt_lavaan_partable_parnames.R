## File Name: sirt_lavaan_partable_parnames.R
## File Version: 0.02
## File Last Change: 2019-04-26

sirt_lavaan_partable_parnames <- function(partable)
{
    res <- paste0( partable$lhs, partable$op, partable$rhs )
    return(res)
}
