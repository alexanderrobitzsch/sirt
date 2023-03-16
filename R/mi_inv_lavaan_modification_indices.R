## File Name: mi_inv_lavaan_modification_indices.R
## File Version: 0.05
## File Last Change: 2022-04-17

mi_inv_lavaan_modification_indices <- function(mod, op=c("~1","=~"))
{
    requireNamespace("lavaan")
    res <- lavaan::modificationIndices(object=mod, free.remove=FALSE,
                    op=op, sort=TRUE)
    return(res)
}
