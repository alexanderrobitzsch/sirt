## File Name: mgsem_partable2coef.R
## File Version: 0.05

mgsem_partable2coef <- function(partable)
{
    dfr <- partable
    ind <- which(dfr$unique==1)
    coef <- dfr[ ind,'est']
    names(coef) <- dfr[ind,'name']
    return(coef)
}
