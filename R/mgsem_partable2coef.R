## File Name: mgsem_partable2coef.R
## File Version: 0.04
## File Last Change: 2022-02-25

mgsem_partable2coef <- function(partable)
{
    dfr <- partable
    ind <- which(dfr$unique==1)
    coef <- dfr[ ind,"est"]
    names(coef) <- dfr[ind,"name"]
    return(coef)
}
