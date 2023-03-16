## File Name: mgsem_coef2partable.R
## File Version: 0.01
## File Last Change: 2022-01-19

mgsem_coef2partable <- function(coef, partable)
{
    dfr <- partable
    dfr$est <- coef[ dfr$index ]
    return(dfr)
}
