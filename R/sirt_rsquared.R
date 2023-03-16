## File Name: sirt_rsquared.R
## File Version: 0.03
## File Last Change: 2018-12-30

sirt_rsquared <- function(x, expl, na.rm=TRUE)
{
    res <- 1 - sum( ( x - expl)^2, na.rm=na.rm ) / sum( x^2, na.rm=na.rm)
    return(res)
}
