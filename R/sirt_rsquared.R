## File Name: sirt_rsquared.R
## File Version: 0.01

sirt_rsquared <- function(x, expl, na.rm=TRUE)
{
    res <- 1 - sum( ( x - expl)^2, na.rm=na.rm ) / sum( x^2, na.rm=na.rm)
    return(res)
}
