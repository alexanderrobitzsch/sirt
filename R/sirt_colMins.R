## File Name: sirt_colMins.R
## File Version: 0.04
## File Last Change: 2018-12-30

sirt_colMins <- function(x, na.rm=TRUE)
{
    res <- apply(x, 2, min, na.rm=na.rm)
    return(res)
}

sirt_colMin <- sirt_colMins
