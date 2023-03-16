## File Name: sirt_colMedians.R
## File Version: 0.04

sirt_colMedians <- function(x, na.rm=TRUE)
{
    res <- apply(x, 2, stats::median, na.rm=na.rm)
    return(res)
}

sirt_colMedian <- sirt_colMedians
