## File Name: sirt_colMeans.R
## File Version: 0.01

sirt_colMeans <- function(x, na.rm=TRUE)
{
    res <- colMeans(x, na.rm=na.rm)
    return(res)
}
