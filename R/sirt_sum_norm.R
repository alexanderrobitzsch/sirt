## File Name: sirt_sum_norm.R
## File Version: 0.04

sirt_sum_norm <- function(x, na.rm=TRUE)
{
    y <- x/sum(x, na.rm=na.rm)
    return(y)
}
