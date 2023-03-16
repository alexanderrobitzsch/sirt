## File Name: sirt_sum_norm.R
## File Version: 0.04
## File Last Change: 2019-01-02

sirt_sum_norm <- function(x, na.rm=TRUE)
{
    y <- x/sum(x, na.rm=na.rm)
    return(y)
}
