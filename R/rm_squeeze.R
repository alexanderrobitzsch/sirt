## File Name: rm_squeeze.R
## File Version: 0.03
## File Last Change: 2018-12-30

rm_squeeze <- function(x, lower, upper )
{
    x[ x < lower ] <- lower
    x[ x > upper ] <- upper
    return(x)
}
