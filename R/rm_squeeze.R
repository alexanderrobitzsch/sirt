## File Name: rm_squeeze.R
## File Version: 0.03

rm_squeeze <- function(x, lower, upper )
{
    x[ x < lower ] <- lower
    x[ x > upper ] <- upper
    return(x)
}
