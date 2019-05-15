## File Name: sirt_add_increment.R
## File Version: 0.01

sirt_add_increment <- function(x, pos, value)
{
    y <- x
    y[pos] <- x[pos] + value
    return(y)
}
