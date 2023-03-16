## File Name: sirt_add_pos.R
## File Version: 0.01

sirt_add_pos <- function(x, pos, val)
{
    y <- x
    y[pos] <- y[pos] + val
    return(y)
}
