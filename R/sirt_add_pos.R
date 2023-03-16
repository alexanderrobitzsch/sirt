## File Name: sirt_add_pos.R
## File Version: 0.01
## File Last Change: 2019-01-05

sirt_add_pos <- function(x, pos, val)
{
    y <- x
    y[pos] <- y[pos] + val
    return(y)
}
