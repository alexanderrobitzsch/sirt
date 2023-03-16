## File Name: m_est_add_increment.R
## File Version: 0.01
## File Last Change: 2023-03-11

m_est_add_increment <- function(x, pos, h)
{
    y <- x
    y[pos] <- x[pos] + h
    return(y)
}
