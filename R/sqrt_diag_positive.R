## File Name: sqrt_diag_positive.R
## File Version: 0.02

sqrt_diag_positive <- function(x)
{
    y <- diag(x)
    y <- ifelse(y<0,0,y)
    names(y) <- rownames(x)
    res <- sqrt(y)
    return(res)
}
