## File Name: sirt_cormat2.R
## File Version: 0.01

sirt_cormat2 <- function(r)
{
    mat <- matrix(c(1, r, r, 1), 2, 2)
    return(mat)
}
