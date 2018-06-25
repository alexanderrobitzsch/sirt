## File Name: sirt_matrix2.R
## File Version: 0.02

sirt_matrix2 <- function(x, nrow)
{
    matr <- matrix( x, nrow=nrow, ncol=length(x), byrow=TRUE )
    return(matr)
}
