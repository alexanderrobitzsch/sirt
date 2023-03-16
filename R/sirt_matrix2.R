## File Name: sirt_matrix2.R
## File Version: 0.04
## File Last Change: 2018-12-30

sirt_matrix2 <- function(x, nrow)
{
    matr <- matrix( x, nrow=nrow, ncol=length(x), byrow=TRUE )
    return(matr)
}
