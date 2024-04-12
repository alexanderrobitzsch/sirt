## File Name: mgsem_bdiag.R
## File Version: 0.06

mgsem_bdiag <- function(x1, x2)
{
    vars <- c(rownames(x1),rownames(x2))
    n1 <- ncol(x1)
    n2 <- ncol(x2)
    mat <- matrix(0,nrow=n1+n2,ncol=n1+n2)
    rownames(mat) <- colnames(mat)
    mat[1L:n1,1L:n1] <- x1
    mat[n1+1L:n2,n1+1L:n2] <- x2
    return(mat)
}
