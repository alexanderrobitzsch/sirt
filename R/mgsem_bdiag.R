## File Name: mgsem_bdiag.R
## File Version: 0.05

mgsem_bdiag <- function(x1, x2)
{
    vars <- c(rownames(x1),rownames(x2))
    n1 <- ncol(x1)
    n2 <- ncol(x2)
    mat <- matrix(0,nrow=n1+n2,ncol=n1+n2)
    rownames(mat) <- colnames(mat)
    mat[1:n1,1:n1] <- x1
    mat[n1+1:n2,n1+1:n2] <- x2
    return(mat)
}
