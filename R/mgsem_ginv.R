## File Name: mgsem_ginv.R
## File Version: 0.04

mgsem_ginv <- function(X)
{
    requireNamespace('MASS')
    res <- MASS::ginv(X=X)
    rownames(res) <- colnames(res) <- colnames(X)
    return(res)
}
