## File Name: dmlavaan_sandwich_formula.R
## File Version: 0.06
## File Last Change: 2023-03-08

dmlavaan_sandwich_formula <- function(A, B, parnames=NULL)
{
    requireNamespace('MASS')
    B1 <- MASS::ginv(X=B)
    V <- B1 %*% A %*% B1
    if (!is.null(parnames)){
        rownames(V) <- colnames(V) <- parnames
    }
    se_sw <- sqrt( diag(V) )
    #-- output
    res <- list(V=V, se_sw=se_sw)
    return(res)
}
