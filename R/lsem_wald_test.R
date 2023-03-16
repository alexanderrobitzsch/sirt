## File Name: lsem_wald_test.R
## File Version: 0.04
## File Last Change: 2023-03-15

lsem_wald_test <- function(theta, V, A)
{
    requireNamespace('MASS')
    r <- ( A %*% theta )[,1]
    W <- A %*% V %*% t(A)
    W1 <- MASS::ginv( X=W )
    chisq <- ( t(r) %*% W1 %*% r )[1,1]
    df <- nrow(A)
    p <- 1 - stats::pchisq(q=chisq, df=df)

    #-- output
    res <- list(chisq=chisq, df=df, p=p)
    return(res)
}
