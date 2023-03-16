## File Name: noharm_sirt_efa_rotated_solution.R
## File Version: 0.05
## File Last Change: 2018-12-30


noharm_sirt_efa_rotated_solution <- function(res, items, F_dimnames)
{
    L1 <- res$loadings
    I <- nrow(L1)
    D <- ncol(L1)
    m1 <- stats::promax(L1)
    p1 <- matrix( 0, nrow=I, ncol=D)
    for (dd in 1:D){
        p1[,dd] <- m1$loadings[,dd]
    }
    colnames(p1) <- F_dimnames
    rownames(p1) <- items
    res$promax <- p1
    res$factor.cor <- solve( crossprod( m1$rotmat ) )
    rownames(res$factor.cor) <- colnames(res$factor.cor) <- F_dimnames
    # conversion to THETA parametrization
    h1 <- rowSums(p1^2)
    p2 <- p1 / sqrt( max(1 - h1, 1e-4 ) )
    res$promax.theta <- p2
    #--- output
    return(res)
}
