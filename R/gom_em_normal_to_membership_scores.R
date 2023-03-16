## File Name: gom_em_normal_to_membership_scores.R
## File Version: 0.03


gom_em_normal_to_membership_scores <- function(theta_grid, K, TP)
{
    theta0 <- matrix(0, nrow=TP, ncol=K)
    for (kk in 1:(K-1)){
        theta0[,kk]    <- theta_grid[,kk]
    }
    theta0_rowmax <- rowMaxs.sirt(matr=theta0)$maxval
    theta0 <- theta0 - theta0_rowmax
    theta.k <- exp(theta0)
    theta.k <- theta.k / rowSums(theta.k)
    return(theta.k)
}
