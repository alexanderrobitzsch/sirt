## File Name: mgsem_compute_model_implied_moments.R
## File Version: 0.165


mgsem_compute_model_implied_moments <- function(est, is_B=FALSE, calc_Sigma=TRUE,
    calc_Mu=TRUE)
{
    Mu <- NULL
    Sigma <- NULL
    if (is_B){
        requireNamespace("MASS")
        D <- ncol(est$PHI)
        ID <- diag(rep(1,D))
        B <- res$B
        B1 <- MASS::ginv(ID-B)
        LAMB <- est$LAM %*% B1
        Mu <- LAMB %*% est$ALPHA + est$NU
        # Sigma <- LAMB %*% tcrossprod(est$PHI, LAMB) + est$PSI
        if (calc_Sigma){
            Sigma <- sirt_rcpp_mgsem_compute_cov(LAM=LAMB, PHI=est$PHI, PSI=est$PSI)
        }
    } else {
        if (calc_Mu){
            Mu <- est$LAM %*% est$ALPHA + est$NU
        }
        # Sigma <- est$LAM %*% tcrossprod(est$PHI, est$LAM) + est$PSI
        if (calc_Sigma){
            Sigma <- sirt_rcpp_mgsem_compute_cov(LAM=est$LAM, PHI=est$PHI, PSI=est$PSI)
        }
    }
    #*** output
    res <- list(Mu=Mu, Sigma=Sigma)
    return(res)
}
