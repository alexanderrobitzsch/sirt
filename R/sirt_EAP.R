## File Name: sirt_EAP.R
## File Version: 0.02

sirt_EAP <- function(post, theta)
{
    TP <- ncol(post)
    N <- nrow(post)
    D <- ncol(theta)
    EAP <- matrix(NA, nrow=N, ncol=D)
    for (tt in 1L:D){
        thetaM <- sirt_matrix2(theta[,tt], nrow=N)
        EAP[,tt] <- rowSums(thetaM * post)
    }
    return(EAP)
}
