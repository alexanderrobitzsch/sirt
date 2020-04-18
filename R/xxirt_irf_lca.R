## File Name: xxirt_irf_lca.R
## File Version: 0.071

xxirt_irf_lca <- function(par, Theta, ncat)
{
    K <- nrow(Theta)
    P <- matrix( NA, nrow=K, ncol=ncat)
    P[,1] <- 0
    for (hh in 2L:ncat){
        b <- par[ K*(hh-2) + seq(1,K) ]
        P[,hh] <- b
    }
    a1M <- matrix(apply(P, 1, max), nrow=K, ncol=ncat)
    P <- exp(P-a1M)
    P <- P / rowSums(P)
    return(P)
}
