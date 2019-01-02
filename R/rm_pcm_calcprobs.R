## File Name: rm_pcm_calcprobs.R
## File Version: 0.12



#- calculation of probabilities in the partial credit model
rm_pcm_calcprobs <- function( a, b, Qmatrix, theta.k, I, K, TP )
{
    probs <- array( 0, dim=c(I,K+1,TP) )   # categories 0, ..., K
    for (kk in 1:K){
        l0 <- matrix( - b[,kk], nrow=I,ncol=TP)
        l0 <- l0 + TAM::tam_outer( a * Qmatrix[, kk], theta.k )
        probs[,kk+1,] <- l0
    }
    probs <- exp(probs)
    probs1 <- probs[,1,]
    for (kk in 2:(K+1)){
        probs1 <- probs1 + probs[,kk,]
    }
    for (kk in 1:(K+1)){
        probs[,kk,] <- probs[,kk,] / probs1
    }
   return(probs)
}
