## File Name: xxirt_classprobs_lca.R
## File Version: 0.09

xxirt_classprobs_lca <- function(par, Theta, G)
{
    K <- nrow(Theta)
    probs <- matrix(NA, nrow=K, ncol=G)
    for (gg in 1L:G){
        logitprobs_gg <- par[ (K-1)*(gg-1) + ( 1L:(K-1) ) ]
        probs[,gg] <- xxirt_classprobs_lca_compute_probs(logitprobs=logitprobs_gg)
    }
    return(probs)
}
