## File Name: regpolca_postproc_prob_Theta.R
## File Version: 0.07
## File Last Change: 2020-02-24

regpolca_postproc_prob_Theta <- function(probs_Theta)
{
    K <- nrow(probs_Theta)
    G <- ncol(probs_Theta)
    rownames(probs_Theta) <- paste0("Class", 1:K)
    colnames(probs_Theta) <- paste0("Gr", 1:G)
    return(probs_Theta)
}
