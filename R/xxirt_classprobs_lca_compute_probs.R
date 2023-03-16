## File Name: xxirt_classprobs_lca_compute_probs.R
## File Version: 0.04
## File Last Change: 2020-04-27

xxirt_classprobs_lca_compute_probs <- function(logitprobs)
{
    v1 <- c( logitprobs, 0 )
    v1 <- v1 - max(v1)
    l1 <- exp(v1)
    res <- sirt_sum_norm(x=l1)
    return(res)
}
