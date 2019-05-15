## File Name: sirt_logit_to_probs.R
## File Version: 0.01


sirt_logit_to_probs <- function(y)
{
    K1 <- length(y)
    x <- rep(0,K1+1)
    x[1:K1] <- y
    x <- exp(x)
    x <- x / sum(x)
    return(x)
}
