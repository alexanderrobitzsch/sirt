## File Name: sirt_probs_to_logit.R
## File Version: 0.02

sirt_probs_to_logit <- function(y)
{
    K <- length(y)
    y_logit <- log( y[-K] / y[K] )
    return(y_logit)
}
