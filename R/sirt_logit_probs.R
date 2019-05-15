## File Name: sirt_logit_probs.R
## File Version: 0.01

sirt_logit_probs <- function(y)
{
    K <- length(y)
    y_logit <- log( y[-K] / y[K] )
    return(y_logit)
}
