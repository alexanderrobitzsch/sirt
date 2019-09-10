## File Name: squeeze_probs.R
## File Version: 0.01

squeeze_probs <- function(probs, eps)
{
    res <- ( probs + eps ) / ( 1 + 2*eps)
    return(res)
}
