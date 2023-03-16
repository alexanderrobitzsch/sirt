## File Name: sirt_squeeze_probs.R
## File Version: 0.02
## File Last Change: 2019-09-12

sirt_squeeze_probs <- function(probs, eps)
{
    res <- ( probs + eps ) / ( 1 + 2*eps)
    return(res)
}
