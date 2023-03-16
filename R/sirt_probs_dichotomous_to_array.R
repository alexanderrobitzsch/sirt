## File Name: sirt_probs_dichotomous_to_array.R
## File Version: 0.02
## File Last Change: 2019-05-17

sirt_probs_dichotomous_to_array <- function(probs)
{
    I <- nrow(probs)
    TP <- ncol(probs)
    ncat <- 2
    probs1 <- array(0, dim=c(I,ncat,TP) )
    probs1[,2,] <- probs
    probs1[,1,] <- 1 - probs
    return(probs1)
}
