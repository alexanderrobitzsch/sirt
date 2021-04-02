## File Name: rasch_pairwise_zerosum.R
## File Version: 0.02


rasch_pairwise_zerosum <- function(eps)
{
    b1 <- - log(eps)
    b2 <- b1 - mean(b1)
    eps <- exp(-b2)
    #- output
    return(eps)
}
