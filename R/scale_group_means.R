## File Name: scale_group_means.R
## File Version: 0.04

scale_group_means <- function(M, SD, probs=NULL)
{
    if (is.null(probs)){
        probs <- rep(1,length(M))
    }
    probs <- probs/sum(probs)
    M1 <- sum(probs*M)
    SD1 <- sqrt( sum(probs*SD^2)+sum(probs*(M-M1)^2) )
    M_z <- ( M - M1 ) / SD1
    SD_z <- SD / SD1
    #--- output
    res <- list(M1=M1, SD1=SD1, M_z=M_z, SD_z=SD_z)
    return(res)
}
