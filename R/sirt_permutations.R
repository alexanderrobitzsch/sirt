## File Name: sirt_permutations.R
## File Version: 0.01


sirt_permutations <- function(r,v)
{
    NL <- length(v)
    NC <- NL^r
    mat <- matrix(0, nrow=NC, ncol=r)
    hh <- 1
    for (dd in seq(r,1,by=-1)){
        m1 <- rep(v, each=NL^(hh-1) )
        m1 <- rep(m1, NC/length(m1) )
        mat[,dd] <- m1        
        hh <- hh + 1
    }
    return(mat)
}
