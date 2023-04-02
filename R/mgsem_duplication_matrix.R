## File Name: mgsem_duplication_matrix.R
## File Version: 0.11

mgsem_duplication_matrix <- function(x)
{
    NV <- ncol(x)
    dfr1 <- data.frame( index1=rep(1:NV, NV), index2=rep(1:NV, each=NV))
    N1 <- nrow(dfr1)
    dfr2 <- dfr1[ dfr1[,1] >=dfr1[,2], ]
    ND <- nrow(dfr2)
    dupl <- matrix(0, nrow=ND, ncol=N1)
    dd <- 1
    for (dd in 1:ND){
        h1 <- dfr2$index1[dd]
        h2 <- dfr2$index2[dd]
        i1 <- which( ( dfr1$index1==h1 ) & ( dfr1$index2==h2 ) )
        i2 <- which( ( dfr1$index2==h1 ) & ( dfr1$index1==h2 ) )
        if (i1==i2){
            dupl[dd,i1] <- 1
        } else {
            dupl[dd,i1] <- 0.5
            dupl[dd,i2] <- 0.5
        }
    }
    return(dupl)
}
