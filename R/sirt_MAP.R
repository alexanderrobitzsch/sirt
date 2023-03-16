## File Name: sirt_MAP.R
## File Version: 0.01
## File Last Change: 2019-05-11

sirt_MAP <- function(post, theta)
{
    TP <- ncol(post)
    maxval <- post[,1]
    indval <- 1
    for (tt in 2:TP){
        m0 <- maxval
        maxval <- ifelse( post[,tt] > m0, post[,tt], m0 )
        indval <- ifelse( post[,tt] > m0, tt, indval )
    }
    MAP <- theta[ indval, ]
    return(MAP)
}
