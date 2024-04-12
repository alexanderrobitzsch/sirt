## File Name: mgsem_create_index.R
## File Version: 0.04


mgsem_create_index <- function(x, all=TRUE, start=0, symm=FALSE, onlydiag=FALSE)
{
    ND <- prod(dim(x))
    if (start>0){
        v <- matrix( start + ( 1L:ND ) - 1, nrow=dim(x)[1], ncol=dim(x)[2])
        if (symm){
            v <- v + t(v)
            v <- as.vector(v)
            g1 <- match( v, unique(v) )
            v <- matrix( g1 + start - 1, nrow=dim(x)[1], ncol=dim(x)[2] )
        }
    } else {
        v <- 1+0*x
        if (!all){
            v <- 0*v
        }
    }
    #*** only diagonal entries
    if (onlydiag){
        ND <- dim(x)[1]
        if (start>0){
            v <- diag( start - 1 + 1L:ND )
        } else {
            v <- diag( rep(1, dim(x)[1]))
        }
    }
    return(v)
}
