## File Name: nedelsky.latresp.R
## File Version: 0.10


#---- latent responses for Nedelksy function
nedelsky.latresp <- function(K)
{
    nodes <- c(0,1)
    ndim <- K
    combis <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim), ncol=ndim ))))
    return(combis)
}
