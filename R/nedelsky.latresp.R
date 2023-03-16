## File Name: nedelsky.latresp.R
## File Version: 0.10
## File Last Change: 2018-12-30


#---- latent responses for Nedelksy function
nedelsky.latresp <- function(K)
{
    nodes <- c(0,1)
    ndim <- K
    combis <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim), ncol=ndim ))))
    return(combis)
}
