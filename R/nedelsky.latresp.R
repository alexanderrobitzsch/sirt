## File Name: nedelsky.latresp.R
## File Version: 0.03
## File Last Change: 2017-01-18 11:02:50

##########################################
# latent responses for Nedelksy function
nedelsky.latresp <- function(K){
	nodes <- c(0,1)
	ndim <- K
	combis <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, ndim) , 
            ncol = ndim ) ) ) )
	return(combis)
		}
