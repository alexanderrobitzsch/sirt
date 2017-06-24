

#########################################################################
# simulation of Rasch type models
sim.raschtype <- function( theta , b , alpha1 = 0, alpha2 = 0 , fixed.a = NULL , 
		fixed.c = NULL , fixed.d = NULL )
{ 
    if ( is.null(fixed.a)){ fixed.a <- 1+0*b }
    if ( is.null(fixed.c)){ fixed.c <- 0*b }
    if ( is.null(fixed.d)){ fixed.d <- 1 + 0*b}
    # latent response (subtraction)
    latresp <- outer( theta , b , "-" )
    # include slope simulation
    latresp <- outer( rep(1,length(theta)) , fixed.a ) * latresp 
    # transformed response
    cM <- outer( rep(1,length(theta)) , fixed.c )
    dM <- outer( rep(1,length(theta)) , fixed.d )    
    trlat <- pgenlogis( latresp , alpha1 = alpha1 , alpha2 = alpha2 )
    trlat <- cM + ( dM - cM )*trlat
    expprob <- trlat
    # define response matrix
    dat.resp <- 1 * ( expprob > matrix( stats::runif( nrow(expprob)*ncol(expprob) ) , 
							ncol= ncol(expprob )) )
	I <- length(b)
	pot <- max( 2 , log10(I) + 1 )
	colnames(dat.resp) <- paste( "I" , substring( 10^pot + 1:I,2) , sep="")
	dat.resp <- as.data.frame(dat.resp)
    return(dat.resp)
}
#########################################################################

