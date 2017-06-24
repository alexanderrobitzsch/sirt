
lsem_weighted_cov <- function( x , weights )
{
	x <- as.matrix(x)
	x_resp <- 1 - is.na(x)
	eps <- 1E-100 * max(weights)
	weights_m <- sqrt( weights + eps ) * x_resp
	x[ is.na(x) ] <- 0	
	x_center <- lsem_weighted_mean( x=x , weights=weights )$mean
	XC <- matrix( x_center , nrow=nrow(x) , ncol=ncol(x) , byrow=TRUE )
	x <- x - XC
	weightsN <- crossprod(weights_m)	
	xw <- as.matrix( x * weights_m )
	covw <- crossprod(xw) / weightsN
	Nobs <- mean( weightsN[ ! upper.tri(weightsN) ] )
	res <- list( weightsN = weightsN , cov = covw , Nobs = Nobs )
	return(res)
}