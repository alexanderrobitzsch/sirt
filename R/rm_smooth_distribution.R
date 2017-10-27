## File Name: rm_smooth_distribution.R
## File Version: 0.14

rm_smooth_distribution <- function( theta.k , pi.k , est.mean = FALSE ,
			skillspace="normal" )
{
	m2 <- 0
	if ( est.mean ){
		m2 <- sum( theta.k * pi.k )
	}
	w2 <- sum( theta.k^2 * pi.k ) - m2^2
	sigma <- sqrt(w2)
	if ( skillspace == "normal" ){
		pi.k <- sirt_dnorm_discrete(x=theta.k, mean=m2, sd=sigma)				
	}
	res <- list( mu=m2 , sigma = sigma , pi.k=pi.k)
	return(res)
}

rm.smooth.distribution <- rm_smooth_distribution
