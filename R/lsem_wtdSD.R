## File Name: lsem_wtdSD.R
## File Version: 0.11


lsem_wtdSD <- function( x , w ){
	res1 <- sum( x*w )
	res2 <- sum( x^2*w)
	res <- sqrt( res2 - res1^2 )
	return(res)
}
