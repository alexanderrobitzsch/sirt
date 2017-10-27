## File Name: weighted_stats_extend_wgt.R
## File Version: 0.03

weighted_stats_extend_wgt <- function( wgt , mat )
{
	N1 <- nrow(mat)
	N2 <- ncol(mat)
	if ( is.null(wgt) ){
		wgt <- rep( 1 , N1 )
	}
	if ( is.vector(wgt) ){
		wgt <- matrix( wgt , nrow = N1 , ncol= N2 )	
	}	
	return(wgt)
}
