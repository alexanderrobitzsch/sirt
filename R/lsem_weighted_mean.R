## File Name: lsem_weighted_mean.R
## File Version: 0.03
## File Last Change: 2017-04-19 10:38:33

lsem_weighted_mean <- function( x , weights )
{
	x <- as.matrix(x)
	x_resp <- 1 - is.na(x)
	weights_m <- weights * x_resp
	x[ is.na(x) ] <- 0	
	weightsN <- colSums(weights_m)
	wm <- colSums( x * weights_m ) / ( weightsN - 1 )
	res <- list( weightsN = weightsN , mean=wm )
	return(res)
}
