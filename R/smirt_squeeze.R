## File Name: smirt_squeeze.R
## File Version: 0.03
## File Last Change: 2017-02-28 20:15:31

smirt_squeeze <- function( val , lower , upper, est)
{
	D <- 1
	is_matrix <- FALSE
	if ( is.matrix(est) ){
		D <- ncol(est)
		is_matrix <- TRUE
	}
	if( ! is.matrix(val) ){
		val <- matrix(val, ncol=1)
	}
	est <- matrix( est , ncol=D)
	val0 <- val
	for (dd in 1:D){
		val[,dd] <- ifelse( val[,dd] < lower , lower , val[,dd] )
		val[,dd] <- ifelse( val[,dd] > upper , upper , val[,dd] )
		ind_dd <- which(est[,dd] == 0)
		if ( length(ind_dd) > 0 ){
			val[ ind_dd , dd] <- val0[ ind_dd,dd]
		}
	}
	if ( ! is_matrix){
		val <- val[,1]
	}
	return(val)
}
