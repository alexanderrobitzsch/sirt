## File Name: weighted_colMeans.R
## File Version: 0.04
## File Last Change: 2017-01-18 11:02:55

weighted_colMeans <- function( mat , wgt=NULL){
	wgt <- weighted_stats_extend_wgt( wgt=wgt , mat=mat )
	mat1 <- colSums( mat * wgt , na.rm=TRUE) 
	mat2 <- colSums( wgt , na.rm=TRUE) 
	mat1 <- mat1 / mat2
	return(mat1)
}
