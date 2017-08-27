## File Name: sirtcat.R
## File Version: 0.04
## File Last Change: 2017-03-17 18:05:29

######################################################
sirtcat <- function( label , time0 , active ){
	if (active){
		z0 <- time0
		cat( label , "  " )
		z1 <- Sys.time()
		print(z1-z0)
		z0 <- z1 	
		zout <- z0
	} else {
		zout <- NULL
	}
	return(zout)
}
######################################################	
