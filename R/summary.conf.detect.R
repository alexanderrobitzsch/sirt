## File Name: summary.conf.detect.R
## File Version: 0.04
## File Last Change: 2017-09-20 10:44:55
#*******************************************************
# Summary for conf.detect object
summary.conf.detect <- function( object , digits = 3 , file=NULL , ...){

    # open sink
    sirt_osink( file = file )

	cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
	itemcluster <- object$itemcluster
	IC <- length( unique(itemcluster) )
				
	des1 <- paste0( "Confirmatory DETECT Analysis with " , IC , " Item Clusters" )
	cat(des1,"\n")
    cat(paste("Bandwidth Scale:" , object$bwscale , "\n" ) ) 	

	cat("-----------------------------------------------------------------\n")
	cat("Dimensionality Statistics \n")
	obji <- object$detect.summary
	obji <- round( obji , digits)
	print( obji ) 

	# close sink
    sirt_csink( file = file )		
}
#*******************************************************
