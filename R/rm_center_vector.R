## File Name: rm_center_vector.R
## File Version: 0.04

rm_center_vector <- function( vec, center_type, do_log=FALSE )
{
	# log metric
	if (do_log){
		vec <- log(vec)
	}
	#--- center_type = 1
	if ( center_type == 1){
		RR <- length(vec)
		vec[RR] <- - sum( vec[-RR] )
	}
	#--- center_type = 2
	if ( center_type == 2){
		vec <- vec - mean(vec)
	}
	# reconvert to exp metric
	if (do_log){
		vec <- exp(vec)
	}	
	#--- output
	return(vec)
}
