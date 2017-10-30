## File Name: rm_proc_fixed_values_reference_rater.R
## File Version: 0.03

rm_proc_fixed_values_reference_rater <- function( rater.index1, b.rater.fixed, a.rater.fixed, rater_item_int,
		reference_rater, est.b.rater, est.a.rater )
{	
	if ( ! is.null(reference_rater) ){
		#----- no item-rater-interactions
		if ( ! rater_item_int ){
			RR <- nrow(rater.index1)
			#** b parameters
			if (est.b.rater){
				if ( is.null(b.rater.fixed) ){
					val <- 0
					b.rater.fixed <- rep(NA, RR)
					b.rater.fixed[ rater.index1[ paste(rater.index1$rater) == reference_rater , "rater.id" ] ] <- val
				}
			}
			#** a parameters
			if (est.a.rater){
				if ( is.null(a.rater.fixed) ){
					val <- 1
					a.rater.fixed <- rep(NA, RR)
					a.rater.fixed[ rater.index1[ paste(rater.index1$rater) == reference_rater , "rater.id" ] ] <- val
				}
			}		
			
		}
		#---- item-rater-interactions
		if ( rater_item_int ){
			RR <- nrow(rater.index1)
			ind <- match( reference_rater, paste(rater.index1$rater) )
			#** b parameters
			if (est.b.rater){
				val <- 0
				if ( is.null(b.rater.fixed) ){
					b.rater.fixed <- rep(NA, RR)					
					b.rater.fixed[ ind ] <- val
				}
			}
			#** a parameters
			if (est.a.rater){
				if ( is.null(a.rater.fixed) ){
					val <- 1
					a.rater.fixed <- rep(NA, RR)
					a.rater.fixed[ ind ] <- val
				}
			}		
		}
	}
	#--- output
	res <- list( a.rater.fixed=a.rater.fixed, b.rater.fixed=b.rater.fixed )
	return(res)
}
