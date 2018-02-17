## File Name: mcmc_rename_parameter_names.R
## File Version: 0.04

mcmc_rename_parameter_names <- function( vec, orig, trans)
{
	NO <- length(orig)
	for (oo in 1:NO){
		trans_oo <- mcmc_rename_helper( string=trans[oo] )		
		vec <- gsub( orig[oo] , trans_oo , vec , fixed=TRUE)
	}
	return(vec)
}
