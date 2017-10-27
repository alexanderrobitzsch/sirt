## File Name: rm_facets_print_progress_parameter.R
## File Version: 0.01


rm_facets_print_progress_parameter <- function( parm, parm0, parmlabel, digits_parm )
{
	cat( paste( "    Maximum " , parmlabel , " parameter change = " , 
			paste( round(max(abs(parm0-parm)), digits_parm) , collapse=" " ) , "\n" , sep=""))	
}
